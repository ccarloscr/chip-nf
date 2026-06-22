#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Main Nextflow pipeline for ChIP-seq data processing
// This pipeline includes the following steps:
// 0. Parsing the metadata CSV to prepare input channels
// 1. Quality control with FastQC
// 2. MultiQC report generation
// 3. Read alignment with Bowtie2
// 4. Peak calling with MACS2
// 5. Peak annotation with HOMER
// 6. Visualization: PCA, volcano plot, heatmap, genomic distribution, distance-to-gene


// Include process definitions from separate module files for better organization and reusability
include { validateParameters        } from 'plugin/nf-schema'
include { FASTQC as FASTQC_RAW      } from './modules/qc/fastqc'
include { FASTP                     } from './modules/qc/fastp'
include { FASTQC as FASTQC_TRIMMED  } from './modules/qc/fastqc'
include { ALIGN_READS               } from './modules/alignment/bowtie2'
include { CALL_PEAKS_BROAD          } from './modules/peaks/macs2_broad'
include { CALL_PEAKS_NARROW         } from './modules/peaks/macs2_narrow'
include { DESEQ2_PEAKS              } from './modules/peaks/deseq2_peaks'
include { MERGE_PEAKS               } from './modules/peaks/merge_peaks'
include { HOMER_SETUP               } from './modules/annotation/homer_setup'
include { ANNOTATE_PEAKS            } from './modules/annotation/homer_annot'
include { PCA_PEAKS                 } from './modules/visualization/pca_peaks'
include { VOLCANO_PLOT              } from './modules/visualization/volcano_plot'
include { HEATMAP_PLOT		        } from './modules/visualization/plot_heatmap'
include { GENOMIC_DISTRIBUTION      } from './modules/visualization/genomic_distribution'
include { MULTIQC                   } from './modules/qc/multiqc'


// Sub-workflow to parse the metadata CSV and prepare input channels for the main workflow
workflow PARSE_SAMPLESHEET {
    // The metadata CSV is expected to have columns: sample, fastq_1, fastq_2, control (optional)
    take: csv_file

    // Create a channel that emits tuples of (metadata, fastq_files) for each sample
    main:
    ch_reads = Channel
        .fromPath(csv_file)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id         : row.sample,                                        // Unique sample identifier
                control    : row.control ?: null,                               // Control sample ID if available, else null
                is_control : (row.control == null || row.control.trim() == '')  // Indicates if a sample is control (input)
            ]
            // Prepare the fastq files as a tuple of paths for paired-end data
            def fastq_files = [ file(row.fastq_1, checkIfExists: true),
                                 file(row.fastq_2, checkIfExists: true) ]
            return [ meta, fastq_files ]
        }

    // Emit the channel for use in downstream processes
    emit:
    reads = ch_reads
}



// Main workflow that orchestrates the entire ChIP-seq processing pipeline
workflow {

    // --- Schema validation ---
    validateParameters()

    // --- Parse samplesheet ---
    PARSE_SAMPLESHEET(params.metadata)


    // --- Quality control with FastQC on raw reads ---
    FASTQC_RAW(PARSE_SAMPLESHEET.out.reads)

    // --- Adapter trimming and quality filtering with fastp ---
    FASTP(PARSE_SAMPLESHEET.out.reads)

    // --- Quality control with FastQC on trimmed reads ---
    FASTQC_TRIMMED(FASTP.out.reads)


    
    // --- Alignment ---
    
    // 1. Load channel to bowtie2 index
    index_ch        = Channel.value(params.bowtie2_index)

    // 2. Bowtie2 alignment to the reference genome, producing sorted BAM files
    ALIGN_READS(FASTP.out.reads, index_ch)


    // --- Prepare input for peak calling ---
    // Define ChIP BAMs by filtering the ALIGN_READS output for non-control samples
    chip_bams = ALIGN_READS.out.bam
        .filter { meta, bam, bai -> !meta.is_control }

    // Define control BAMs by filtering the ALIGN_READS output for control samples
    ctrl_bams = ALIGN_READS.out.bam
        .filter { meta, bam, bai -> meta.is_control }
        .map    { meta, bam, bai -> [ meta.id, bam, bai ] }  // key by id for joining

    // Join ChIP and control BAMs using the metadata to prepare input for peak calling
    peak_input_ch = chip_bams
        .map    { meta, bam, bai -> [ meta.control, meta, bam, bai ] }      // key by control id
        .combine(ctrl_bams, by: 0)                                          // join on control id
        .map    { ctrl_id, meta, chip_bam, chip_bai, ctrl_bam, ctrl_bai ->
                  [ meta, chip_bam, chip_bai, ctrl_bam, ctrl_bai ]
                }


    // --- Peak calling with MACS2 using the prepared input channel ---
    if (params.peak_mode == 'broad') {
        CALL_PEAKS_BROAD(peak_input_ch)
        peaks_ch = CALL_PEAKS_BROAD.out.peaks
    }
    else if (params.peak_mode == 'narrow') {
        CALL_PEAKS_NARROW(peak_input_ch)
        peaks_ch = CALL_PEAKS_NARROW.out.peaks
    }
    else {
        error "Invalid peak_mode: ${params.peak_mode}. Use 'narrow' or 'broad'."
    }


    // --- HOMER genome setup (runs once; cached by Nextflow on -resume) ---
    HOMER_SETUP()

    // --- Annotate peaks with HOMER using the output from MACS2 ---
    ANNOTATE_PEAKS(peaks_ch, file(params.gtf_file), HOMER_SETUP.out.homer_data_dir)


    // --- Generation of a master consensus peak set for downstream processes ---

    // Collect all peak files into a single list, stripping the meta map
    all_peaks_ch = peaks_ch
        .map { meta, peaks -> peaks }
        .collect()

    // Run the process to get a consensus list of peaks
    MERGE_PEAKS(all_peaks_ch) 


    // --- Principal Component Analysis (PCA) ---

    // 1. Separate and collect all bam, bai files and sample IDs (excluding inputs).
    //    Sort by meta.id so BAMs and IDs are in the same deterministic order,
    //    which is critical for PCA_PEAKS to name tag directories correctly.
    all_chip_bams_ch    = chip_bams.toSortedList { a, b -> a[0].id <=> b[0].id }
                            .flatMap { it }
                            .map { meta, bam, bai -> bam }.collect()
    all_chip_bais_ch    = chip_bams.toSortedList { a, b -> a[0].id <=> b[0].id }
                            .flatMap { it }
                            .map { meta, bam, bai -> bai }.collect()
    chip_sample_ids_ch  = chip_bams.toSortedList { a, b -> a[0].id <=> b[0].id }
                            .flatMap { it }
                            .map { meta, bam, bai -> meta.id }.collect()

    // 2. Run the PCA process using the consensus merged peaks and BAM alignment data
    PCA_PEAKS(
        MERGE_PEAKS.out,
        all_chip_bams_ch,
        all_chip_bais_ch,
        chip_sample_ids_ch,
        file(params.gtf_file),
        HOMER_SETUP.out.homer_data_dir
    )


    // --- Heatmap ---

    HEATMAP_PLOT(
        MERGE_PEAKS.out,
        all_chip_bams_ch,
        all_chip_bais_ch,
        chip_sample_ids_ch
    )



    // -- Differential peak analysis and volcano plot ---
 
    // 1. DESeq2 differential analysis on the consensus peak count matrix
    DESEQ2_PEAKS(PCA_PEAKS.out.matrix, chip_sample_ids_ch)

    // 2. Volcano plot of significant differential peaks between conditions
    VOLCANO_PLOT(DESEQ2_PEAKS.out.results)


    // --- Genomic Distribution: stacked bar plot showing mean % of all replicates per condition
   
    // 1. Extract the condition for each replicate
    all_annotations_ch = ANNOTATE_PEAKS.out.annotated
        .map { meta, tsv -> 
            // Extracting the condition from the ID: WT_rep1 --> WT
            def condition = meta.id.split('_')[0] 
            // Define a new meta containing the condition
            def meta_cond = meta + [condition: condition]
            return [ meta_cond, tsv ]
        }

    // 2. Define a channel for collected conditions, and for collectes tsv files
    conditions_ch = all_annotations_ch.map { meta_cond, tsv -> meta_cond.condition }.collect()
    tsv_files_ch  = all_annotations_ch.map { cond, tsv -> tsv }.collect()

    // 3. Proportion of peaks in promoters vs intergenic regions (from HOMER annotation)
    GENOMIC_DISTRIBUTION(conditions_ch, tsv_files_ch)


    // --- Aggregation of QC reports with MultiQC ---
    MULTIQC(
        FASTQC_RAW.out.reports.map { meta, reports -> reports }.collect(),
        FASTQC_TRIMMED.out.reports.map { meta, reports -> reports }.collect(),
        FASTP.out.json.map { meta, json -> json }.collect(),
        ALIGN_READS.out.log.collect(),
        ALIGN_READS.out.flagstat.collect(),
        ALIGN_READS.out.markdup_stats.collect(),
        params.peak_mode == 'broad' ? CALL_PEAKS_BROAD.out.xls.collect() : CALL_PEAKS_NARROW.out.xls.collect()
    )


}
