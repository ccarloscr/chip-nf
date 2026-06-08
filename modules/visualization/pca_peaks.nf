// PCA of ChIP-seq signal across all replicates
//
// Strategy:
//   1. Build per-sample HOMER tag directories from the ChIP BAMs.
//   2. Quantify signal at consensus peaks across all tag directories
//      (annotatePeaks.pl -d) → raw_matrix.txt
//   3. Clean the header so column names are bare sample IDs → homer_matrix.txt
//   4. PCA in R (bin/pca_merged_peaks.R)
//
// Receives homer_data_dir from HOMER_SETUP — no network access on compute nodes.

process PCA_PEAKS {
    tag "consensus"
    label 'process_high'
    publishDir "${params.output_dir}/Visualization/PCA", mode: 'copy'
    conda "bioconda::homer=4.11 conda-forge::perl conda-forge::r-base=4.3 conda-forge::r-ggplot2 conda-forge::r-ggrepel"

    input:
    path consensus_peaks
    path bam_files, stageAs: "bams/*.bam"
    path bai_files, stageAs: "bams/*.bai"
    val  sample_ids           // ordered list matching bam_files order
    path gtf_file 
    path homer_data_dir   // pre-built genome data dir from HOMER_SETUP

    output:
    path "pca_peaks.png", emit: plot
    path "homer_matrix.txt", emit: matrix

    script:
    // Build a bash array literal of sample IDs in the same order as the staged BAMs
    def id_array = sample_ids.collect { "\"${it}\"" }.join(' ')
    """
    export PATH=\$CONDA_PREFIX/share/homer/bin:\$PATH

    # 1. Point HOMER at the pre-built data directory
    export HOMER_DATA_DIR=\$(realpath ${homer_data_dir})

    # 2. Build a HOMER tag directory for each BAM, using the real sample ID as the name
    mkdir -p tagdirs
    bam_array=(bams/*.bam)
    id_array=(${id_array})
    for i in "\${!bam_array[@]}"; do
        makeTagDirectory tagdirs/\${id_array[\$i]} \${bam_array[\$i]}
    done

    # 3. Annotate consensus peaks and quantify across all tag directories
    annotatePeaks.pl ${consensus_peaks} ${params.genome} \\
        -gtf ${gtf_file} \\
        -cpu ${task.cpus} \\
        -d tagdirs/* \\
        -raw \\
        > raw_matrix.txt

    # 4. Clean header: HOMER writes "tagdirs/<sample> Tag Count" per sample column;
    #    use awk to strip path prefix and suffix field-by-field.
    sed '1s|tagdirs/||g; 1s| Tag Count||g' raw_matrix.txt > homer_matrix.txt

    # 5. PCA plot
    Rscript ${projectDir}/bin/pca_merged_peaks.R homer_matrix.txt
    """
}
