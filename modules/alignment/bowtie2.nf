// Bowtie2 alignment process for ChIP-seq data
// Post-alignment steps (all within one process to avoid redundancy and intermediate files):
//   1. bowtie2       – paired-end alignment
//   2. samtools view – MAPQ filter + BAM conversion
//   3. samtools sort – coordinate sort (required by markdup)
//   4. samtools markdup – mark PCR duplicates (reads remain; marked in FLAG)
//   5. samtools index   – index the final BAM
//   6. samtools flagstat – alignment + duplicate stats for MultiQC
//
// Output BAMs are duplicate-MARKED (not removed), so MACS2 --keep-dup auto
// can still apply its own statistical duplicate handling, while MultiQC picks
// up the per-sample duplicate rate from the flagstat file.

process ALIGN_READS {
    // Tag the process with the sample ID for easier tracking in logs and reports
    tag "${meta.id}"
    label 'process_high'

    // Publish aligned BAM files to a specific directory for alignment results
    publishDir "${params.output_dir}/Aligned", mode: 'copy'

    // Use a specific Conda environment that includes Bowtie2 and Samtools for alignment and BAM processing
    conda "bioconda::bowtie2=2.5.1 bioconda::samtools=1.18"

    // Input: a tuple of metadata and a list of FASTQ file paths, along with the Bowtie2 index path
    input:
    tuple val(meta), path(fastq_files)
    val bowtie2_index

    // The published BAM is the duplicate-marked, coordinate-sorted file.
    // flagstat now contains duplicate stats readable by MultiQC.
    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"),   emit: bam
    path "${meta.id}_bowtie2.log",                                         emit: log
    path "${meta.id}_flagstat.txt",                                        emit: flagstat
    path "${meta.id}_markdup_stats.txt",                                   emit: markdup_stats

    // Script to run Bowtie2 for alignment, followed by Samtools for sorting and indexing the BAM file
    script:
    def (r1, r2) = fastq_files
    """
    # 1. Align, filter by MAPQ, name-sort (required by fixmate)
    bowtie2 \\
        --very-sensitive \\
        --no-mixed \\
        --no-discordant \\
        -X ${params.max_insert_size} \\
        -x ${bowtie2_index} \\
        -1 ${r1} \\
        -2 ${r2} \\
        -p ${task.cpus} \\
        2> ${meta.id}_bowtie2.log \\
    | samtools view -b -q ${params.min_mapq} - \\
    | samtools sort -n -@ ${task.cpus} -o ${meta.id}_namesorted.bam

    # 2. Run fixmate to add MC and ms tags (required by markdup)
    samtools fixmate \\
        -m \\
        -@ ${task.cpus} \\
        ${meta.id}_namesorted.bam \\
        ${meta.id}_fixmate.bam

    # 3. Coordinate-sort (required by markdup)
    samtools sort -@ ${task.cpus} -o ${meta.id}_sorted.bam ${meta.id}_fixmate.bam

    # 4. Mark PCR duplicates
    samtools markdup \\
        -@ ${task.cpus} \\
        -s \\
        -f ${meta.id}_markdup_stats.txt \\
        ${meta.id}_sorted.bam \\
        ${meta.id}.bam

    # 5. Index and flagstat
    samtools index -@ ${task.cpus} ${meta.id}.bam
    samtools flagstat ${meta.id}.bam > ${meta.id}_flagstat.txt

    # 6. Clean up intermediates
    rm ${meta.id}_namesorted.bam ${meta.id}_fixmate.bam ${meta.id}_sorted.bam
    """
}
