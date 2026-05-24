// FastQC process definition for quality control of sequencing reads
process FASTQC {
    // Tag the process with the sample ID for easier tracking in logs and reports
    tag "${meta.id}"
    label 'process_high'

    // Publish FastQC reports to a specific directory for QC results
    publishDir "${params.output_dir}/QC/FastQC", mode: 'copy'

    // Use a specific Conda environment for FastQC
    conda "bioconda::fastqc=0.12.1"

    // Input: a tuple of metadata and a list of FASTQ file paths
    input:
    tuple val(meta), path(fastq_files)

    // Output: a tuple of metadata and the generated FastQC report files (zip and html)
    output:
    tuple val(meta), path("*.{zip,html}"), emit: reports

    // Script to run FastQC on the input FASTQ files using the specified number of threads
    script:
    """
    fastqc --threads ${task.cpus} ${fastq_files}
    """
}