// Fastp process definition for adapter trimming and quality filtering of sequencing reads
process FASTP {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.output_dir}/QC/Fastp", mode: 'copy', pattern: "*.{json,html}"
    publishDir "${params.output_dir}/trimmed",  mode: 'copy', pattern: "*.fastq.gz"

    conda "bioconda::fastp=0.23.4"

    input:
    tuple val(meta), path(fastq_files)

    // trimmed reads (paired-end), JSON report (for MultiQC), HTML report (human inspection)
    output:
    tuple val(meta), path("${meta.id}_R{1,2}_trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}_fastp.json"),               emit: json
    tuple val(meta), path("${meta.id}_fastp.html"),               emit: html

    script:
    """
    fastp \\
        --in1  ${fastq_files[0]} \\
        --in2  ${fastq_files[1]} \\
        --out1 ${meta.id}_R1_trimmed.fastq.gz \\
        --out2 ${meta.id}_R2_trimmed.fastq.gz \\
        --json ${meta.id}_fastp.json \\
        --html ${meta.id}_fastp.html \\
        --detect_adapter_for_pe \\
        --cut_front \\
        --cut_tail \\
        --cut_window_size 4 \\
        --cut_mean_quality  ${params.fastp_cut_mean_quality} \\
        --length_required   ${params.fastp_min_read_length} \\
        --thread ${task.cpus}
    """
}
