// MACS2 broas peak calling process for ChIP-seq data
process CALL_PEAKS_BROAD {
    // Tag the process with the sample ID for easier tracking in logs and reports
    tag "${meta.id}"
    label 'process_low'

    // Publish MACS2 peak files to a specific directory for peak calling results
    publishDir "${params.output_dir}/Peaks", mode: 'copy'

    // Use a specific Conda environment that includes MACS2 for peak calling
    conda "bioconda::macs2=2.2.9.1"

    // Input: a tuple of sample ID, ChIP BAM and its index, control BAM and its index, along with the extension size for peak calling
    input:
    // ChIP BAM + its control (Input) BAM, joined by sample logic in main.nf
    tuple val(meta), path(chip_bam), path(chip_bai), path(ctrl_bam), path(ctrl_bai)

    output:
    tuple val(meta), path("${meta.id}_peaks.broadPeak"),    emit: peaks
    path "${meta.id}_peaks.xls",                            emit: xls

    script:
    """
    macs2 callpeak \
        -t ${chip_bam} \
        -c ${ctrl_bam} \
        -f BAMPE \
        -g ${params.genome_size} \
        -n ${meta.id} \
        --broad \
        --broad-cutoff ${params.broad_cutoff} \
        -q ${params.min_qval} \
        --outdir .
    """
}