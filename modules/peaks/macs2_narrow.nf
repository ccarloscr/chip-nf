// MACS2 narrow peak calling process for ChIP-seq data
process CALL_PEAKS_NARROW {
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
    tuple val(meta), path("${meta.id}_peaks.narrowPeak"),   emit: peaks
    tuple val(meta), path("${meta.id}_summits.bed"),        emit: summits
    path "${meta.id}_peaks.xls",                            emit: xls

    script:
    """
    macs2 callpeak \
        -t ${chip_bam} \
        -c ${ctrl_bam} \
        -f BAMPE \
        -g ${params.genome_size} \
        -n ${meta.id} \
        -q ${params.min_qval} \
        --keep-dup auto \
        --call-summits \
        --outdir .
    """
}