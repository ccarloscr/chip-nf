// HOMER peak annotation process for ChIP-seq data
process ANNOTATE_PEAKS {
    // Tag the process with the sample ID for easier tracking in logs and reports
    tag "${meta.id}"
    label 'process_high'

    // Publish annotated peak files to a specific directory for annotation results
    publishDir "${params.output_dir}/Annotated", mode: 'copy'

    // Use a specific Conda environment that includes HOMER for peak annotation
    conda "bioconda::homer=4.11 conda-forge::perl"

    //
    input:
    tuple val(meta), path(peaks)
    path gtf_file       // Pass the GTF annotation file
    path homer_data_dir     // pre-built genome data dir from HOMER_SETUP

    // Output annotated peaks and annotation stats
    output:
    tuple val(meta), path("${meta.id}_annotated.txt"), emit: annotated
    path("${meta.id}_annotation_stats.txt")          , emit: stats

    // Run the HOMER annotatePeaks.pl script to annotate the peaks with genomic features
    script:
    """
    export PATH=\$CONDA_PREFIX/share/homer/bin:\$PATH

    # Point HOMER at the pre-built data directory instead of the conda default,
    # so no download is attempted on the compute node.
    export HOMER_DATA_DIR=\$(realpath ${homer_data_dir})

    annotatePeaks.pl ${peaks} ${params.genome} \\
        -gtf ${gtf_file} \\
        -annStats ${meta.id}_annotation_stats.txt \\
        -raw \\
        > ${meta.id}_annotated.txt
    """
}