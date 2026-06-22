// DESEq2 process for differential analysis of peaks between conditions prior to Volcano Plot
process DESEQ2_PEAKS {
    tag "consensus"
    conda "conda-forge::r-base=4.3 bioconda::bioconductor-deseq2 conda-forge::r-ggplot2 conda-forge::r-ggrepel"
    publishDir "${params.output_dir}/DESeq2", mode: 'copy'

    input:
    path homer_matrix
    val  sample_ids      // list like [WT_rep1, WT_rep2, PH_rep1, PH_rep2]

    output:
    path "deseq2_results.tsv", emit: results

    script:
    def r_samples = sample_ids.join(",")
    """
    deseq2_peaks.R ${homer_matrix} "${r_samples}" ${params.reference_condition} ${params.treatment_condition}
    """
}
