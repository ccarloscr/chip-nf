process MERGE_PEAKS {
    label 'process_low'
    conda "bioconda::bedtools=2.31.1"

    input:
    path all_peaks_ch

    output:
    path "consensus_peaks.bed"

    script:
    """
    # Sort and merge all overlapping peaks across all samples into a master list
    cat ${all_peaks_ch} \
        | awk 'BEGIN{OFS="\t"} !/^#/ {print \$1,\$2,\$3}' \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i stdin \
        > consensus_peaks.bed
    """
}