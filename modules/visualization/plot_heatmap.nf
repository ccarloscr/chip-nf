// deeptools heatmap process for ChIP-seq data
// Computes a signal matrix over consensus peaks with bamCompare/multiBamSummary
// and visualises it as a clustered heatmap using computeMatrix + plotHeatmap.
//
// Strategy:
//   1. bamCoverage  → normalised bigWig per ChIP sample (RPKM, duplicate-filtered)
//   2. computeMatrix → score matrix centred on consensus peaks
//   3. plotHeatmap  → clustered heatmap comparing signal across conditions

process HEATMAP_PLOT {
    tag "consensus_peaks"
    label 'process_high'
    publishDir "${params.output_dir}/Visualization/Heatmap", mode: 'copy'

    conda "bioconda::deeptools=3.5.4 conda-forge::python=3.10 conda-forge::numpy=1.26"

    // ------------------------------------------------------------------ //
    // Inputs
    // ------------------------------------------------------------------ //
    input:
    // Consensus peak set produced by MERGE_PEAKS
    path consensus_peaks

    // All ChIP BAM files (non-control) and their indices, collected into lists
    path bam_files, stageAs: "bams/*.bam"
    path bai_files, stageAs: "bams/*.bai"

    // Sample IDs in the same order as the BAM files (used for labels)
    val  sample_ids

    // ------------------------------------------------------------------ //
    // Outputs
    // ------------------------------------------------------------------ //
    output:
    path "heatmap_consensus_peaks.png",     emit: heatmap
    path "heatmap_consensus_peaks.pdf",     emit: heatmap_pdf
    path "score_matrix.gz",                 emit: matrix      // re-usable for plotProfile etc.

    // ------------------------------------------------------------------ //
    // Script
    // ------------------------------------------------------------------ //
    script:
    // Build a space-separated string of labels from the sample_ids list
    def labels = sample_ids instanceof List
                    ? sample_ids.join(' ')
                    : sample_ids.toString()

    """
    # ── 1. Generate a normalised bigWig for every ChIP BAM ──────────────
    mkdir -p bigwigs

    for bam in bams/*.bam; do
        sample=\$(basename \$bam .bam)
        bamCoverage \\
            --bam          \$bam \\
            --outFileName  bigwigs/\${sample}.bw \\
            --outFileFormat bigwig \\
            --normalizeUsing RPKM \\
            --ignoreDuplicates \\
            --numberOfProcessors ${task.cpus} \\
            --extendReads
    done

    # Collect bigWigs in the same order as the input BAMs
    # (preserves the label order defined by sample_ids)
    bigwig_list=\$(ls bams/*.bam \\
        | xargs -I{} basename {} .bam \\
        | xargs -I{} echo bigwigs/{}.bw)

    # ── 2. Build score matrix centred on consensus peaks ────────────────
    computeMatrix reference-point \\
        --referencePoint  center \\
        --regionsFileName ${consensus_peaks} \\
        --scoreFileName   \$bigwig_list \\
        --samplesLabel    ${labels} \\
        --beforeRegionStartLength 5000 \\
        --afterRegionStartLength  5000 \\
        --skipZeros \\
        --numberOfProcessors ${task.cpus} \\
        --outFileName score_matrix.gz

    # ── 3. Plot the heatmap ──────────────────────────────────────────────
    plotHeatmap \\
        --matrixFile        score_matrix.gz \\
        --outFileName       heatmap_consensus_peaks.png \\
        --outFileSortedRegions sorted_regions.bed \\
        --colorMap          RdBu_r \\
        --whatToShow        'heatmap and colorbar' \\
        --kmeans            1 \\
        --heatmapHeight     15 \\
        --heatmapWidth      5 \\
        --xAxisLabel        "Distance from peak centre (bp)" \\
        --regionsLabel      "Consensus peaks" \\
        --plotTitle         "ChIP-seq signal at consensus peaks"

    # PDF copy for publication-quality figures
    plotHeatmap \\
        --matrixFile        score_matrix.gz \\
        --outFileName       heatmap_consensus_peaks.pdf \\
        --colorMap          RdBu_r \\
        --whatToShow        'heatmap and colorbar' \\
        --kmeans            1 \\
        --heatmapHeight     15 \\
        --heatmapWidth      5 \\
        --xAxisLabel        "Distance from peak centre (bp)" \\
        --regionsLabel      "Consensus peaks" \\
        --plotTitle         "ChIP-seq signal at consensus peaks"
    """
}
