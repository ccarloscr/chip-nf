// Volcano plot for differential ChIP-seq peak analysis.
// Input : deseq2_results.tsv from DESEQ2_PEAKS
//         Expected columns: peak, log2FoldChange, padj, contrast_ref, contrast_trt
//         Optional column : gene (nearest gene name; falls back to peak ID if absent)
// Output: volcano_plot.png
process VOLCANO_PLOT {
    tag "consensus"
    label 'process_low'

    publishDir "${params.output_dir}/Visualization/Volcano", mode: 'copy'

    conda "conda-forge::r-base=4.3 conda-forge::r-ggplot2 conda-forge::r-ggrepel conda-forge::r-dplyr"

    input:
    path deseq2_results     // deseq2_results.tsv from DESEQ2_PEAKS.out.results
 
    output:
    path "volcano_plot.png", emit: plot

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(ggplot2)
        library(ggrepel)
        library(dplyr)
    })

    # -------------------------------------------------------------------------
    # 1. Load DESeq2 results
    # -------------------------------------------------------------------------
    res_df <- read.table("${deseq2_results}", sep = "\\t", header = TRUE,
                         quote = "", stringsAsFactors = FALSE)

    # Keep only peaks with LFC and padj estimates
    plot_df <- res_df[!is.na(res_df\$log2FoldChange) & !is.na(res_df\$padj), ]

    # -------------------------------------------------------------------------
    # 2. Label: use gene column if present and non-empty, otherwise peak ID
    # -------------------------------------------------------------------------
    if ("gene_id" %in% colnames(plot_df) ) {
        plot_df\$display_label <- ifelse(
            !is.na(plot_df\$gene_id) & nchar(trimws(plot_df\$gene_id)) > 0,
            plot_df\$gene_id,
            plot_df\$peak
        )
    } else {
        plot_df\$display_label <- plot_df\$peak
    }

    # -------------------------------------------------------------------------
    # 3. Classify significance
    # -------------------------------------------------------------------------
    plot_df\$sig_class <- case_when(
        plot_df\$padj < 0.05 & plot_df\$log2FoldChange >=  1 ~ "Up",
        plot_df\$padj < 0.05 & plot_df\$log2FoldChange <= -1 ~ "Down",
        TRUE                                                  ~ "NS"
    )
    plot_df\$sig_class <- factor(plot_df\$sig_class, levels = c("Up", "Down", "NS"))

    # -------------------------------------------------------------------------
    # 4. Condition names from contrast columns written by DESEQ2_PEAKS
    # -------------------------------------------------------------------------
    cond_a <- unique(plot_df\$contrast_ref)[1]   # reference (denominator)
    cond_b <- unique(plot_df\$contrast_trt)[1]   # treatment (numerator)

    # -------------------------------------------------------------------------
    # 5. Top N peaks to label: lowest padj among significant, ties by |log2FC|
    # -------------------------------------------------------------------------
    top_peaks <- plot_df %>%
        filter(sig_class != "NS") %>%
        arrange(padj, desc(abs(log2FoldChange))) %>%
        slice_head(n = ${params.volcano_top})

    # -------------------------------------------------------------------------
    # 6. Clip extreme LFC values for display (real values stay in the TSV)
    # -------------------------------------------------------------------------
    lfc_limit          <- max(ceiling(quantile(abs(plot_df\$log2FoldChange), 0.99)), 4)
    plot_df\$lfc_plot   <- pmax(pmin(plot_df\$log2FoldChange,  lfc_limit), -lfc_limit)
    top_peaks\$lfc_plot <- pmax(pmin(top_peaks\$log2FoldChange, lfc_limit), -lfc_limit)

    # -------------------------------------------------------------------------
    # 7. Plot
    # -------------------------------------------------------------------------
    n_up   <- sum(plot_df\$sig_class == "Up")
    n_down <- sum(plot_df\$sig_class == "Down")
    subtitle_text <- paste0(
        cond_b, " vs ", cond_a,
        "  |  Up: ", n_up, "  Down: ", n_down,
        "  (padj < 0.05, |log2FC| \u2265 1)"
    )

    p <- ggplot(plot_df, aes(x = lfc_plot, y = -log10(padj + 1e-300), colour = sig_class)) +
        geom_point(alpha = 0.55, size = 1.2) +
        scale_colour_manual(
            values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey70"),
            labels = c(Up   = paste0("Up in ", cond_b),
                       Down = paste0("Down in ", cond_b),
                       NS   = "NS"),
            name = NULL
        ) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black", linewidth = 0.35) +
        geom_vline(xintercept = c(-1, 1),     linetype = "dashed", colour = "black", linewidth = 0.35) +
        geom_text_repel(
            data           = top_peaks,
            aes(x = lfc_plot, y = -log10(padj + 1e-300), label = display_label),
            colour         = "black",
            size           = 2.5,
            max.overlaps   = 40,
            segment.colour = "grey40",
            segment.size   = 0.3,
            box.padding    = 0.3
        ) +
        labs(
            title    = "Differential Peak Analysis",
            subtitle = subtitle_text,
            x        = paste0("log2 Fold Change  (", cond_b, " / ", cond_a, ")"),
            y        = expression(-log[10](p[adj]))
        ) +
        theme_bw(base_size = 11) +
        theme(legend.position  = "bottom",
              panel.grid.minor = element_blank())

    ggsave("volcano_plot.png", plot = p, width = 7, height = 6, dpi = 180)
    """
}
