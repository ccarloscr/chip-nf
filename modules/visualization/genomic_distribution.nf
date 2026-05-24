// Genomic distribution of peaks: proportion in promoters, intergenic, introns, exons, etc.
// Parses the "Annotation" column from HOMER's annotatePeaks.pl output.
// Input : HOMER annotated peaks TSV (one per sample).
// Output: stacked bar / pie chart PNG per sample.
process GENOMIC_DISTRIBUTION {
    tag "consensus"
    label 'process_low'

    publishDir "${params.output_dir}/Visualization/GenomicDistribution", mode: 'copy'

    conda "conda-forge::r-base=4.3 conda-forge::r-ggplot2 conda-forge::r-dplyr"

    input:
    val conditions
    path tsv_files

    output:
    path "conditions_genomic_distribution.png", emit: plot  
    path "conditions_genomic_distribution.tsv", emit: table

    script:
    // Convert nextflow listed inputs into R vectors
    def r_conditions = 'c("' + conditions.join('", "') + '")'
    def r_files      = 'c("' + tsv_files.join('", "') + '")'
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(ggplot2)
        library(dplyr)
    })

    # Vectors transferred from nextflow
    conds <- ${r_conditions}
    files <- ${r_files}

    # Initialize empty df to save data from all replicates
    df_all_samples <- data.frame()

    # Loop to read each file and identify the condition
    for (i in 1:length(files)) {
        df <- read.table(files[i], sep = "\\t", header = TRUE, quote = "",
                         comment.char = "", stringsAsFactors = FALSE)
        colnames(df)[1] <- "PeakID"

        # Identify the annotation column from homer
        annot_col <- grep("^Annotation\$", colnames(df), value = TRUE)
        if (length(annot_col) == 0) annot_col <- grep("Annotation", colnames(df), value = TRUE)[1]
        
        # Simplify categories
        df\$category <- df[[annot_col]]
        df\$category <- ifelse(grepl("promoter",    df\$category, ignore.case = TRUE), "Promoter",    df\$category)
        df\$category <- ifelse(grepl("intron",      df\$category, ignore.case = TRUE), "Intron",      df\$category)
        df\$category <- ifelse(grepl("exon",        df\$category, ignore.case = TRUE), "Exon",        df\$category)
        df\$category <- ifelse(grepl("intergenic",  df\$category, ignore.case = TRUE), "Intergenic",  df\$category)
        df\$category <- ifelse(grepl("TTS",         df\$category, ignore.case = TRUE), "TTS",         df\$category)
        df\$category <- ifelse(grepl("non-coding",  df\$category, ignore.case = TRUE), "Non-coding",  df\$category)
        df\$category[!df\$category %in% c("Promoter","Intron","Exon","Intergenic","TTS","Non-coding")] <- "Other"

        # Calculate percentages for the replicate
        tbl_sample <- df %>%
            count(category) %>%
            mutate(pct = 100 * n / sum(n),
                   condition = conds[i],      # Save condition
                   file_name = files[i])      # Save replicate name
        
        df_all_samples <- rbind(df_all_samples, tbl_sample)
    }

    # Group per condition and category --> mean percentages of all replicates
    tbl_final <- df_all_samples %>%
        group_by(condition, category) %>%
        summarise(mean_pct = mean(pct), .groups = 'drop') %>%
        mutate(category = factor(category, levels = rev(sort(unique(category)))))

    # Save table with mean percentages
    write.table(tbl_final, "conditions_genomic_distribution.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)


    # Desired order of conditions
    tbl_final\$condition <- factor(tbl_final\$condition, levels = c("WT", "PH"))

    # Define the order of categories for the stacked plot
    desired_order <- c("Promoter", "Exon", "Intron", "TTS", "Non-coding", "Intergenic", "Other")
    tbl_final\$category <- factor(tbl_final\$category, levels = desired_order)


    # Define plot colours
    colours <- c(
        Promoter    = "#E41A1C",
        Intron      = "#377EB8",
        Exon        = "#4DAF4A",
        Intergenic  = "#FF7F00",
        TTS         = "#984EA3",
        `Non-coding`= "#A65628",
        Other       = "#CCCCCC"
    )

    # Stacked bar plot
    p <- ggplot(tbl_final, aes(x = condition, y = mean_pct, fill = category)) +
        geom_bar(stat = "identity", width = 0.6, colour = "white") +
        scale_fill_manual(values = colours, name = "Region") +

        # Add labels if percentage is > 5%
        geom_text(aes(label = ifelse(mean_pct > 5, paste0(round(mean_pct, 1), "%"), "")),
                  position = position_stack(vjust = 0.5), size = 3) +
        labs(title = "Genomic Distribution of Peaks by Condition",
             subtitle = "Average percentage across biological replicates",
             x = "Condition",
             y = "Mean Percentage (%)") +
        theme_minimal(base_size = 11) +
        theme(legend.position = "right",
              panel.grid.major.x = element_blank())

    ggsave("conditions_genomic_distribution.png", plot = p, width = 7, height = 5, dpi = 180)
    """
}
