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
    def r_samples = 'c("' + sample_ids.join('", "') + '")'
    """
    #!/usr/bin/env Rscript
    library(DESeq2)

    mat <- read.table("${homer_matrix}", sep="\t", header=TRUE,
                      quote="", comment.char="", check.names=FALSE)
    rownames(mat) <- mat[,1]
    counts <- as.matrix(mat[, 20:ncol(mat)])
    colnames(counts) <- sub("^.*/", "", colnames(counts))
    colnames(counts) <- sub(" .*",  "", colnames(counts))

    sample_ids <- ${r_samples}
    counts     <- counts[, sample_ids, drop = FALSE]

    condition <- sub("_.*", "", sample_ids)
    coldata   <- data.frame(condition = factor(condition), row.names = sample_ids)

    dds <- DESeqDataSetFromMatrix(round(counts), colData = coldata, design = ~condition)
    dds <- DESeq(dds)
    
    res_df <- as.data.frame(results(dds))
    res_df\$peak         <- rownames(res_df)
    cond_levels          <- levels(coldata\$condition)   # Factors are set alphabetically!
    res_df\$contrast_ref <- cond_levels[2]               # [2] alphabetically is "WT"! --> Reference
    res_df\$contrast_trt <- cond_levels[1]               # [1] alphabetically is "PH"!

    # Extract the gene mapping from the homer_matrix and add to the output table:
    gene_mapping <- mat[, 11]
    names(gene_mapping) <- mat[, 1]
    res_df\$gene_id <- gene_mapping[res_df\$peak]

    write.table(res_df, "deseq2_results.tsv", sep="\t", quote=FALSE, row.names=FALSE)
    """
}