#!/usr/bin/env Rscript
library(DESeq2)

args              <- commandArgs(trailingOnly = TRUE)
homer_matrix      <- args[1]
sample_ids        <- strsplit(args[2], ",")[[1]]
reference_cond    <- args[3]
treatment_cond    <- args[4]

mat <- read.table(homer_matrix, sep="\t", header=TRUE,
                    quote="", comment.char="", check.names=FALSE)

rownames(mat) <- mat[,1]

counts <- as.matrix(mat[, 20:ncol(mat)])
colnames(counts) <- sub("^.*/", "", colnames(counts))
colnames(counts) <- sub(" .*",  "", colnames(counts))
counts    <- counts[, sample_ids, drop = FALSE]

condition <- sub("_.*", "", sample_ids)
coldata   <- data.frame(condition = factor(condition), row.names = sample_ids)

dds <- DESeqDataSetFromMatrix(round(counts), colData = coldata, design = ~condition)
dds <- DESeq(dds)

res_df <- as.data.frame(results(dds))
res_df$peak         <- rownames(res_df)
res_df$contrast_ref <- reference_cond
res_df$contrast_trt <- treatment_cond

# Extract the gene mapping from the homer_matrix and add to the output table:
gene_mapping <- mat[, 11]
names(gene_mapping) <- mat[, 1]
res_df$gene_id <- gene_mapping[res_df$peak]

write.table(res_df, "deseq2_results.tsv", sep="\t", quote=FALSE, row.names=FALSE)