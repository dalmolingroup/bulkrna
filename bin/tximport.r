#!/usr/bin/env Rscript

library(tximport)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: tximport.r <kallisto_results> <tx2gene>", call. = FALSE)
}

kallisto_results <- args[1]

tx2gene <- readRDS(args[2])

files <- list.files(path = kallisto_results, pattern = "tsv", recursive = TRUE, full.names = TRUE)

names(files) <- sapply(strsplit(files, "\\/"), "[[", 2)

txi_gene <- tximport(
    files = files,
    type = "kallisto",
    tx2gene = tx2gene,
    ignoreTxVersion = TRUE,
    ignoreAfterBar = TRUE,
    dropInfReps=TRUE
)

txi_tx <- tximport(
    files = files,
    type = "kallisto",
    txOut = TRUE,
    ignoreTxVersion = TRUE,
    ignoreAfterBar = TRUE,
    dropInfReps=TRUE
)

txi_tx_scaled <- tximport(
    files = files,
    type = "kallisto",
    txOut = TRUE,
    countsFromAbundance = "scaledTPM",
    ignoreTxVersion = TRUE,
    ignoreAfterBar = TRUE,
    dropInfReps=TRUE
)

write.table(txi_gene$abundance,  file="gene_tpm.tsv", sep="\t", quote=FALSE)
write.table(txi_gene$counts,  file="gene_counts.tsv", sep="\t", quote=FALSE)
write.table(txi_tx$abundance,  file="transcript_tpm.tsv", sep="\t", quote=FALSE)
write.table(txi_tx$counts,  file="transcript_counts.tsv", sep="\t", quote=FALSE)
write.table(txi_tx_scaled$abundance,  file="transcript_scaled_tpm.tsv", sep="\t", quote=FALSE)
write.table(txi_tx_scaled$counts,  file="transcript_scaled_counts.tsv", sep="\t", quote=FALSE)
