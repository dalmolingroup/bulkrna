#!/usr/bin/env Rscript

library(GenomicFeatures)
library(AnnotationDbi)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: gen_tax2gene.r <gtf_file> <output_name>", call. = FALSE)
}

gtf_file <- args[1]

output_name <- args[2]

txdb_filename <- "gtf.sqlite"

txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
saveDb(txdb, txdb_filename)

txdb <- loadDb(txdb_filename)

txdf <- AnnotationDbi::select(
    txdb, keys(txdb, "GENEID"),
    "TXNAME", "GENEID"
)
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

tx2gene <- data.frame(
    tx = txdf$TXNAME,
    gene = txdf$GENEID
)

saveRDS(tx2gene, file = output_name)
