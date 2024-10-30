# __author__ = "Leandro C. Hermida"
# __email__ = "hermidalc@pitt.edu"
# __license__ = "MIT"

suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(rtracklayer)
})

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


gtf_file <- snakemake@input[[1]]
cat("Loading", gtf_file, "\n")
gtf <- import(gtf_file)

gene_annots <- unique(data.frame(
    ID_REF = gtf$gene_id, Symbol = gtf$gene_name, stringsAsFactors = FALSE
))

gene_lengths <- sum(width(reduce(
    exonsBy(makeTxDbFromGFF(gtf_file, format = "gtf"), by = "gene")
)))
gene_length_df <- data.frame(
    names(gene_lengths), gene_lengths,
    stringsAsFactors = FALSE
)
colnames(gene_length_df) <- c("ID_REF", snakemake@params[["length_col"]])

gene_annots <- merge(gene_annots, gene_length_df, by = "ID_REF")

gene_annots <- gene_annots[order(gene_annots$ID_REF), , drop = FALSE]

cat("Writing", snakemake@output[[1]], "\n")
write.table(
    gene_annots,
    file = snakemake@output[[1]], quote = FALSE, row.names = FALSE,
    sep = "\t"
)

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()
