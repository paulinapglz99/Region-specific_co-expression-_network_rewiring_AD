# Subset a residualized expression matrix to one region and one diagnosis.
#
# Input matrix: genes in rows ('feature' column), samples (specimenID) in columns.
# Output: same format, only the samples listed in the selected-samples table
# (one sample per donor, see select_samples.R).

library(data.table)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

samples <- fread(snakemake@input[["samples"]])$specimenID

# Read only the gene column and the selected samples
counts <- fread(snakemake@input[["counts"]], select = c("feature", samples))

fwrite(counts, snakemake@output[[1]], sep = "\t")

message("region=", snakemake@wildcards[["region"]], " diagnosis=", snakemake@wildcards[["dx"]])
message("samples=", length(samples), " genes=", nrow(counts))
