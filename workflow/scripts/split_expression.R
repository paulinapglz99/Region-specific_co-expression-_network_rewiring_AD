# Subset a residualized expression matrix to one region and one diagnosis.
#
# Input matrix: genes in rows ('feature' column), samples (specimenID) in columns.
# Output: same format, only the samples listed in the selected-samples table
# (one sample per donor, see select_samples.R).

library(dplyr)
library(vroom)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

samples <- vroom(snakemake@input[["samples"]]) |>
  pull(specimenID)

# Read only the gene column and the selected samples
counts <- vroom(
  snakemake@input[["counts"]],
  col_select = c(feature, all_of(samples))
)

vroom_write(counts, snakemake@output[[1]])

message("region=", snakemake@wildcards[["region"]], " diagnosis=", snakemake@wildcards[["dx"]])
message("samples=", length(samples), " genes=", nrow(counts))
