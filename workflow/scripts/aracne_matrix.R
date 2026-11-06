# Expression matrix and regulator list for one ARACNe-AP network.
#
# Input: the expression matrix of the network's region x diagnosis group
# (results/expression/, all donors of the group) and the subsample table.
# Matrix: genes in rows (first column 'gene'), only the samples of this network
# in columns, tab-separated (ARACNe-AP input format).
# Regulators: every gene of the matrix (the full gene universe).

library(dplyr)
library(vroom)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

net <- snakemake@wildcards[["network"]]
subsamples <- vroom(snakemake@input[["subsamples"]])

if ("region" %in% names(subsamples)) {
  stop("Paired tables (with a 'region' column) are not supported yet")
}

samples <- subsamples |>
  filter(network == net) |>
  pull(specimenID)

if (length(samples) == 0) {
  stop("Network ", net, " not found in ", snakemake@input[["subsamples"]])
}

# Columns of the group matrix = the samples of this network
expression <- vroom(
  snakemake@input[["expression"]],
  col_select = c(feature, all_of(samples))
) |>
  rename(gene = feature)

vroom_write(expression, snakemake@output[["matrix"]])
vroom_write_lines(expression$gene, snakemake@output[["regulators"]])

message("network=", net, " | group matrix=", snakemake@input[["expression"]])
message("samples=", length(samples), " genes=", nrow(expression))
