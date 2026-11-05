# Design 3 (reproducibility floor): two disjoint subsamples of size n from the
# same group ("rep01_A", "rep01_B", ...). No donor is in both A and B.
# The female proportion is the one of the group itself.

source(snakemake@input[["utils"]])

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

pool <- fread(snakemake@input[["samples"]])
n <- as.integer(snakemake@wildcards[["n"]])
replicates <- snakemake@params[["replicates"]]

if (2 * n > nrow(pool)) {
  stop("Two disjoint subsamples of n = ", n, " need ", 2 * n, " donors; the group has ", nrow(pool))
}

set_seed_for(
  snakemake@params[["seed"]],
  paste("curve", snakemake@wildcards[["region"]], snakemake@wildcards[["dx"]], n)
)
female_prop <- mean(pool$sex == "female")

networks <- list()
for (r in seq_len(replicates)) {
  halves <- sample_two_disjoint(pool, n, female_prop)
  networks[[paste0(rep_label(r), "_A")]] <- halves$A
  networks[[paste0(rep_label(r), "_B")]] <- halves$B
}

subsamples <- rbindlist(networks, idcol = "network")
subsamples[, diagnosis := snakemake@wildcards[["dx"]]]
fwrite(
  subsamples[, .(network, diagnosis, specimenID, individualID, sex, RIN)],
  snakemake@output[[1]],
  sep = "\t"
)

message("n=", n, " x 2 | group=", nrow(pool), " (", round(100 * female_prop, 1), "% female)")
print(summarise_networks(subsamples))
