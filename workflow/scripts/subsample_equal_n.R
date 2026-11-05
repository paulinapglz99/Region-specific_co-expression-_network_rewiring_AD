# Design 2 (robustness): AD and CT both subsampled to the same n,
# with the female proportion of CT ("AD_rep01", "CT_rep01", ...).

source(snakemake@input[["utils"]])

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

ad <- fread(snakemake@input[["ad"]])
ct <- fread(snakemake@input[["ct"]])
n <- snakemake@params[["n"]]
replicates <- snakemake@params[["replicates"]]

if (n > min(nrow(ad), nrow(ct))) {
  stop("n = ", n, " is larger than the AD (", nrow(ad), ") or CT (", nrow(ct), ") group")
}

set_seed_for(
  snakemake@params[["seed"]],
  paste("equal_n", snakemake@wildcards[["set"]], snakemake@wildcards[["region"]])
)
female_prop <- mean(ct$sex == "female")

networks <- list()
for (r in seq_len(replicates)) {
  networks[[paste0("AD_", rep_label(r))]] <- sample_by_sex(ad, n, female_prop)
  networks[[paste0("CT_", rep_label(r))]] <- sample_by_sex(ct, n, female_prop)
}

subsamples <- rbindlist(networks, idcol = "network")
subsamples[, diagnosis := sub("_.*", "", network)]
fwrite(
  subsamples[, .(network, diagnosis, specimenID, individualID, sex, RIN)],
  snakemake@output[[1]],
  sep = "\t"
)

message("n=", n, " | target female=", round(100 * female_prop, 1), "% (CT)",
        " | pools AD=", nrow(ad), " CT=", nrow(ct))
print(summarise_networks(subsamples))
