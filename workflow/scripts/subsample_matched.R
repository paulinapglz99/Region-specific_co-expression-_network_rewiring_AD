# Design 1 (main): AD vs CT with the same n.
#
# CT is used complete (one network, "CT_all").
# AD is subsampled to the n of CT, with the female proportion of CT ("AD_rep01", ...).

source(snakemake@input[["utils"]])

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

ad <- fread(snakemake@input[["ad"]])
ct <- fread(snakemake@input[["ct"]])
replicates <- snakemake@params[["replicates"]]

set_seed_for(snakemake@params[["seed"]], paste("matched", snakemake@wildcards[["region"]]))
female_prop <- mean(ct$sex == "female")

networks <- list(CT_all = ct)
for (r in seq_len(replicates)) {
  networks[[paste0("AD_", rep_label(r))]] <- sample_by_sex(ad, nrow(ct), female_prop)
}

subsamples <- rbindlist(networks, idcol = "network")
subsamples[, diagnosis := sub("_.*", "", network)]
fwrite(
  subsamples[, .(network, diagnosis, specimenID, individualID, sex, RIN)],
  snakemake@output[[1]],
  sep = "\t"
)

message("CT n=", nrow(ct), " (", round(100 * female_prop, 1), "% female) | AD pool=", nrow(ad))
print(summarise_networks(subsamples))
