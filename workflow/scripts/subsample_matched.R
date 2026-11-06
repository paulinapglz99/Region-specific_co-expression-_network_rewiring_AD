# Design 1 (main): AD vs CT with the same n.
#
# CT is used complete (one network, "CT_all").
# AD is subsampled to the n of CT, with the female proportion of CT ("AD_rep01", ...).

source(snakemake@input[["utils"]])

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

ad <- read_samples(snakemake@input[["ad"]])
ct <- read_samples(snakemake@input[["ct"]])
reps <- rep_label(seq_len(snakemake@params[["replicates"]]))

set_seed_for(snakemake@params[["seed"]], paste("matched", snakemake@wildcards[["region"]]))
female_prop <- mean(ct$sex == "female")

# One AD subsample per replicate
ad_networks <- reps |>
  map(\(rep) sample_by_sex(ad, nrow(ct), female_prop)) |>
  set_names(paste0("AD_", reps))

subsamples <- c(list(CT_all = ct), ad_networks) |>
  bind_rows(.id = "network") |>
  mutate(diagnosis = sub("_.*", "", network)) |>
  select(network, diagnosis, specimenID, individualID, sex, RIN)

vroom_write(subsamples, snakemake@output[[1]])

message("CT n=", nrow(ct), " (", round(100 * female_prop, 1), "% female) | AD pool=", nrow(ad))
print(summarise_networks(subsamples), n = Inf)
