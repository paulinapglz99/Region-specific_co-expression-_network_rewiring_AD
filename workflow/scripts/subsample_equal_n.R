# Design 2 (robustness): AD and CT both subsampled to the same n,
# with the female proportion of CT ("AD_rep01", "CT_rep01", ...).

source(snakemake@input[["utils"]])

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

ad <- read_samples(snakemake@input[["ad"]])
ct <- read_samples(snakemake@input[["ct"]])
n <- snakemake@params[["n"]]
reps <- rep_label(seq_len(snakemake@params[["replicates"]]))

if (n > min(nrow(ad), nrow(ct))) {
  stop("n = ", n, " is larger than the AD (", nrow(ad), ") or CT (", nrow(ct), ") group")
}

set_seed_for(
  snakemake@params[["seed"]],
  paste("equal_n", snakemake@wildcards[["region"]])
)
female_prop <- mean(ct$sex == "female")

# One AD and one CT subsample per replicate
subsamples <- reps |>
  map(\(rep) list(
    AD = sample_by_sex(ad, n, female_prop),
    CT = sample_by_sex(ct, n, female_prop)
  )) |>
  set_names(reps) |>
  map(\(pair) bind_rows(pair, .id = "diagnosis")) |>
  bind_rows(.id = "replicate") |>
  mutate(network = paste0(diagnosis, "_", replicate)) |>
  select(network, diagnosis, specimenID, individualID, sex, RIN)

vroom_write(subsamples, snakemake@output[[1]])

message(
  "n=", n, " | target female=", round(100 * female_prop, 1), "% (CT)",
  " | pools AD=", nrow(ad), " CT=", nrow(ct)
)
print(summarise_networks(subsamples), n = Inf)
