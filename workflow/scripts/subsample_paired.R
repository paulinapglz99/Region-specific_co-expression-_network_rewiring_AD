# Paired design (not in the current design): paired regions, same donors in both networks.
#
# Only donors sampled in both regions are used.
# CT: all shared CT donors ("CT_all").
# AD: shared AD donors subsampled to the n of CT, with the female proportion of CT
#     ("AD_rep01", ...). The same AD donors are used in both regions.

source(snakemake@input[["utils"]])

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

region_a <- snakemake@wildcards[["region_a"]]
region_b <- snakemake@wildcards[["region_b"]]
a_ad <- read_samples(snakemake@input[["a_ad"]])
a_ct <- read_samples(snakemake@input[["a_ct"]])
b_ad <- read_samples(snakemake@input[["b_ad"]])
b_ct <- read_samples(snakemake@input[["b_ct"]])
reps <- rep_label(seq_len(snakemake@params[["replicates"]]))

# Donors present in both regions (one row per donor, taken from region A)
shared_ct <- a_ct |> filter(individualID %in% b_ct$individualID)
shared_ad <- a_ad |> filter(individualID %in% b_ad$individualID)

set_seed_for(snakemake@params[["seed"]], paste("paired", region_a, region_b))
female_prop <- mean(shared_ct$sex == "female")

# Donor IDs of every network: all shared CT + one AD draw per replicate
ad_donors <- reps |>
  map(\(rep) sample_by_sex(shared_ad, nrow(shared_ct), female_prop) |> pull(individualID)) |>
  set_names(paste0("AD_", reps))
donors <- c(list(CT_all = shared_ct$individualID), ad_donors)

# For each network, the specimen of each of its donors in one region
region_samples <- function(region, ad, ct) {
  donors |>
    imap(\(ids, net) {
      group <- if (startsWith(net, "AD")) ad else ct
      group |>
        filter(individualID %in% ids) |>
        mutate(network = net, region = region)
    }) |>
    bind_rows()
}

subsamples <- bind_rows(
  region_samples(region_a, a_ad, a_ct),
  region_samples(region_b, b_ad, b_ct)
) |>
  mutate(diagnosis = sub("_.*", "", network)) |>
  select(network, region, diagnosis, specimenID, individualID, sex, RIN)

vroom_write(subsamples, snakemake@output[[1]])

message(
  region_a, " + ", region_b, " | shared CT=", nrow(shared_ct),
  " (", round(100 * female_prop, 1), "% female) | shared AD=", nrow(shared_ad)
)
print(summarise_networks(subsamples, region), n = Inf)
