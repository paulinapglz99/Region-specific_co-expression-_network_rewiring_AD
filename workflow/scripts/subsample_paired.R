# Design 4 (sensitivity): paired regions, same donors in both networks.
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
samples <- list(
  a_ad = fread(snakemake@input[["a_ad"]]), a_ct = fread(snakemake@input[["a_ct"]]),
  b_ad = fread(snakemake@input[["b_ad"]]), b_ct = fread(snakemake@input[["b_ct"]])
)
replicates <- snakemake@params[["replicates"]]

# Donors present in both regions (one row per donor, from region A)
shared_ct <- samples$a_ct[individualID %in% samples$b_ct$individualID]
shared_ad <- samples$a_ad[individualID %in% samples$b_ad$individualID]

set_seed_for(snakemake@params[["seed"]], paste("paired", region_a, region_b))
female_prop <- mean(shared_ct$sex == "female")

donors <- list(CT_all = shared_ct$individualID)
for (r in seq_len(replicates)) {
  picked <- sample_by_sex(shared_ad, nrow(shared_ct), female_prop)
  donors[[paste0("AD_", rep_label(r))]] <- picked$individualID
}

# For each network and region, the specimen of each donor in that region
region_samples <- function(region, ad, ct) {
  rbindlist(lapply(names(donors), function(net) {
    group <- if (startsWith(net, "AD")) ad else ct
    group[individualID %in% donors[[net]]][, `:=`(network = net, region = region)]
  }))
}
subsamples <- rbind(
  region_samples(region_a, samples$a_ad, samples$a_ct),
  region_samples(region_b, samples$b_ad, samples$b_ct)
)
subsamples[, diagnosis := sub("_.*", "", network)]
fwrite(
  subsamples[, .(network, region, diagnosis, specimenID, individualID, sex, RIN)],
  snakemake@output[[1]],
  sep = "\t"
)

message(region_a, " + ", region_b, " | shared CT=", nrow(shared_ct),
        " (", round(100 * female_prop, 1), "% female) | shared AD=", nrow(shared_ad))
print(subsamples[, .(n = .N, pct_female = round(100 * mean(sex == "female"), 1)),
                 by = .(network, region)])
