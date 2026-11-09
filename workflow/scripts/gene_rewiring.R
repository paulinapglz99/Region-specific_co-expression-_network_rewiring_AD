# Analysis 4, step 1: rewiring of each gene, AD_repk vs CT_all, in one region.
#   neighbor_change = 1 - Jaccard of its neighbors in AD and in CT
#   degree_change   = log2((degree AD + 1) / (degree CT + 1))
# Universe: every gene of the DE table of the region (degree 0 if absent from a network).
# Outputs: one row per gene and replicate, and one row per gene (mean over replicates).

source(snakemake@input[["network_utils"]])
source(snakemake@input[["de_utils"]])
library(tidyr)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

region <- snakemake@wildcards[["region"]]
de <- read_de(snakemake@input[["de"]], region, snakemake@params[["adj_p"]], snakemake@params[["fold_change"]])

files <- unlist(snakemake@input[["networks"]])
networks <- set_names(map(files, read_network), basename(dirname(files)))

# Each edge seen from both ends: one row per gene and neighbor
neighbor_table <- function(edges) {
  bind_rows(transmute(edges, gene = a, neighbor = b), transmute(edges, gene = b, neighbor = a))
}

ct <- neighbor_table(networks[["CT_all"]])
degree_ct <- count(ct, gene, name = "degree_ct")

rewiring_one <- function(edges, network) {
  ad <- neighbor_table(edges)
  shared <- inner_join(ad, ct, by = c("gene", "neighbor")) |> count(gene, name = "shared")
  tibble(gene = de$gene) |>
    left_join(count(ad, gene, name = "degree_ad"), by = "gene") |>
    left_join(degree_ct, by = "gene") |>
    left_join(shared, by = "gene") |>
    replace_na(list(degree_ad = 0L, degree_ct = 0L, shared = 0L)) |>
    mutate(
      replicate = sub("^AD_", "", network),
      union = degree_ad + degree_ct - shared,
      neighbor_change = if_else(union > 0, 1 - shared / union, NA_real_),
      degree_change = log2((degree_ad + 1) / (degree_ct + 1)),
      .before = 1
    )
}

ad_names <- grep("^AD_", names(networks), value = TRUE)
replicates <- map(ad_names, \(name) rewiring_one(networks[[name]], name)) |> list_rbind()

genes <- replicates |>
  group_by(gene) |>
  summarise(
    degree_ct = first(degree_ct),
    degree_ad = mean(degree_ad),
    neighbor_change = mean(neighbor_change, na.rm = TRUE),
    degree_change = mean(degree_change)
  ) |>
  mutate(neighbor_change = if_else(is.nan(neighbor_change), NA_real_, neighbor_change)) |>
  left_join(de, by = "gene")

vroom_write(replicates, snakemake@output[["replicates"]])
vroom_write(genes, snakemake@output[["genes"]])

message(
  region, ": ", nrow(genes), " genes | ", length(ad_names), " AD replicates | ",
  sum(genes$degree_ct > 0), " genes with edges in CT"
)
