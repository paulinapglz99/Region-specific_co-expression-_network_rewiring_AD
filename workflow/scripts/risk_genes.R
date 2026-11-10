# Analysis 7: are the AD risk genes among the genes that change in the network?
#   1. Overlap tests (de_utils.R), per region, BH across regions and gene sets:
#        risk x rewired genes (analysis 4)  universe: genes of analysis 4,
#                                           strata: AveExpr x mean degree
#        risk x gained hubs (analysis 5)    universe: genes of analysis 5,
#        risk x lost hubs (analysis 5)      strata: AveExpr
#   2. Table: each risk gene in each region, with its rewiring percentile,
#      hub type and DE status.
# Supplementary analysis: with a few tens of risk genes the tests have little power.

source(snakemake@input[["de_utils"]])
library(tidyr)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

n_perm <- snakemake@params[["n_perm"]]
bins <- snakemake@params[["strata"]]
set.seed(snakemake@params[["seed"]])

risk <- vroom(snakemake@input[["risk"]], col_types = cols(.default = "c"))
rewiring <- vroom(snakemake@input[["rewiring"]], col_types = cols(region = "c", gene = "c", symbol = "c", de_class = "c", group = "c")) |>
  mutate(is_risk = gene %in% risk$ensembl_gene_id)
hubs <- vroom(snakemake@input[["hubs"]], col_types = cols(region = "c", gene = "c", symbol = "c", de_class = "c", hub_type = "c")) |>
  mutate(is_risk = gene %in% risk$ensembl_gene_id)


# ---- 1. Overlap tests -------------------------------------------------------

test_set <- function(df, in_set, strata_vars, gene_set) {
  df |>
    mutate(in_set = {{ in_set }}) |>
    group_by(region) |>
    group_modify(\(d, key) {
      d <- add_strata(d, strata_vars, bins)
      overlap_test(d$in_set, d$is_risk, d$stratum, n_perm)
    }) |>
    ungroup() |>
    mutate(gene_set = gene_set, .after = region)
}

tests <- bind_rows(
  test_set(rewiring, rewired, c("AveExpr", "mean_degree"), "rewired"),
  test_set(hubs, hub_type == "gained", "AveExpr", "gained hubs"),
  test_set(hubs, hub_type == "lost", "AveExpr", "lost hubs")
) |>
  mutate(p_adj = p.adjust(p_perm, method = "BH"))


# ---- 2. Table of the risk genes ---------------------------------------------

table <- risk |>
  rename(gene = ensembl_gene_id) |>
  cross_join(tibble(region = unique(hubs$region))) |>
  left_join(
    rewiring |>
      mutate(rewiring_percentile = 100 * percent_rank(rewiring_adj), .by = region) |>
      select(region, gene, rewiring_percentile, rewired),
    by = c("region", "gene")
  ) |>
  left_join(select(hubs, region, gene, hub_type, is_de, de_class, logFC, adj_p), by = c("region", "gene")) |>
  mutate(
    in_networks = !is.na(hub_type), # >= 1 edge in a network of the region
    passes_degree_filter = !is.na(rewired) # in analysis 4
  ) |>
  relocate(region, .after = tier)

vroom_write(tests, snakemake@output[["tests"]])
vroom_write(table, snakemake@output[["table"]])

message(nrow(risk), " risk genes")
print(summarise(table, in_networks = sum(in_networks), in_analysis_4 = sum(passes_degree_filter), rewired = sum(rewired, na.rm = TRUE), .by = region))
print(mutate(tests, across(where(is.double), \(v) signif(v, 3))), n = Inf, width = Inf)
