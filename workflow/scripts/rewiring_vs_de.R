# Analysis 4: are the rewired genes the DE genes? (all regions)
#   1. Degree filter: genes with >= min_degree neighbors in AD (mean) or in CT
#   2. Adjusted rewiring = residual of neighbor_change ~ log(mean degree) + AveExpr
#   3. Rewired genes = top fraction by adjusted rewiring (main: 10%)
#   4. Overlap test rewired x DE (de_utils.R): odds ratio, with DE labels
#      shuffled within strata of AveExpr x mean degree. BH across regions.
#      Sensitivity: the same test with other top fractions (5% and 20%)
#   5. Four groups of genes: rewired and DE, rewired only, DE only, no change
#   6. Stability, from the adjusted rewiring of each AD replicate:
#        - odds ratio rewired x DE of each replicate (estimate only)
#        - robust rewired gene = in the top fraction in >= min_replicates of
#          the replicates; table of these genes, ordered by adjusted rewiring
# Input: the gene tables of gene_rewiring.R, one per region.

source(snakemake@input[["de_utils"]])
library(tidyr)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

min_degree <- snakemake@params[["min_degree"]]
top_rewired <- snakemake@params[["top_rewired"]]
thresholds <- unique(c(top_rewired, unlist(snakemake@params[["sensitivity"]])))
min_replicates <- snakemake@params[["min_replicates"]]
n_perm <- snakemake@params[["n_perm"]]
bins <- snakemake@params[["strata"]]
set.seed(snakemake@params[["seed"]])

# Files are named {region}.tsv and {region}_replicates.tsv
region_of <- \(files) sub("(_replicates)?\\.tsv$", "", basename(files))
read_tables <- function(files, ...) {
  set_names(files, region_of(files)) |>
    map(\(f) vroom(f, ...)) |>
    list_rbind(names_to = "region")
}
all_genes <- read_tables(unlist(snakemake@input[["genes"]]), col_types = cols(gene = "c", symbol = "c", de_class = "c"))
replicates <- read_tables(unlist(snakemake@input[["replicates"]]), col_types = cols(gene = "c", replicate = "c"))

# Adjusted rewiring: what is left of neighbor_change after degree and expression
adjust_rewiring <- function(df) {
  model <- lm(neighbor_change ~ log(mean_degree) + AveExpr, data = df)
  df |> mutate(rewiring_adj = residuals(model))
}
is_top <- \(x, fraction) x >= quantile(x, 1 - fraction)


# ---- 1-3. Filter, adjusted rewiring and rewired genes ------------------------

genes <- all_genes |>
  filter(degree_ad >= min_degree | degree_ct >= min_degree, !is.na(AveExpr)) |>
  mutate(mean_degree = (degree_ad + degree_ct) / 2) |>
  group_by(region) |>
  group_modify(\(df, key) {
    df |>
      adjust_rewiring() |>
      mutate(rewired = is_top(rewiring_adj, top_rewired)) |>
      add_strata(c("AveExpr", "mean_degree"), bins)
  }) |>
  ungroup()

genes |>
  summarise(genes = n(), de = sum(is_de), .by = region) |>
  left_join(count(all_genes, region, name = "universe"), by = "region") |>
  print()


# ---- 4. Overlap test, main threshold and sensitivity -----------------------

tests <- expand_grid(region = unique(genes$region), top_fraction = thresholds) |>
  pmap(\(region, top_fraction) {
    df <- genes[genes$region == region, ]
    overlap_test(is_top(df$rewiring_adj, top_fraction), df$is_de, df$stratum, n_perm) |>
      mutate(region = region, top_fraction = top_fraction, .before = 1)
  }) |>
  list_rbind() |>
  mutate(p_adj = p.adjust(p_perm, method = "BH"), .by = top_fraction) |>
  mutate(main = top_fraction == top_rewired, .after = top_fraction)


# ---- 5. Four groups of genes ------------------------------------------------

genes <- genes |>
  mutate(group = case_when(
    rewired & is_de ~ "rewired and DE",
    rewired ~ "rewired only",
    is_de ~ "DE only",
    TRUE ~ "no change"
  ))


# ---- 6. Stability: adjusted rewiring in each replicate ----------------------

# Same genes, model and top fraction as above, in each AD replicate alone
replicate_genes <- replicates |>
  semi_join(genes, by = c("region", "gene")) |>
  left_join(select(genes, region, gene, AveExpr, is_de), by = c("region", "gene")) |>
  mutate(mean_degree = (degree_ad + degree_ct) / 2) |>
  filter(!is.na(neighbor_change), mean_degree > 0) |>
  group_by(region, replicate) |>
  group_modify(\(df, key) {
    df |>
      adjust_rewiring() |>
      mutate(top = is_top(rewiring_adj, top_rewired))
  }) |>
  ungroup()

per_replicate <- replicate_genes |>
  summarise(odds_ratio = odds_ratio(top, is_de), .by = c(region, replicate))

# Robust rewired genes: in the top fraction in >= min_replicates replicates
genes <- genes |>
  left_join(
    summarise(replicate_genes, rewired_replicates = sum(top), .by = c(region, gene)),
    by = c("region", "gene")
  ) |>
  mutate(
    rewired_replicates = replace_na(rewired_replicates, 0L),
    robust_rewired = rewired_replicates >= min_replicates
  )

robust <- genes |>
  filter(robust_rewired) |>
  add_count(gene, name = "n_regions") |> # regions where the gene is robust rewired
  select(
    region, gene, symbol, rewiring_adj, neighbor_change, rewired_replicates,
    degree_ct, degree_ad, is_de, de_class, logFC, adj_p, n_regions
  ) |>
  arrange(region, desc(rewiring_adj))


vroom_write(select(genes, -stratum), snakemake@output[["genes"]])
vroom_write(tests, snakemake@output[["tests"]])
vroom_write(per_replicate, snakemake@output[["replicates"]])
vroom_write(robust, snakemake@output[["robust"]])

print(mutate(tests, across(where(is.double), \(v) signif(v, 3))), n = Inf, width = Inf)
print(count(genes, region, group) |> pivot_wider(names_from = group, values_from = n))
print(summarise(per_replicate, min_or = min(odds_ratio), max_or = max(odds_ratio), .by = region))
message("Robust rewired genes (top ", 100 * top_rewired, "% in >= ", min_replicates, " replicates):")
print(genes |> summarise(rewired = sum(rewired), robust = sum(robust_rewired), robust_and_rewired = sum(rewired & robust_rewired), .by = region))
