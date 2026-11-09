# Analysis 4, steps 2-6: are the most rewired genes the DE genes? (one region)
#   2. Degree filter: genes with >= min_degree neighbors in AD (mean) or in CT
#   3. Adjusted rewiring = residual of neighbor_change ~ log(mean degree) + AveExpr
#   4. Tests, with a null that shuffles the DE labels within strata of
#      AveExpr x mean degree (deciles):
#        de_effect   mean adjusted rewiring of DE genes - of non-DE genes (+ 95% CI of lm)
#        overlap     log odds ratio of "rewired" (top fraction) vs DE
#        spearman    Spearman rho of adjusted rewiring vs |t|
#        direction   mean degree_change of DE up - of DE down
#   5. Four groups of genes: rewired and DE, rewired only, DE only, no change
#   6. Robustness: de_effect and overlap in each replicate

source(snakemake@input[["de_utils"]])
library(tidyr)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

region <- snakemake@wildcards[["region"]]
min_degree <- snakemake@params[["min_degree"]]
top_rewired <- snakemake@params[["top_rewired"]]
n_perm <- snakemake@params[["n_perm"]]
bins <- snakemake@params[["strata"]]
set.seed(snakemake@params[["seed"]])

all_genes <- vroom(snakemake@input[["genes"]], show_col_types = FALSE)
replicates <- vroom(snakemake@input[["replicates"]], col_types = cols(gene = "c", replicate = "c"))

# Adjusted rewiring: what is left of neighbor_change after degree and expression
adjust_rewiring <- function(df) {
  model <- lm(neighbor_change ~ log(mean_degree) + AveExpr, data = df)
  df |> mutate(rewiring_adj = residuals(model))
}
flag_rewired <- function(df) {
  df |> mutate(rewired = rewiring_adj >= quantile(rewiring_adj, 1 - top_rewired))
}


# ---- 2-3. Filter and adjusted rewiring --------------------------------------

genes <- all_genes |>
  filter(degree_ad >= min_degree | degree_ct >= min_degree, !is.na(AveExpr)) |>
  mutate(mean_degree = (degree_ad + degree_ct) / 2) |>
  adjust_rewiring() |>
  flag_rewired() |>
  add_strata(c("AveExpr", "mean_degree"), bins)

message(
  region, ": ", nrow(genes), " of ", nrow(all_genes), " genes pass the degree filter (>= ", min_degree,
  ") | DE: ", sum(genes$is_de), " (up ", sum(genes$de_class == "DE up"), ", down ", sum(genes$de_class == "DE down"), ")"
)


# ---- 4. Tests -------------------------------------------------------------

statistics <- function(df) {
  c(
    de_effect = mean(df$rewiring_adj[df$is_de]) - mean(df$rewiring_adj[!df$is_de]),
    overlap = log_odds_ratio(df$rewired, df$is_de),
    spearman = cor(df$rewiring_adj, abs(df$t), method = "spearman"),
    direction = mean(df$degree_change[df$de_class == "DE up"]) - mean(df$degree_change[df$de_class == "DE down"])
  )
}

observed <- statistics(genes)

# Null: the DE columns move together to other genes of the same stratum
de_columns <- c("is_de", "de_class", "t")
null <- map(seq_len(n_perm), \(i) {
  shuffled <- genes
  shuffled[de_columns] <- genes[shuffle_within(genes$stratum), de_columns]
  statistics(shuffled)
}) |>
  list_transpose()

ci <- confint(lm(rewiring_adj ~ is_de, data = genes))["is_deTRUE", ]

tests <- tibble(
  region = region,
  test = names(observed),
  estimate = unname(observed),
  ci_low = c(ci[[1]], NA, NA, NA),
  ci_high = c(ci[[2]], NA, NA, NA),
  null_mean = map_dbl(names(observed), \(s) mean(null[[s]], na.rm = TRUE)),
  p_perm = map_dbl(names(observed), \(s) perm_p(observed[[s]], null[[s]])),
  n_genes = nrow(genes),
  n_de = sum(genes$is_de),
  n_perm = n_perm
)

direction <- genes |>
  group_by(region = region, de_class) |>
  summarise(n = n(), mean_degree_change = mean(degree_change), .groups = "drop")


# ---- 5. Four groups of genes ------------------------------------------------

genes <- genes |>
  mutate(group = case_when(
    rewired & is_de ~ "rewired and DE",
    rewired ~ "rewired only",
    is_de ~ "DE only",
    TRUE ~ "no change"
  ))


# ---- 6. Robustness: the tests in each replicate --------------------------

per_replicate <- replicates |>
  semi_join(genes, by = "gene") |>
  left_join(select(genes, gene, AveExpr, is_de), by = "gene") |>
  mutate(mean_degree = (degree_ad + degree_ct) / 2) |>
  filter(!is.na(neighbor_change)) |>
  group_by(replicate) |>
  group_modify(\(df, key) {
    df <- df |> adjust_rewiring() |> flag_rewired()
    tibble(
      de_effect = mean(df$rewiring_adj[df$is_de]) - mean(df$rewiring_adj[!df$is_de]),
      overlap = log_odds_ratio(df$rewired, df$is_de)
    )
  }) |>
  ungroup() |>
  mutate(region = region, .before = 1)


vroom_write(mutate(genes, region = region, .before = 1), snakemake@output[["genes"]])
vroom_write(tests, snakemake@output[["tests"]])
vroom_write(direction, snakemake@output[["direction"]])
vroom_write(per_replicate, snakemake@output[["replicates"]])

print(mutate(tests, across(where(is.double), \(v) signif(v, 3))))
print(count(genes, group))
