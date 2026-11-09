# Helpers for the analyses that relate networks and differential expression
# (analyses 4, 5 and 7): DE table of one region and stratified permutation null.

library(dplyr)
library(purrr)
library(vroom)

# Sage's DE results of one region (comparison "AD_{region} - CT_{region}").
# DE gene = adj.P.Val < max_adj_p and |FC| > fold_change
read_de <- function(file, region, max_adj_p, fold_change) {
  vroom(file, col_types = cols(Comparison = "c", ensembl_gene_id = "c", hgnc_symbol = "c")) |>
    filter(Comparison == sprintf("AD_%s - CT_%s", region, region)) |>
    transmute(gene = ensembl_gene_id, symbol = hgnc_symbol, logFC, t, AveExpr, adj_p = adj.P.Val) |>
    distinct(gene, .keep_all = TRUE) |>
    mutate(
      is_de = adj_p < max_adj_p & abs(logFC) > log2(fold_change),
      de_class = case_when(!is_de ~ "not DE", logFC > 0 ~ "DE up", TRUE ~ "DE down")
    )
}

# Stratum of each gene: bins (deciles by default) of one or more variables
add_strata <- function(df, vars, bins) {
  df |>
    mutate(across(all_of(vars), \(v) ntile(v, bins), .names = "bin_{.col}")) |>
    tidyr::unite("stratum", starts_with("bin_"), sep = "_", remove = TRUE)
}

# Row order with the rows shuffled only within each stratum
shuffle_within <- function(stratum) {
  tibble(i = seq_along(stratum), stratum = stratum) |>
    mutate(j = i[sample.int(n())], .by = stratum) |>
    pull(j)
}

# Two-sided permutation p-value: how extreme the observed value is in the null
perm_p <- function(observed, null) {
  null <- null[!is.na(null)]
  upper <- (1 + sum(null >= observed)) / (1 + length(null))
  lower <- (1 + sum(null <= observed)) / (1 + length(null))
  min(1, 2 * min(upper, lower))
}

# Log odds ratio of two logical vectors (0.5 added to every cell, so it is
# finite when a cell is empty)
log_odds_ratio <- function(x, y) {
  log(((sum(x & y) + 0.5) * (sum(!x & !y) + 0.5)) / ((sum(x & !y) + 0.5) * (sum(!x & y) + 0.5)))
}
