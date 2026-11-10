# Helpers for the analyses that relate networks and differential expression
# (analyses 4, 5 and 7): DE table of one region and the overlap test that all
# three analyses share (odds ratio + stratified permutation null).

library(dplyr)
library(purrr)
library(vroom)

# DE file of one study among the DE files of all studies
# ("{study}_Differential_Expression_(diagnosis).tsv")
de_file_of <- function(files, study) {
  files[startsWith(basename(files), paste0(study, "_"))]
}

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

# Function that returns a row order with the rows shuffled only within each
# stratum. The strata are split once, because the shuffle runs once per
# permutation (thousands of times per test).
shuffler_within <- function(stratum) {
  groups <- unname(split(seq_along(stratum), stratum))
  positions <- unlist(groups, use.names = FALSE)
  function() {
    order <- integer(length(stratum))
    order[positions] <- unlist(map(groups, \(g) g[sample.int(length(g))]), use.names = FALSE)
    order
  }
}

# Two-sided permutation p-value: how extreme the observed value is in the null
perm_p <- function(observed, null) {
  null <- null[!is.na(null)]
  upper <- (1 + sum(null >= observed)) / (1 + length(null))
  lower <- (1 + sum(null <= observed)) / (1 + length(null))
  min(1, 2 * min(upper, lower))
}

# Odds ratio of the 2 x 2 table of two logical vectors (0 or Inf when a cell is empty)
odds_ratio <- function(x, y) {
  (sum(x & y) * sum(!x & !y)) / (sum(x & !y) * sum(!x & y))
}

# Overlap test shared by analyses 4, 5 and 7: does gene set `in_set` contain
# more genes with `label` than expected among genes of the same stratum?
#   null:      `label` shuffled among genes of the same stratum, n_perm times
#   statistic: genes in both. The shuffle keeps the size of the set and the
#              number of labels, so the odds ratio only depends on this count
#   expected:  mean of the genes in both in the null
#   effect:    odds ratio of the 2 x 2 table in_set x label
# All three vectors have one element per gene of the universe.
overlap_test <- function(in_set, label, stratum, n_perm) {
  observed <- sum(in_set & label)
  shuffle <- shuffler_within(stratum)
  null <- map_int(seq_len(n_perm), \(i) sum(in_set & label[shuffle()]))
  tibble(
    n_universe = length(in_set),
    n_set = sum(in_set),
    n_label = sum(label),
    n_both = observed,
    expected_both = mean(null),
    odds_ratio = odds_ratio(in_set, label),
    p_perm = if (sum(in_set) == 0 || sum(label) == 0) NA_real_ else perm_p(observed, null),
    n_perm = n_perm
  )
}
