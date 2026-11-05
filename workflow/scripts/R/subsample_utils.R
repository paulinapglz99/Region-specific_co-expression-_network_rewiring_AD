# Helpers shared by the subsampling scripts (subsample_*.R).
# Every pool is a data.table with one row per donor and a 'sex' column.

library(data.table)

# Reproducible seed that is different for each file: base seed + a number from a label
set_seed_for <- function(base_seed, label) {
  set.seed(base_seed + sum(utf8ToInt(label)))
}

# "rep01", "rep02", ...
rep_label <- function(r) sprintf("rep%02d", r)

# How many females and males to draw for a subsample of size n with a given
# female proportion. If one sex runs out, the other sex fills the gap.
# k = number of subsamples drawn at once without repeating donors (2 = disjoint pair).
sex_counts <- function(pool, n, female_prop, k = 1) {
  avail_f <- sum(pool$sex == "female") %/% k
  avail_m <- sum(pool$sex == "male") %/% k

  n_f <- min(round(n * female_prop), avail_f)
  n_m <- n - n_f
  if (n_m > avail_m) {
    n_m <- avail_m
    n_f <- n - n_m
  }
  if (n_f > avail_f) stop("Not enough donors to draw ", k, " x n = ", n)

  c(female = n_f, male = n_m)
}

# Draw n donors keeping the female proportion (as close as the pool allows)
sample_by_sex <- function(pool, n, female_prop) {
  counts <- sex_counts(pool, n, female_prop)
  rbind(
    pool[sex == "female"][sample(.N, counts[["female"]])],
    pool[sex == "male"][sample(.N, counts[["male"]])]
  )
}

# Draw two subsamples of n donors each, with no donor in common
sample_two_disjoint <- function(pool, n, female_prop) {
  counts <- sex_counts(pool, n, female_prop, k = 2)
  n_f <- counts[["female"]]
  n_m <- counts[["male"]]
  females <- pool[sex == "female"][sample(.N, 2 * n_f)]
  males <- pool[sex == "male"][sample(.N, 2 * n_m)]
  list(
    A = rbind(females[seq_len(n_f)], males[seq_len(n_m)]),
    B = rbind(females[n_f + seq_len(n_f)], males[n_m + seq_len(n_m)])
  )
}

# Size and % female of every network, for the log
summarise_networks <- function(subsamples) {
  subsamples[, .(n = .N, pct_female = round(100 * mean(sex == "female"), 1)), by = network]
}
