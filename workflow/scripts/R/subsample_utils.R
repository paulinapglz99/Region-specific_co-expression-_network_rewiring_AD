# Helpers shared by the subsampling scripts (subsample_*.R).
# Every pool is a tibble with one row per donor and a 'sex' column.

library(dplyr)
library(purrr)
library(vroom)

# Read a selected-samples table (results/samples/...).
# IDs as text: Mayo donor IDs look like numbers.
read_samples <- function(file) {
  vroom(file, col_types = cols(
    specimenID = "c", individualID = "c", age_death = "c", apoe4_allele = "c"
  ))
}

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
  bind_rows(
    pool |> filter(sex == "female") |> slice_sample(n = counts[["female"]]),
    pool |> filter(sex == "male") |> slice_sample(n = counts[["male"]])
  )
}

# Draw two subsamples (A and B) of n donors each, with no donor in common
sample_two_disjoint <- function(pool, n, female_prop) {
  counts <- sex_counts(pool, n, female_prop, k = 2)
  n_f <- counts[["female"]]
  n_m <- counts[["male"]]

  # Draw twice the donors needed, then give the first half to A and the rest to B
  females <- pool |> filter(sex == "female") |> slice_sample(n = 2 * n_f)
  males <- pool |> filter(sex == "male") |> slice_sample(n = 2 * n_m)

  list(
    A = bind_rows(slice_head(females, n = n_f), slice_head(males, n = n_m)),
    B = bind_rows(slice_tail(females, n = n_f), slice_tail(males, n = n_m))
  )
}

# Size and % female of every network, for the log
summarise_networks <- function(subsamples, ...) {
  subsamples |>
    group_by(network, ...) |>
    summarise(n = n(), pct_female = round(100 * mean(sex == "female"), 1), .groups = "drop")
}
