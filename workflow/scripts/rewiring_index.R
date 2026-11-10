# Reference curve and rewiring index (analyses 1-3).
#   1. J_ref(n): J between two disjoint halves of one group, as a function of n.
#      Model: J_ref(n) = a + b * log(n), fitted on the calibration group
#      (the curve entry with several n); the other curve groups validate it.
#   2. Rewiring index of each AD vs CT pair = J(AD, CT) / J_ref(n).
#   3. Analysis 2: order of the regions by index in each equal_n replicate,
#      compared with the order of the matched design (1 = most rewiring).
#   4. Analysis 6, prediction registered on 2026-11-06: TCX more rewired than CBE.
#      Holds if the mean index of TCX < that of CBE and the Wilcoxon test of
#      their 10 matched replicates gives p < 0.05. Also reported with the direct J.

library(dplyr)
library(purrr)
library(vroom)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

calibration <- snakemake@params[["calibration"]] # e.g. "DLPFC_AD"
equal_n_regions <- unlist(snakemake@params[["equal_n_regions"]])

metrics <- unlist(snakemake@input[["metrics"]]) |>
  map(\(f) vroom(f, col_types = cols(
    design = "c", table = "c", replicate = "c", net_x = "c", net_y = "c", .default = "d"
  ))) |>
  list_rbind()


# ---- 1. Reference curve ---------------------------------------------------------

curve <- metrics |>
  filter(design == "curve") |>
  mutate(
    group = sub("_n\\d+$", "", table),
    n = n_x,
    role = if_else(group == calibration, "calibration", "validation")
  )

fit <- lm(jaccard ~ log(n), data = filter(curve, role == "calibration"))
j_ref <- \(n) unname(predict(fit, newdata = tibble(n = n)))

curve <- curve |>
  mutate(j_ref = j_ref(n)) |>
  select(table, group, role, replicate, n, jaccard, j_ref)

message("J_ref(n) = ", round(coef(fit)[[1]], 4), " + ", round(coef(fit)[[2]], 4), " * log(n)")
message("Mean J(A, B) per curve table vs the fitted J_ref(n):")
print(
  curve |>
    group_by(table, role, n) |>
    summarise(mean_J = mean(jaccard), j_ref = first(j_ref), .groups = "drop") |>
    mutate(across(c(mean_J, j_ref), \(v) round(v, 3))),
  n = Inf
)


# ---- 2. Rewiring index ------------------------------------------------------------

index <- metrics |>
  filter(design != "curve") |>
  rename(region = table) |>
  mutate(n = n_x, j_ref = j_ref(n), index = jaccard / j_ref)

message("\nMean rewiring index per design and region:")
print(
  index |>
    group_by(design, region, n) |>
    summarise(mean_J = mean(jaccard), mean_index = mean(index), sd_index = sd(index), .groups = "drop") |>
    mutate(across(c(mean_J, mean_index, sd_index), \(v) round(v, 3))),
  n = Inf
)


# ---- 3. Region order (analysis 2) ------------------------------------------

# Order of the matched design: lowest mean index first (most rewiring)
reference <- index |>
  filter(design == "matched", region %in% equal_n_regions) |>
  group_by(region) |>
  summarise(index = mean(index)) |>
  arrange(index) |>
  pull(region)

region_order <- index |>
  filter(design == "equal_n") |>
  arrange(replicate, index) |>
  summarise(
    order = paste(region, collapse = " > "),
    kendall_tau = cor(seq_along(reference), match(reference, region), method = "kendall"),
    .by = replicate
  ) |>
  mutate(
    reference_order = paste(reference, collapse = " > "),
    same_order = order == reference_order
  )

message("\nOrder of the matched design: ", paste(reference, collapse = " > "))
message(
  "equal_n replicates with the same order: ", sum(region_order$same_order), " of ", nrow(region_order),
  " | mean Kendall tau: ", round(mean(region_order$kendall_tau), 3)
)


# ---- 4. Registered prediction (analysis 6) ---------------------------------

more <- snakemake@params[["more_rewired"]] # TCX
less <- snakemake@params[["less_rewired"]] # CBE
max_p <- snakemake@params[["max_p"]]

prediction <- c("index", "jaccard") |>
  map(\(metric) {
    x <- index[[metric]][index$design == "matched" & index$region == more]
    y <- index[[metric]][index$design == "matched" & index$region == less]
    test <- wilcox.test(x, y, exact = FALSE)
    tibble(
      metric = metric,
      more_rewired = more, less_rewired = less,
      mean_more = mean(x), mean_less = mean(y),
      n_more = length(x), n_less = length(y),
      wilcoxon_W = unname(test$statistic), p = test$p.value,
      # Lower index (or J) = AD and CT less alike = more rewiring
      holds = mean(x) < mean(y) & test$p.value < max_p
    )
  }) |>
  list_rbind()

message("\nRegistered prediction (", more, " more rewired than ", less, "):")
print(mutate(prediction, across(where(is.double), \(v) signif(v, 3))), width = Inf)

vroom_write(index, snakemake@output[["index"]])
vroom_write(curve, snakemake@output[["curve"]])
vroom_write(region_order, snakemake@output[["order"]])
vroom_write(prediction, snakemake@output[["prediction"]])
