# Figure of analyses 1-3: network comparisons between AD and CT.
#   a. Edge Jaccard AD vs CT (matched n)
#   b. Reference curve J_ref(n)
#   c. Rewiring index (matched n)
#   d. Rewiring index with n = 90 in every network (equal_n)
#   e. Spearman rho of degree and NMI of communities (matched n)
#   f. Mean degree of the AD and CT networks (matched n)

library(dplyr)
library(tidyr)
library(purrr)
library(vroom)
library(patchwork)
source(snakemake@input[["theme"]])

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

study_of <- unlist(snakemake@params[["study_of"]]) # region -> study

index <- vroom(snakemake@input[["index"]], col_types = cols(design = "c", region = "c", replicate = "c", net_x = "c", net_y = "c", .default = "d")) |>
  mutate(study = study_of[region]) |>
  set_pub_levels()
curve <- vroom(snakemake@input[["curve"]], col_types = cols(table = "c", group = "c", role = "c", replicate = "c", .default = "d"))
region_order <- vroom(snakemake@input[["order"]], col_types = cols(replicate = "c", order = "c", reference_order = "c", same_order = "l", .default = "d"))

matched <- filter(index, design == "matched")
equal_n <- filter(index, design == "equal_n")


# ---- Shared pieces -------------------------------------------------------------

# 10 replicates as small points + mean and SD in black
replicate_points <- function(data, y) {
  summary <- data |>
    group_by(study, region) |>
    summarise(mean = mean({{ y }}), sd = sd({{ y }}), .groups = "drop")
  list(
    geom_point(
      data = data, aes(x = region, y = {{ y }}, colour = region),
      position = position_jitter(width = 0.15, height = 0, seed = 1), size = 0.7, alpha = 0.8
    ),
    geom_errorbar(
      data = summary, aes(x = region, ymin = mean - sd, ymax = mean + sd),
      width = 0.25, linewidth = PUB_LINE_WIDTH, inherit.aes = FALSE
    ),
    geom_point(data = summary, aes(x = region, y = mean), size = 1.2, inherit.aes = FALSE),
    scale_colour_manual(values = COL_REGION, guide = "none"),
    facet_grid(. ~ study, scales = "free_x", space = "free_x"),
    labs(x = NULL)
  )
}


# ---- Panels ----------------------------------------------------------------------

p_jaccard <- ggplot() +
  replicate_points(matched, jaccard) +
  labs(title = "AD vs CT, same n", y = "Edge Jaccard (J)") +
  theme_pub()

curve_means <- curve |>
  group_by(table, group, role, n) |>
  summarise(jaccard = mean(jaccard), .groups = "drop") |>
  # Validation labels: the smallest n above-left of its point, the rest below-right,
  # so the labels of close points do not overlap
  arrange(role, n) |>
  group_by(role) |>
  mutate(
    label_hjust = if_else(row_number() == 1, 1.2, -0.2),
    label_vjust = if_else(row_number() == 1, -0.6, 1.6)
  ) |>
  ungroup()
p_curve <- ggplot(curve, aes(x = n, y = jaccard)) +
  geom_smooth(
    data = filter(curve, role == "calibration"),
    method = "lm", formula = y ~ log(x), # same model as rewiring_index.R
    se = FALSE, fullrange = TRUE, colour = INK_SECONDARY, linewidth = 0.4
  ) +
  geom_point(aes(shape = role), size = 0.7, alpha = 0.6, colour = INK_SECONDARY) +
  geom_point(data = curve_means, aes(shape = role), size = 1.6, fill = "white") +
  geom_text(
    data = filter(curve_means, role == "validation"),
    aes(label = sub("_", " ", group), hjust = label_hjust, vjust = label_vjust), size = PUB_LABEL_SIZE
  ) +
  scale_shape_manual(values = c(calibration = 16, validation = 21), name = NULL) +
  scale_x_continuous(limits = c(55, 145), breaks = c(60, 80, 100, 120, 140)) +
  labs(title = "Reference J: two halves of one group", x = "n (samples per network)", y = expression(J[ref](n))) +
  theme_pub() +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1))

p_index <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = PUB_LINE_WIDTH, colour = COL_NEUTRAL) +
  replicate_points(matched, index) +
  labs(title = "Rewiring index, same n", y = expression(J(AD, CT) / J[ref](n))) +
  theme_pub()

same <- sum(region_order$same_order)
p_equal_n <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = PUB_LINE_WIDTH, colour = COL_NEUTRAL) +
  replicate_points(equal_n, index) +
  labs(
    title = "Rewiring index, n = 90 in every network",
    subtitle = sprintf(
      "Same region order as in c: %d of %d replicates (mean Kendall τ = %.2f)",
      same, nrow(region_order), mean(region_order$kendall_tau)
    ),
    y = expression(J(AD, CT) / J[ref](90))
  ) +
  theme_pub() +
  theme(plot.subtitle = element_text(size = PUB_BASE_SIZE - 1, colour = INK_SECONDARY))

structure_long <- matched |>
  select(study, region, replicate, spearman_degree, nmi) |>
  pivot_longer(c(spearman_degree, nmi), names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric, levels = c("spearman_degree", "nmi"), labels = c("Degree (Spearman ρ)", "Communities (NMI)")))
p_structure <- ggplot(structure_long, aes(x = region, y = value, colour = region)) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1), size = 0.7, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", colour = INK, size = 1.2) +
  scale_colour_manual(values = COL_REGION, guide = "none") +
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "AD vs CT: degree and communities, same n", x = NULL, y = "Similarity") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

degree_long <- bind_rows(
  matched |> transmute(study, region, dx = "AD", mean_degree = mean_degree_x),
  matched |> distinct(study, region, mean_degree_y) |> transmute(study, region, dx = "CT", mean_degree = mean_degree_y)
) |>
  mutate(dx = factor(dx, levels = DX_ORDER))
p_degree <- ggplot(degree_long, aes(x = region, y = mean_degree, colour = dx)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6, seed = 1), size = 0.8) +
  scale_colour_dx() +
  facet_grid(. ~ study, scales = "free_x", space = "free_x") +
  labs(title = "Network density, same n", x = NULL, y = "Mean degree") +
  theme_pub()


# ---- Figure ---------------------------------------------------------------------

figure <- (p_jaccard | p_curve) / (p_index | p_equal_n) / (p_structure | p_degree) +
  plot_annotation(tag_levels = "a")

save_pub(figure, snakemake@output[["pdf"]], width = "double", height = 190)
message("Saved ", snakemake@output[["pdf"]])
