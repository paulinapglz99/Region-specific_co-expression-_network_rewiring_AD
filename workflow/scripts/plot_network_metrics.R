# Figures of analyses 1-3 and 6: network comparisons between AD and CT.
# Each figure goes to its own folder, with one file per panel (save_figure).
# Figure 1:
#   a. Schematic made by hand (figure1/1a.pdf): not made here. The whole
#      figure leaves an empty slot for it
#   b. Reference curve J_ref(n)
#   c. Rewiring index (matched n), with the registered prediction TCX vs CBE
#   d. Spearman rho of degree and NMI of communities (matched n)
#   e. Mean degree of the AD and CT networks (matched n)
# Supplementary figure S3:
#   a. Edge Jaccard AD vs CT (matched n), with the prediction on the direct J
#   b. Rewiring index with n = 90 in every network (equal_n)

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
prediction <- vroom(snakemake@input[["prediction"]], col_types = cols(metric = "c", more_rewired = "c", less_rewired = "c", holds = "l", .default = "d"))

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


# Registered prediction (analysis 6): bracket between the two regions with the
# Wilcoxon p, drawn above the highest replicate of `y`
prediction_bracket <- function(data, y, metric) {
  pred <- prediction[prediction$metric == metric, ]
  top <- max(data[[y]][data$region %in% c(pred$more_rewired, pred$less_rewired)])
  height <- top + 0.08 * diff(range(data[[y]]))
  bracket <- tibble(
    study = factor(study_of[[pred$more_rewired]], levels = STUDY_ORDER),
    x = pred$more_rewired, xend = pred$less_rewired, y = height,
    label = sprintf("Registered: %s > %s\nWilcoxon p = %.2g", pred$more_rewired, pred$less_rewired, pred$p)
  )
  list(
    geom_segment(data = bracket, aes(x = x, xend = xend, y = y, yend = y), linewidth = PUB_LINE_WIDTH, inherit.aes = FALSE),
    geom_text(
      data = bracket, aes(x = 1.5, y = y, label = label), vjust = -0.3,
      size = PUB_LABEL_SIZE, lineheight = 0.9, inherit.aes = FALSE
    ),
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3)))
  )
}


# ---- Panels ----------------------------------------------------------------------

p_jaccard <- ggplot() +
  replicate_points(matched, jaccard) +
  prediction_bracket(matched, "jaccard", "jaccard") +
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
  prediction_bracket(matched, "index", "index") +
  labs(title = "Rewiring index, same n", y = expression(J(AD, CT) / J[ref](n))) +
  theme_pub()

same <- sum(region_order$same_order)
p_equal_n <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = PUB_LINE_WIDTH, colour = COL_NEUTRAL) +
  replicate_points(equal_n, index) +
  labs(
    title = "Rewiring index, n = 90 in every network",
    subtitle = sprintf(
      "Same region order as with matched n: %d of %d replicates (mean Kendall τ = %.2f)",
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

# Empty slot for panel a (hand-made schematic, figure1/1a.pdf)
p_schematic <- ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = NA, colour = COL_NEUTRAL, linetype = "dashed") +
  annotate("text", x = 0.5, y = 0.5, label = "Panel a: hand-made schematic\n(figure1/1a.pdf)", size = PUB_LABEL_SIZE, colour = INK_SECONDARY) +
  theme_void(base_size = PUB_BASE_SIZE)

# Rows: schematic | reference curve; rewiring index (main panel, full width);
# structure | density
figure <- (p_schematic | p_curve) / p_index / (p_structure | p_degree) +
  plot_layout(heights = c(1, 0.9, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = PUB_TAG_SIZE, face = "bold"))
save_figure(
  figure, list(b = p_curve, c = p_index, d = p_structure, e = p_degree),
  dir = dirname(snakemake@output[["figure"]]), id = "1", height = 200
)

supplementary <- (p_jaccard | p_equal_n) +
  plot_annotation(tag_levels = "a")
save_figure(
  supplementary, list(a = p_jaccard, b = p_equal_n),
  dir = dirname(snakemake@output[["supplementary"]]), id = "S3", height = 70
)
