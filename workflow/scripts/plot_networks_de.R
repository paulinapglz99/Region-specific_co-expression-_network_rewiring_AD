# Figure 2 (analyses 4 and 5): the network finds genes that DE does not see.
#   a. Rewired genes (top 10%) of each region: fraction that is DE ("rewired
#      and DE") and not DE ("rewired only"), with the fraction expected by
#      chance (stratified null) and the odds ratio with its BH-adjusted p
#   b. GO terms enriched among the "rewired only" genes (top 5 per region)
#   c. Hubs: fraction of DE genes among shared, gained and lost hubs, observed
#      vs expected (null with the same expression)
# Figure S4 (robustness of analysis 4):
#   a. Odds ratio rewired x DE with the top 5%, 10% and 20% as "rewired"
#   b. Odds ratio of each AD replicate, next to the one of the mean (top 10%)

library(dplyr)
library(tidyr)
library(purrr)
library(vroom)
library(patchwork)
source(snakemake@input[["theme"]])

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

study_of <- unlist(snakemake@params[["study_of"]])
max_adj_p <- snakemake@params[["max_adj_p"]]

add_study <- \(df) df |> mutate(study = study_of[region]) |> set_pub_levels()

all_rewiring_tests <- vroom(snakemake@input[["rewiring_tests"]], col_types = cols(region = "c", main = "l")) |>
  add_study()
rewiring_tests <- filter(all_rewiring_tests, main)
rewiring_replicates <- vroom(snakemake@input[["rewiring_replicates"]], col_types = cols(region = "c", replicate = "c")) |>
  add_study()
go <- vroom(snakemake@input[["go"]], col_types = cols(region = "c", term = "c", description = "c", .default = "?")) |>
  add_study()
hubs_tests <- vroom(snakemake@input[["hubs_tests"]], col_types = cols(region = "c", hub_type = "c")) |>
  add_study() |>
  mutate(hub_type = factor(hub_type, levels = c("shared", "gained", "lost"), labels = c("Shared", "Gained in AD", "Lost in AD")))

COL_GROUP <- c("rewired and DE" = NPG[["red"]], "rewired only" = NPG[["blue"]])
p_label <- \(or, p) sprintf("OR %.2f\nq = %.2g", or, p)


# ---- a. Rewired genes: DE or not ---------------------------------------------

bars <- rewiring_tests |>
  transmute(study, region, `rewired and DE` = n_both / n_set, `rewired only` = 1 - n_both / n_set) |>
  pivot_longer(c(`rewired and DE`, `rewired only`), names_to = "group", values_to = "fraction") |>
  mutate(group = factor(group, levels = rev(names(COL_GROUP))))

p_groups <- ggplot(bars, aes(x = region, y = fraction, fill = group)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    data = rewiring_tests, aes(x = region, ymin = expected_both / n_set, ymax = expected_both / n_set),
    width = 0.8, linewidth = 0.5, colour = INK, inherit.aes = FALSE
  ) +
  geom_text(
    data = rewiring_tests, aes(x = region, y = 1, label = p_label(odds_ratio, p_adj)),
    vjust = -0.2, size = PUB_LABEL_SIZE, lineheight = 0.9, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = COL_GROUP, breaks = names(COL_GROUP), name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.2)), breaks = seq(0, 1, 0.25)) +
  facet_grid(. ~ study, scales = "free_x", space = "free_x") +
  labs(
    title = "Rewired genes (top 10%): DE or not",
    subtitle = "Black line: DE fraction expected by chance (same expression and degree)",
    x = NULL, y = "Rewired genes"
  ) +
  theme_pub() +
  theme(plot.subtitle = element_text(size = PUB_BASE_SIZE - 1, colour = INK_SECONDARY))


# ---- b. GO terms of the "rewired only" genes --------------------------------

top_terms <- go |>
  filter(p_adj < max_adj_p) |>
  slice_min(p_adj, n = 5, with_ties = FALSE, by = region)

p_go <- if (nrow(top_terms) == 0) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = sprintf("No GO term with q < %.2g\nin any region", max_adj_p), size = PUB_LABEL_SIZE) +
    labs(title = "GO terms of the \"rewired only\" genes") +
    theme_void(base_size = PUB_BASE_SIZE) +
    theme(plot.title = element_text(size = PUB_TITLE_SIZE, face = "bold"))
} else {
  # Every region's value for the terms that are in the top 5 of some region
  shown <- go |>
    filter(term %in% top_terms$term) |>
    mutate(description = factor(description, levels = rev(unique(top_terms$description[order(top_terms$region, top_terms$p_adj)]))))
  ggplot(shown, aes(x = region, y = description, size = count, colour = -log10(p_adj))) +
    geom_point() +
    geom_point(data = filter(shown, p_adj >= max_adj_p), shape = 1, colour = COL_NEUTRAL) +
    scale_colour_gradient(low = "grey80", high = NPG[["blue"]], name = expression(-log[10](q))) +
    scale_size_area(max_size = 3, name = "Genes") +
    scale_x_discrete(drop = FALSE) +
    labs(title = "GO terms of the \"rewired only\" genes", x = NULL, y = NULL) +
    theme_pub() +
    theme(legend.position = "right", legend.direction = "vertical")
}


# ---- c. Hubs: fraction DE, observed vs expected ---------------------------

hub_points <- hubs_tests |>
  mutate(observed = n_both / n_set, expected = expected_both / n_set) |>
  filter(n_set > 0)

p_hubs <- ggplot(hub_points, aes(x = region)) +
  geom_linerange(aes(ymin = expected, ymax = observed), linewidth = PUB_LINE_WIDTH, colour = COL_NEUTRAL) +
  geom_point(aes(y = expected), shape = 21, fill = "white", size = 1.4) +
  geom_point(aes(y = observed, colour = region), size = 1.6) +
  geom_text(aes(y = pmax(observed, expected), label = paste0("n = ", n_set, if_else(p_adj < 0.05, "\n*", ""))),
            vjust = -0.4, size = PUB_LABEL_SIZE, lineheight = 0.8) +
  scale_colour_manual(values = COL_REGION, guide = "none") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.05, 0.25))) +
  facet_grid(. ~ hub_type) +
  labs(
    title = "Hubs: fraction of DE genes",
    subtitle = "Filled: observed · open: expected by chance (same expression) · * q < 0.05",
    x = NULL, y = "DE genes among the hubs"
  ) +
  theme_pub() +
  theme(
    plot.subtitle = element_text(size = PUB_BASE_SIZE - 1, colour = INK_SECONDARY),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# ---- Figure ---------------------------------------------------------------------

figure <- (p_groups | p_go) / p_hubs +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "a")

save_figure(
  figure, list(a = p_groups, b = p_go, c = p_hubs),
  dir = dirname(snakemake@output[["figure"]]), id = "2", height = 150,
  panel_width = c(a = "single", b = "single", c = "double"),
  panel_height = c(a = 85, b = 85, c = 70)
)


# ---- Figure S4: robustness of analysis 4 ----------------------------------

or_axis <- list(
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = PUB_LINE_WIDTH, colour = COL_NEUTRAL),
  scale_y_log10(breaks = c(0.25, 0.5, 0.7, 1, 1.5, 2, 3, 4)),
  labs(y = "Odds ratio, rewired × DE (log)")
)

# a. Same test with other top fractions
thresholds <- all_rewiring_tests |>
  mutate(
    threshold = factor(paste0(100 * top_fraction, "%"), levels = paste0(100 * sort(unique(top_fraction)), "%")),
    significant = if_else(p_adj < 0.05, "q < 0.05", "q ≥ 0.05")
  )
p_thresholds <- ggplot(thresholds, aes(x = threshold, y = odds_ratio, group = region, colour = region)) +
  or_axis +
  geom_line(linewidth = 0.4) +
  geom_point(aes(shape = significant), size = 1.6, fill = "white") +
  scale_colour_manual(values = COL_REGION, name = NULL) +
  scale_shape_manual(values = c("q < 0.05" = 16, "q ≥ 0.05" = 21), name = NULL) +
  labs(title = "Other definitions of \"rewired\"", x = "Rewired genes = top ...") +
  theme_pub()

# b. Each AD replicate alone vs the mean of the 10 (main test)
p_replicates <- ggplot(rewiring_replicates, aes(x = region, y = odds_ratio)) +
  or_axis +
  geom_point(
    aes(colour = region),
    position = position_jitter(width = 0.15, height = 0, seed = 1), size = 0.8, alpha = 0.8
  ) +
  geom_point(data = rewiring_tests, shape = 23, size = 2, fill = "white", colour = INK) +
  scale_colour_manual(values = COL_REGION, guide = "none") +
  facet_grid(. ~ study, scales = "free_x", space = "free_x") +
  labs(
    title = "Each AD replicate alone (top 10%)",
    subtitle = "Dots: one replicate · diamond: mean of the 10 (main test)",
    x = NULL
  ) +
  theme_pub() +
  theme(plot.subtitle = element_text(size = PUB_BASE_SIZE - 1, colour = INK_SECONDARY))

supplementary <- (p_thresholds | p_replicates) +
  plot_annotation(tag_levels = "a")
save_figure(
  supplementary, list(a = p_thresholds, b = p_replicates),
  dir = dirname(snakemake@output[["supplementary"]]), id = "S4", height = 75
)
