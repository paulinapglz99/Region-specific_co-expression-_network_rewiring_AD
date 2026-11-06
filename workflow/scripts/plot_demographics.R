# Demographics of every region x diagnosis group before building networks.
#
# Figure: n individuals, sex, age at death, APOE4 alleles, RIN, PMI and
# how many donors were sampled in more than one region.
# Table: summary per group + AD vs CT tests within each region.

source(snakemake@input[["theme"]])
library(dplyr)
library(tidyr)
library(purrr)
library(vroom)
library(patchwork)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

regions <- unlist(snakemake@params[["regions"]])
diagnoses <- unlist(snakemake@params[["diagnoses"]])


# ---- Read covariates of all studies -----------------------------------------

read_covariates <- function(file) {
  # IDs as text: Mayo donor IDs look like numbers
  vroom(file, col_types = cols(
    specimenID = "c", individualID = "c", age_death = "c", apoe4_allele = "c"
  )) |>
    mutate(study = sub("_Covariates_ages_censored.tsv$", "", basename(file))) |>
    select(study, specimenID, individualID, region = tissue, diagnosis,
           sex, age_death, apoe4_allele, RIN, pmi)
}

meta <- snakemake@input[["metadata"]] |>
  map(read_covariates) |>
  bind_rows() |>
  filter(region %in% regions, diagnosis %in% diagnoses) |>
  mutate(
    # Ages above 89 are censored as "90+" (counted as 90 for the median)
    age = as.numeric(if_else(age_death == "90+", "90", age_death)),
    age_group = case_when(age < 80 ~ "<80", age < 90 ~ "80-89", .default = "90+"),
    age_group = factor(age_group, levels = AGE_ORDER),
    apoe4_allele = factor(apoe4_allele, levels = APOE_ORDER)
  ) |>
  set_pub_levels()


# ---- Panels -----------------------------------------------------------------

counts <- meta |>
  group_by(study, region, diagnosis) |>
  summarise(n = n_distinct(individualID), .groups = "drop")

p_n <- ggplot(counts, aes(diagnosis, n, fill = diagnosis)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.4, size = PUB_LABEL_SIZE) +
  facet_regions() +
  scale_x_dx() +
  scale_fill_dx(guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Individuals", x = NULL, y = "n") +
  theme_pub()

# Stacked 100% bars of one categorical variable
plot_proportion <- function(var, fill_scale, title) {
  props <- meta |>
    count(study, region, diagnosis, .data[[var]]) |>
    group_by(study, region, diagnosis) |>
    mutate(prop = n / sum(n)) |>
    ungroup()

  ggplot(props, aes(diagnosis, prop, fill = .data[[var]])) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    facet_regions() +
    scale_x_dx() +
    fill_scale +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(title = title, x = NULL, y = "Individuals") +
    theme_pub()
}

p_sex <- plot_proportion("sex", scale_fill_sex(), "Sex")
p_age <- plot_proportion("age_group", scale_fill_age(), "Age at death")
p_apoe <- plot_proportion("apoe4_allele", scale_fill_apoe(), "APOE ε4")

# Distribution of one continuous variable
plot_box <- function(var, title) {
  ggplot(meta, aes(diagnosis, .data[[var]])) +
    geom_point(
      aes(colour = diagnosis),
      position = position_jitter(width = 0.15, seed = 1),
      size = 0.4, alpha = 0.5
    ) +
    geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, linewidth = 0.3) +
    facet_regions() +
    scale_x_dx() +
    scale_colour_dx(guide = "none") +
    labs(title = title, x = NULL, y = NULL) +
    theme_pub()
}

p_rin <- plot_box("RIN", "RNA integrity (RIN)")
p_pmi <- plot_box("pmi", "Post-mortem interval (h)")

# UpSet plot of the regions sampled per donor, one per study:
# top = donors in each combination, bottom = which regions form the combination
plot_upset <- function(s, title = NULL) {
  study_meta <- meta |> filter(study == s)

  # Combination of regions of each donor, e.g. "DLPFC + PCC"
  combos <- study_meta |>
    distinct(individualID, region) |>
    arrange(region) |>
    group_by(individualID) |>
    summarise(combo = paste(region, collapse = " + ")) |>
    count(combo) |>
    arrange(desc(n)) |>
    mutate(combo = factor(combo, levels = combo))

  study_regions <- REGION_ORDER[REGION_ORDER %in% study_meta$region]
  dots <- expand_grid(combo = combos$combo, region = study_regions) |>
    mutate(
      present = map2_lgl(as.character(combo), region, \(c, r) r %in% strsplit(c, " + ", fixed = TRUE)[[1]]),
      region = factor(region, levels = rev(study_regions))
    )

  bars <- ggplot(combos, aes(combo, n)) +
    geom_col(width = 0.7, fill = COL_NEUTRAL) +
    geom_text(aes(label = n), vjust = -0.4, size = PUB_LABEL_SIZE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = title, subtitle = s, x = NULL, y = "Donors (n)") +
    theme_pub() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  matrix <- ggplot(dots, aes(combo, region)) +
    geom_line(data = filter(dots, present), aes(group = combo), linewidth = 0.5) +
    geom_point(aes(colour = present), size = 2.2) +
    scale_colour_manual(values = c(`TRUE` = INK, `FALSE` = "grey85"), guide = "none") +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank())

  bars / matrix + plot_layout(heights = c(3, 1.2))
}

# Widths follow the number of region combinations of each study
p_overlap <- wrap_elements(
  full = (plot_upset("ROSMAP", "Regions sampled per donor") | plot_upset("Mayo", " ")) +
    plot_layout(widths = c(7, 3))
)


# ---- Figure -----------------------------------------------------------------

fig <- (p_n | p_sex) / (p_age | p_apoe) / (p_rin | p_pmi) / p_overlap +
  plot_layout(heights = c(1, 1, 1, 1.1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = PUB_TAG_SIZE, face = "bold"))

save_pub(fig, snakemake@output[["pdf"]], width = "double", height = FIG_MAX_HEIGHT)


# ---- Table ------------------------------------------------------------------

summary_tab <- meta |>
  group_by(study, region, diagnosis) |>
  summarise(
    n_samples = n(),
    n_individuals = n_distinct(individualID),
    pct_female = round(100 * mean(sex == "female"), 1),
    median_age = median(age),
    pct_age_90plus = round(100 * mean(age_group == "90+"), 1),
    pct_apoe4_carrier = round(100 * mean(apoe4_allele != "0"), 1),
    median_RIN = median(RIN),
    median_pmi = median(pmi),
    .groups = "drop"
  )

# AD vs CT within each region
tests <- meta |>
  group_by(study, region) |>
  summarise(
    p_sex = fisher.test(table(diagnosis, sex))$p.value,
    p_age_group = fisher.test(table(diagnosis, age_group))$p.value,
    p_apoe4 = fisher.test(table(diagnosis, apoe4_allele))$p.value,
    p_RIN = wilcox.test(RIN[diagnosis == "AD"], RIN[diagnosis == "CT"])$p.value,
    p_pmi = wilcox.test(pmi[diagnosis == "AD"], pmi[diagnosis == "CT"])$p.value,
    .groups = "drop"
  ) |>
  mutate(across(starts_with("p_"), \(p) signif(p, 3)))

table_out <- summary_tab |>
  left_join(tests, by = c("study", "region")) |>
  arrange(study, region, diagnosis)

vroom_write(table_out, snakemake@output[["table"]])

message(
  "Groups: ", nrow(summary_tab), " | samples: ", nrow(meta),
  " | individuals: ", n_distinct(meta$individualID)
)
print(table_out, n = Inf, width = Inf)
