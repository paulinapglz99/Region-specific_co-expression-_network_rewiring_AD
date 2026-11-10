# Shared style for every figure of the paper.
#
# Change fonts, sizes, colors, group order or figure sizes HERE and all plots
# of the workflow follow. Plot scripts load it with:
#   source(snakemake@input[["theme"]])
# (the theme is an input of each plot rule, so editing it re-makes the plots)

library(dplyr)
library(ggplot2)


# ---- Fonts and sizes --------------------------------------------------------

PUB_FONT <- "sans" # e.g. "Helvetica" or "Arial" if the journal asks for it
PUB_BASE_SIZE <- 7 # pt; most journals want 5-8 pt at final print size
PUB_TITLE_SIZE <- 8 # pt; panel titles
PUB_TAG_SIZE <- 9 # pt; panel letters (a, b, c...)
PUB_LINE_WIDTH <- 0.3 # mm; axis lines and ticks
PUB_LABEL_SIZE <- (PUB_BASE_SIZE - 1) / .pt # geom_text() size (in mm)

INK <- "black"
INK_SECONDARY <- "grey30"
COL_NEUTRAL <- "grey60"


# ---- Figure sizes (mm, Nature column widths) --------------------------------

FIG_WIDTH <- c(single = 89, onehalf = 120, double = 183)
FIG_MAX_HEIGHT <- 247
FIG_DPI <- 600


# ---- Order and labels of variables ------------------------------------------

STUDY_ORDER <- c("ROSMAP", "Mayo")
REGION_ORDER <- c("DLPFC", "ACC", "PCC", "TCX", "CBE")
DX_ORDER <- c("CT", "AD")
DX_LABELS <- c(CT = "Control", AD = "AD") # legends
DX_AXIS_LABELS <- c(CT = "CT", AD = "AD") # short, for narrow axes
SEX_ORDER <- c("female", "male")
SEX_LABELS <- c(female = "Female", male = "Male")
AGE_ORDER <- c("<80", "80-89", "90+")
APOE_ORDER <- c("0", "1", "2")


# ---- Colors: Nature Publishing Group palette (same hex as ggsci::pal_npg) ---

NPG <- c(
  red = "#E64B35", cyan = "#4DBBD5", green = "#00A087", blue = "#3C5488",
  salmon = "#F39B7F", lavender = "#8491B4", mint = "#91D1C2",
  darkred = "#DC0000", brown = "#7E6148", tan = "#B09C85"
)

# Pairs chosen to stay distinguishable with color blindness
COL_DX <- c(CT = NPG[["blue"]], AD = NPG[["red"]])
COL_SEX <- c(female = NPG[["green"]], male = "#D9772B") # orange (not in NPG)

# Ordered categories: light -> dark tints of one NPG hue
COL_AGE <- c("#9AA5C8", "#6676A8", "#3C5488") # NPG blue
COL_APOE <- c("#66C2A8", "#00957A", "#00513F") # NPG green

# Regions: 5 colors cannot all be told apart with color blindness,
# so always label regions too (axis, facet or text), never color alone
COL_REGION <- c(
  DLPFC = NPG[["cyan"]], ACC = NPG[["green"]], PCC = NPG[["lavender"]],
  TCX = NPG[["salmon"]], CBE = NPG[["brown"]]
)


# ---- Theme -------------------------------------------------------------------

theme_pub <- function(base_size = PUB_BASE_SIZE, base_family = PUB_FONT) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(colour = INK),
      axis.text = element_text(size = base_size, colour = INK_SECONDARY),
      axis.title = element_text(size = base_size),
      axis.line = element_line(linewidth = PUB_LINE_WIDTH),
      axis.ticks = element_line(linewidth = PUB_LINE_WIDTH),
      plot.title = element_text(size = PUB_TITLE_SIZE, face = "bold"),
      plot.tag = element_text(size = PUB_TAG_SIZE, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(3, "mm"),
      legend.position = "top",
      legend.justification = "left",
      panel.spacing = unit(2, "mm"),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )
}


# ---- Scales and facets ------------------------------------------------------

scale_fill_dx <- function(...) {
  scale_fill_manual(values = COL_DX, labels = DX_LABELS, name = "Diagnosis", ...)
}
scale_colour_dx <- function(...) {
  scale_colour_manual(values = COL_DX, labels = DX_LABELS, name = "Diagnosis", ...)
}
scale_fill_sex <- function(...) {
  scale_fill_manual(values = COL_SEX, labels = SEX_LABELS, name = "Sex", ...)
}
scale_fill_age <- function(...) {
  scale_fill_manual(values = COL_AGE, name = "Years", drop = FALSE, ...)
}
scale_fill_apoe <- function(...) {
  scale_fill_manual(values = COL_APOE, name = "Alleles", drop = FALSE, ...)
}
scale_x_dx <- function(...) scale_x_discrete(labels = DX_AXIS_LABELS, ...)

# One column of panels per region, grouped by study
facet_regions <- function() {
  facet_grid(. ~ study + region, scales = "free_x", space = "free_x")
}


# ---- Factor order ------------------------------------------------------------

# Apply the order above to the usual columns (only those present in the table)
set_pub_levels <- function(df) {
  df |>
    mutate(
      across(any_of("study"), \(x) factor(x, levels = STUDY_ORDER)),
      across(any_of("region"), \(x) factor(x, levels = REGION_ORDER)),
      across(any_of("diagnosis"), \(x) factor(x, levels = DX_ORDER)),
      across(any_of("sex"), \(x) factor(x, levels = SEX_ORDER))
    )
}


# ---- Saving ------------------------------------------------------------------

# Always saves both <file>.pdf and <file>.jpeg.
# `file` may end in .pdf, .jpeg or have no extension.
# `width`: "single", "onehalf", "double" or a number in mm.
save_pub <- function(plot, file, width = "double", height = 120) {
  if (is.character(width)) width <- FIG_WIDTH[[width]]
  base <- sub("\\.(pdf|jpe?g)$", "", file)

  ggsave(
    paste0(base, ".pdf"), plot,
    width = width, height = height, units = "mm", device = cairo_pdf
  )
  ggsave(
    paste0(base, ".jpeg"), plot,
    width = width, height = height, units = "mm", dpi = FIG_DPI,
    device = "jpeg", bg = "white", quality = 95
  )
  invisible(paste0(base, c(".pdf", ".jpeg")))
}

# Saves a whole figure and each of its panels in one folder:
#   <dir>/figure<id>.pdf  and  <dir>/<id><letter>.pdf  (+ .jpeg), e.g.
#   results/figures/figure1/figure1.pdf, 1a.pdf, 1b.pdf ...
# `panels`: named list, names = panel letters (same order as in the figure).
# `panel_width`, `panel_height`: one value for all panels, or a named vector
# with one value per letter.
save_figure <- function(figure, panels, dir, id, width = "double", height = 120,
                        panel_width = "single", panel_height = 70) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  save_pub(figure, file.path(dir, paste0("figure", id)), width = width, height = height)
  size_of <- \(x, letter) if (length(x) > 1) x[[letter]] else x
  for (letter in names(panels)) {
    save_pub(
      panels[[letter]], file.path(dir, paste0(id, letter)),
      width = size_of(panel_width, letter), height = size_of(panel_height, letter)
    )
  }
  message("Saved figure ", id, " and panels ", paste(names(panels), collapse = ", "), " in ", dir)
}
