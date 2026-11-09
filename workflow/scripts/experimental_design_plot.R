# Diagram of the analyses in docs/experimental_design.md.
#
# Not part of the Snakemake workflow: run it by hand, from the repository root,
# with the plots environment (workflow/envs/plots_r.yaml):
#   Rscript workflow/scripts/experimental_design_plot.R
#
# Output: docs/figures/experimental_design.pdf and .jpeg

library(tibble)
source("workflow/scripts/R/theme_pub.R")

out_file <- "docs/figures/experimental_design"
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)


# ---- Canvas (mm; the plot is drawn at its final size) -----------------------

WIDTH <- FIG_WIDTH[["double"]] # 183 mm
HEIGHT <- 104

# Colors of each block (pipeline steps in grey)
COL_BLOCK <- c(
  pipeline = "grey55",
  A = NPG[["cyan"]],
  B = NPG[["green"]],
  C = NPG[["salmon"]]
)


# ---- Boxes -------------------------------------------------------------------

# One row per box: position (mm), block, bold title and body text
boxes <- tribble(
  ~id, ~xmin, ~xmax, ~ymin, ~ymax, ~block, ~title, ~body,
  # Data and networks
  "data", 10, 63, 88, 102, "pipeline", "Data",
  "AMP-AD residualized counts\nROSMAP: DLPFC, ACC, PCC\nMayo: TCX, CBE",
  "donor", 69, 122, 88, 102, "pipeline", "One sample per donor",
  "Highest RIN\nAD and CT group per region",
  "aracne", 128, 181, 88, 102, "pipeline", "Networks",
  "Sex-stratified subsampling\nARACNe-AP: 100 bootstraps, p < 0.05\n175 networks",
  # Block A: rewiring between AD and CT
  "a1", 10, 63, 68, 82, "A", "1 · Matched n (main)",
  "CT complete vs AD subsampled\nto the n of CT\n5 regions × 10 replicates: 55 networks",
  "a2", 69, 122, 68, 82, "A", "2 · Equal n (robustness)",
  "AD and CT subsampled to n = 90\n3 ROSMAP regions × 10 replicates\n60 networks",
  "a3", 128, 181, 68, 82, "A", "3 · Reference J",
  "Two disjoint halves of one group\nJ_ref(n), n = 60 to 140, 5 replicates\n60 networks",
  "index", 10, 181, 54, 63, "A", "Rewiring index = J(AD, CT) / J_ref(n)",
  "Metrics: edge Jaccard (J), Spearman ρ of degree, NMI of Leiden communities, mean degree",
  # Block B: networks and differential expression
  "b4", 10, 93, 27, 49, "B", "4 · Rewiring vs DE",
  "Per-gene rewiring, adjusted for degree and expression\nEgo networks of the top 25 up and 25 down DE genes\nDE: Sage, adj. P < 0.05 and |FC| > 1.2\nInput: networks of 1",
  "b5", 98, 181, 27, 49, "B", "5 · Hubs vs DE",
  "Hub: top 1% PageRank, in ≥ 8 of 10 AD replicates\nShared, gained and lost hubs\nRegion-exclusive hubs vs binomial and permutation nulls\nInput: networks of 1",
  # Block C: biological validation
  "c6", 10, 93, 2, 22, "C", "6 · Regional vulnerability (pre-registered)",
  "Prediction: TCX more rewired than CBE\nRewiring index and direct J\nWilcoxon p < 0.05 over 10 replicates\nInput: 1 and 3",
  "c7", 98, 181, 2, 22, "C", "7 · AD risk genes",
  "Prioritized genes of Bellenguez et al. (2022)\nEnrichment among rewired genes\nand among gained and lost hubs\nInput: 4 and 5"
)

# Block labels: a colored strip at the left of each block
strips <- tribble(
  ~xmin, ~xmax, ~ymin, ~ymax, ~block, ~label,
  2, 7, 54, 82, "A", "A · Rewiring AD vs CT",
  2, 7, 27, 49, "B", "B · Networks and DE",
  2, 7, 2, 22, "C", "C · Validation"
)


# ---- Arrows ------------------------------------------------------------------

arrow_head <- arrow(length = unit(1.2, "mm"), type = "closed")

# Arrows between boxes (x, y = start; xend, yend = end)
arrows <- tribble(
  ~x, ~y, ~xend, ~yend,
  63, 95, 69, 95, # data -> one sample per donor
  122, 95, 128, 95, # one sample per donor -> networks
  36.5, 85, 36.5, 82, # networks -> designs 1, 2 and 3
  95.5, 85, 95.5, 82,
  154.5, 85, 154.5, 82,
  36.5, 68, 36.5, 63, # designs -> rewiring index
  95.5, 68, 95.5, 63,
  154.5, 68, 154.5, 63
)

# Line without head: from the networks box down to the bus over the designs
bus <- tribble(
  ~x, ~y, ~xend, ~yend,
  154.5, 88, 154.5, 85,
  36.5, 85, 154.5, 85
)


# ---- Plot --------------------------------------------------------------------

text_size <- PUB_LABEL_SIZE

p <- ggplot() +
  geom_rect(
    data = boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = block),
    fill = "white", linewidth = PUB_LINE_WIDTH * 1.5
  ) +
  geom_rect(
    data = strips,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = block)
  ) +
  geom_text(
    data = strips,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    angle = 90, size = text_size, fontface = "bold", colour = "white",
    family = PUB_FONT
  ) +
  geom_text(
    data = boxes,
    aes(x = xmin + 2, y = ymax - 2, label = title),
    hjust = 0, vjust = 1, size = text_size + 0.3, fontface = "bold",
    family = PUB_FONT, colour = INK
  ) +
  geom_text(
    data = boxes,
    aes(x = xmin + 2, y = ymax - 6, label = body),
    hjust = 0, vjust = 1, size = text_size, lineheight = 1,
    family = PUB_FONT, colour = INK_SECONDARY
  ) +
  geom_segment(
    data = bus, aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = PUB_LINE_WIDTH * 1.5, colour = INK_SECONDARY
  ) +
  geom_segment(
    data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = PUB_LINE_WIDTH * 1.5, colour = INK_SECONDARY,
    arrow = arrow_head
  ) +
  scale_colour_manual(values = COL_BLOCK, guide = "none") +
  scale_fill_manual(values = COL_BLOCK, guide = "none") +
  scale_x_continuous(limits = c(0, WIDTH), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, HEIGHT), expand = c(0, 0)) +
  coord_fixed() +
  theme_pub() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

save_pub(p, out_file, width = WIDTH, height = HEIGHT)
message("Saved ", out_file, ".pdf and .jpeg")
