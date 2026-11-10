# Didactic diagram of block B, analyses 4 and 5 (docs/experimental_design.md),
# and of the overlap test they share with analysis 7, to explain the design to
# people outside the project. Same style as analyses_1_3_explainer_plot.R.
#
# Not part of the Snakemake workflow: run it by hand, from the repository root,
# with the plots environment (workflow/envs/plots_r.yaml):
#   Rscript workflow/scripts/analyses_4_5_explainer_plot.R
#
# Output: docs/figures/analyses_4_5_explainer.pdf and .jpeg
# Labels are in Spanish. The values in the small charts are illustrative, not results.

library(dplyr)
library(purrr)
library(tibble)
source("workflow/scripts/R/theme_pub.R")

out_file <- "docs/figures/analyses_4_5_explainer"
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)

WIDTH <- FIG_WIDTH[["double"]] # 183 mm
HEIGHT <- 246

# Bottom of each panel (mm); each panel is 80 mm tall
PANEL_Y <- c(a = 166, b = 84, c = 2)

COL_CT <- COL_DX[["CT"]]
COL_AD <- COL_DX[["AD"]]
COL_DE <- NPG[["green"]]
COL_SHARED <- INK
FADED <- "grey85"


# ---- Drawing helpers -----------------------------------------------------------
# Every helper adds rows to one of these layers; all coordinates are in mm.

layers <- list(text = list(), rect = list(), seg = list(), point = list())
add <- function(type, df) layers[[type]] <<- c(layers[[type]], list(df))

label <- function(x, y, text, size = 6.5, face = "plain", colour = INK,
                  hjust = 0, vjust = 0.5, angle = 0) {
  add("text", tibble(x, y, text = as.character(text), size = size / .pt, face, colour, hjust, vjust, angle))
}

box <- function(xmin, xmax, ymin, ymax, fill = NA, colour = "grey70", linetype = "solid") {
  add("rect", tibble(xmin, xmax, ymin, ymax, fill, colour, linetype))
}

line <- function(x, y, xend, yend, colour = INK_SECONDARY, head = TRUE,
                 linetype = "solid", linewidth = 0.45) {
  add("seg", tibble(x, y, xend, yend, colour, head, linetype, linewidth))
}

dot <- function(x, y, colour, size = 1.6, shape = 16, fill = "white", stroke = 0.4) {
  add("point", tibble(x, y, colour, size, shape, fill, stroke))
}

title <- function(y0, number, text, question) {
  label(2, y0 + 77, paste0(number, " · ", text), size = 8.5, face = "bold")
  label(2, y0 + 72.5, question, size = 6.5, colour = INK_SECONDARY)
}

step <- function(x, y, text) label(x, y, text, size = 6.5, face = "bold")
note <- function(x, y, text, ...) label(x, y, text, size = 5.5, colour = INK_SECONDARY, vjust = 1, ...)

# Ego network: one gene in the center and its neighbors on a circle.
# `colours`: one color per neighbor
ego <- function(x0, y0, r, colours, start = 90) {
  angle <- (start + seq(0, 360, length.out = length(colours) + 1)[-1]) * pi / 180
  x <- x0 + r * cos(angle)
  y <- y0 + r * sin(angle)
  line(x0, y0, x, y, colour = colours, head = FALSE, linewidth = 0.4)
  dot(x, y, colour = colours, size = 1.5)
  dot(x0, y0, colour = INK, size = 2.6, shape = 21, fill = "white", stroke = 0.8)
}

# Axes of a small chart: x from x0 to x1, y from y0 to y1 (mm)
axes <- function(x0, x1, y0, y1) {
  line(x0, y0, x1, y0, colour = INK, head = FALSE, linewidth = PUB_LINE_WIDTH)
  line(x0, y0, x0, y1, colour = INK, head = FALSE, linewidth = PUB_LINE_WIDTH)
}


# ---- Panel a: analysis 4 ---------------------------------------------------------

y0 <- PANEL_Y[["a"]]
title(y0, "4", "¿Los genes recableados son DE?",
      "Un gen puede cambiar su nivel de expresión (DE), cambiar con quién está conectado (rewiring), o ambas cosas. ¿Son los mismos genes?")

# Step 1: neighbor change of one gene (same example as the design document)
step(3, y0 + 64, "1. Cambio de vecinos de un gen")
ego(13, y0 + 46, 8, c(rep(COL_SHARED, 4), rep(COL_CT, 4)))
ego(36, y0 + 46, 8, c(rep(COL_SHARED, 4), rep(COL_AD, 6)))
label(13, y0 + 35, "en CT: 8 vecinos", size = 6, hjust = 0.5, face = "bold", colour = COL_CT)
label(36, y0 + 35, "en AD: 10 vecinos", size = 6, hjust = 0.5, face = "bold", colour = COL_AD)
dot(4, y0 + 30, colour = COL_SHARED, size = 1.5)
label(5.5, y0 + 30, "vecino en las dos redes (4)", size = 5.5)
label(3, y0 + 23.5, "cambio de vecinos = 1 − 4 / 14 = 0.71", size = 6, face = "bold")
note(3, y0 + 20.5, "14 = vecinos distintos en total\n0 = mismos vecinos · 1 = ninguno en común")
note(3, y0 + 12, "Se calcula en cada réplica AD contra la red CT\ny se promedia. Solo genes con ≥ 5 vecinos.")

line(47, y0 + 46, 52, y0 + 46)

# Step 2: adjusted rewiring = observed - expected for the degree (and expression)
step(54, y0 + 64, "2. Rewiring ajustado")
expected <- \(d) 0.983 - 0.130 * log(d) # 0.75 with 6 neighbors, 0.45 with 60
d_to_x <- \(d) 60 + (log(d) - log(5)) / (log(150) - log(5)) * 36
v_to_y <- \(v) y0 + 18 + v * 38
axes(60, 97, y0 + 18, y0 + 57)
set.seed(4)
cloud <- tibble(d = exp(runif(70, log(5), log(150)))) |>
  mutate(v = pmin(pmax(expected(d) + rnorm(n(), 0, 0.07), 0.02), 0.98))
dot(d_to_x(cloud$d), v_to_y(cloud$v), colour = "grey70", size = 0.7)
curve_d <- exp(seq(log(5), log(150), length.out = 30))
line(d_to_x(head(curve_d, -1)), v_to_y(expected(head(curve_d, -1))), d_to_x(curve_d[-1]), v_to_y(expected(curve_d[-1])),
     colour = INK, head = FALSE, linewidth = 0.5)
# Genes A and B: same observed change (0.70), different degree
walk2(c(6, 60), c("A", "B"), \(d, gene) {
  line(d_to_x(d), v_to_y(expected(d)), d_to_x(d), v_to_y(0.70), colour = COL_AD, head = FALSE, linewidth = 0.5)
  dot(d_to_x(d), v_to_y(0.70), colour = COL_AD, size = 1.8)
  label(d_to_x(d) + 1.2, v_to_y(0.70) + 2.5, gene, size = 6, face = "bold", colour = COL_AD)
})
label(d_to_x(6) + 1.5, v_to_y(0.61), "−0.05", size = 5.5, colour = COL_AD)
label(d_to_x(60) + 1.5, v_to_y(0.58), "+0.25", size = 5.5, colour = COL_AD)
label(78.5, y0 + 14.5, "vecinos (escala log)", size = 5.5, hjust = 0.5)
label(57, v_to_y(0.5), "cambio de vecinos", size = 5.5, hjust = 0.5, angle = 90)
label(d_to_x(9), v_to_y(expected(9)) + 4.5, "esperado para su grado", size = 5.5, face = "bold")
note(54, y0 + 10.5, "A y B cambian igual (0.70), pero B tiene 60 vecinos:\npara su grado se esperaba 0.45 → B se recablea más")

line(100, y0 + 46, 105, y0 + 46)

# Step 3: rewired genes = top 10%
step(106, y0 + 64, "3. Recableados = top 10%")
h_x <- \(z) 108 + (z + 3) / 6 * 30
bins <- seq(-3, 2.75, by = 0.25)
walk(bins, \(z) {
  height <- 32 * exp(-(z + 0.125)^2 / 2)
  box(h_x(z), h_x(z + 0.25), y0 + 26, y0 + 26 + height,
      fill = if (z >= 1.25) COL_AD else "grey80", colour = "white")
})
line(108, y0 + 26, 138, y0 + 26, colour = INK, head = FALSE, linewidth = PUB_LINE_WIDTH)
line(h_x(1.28), y0 + 24, h_x(1.28), y0 + 52, colour = INK_SECONDARY, head = FALSE, linetype = "dashed", linewidth = 0.3)
label(h_x(1.28) + 1, y0 + 50, "10%", size = 5.5, face = "bold", colour = COL_AD)
label(123, y0 + 22.5, "rewiring ajustado", size = 5.5, hjust = 0.5)
box(106, 140, y0 + 3, y0 + 18, fill = "grey96", colour = "grey80")
label(107.5, y0 + 15.5, "Recableado robusto", size = 5.5, face = "bold")
note(107.5, y0 + 12.5, "en el top 10% en ≥ 8 de las 10\nréplicas AD: la lista de los genes\ncon más rewiring")

line(140.5, y0 + 46, 144, y0 + 46)

# Step 4: rewired x DE
step(145, y0 + 64, "4. Recableados × DE")
cell <- function(col, row, text, fill = "white") {
  xmin <- 155 + 13 * (col - 1)
  ymax <- y0 + 56 - 11 * (row - 1)
  box(xmin, xmin + 13, ymax - 11, ymax, fill = fill, colour = "grey60")
  label(xmin + 6.5, ymax - 5.5, text, size = 5, hjust = 0.5, colour = INK)
}
label(161.5, y0 + 58.5, "DE", size = 6, hjust = 0.5, face = "bold", colour = COL_DE)
label(174.5, y0 + 58.5, "no DE", size = 6, hjust = 0.5, face = "bold")
label(154, y0 + 50.5, "recab.", size = 6, hjust = 1, face = "bold", colour = COL_AD)
label(154, y0 + 39.5, "no recab.", size = 6, hjust = 1, face = "bold")
cell(1, 1, "nivel y\nconexiones")
cell(2, 1, "solo la red\nlo detecta", fill = "#FBE0DA")
cell(1, 2, "solo el\nnivel")
cell(2, 2, "sin\ncambio")
note(145, y0 + 29, "¿Hay más DE entre los recableados\nde lo esperado? OR y p con la\nprueba común (panel c)")
note(145, y0 + 16, "Los de \"solo la red\": enriquecimiento\nde funciones (GO)")


# ---- Panel b: analysis 5 ---------------------------------------------------------

y0 <- PANEL_Y[["b"]]
title(y0, "5", "¿Los hubs que cambian son DE?",
      "¿Los genes que ganan o pierden centralidad en la red de AD son genes que cambian de expresión?")

# Step 1: hubs of one network
step(3, y0 + 64, "1. Hubs de cada red")
nodes <- tibble(
  x = c(24, 12, 16, 33, 36, 28, 15, 20, 31, 40, 8),
  y = y0 + c(44, 52, 34, 54, 40, 30, 43, 57, 47, 50, 40)
)
edges <- rbind(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(1, 6), c(1, 7), c(1, 8), c(1, 9), c(5, 10), c(2, 11), c(3, 7), c(4, 9), c(9, 10))
line(nodes$x[edges[, 1]], nodes$y[edges[, 1]], nodes$x[edges[, 2]], nodes$y[edges[, 2]], colour = COL_CT, head = FALSE, linewidth = 0.4)
dot(nodes$x[-1], nodes$y[-1], colour = COL_CT, size = 1.5)
dot(nodes$x[1], nodes$y[1], colour = COL_CT, size = 3.6)
dot(nodes$x[1], nodes$y[1], colour = INK, size = 5.5, shape = 21, fill = NA, stroke = 0.6)
label(nodes$x[1] + 3.5, nodes$y[1] - 3.5, "hub", size = 6, face = "bold")
note(3, y0 + 25, "PageRank: un gen es central si se\nconecta con otros genes centrales")
note(3, y0 + 15.5, "hub = top 1% de los genes de la red\n(~190 genes), sin pesos")

line(47, y0 + 44, 52, y0 + 44)

# Step 2: hub in CT and in the 10 AD replicates -> type
step(54, y0 + 64, "2. ¿Es hub en CT y en las 10 réplicas AD?")
label(66, y0 + 57, "CT", size = 6, hjust = 0.5, face = "bold", colour = COL_CT)
label(89, y0 + 57, "10 réplicas AD", size = 6, hjust = 0.5, face = "bold", colour = COL_AD)
label(113, y0 + 57, "tipo", size = 6, face = "bold")
genes <- tribble(
  ~gene, ~ct, ~ad, ~type, ~colour,
  "gen 1", TRUE, 9, "compartido", INK,
  "gen 2", FALSE, 10, "ganado en AD", COL_AD,
  "gen 3", TRUE, 1, "perdido en AD", COL_CT,
  "gen 4", TRUE, 5, "no se clasifica", INK_SECONDARY
)
pwalk(list(seq_len(nrow(genes)), genes$gene, genes$ct, genes$ad, genes$type, genes$colour), \(i, gene, ct, ad, type, colour) {
  y <- y0 + 51 - 7.5 * (i - 1)
  label(55, y, gene, size = 6)
  box(64.5, 67.5, y - 1.5, y + 1.5, fill = if (ct) COL_CT else "white", colour = COL_CT)
  walk(1:10, \(k) {
    x <- 72 + 3.6 * (k - 1)
    box(x, x + 3, y - 1.5, y + 1.5, fill = if (k <= ad) COL_AD else "white", colour = COL_AD)
  })
  label(113, y, type, size = 6, face = "bold", colour = colour)
})
box(64.5, 67.5, y0 + 18.5, y0 + 21.5, fill = INK_SECONDARY, colour = INK_SECONDARY)
label(69, y0 + 20, "= hub en esa red", size = 5.5, colour = INK_SECONDARY)
note(55, y0 + 15, "hub en AD: en ≥ 8 de las 10 réplicas · no hub en AD: en ≤ 2\nentre 3 y 7 no se sabe si se pierde: no se clasifica")

line(128.5, y0 + 39.75, 133, y0 + 39.75)

# Step 3: fraction of DE genes among each hub type, observed vs expected
step(134, y0 + 64, "3. ¿Son DE más de lo esperado?")
pct_to_y <- \(p) y0 + 22 + p / 50 * 34
axes(142, 180, y0 + 22, y0 + 58)
walk(c(0, 25, 50), \(p) {
  line(141.2, pct_to_y(p), 142, pct_to_y(p), colour = INK, head = FALSE, linewidth = PUB_LINE_WIDTH)
  label(140.6, pct_to_y(p), paste0(p, "%"), size = 5, hjust = 1, colour = INK_SECONDARY)
})
label(135, pct_to_y(25), "hubs que son DE", size = 5.5, hjust = 0.5, angle = 90)
types <- tibble(
  type = c("compartido", "ganado", "perdido"), x = c(150, 162, 174),
  observed = c(20, 42, 33), expected = c(18, 20, 19), colour = c(INK, COL_AD, COL_CT)
)
pwalk(types, \(type, x, observed, expected, colour) {
  line(x, pct_to_y(expected), x, pct_to_y(observed), colour = COL_NEUTRAL, head = FALSE, linewidth = 0.4)
  dot(x, pct_to_y(expected), colour = INK, size = 1.6, shape = 21, fill = "white")
  dot(x, pct_to_y(observed), colour = colour, size = 2)
  label(x, y0 + 19, type, size = 5.5, hjust = 0.5)
})
dot(143.5, y0 + 13, colour = INK, size = 1.8)
label(145, y0 + 13, "observado", size = 5.5)
dot(160, y0 + 13, colour = INK, size = 1.6, shape = 21, fill = "white")
label(161.5, y0 + 13, "esperado por azar", size = 5.5)
note(134, y0 + 9, "esperado: genes con la misma expresión (panel c)")
note(134, y0 + 5.5, "valores ilustrativos")


# ---- Panel c: the shared overlap test -------------------------------------------

y0 <- PANEL_Y[["c"]]
label(2, y0 + 77, "La prueba común de los análisis 4, 5 y 7", size = 8.5, face = "bold")
label(2, y0 + 72.5, "¿Un conjunto de genes (por ejemplo, los recableados) contiene más genes DE de lo esperado por azar?",
      size = 6.5, colour = INK_SECONDARY)

# Genes in 3 expression strata; DE genes in green, genes of the set with a black ring
strata <- c("alta", "media", "baja")
in_set <- list(alta = c(1, 2, 5), media = c(2, 4), baja = c(1))
de_observed <- list(alta = c(1, 2, 3, 6), media = c(2, 7), baja = c(5))
de_shuffled <- list(alta = c(2, 4, 7, 8), media = c(1, 5), baja = c(3))

genes_by_stratum <- function(x0, de) {
  iwalk(strata, \(s, i) {
    y <- y0 + 54 - 12 * (i - 1)
    box(x0, x0 + 44, y - 4.5, y + 4.5, fill = "grey96", colour = "grey80")
    label(x0 + 1.5, y + 2.6, paste("expresión", s), size = 5, colour = INK_SECONDARY)
    walk(1:8, \(k) {
      x <- x0 + 4.5 + 5 * (k - 1)
      ring <- k %in% in_set[[s]]
      dot(x, y - 1, colour = if (ring) INK else "grey60", size = 2.6, shape = 21,
          fill = if (k %in% de[[s]]) COL_DE else "white", stroke = if (ring) 1.1 else 0.4)
    })
  })
}
count_both <- \(de) sum(map_int(strata, \(s) length(intersect(in_set[[s]], de[[s]]))))

step(3, y0 + 64, "1. Lo observado")
genes_by_stratum(3, de_observed)
label(3, y0 + 22, paste0("genes del conjunto que son DE: ", count_both(de_observed)), size = 6, face = "bold")
dot(4.5, y0 + 15.5, colour = "grey60", size = 2.4, shape = 21, fill = COL_DE)
label(6.5, y0 + 15.5, "DE", size = 5.5)
dot(15, y0 + 15.5, colour = INK, size = 2.4, shape = 21, fill = "white", stroke = 1.1)
label(17, y0 + 15.5, "en el conjunto (ej. recableado)", size = 5.5)
note(3, y0 + 11, "Los genes se agrupan por expresión media (deciles;\ntambién por grado cuando el conjunto son recableados)")

line(48.5, y0 + 42, 55.5, y0 + 42)

step(57, y0 + 64, "2. Se barajan las etiquetas DE")
genes_by_stratum(57, de_shuffled)
label(57, y0 + 22, paste0("en esta permutación: ", count_both(de_shuffled)), size = 6, face = "bold")
label(57, y0 + 16, "× 1000 permutaciones", size = 6, face = "bold", colour = COL_DE)
note(57, y0 + 11, "Cada estrato conserva sus genes DE: los DE al azar\ntienen la misma expresión que los observados")

line(102.5, y0 + 42, 109.5, y0 + 42)

# Null distribution of the overlap and the observed value (illustrative)
step(111, y0 + 64, "3. ¿Qué tan raro es lo observado?")
null_counts <- c(`0` = 140, `1` = 330, `2` = 300, `3` = 160, `4` = 55, `5` = 15)
k_to_x <- \(k) 118 + k * 6
n_to_y <- \(n) y0 + 24 + n / 330 * 30
axes(116, 155, y0 + 24, y0 + 56)
iwalk(null_counts, \(n, k) {
  k <- as.numeric(k)
  box(k_to_x(k) - 2.6, k_to_x(k) + 2.6, n_to_y(0), n_to_y(n),
      fill = if (k >= 3) "#B3E3D9" else "grey80", colour = "white")
  label(k_to_x(k), y0 + 21.5, k, size = 5, hjust = 0.5, colour = INK_SECONDARY)
})
line(k_to_x(3) - 3, y0 + 23, k_to_x(3) - 3, y0 + 58, colour = INK, head = FALSE, linetype = "dashed", linewidth = 0.4)
label(k_to_x(3) - 2, y0 + 57, "observado", size = 5.5, face = "bold")
label(135.5, y0 + 18, "genes del conjunto que son DE", size = 5.5, hjust = 0.5)
label(112.5, n_to_y(165), "permutaciones", size = 5.5, hjust = 0.5, angle = 90)
note(111, y0 + 13.5, "p = fracción de permutaciones tan extremas\ncomo lo observado · esperado = promedio")
note(111, y0 + 6, "valores ilustrativos")

box(158, 181, y0 + 4, y0 + 66, fill = "grey96", colour = "grey80")
label(160, y0 + 62, "Resultado", size = 6.5, face = "bold")
note(160, y0 + 58.5, "OR de la tabla 2×2\n(conjunto × DE)\n\np por permutación\n\nBH entre regiones → q")
label(160, y0 + 34, "¿Por qué no Fisher?", size = 6, face = "bold")
note(160, y0 + 30.5, "Los genes conectados\nen la red no son\nindependientes, y\nFisher no controla\nla expresión")


# ---- Panel separators -------------------------------------------------------

walk(PANEL_Y[c("a", "b")] - 1, \(y) line(2, y, 181, y, colour = "grey80", head = FALSE, linewidth = 0.3))


# ---- Plot ----------------------------------------------------------------------

L <- map(layers, list_rbind)

p <- ggplot() +
  geom_rect(
    data = L$rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, colour = colour, linetype = linetype),
    linewidth = 0.3
  ) +
  geom_segment(
    data = filter(L$seg, !head),
    aes(x = x, y = y, xend = xend, yend = yend, colour = colour, linetype = linetype, linewidth = linewidth)
  ) +
  geom_segment(
    data = filter(L$seg, head),
    aes(x = x, y = y, xend = xend, yend = yend, colour = colour, linetype = linetype, linewidth = linewidth),
    arrow = arrow(length = unit(1.2, "mm"), type = "closed")
  ) +
  geom_point(data = L$point, aes(x = x, y = y, colour = colour, size = size, shape = shape, fill = fill, stroke = stroke)) +
  geom_text(
    data = L$text,
    aes(x = x, y = y, label = text, size = size, fontface = face, colour = colour, hjust = hjust, vjust = vjust, angle = angle),
    family = PUB_FONT, lineheight = 0.95
  ) +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_size_identity() +
  scale_shape_identity() +
  scale_linetype_identity() +
  scale_linewidth_identity() +
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
