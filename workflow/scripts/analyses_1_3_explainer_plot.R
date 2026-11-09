# Didactic diagram of analyses 1, 2 and 3 (docs/experimental_design.md),
# to explain the design to people outside the project.
#
# Not part of the Snakemake workflow: run it by hand, from the repository root,
# with the plots environment (workflow/envs/plots_r.yaml):
#   Rscript workflow/scripts/analyses_1_3_explainer_plot.R
#
# Output: docs/figures/analyses_1_3_explainer.pdf and .jpeg
# Labels are in Spanish. The values in the small charts are illustrative, not results.

library(dplyr)
library(purrr)
library(tibble)
source("workflow/scripts/R/theme_pub.R")

out_file <- "docs/figures/analyses_1_3_explainer"
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)

WIDTH <- FIG_WIDTH[["double"]] # 183 mm
HEIGHT <- 246

# Bottom of each panel (mm); each panel is 80 mm tall
PANEL_Y <- c(a = 166, b = 84, c = 2)

COL_CT <- COL_DX[["CT"]]
COL_AD <- COL_DX[["AD"]]
COL_HALF <- c(A = NPG[["cyan"]], B = NPG[["tan"]])
FADED <- "grey85"


# ---- Drawing helpers -----------------------------------------------------------
# Every helper adds rows to one of these layers; all coordinates are in mm.

layers <- list(text = list(), rect = list(), seg = list(), poly = list(), point = list())
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

dot <- function(x, y, colour, size = 1.6, shape = 16, fill = "white") {
  add("point", tibble(x, y, colour, size, shape, fill))
}

# One person: round head and half-ellipse body. (x, y) = bottom center; h = height
person <- function(id, x, y, h, fill) {
  t <- seq(0, 2 * pi, length.out = 30)
  b <- seq(0, pi, length.out = 20)
  bind_rows(
    tibble(x = x + 0.2 * h * cos(t), y = y + 0.72 * h + 0.2 * h * sin(t), group = paste(id, "head")),
    tibble(x = x + 0.34 * h * cos(b), y = y + 0.48 * h * sin(b), group = paste(id, "body"))
  ) |>
    mutate(fill = fill, alpha = 1)
}

# A grid of people. `fills` = one color for all, or one color per person
crowd <- function(id, x0, y0, ncol, nrow, fills, step = 4.2, h = 5) {
  k <- ncol * nrow
  fills <- rep_len(fills, k)
  people <- tibble(i = seq_len(k)) |>
    mutate(x = x0 + ((i - 1) %% ncol) * step, y = y0 + ((i - 1) %/% ncol) * (h + 1.6))
  add("poly", pmap(list(people$i, people$x, people$y, fills), \(i, x, y, f) person(paste(id, i), x, y, h, f)) |>
    list_rbind())
}

# A small network: 7 genes in a ring + 1 in the center; `edges` = pairs of gene numbers
network <- function(x0, y0, r, edges, colour) {
  angle <- seq(90, 90 + 360, length.out = 8)[1:7] * pi / 180
  radius <- r * c(1, 0.8, 1, 0.85, 1, 0.8, 0.95)
  nodes <- tibble(x = c(x0 + radius * cos(angle), x0), y = c(y0 + radius * sin(angle), y0))
  line(
    nodes$x[edges[, 1]], nodes$y[edges[, 1]], nodes$x[edges[, 2]], nodes$y[edges[, 2]],
    colour = colour, head = FALSE, linewidth = 0.4
  )
  dot(nodes$x, nodes$y, colour = colour, size = 1.4)
}

circle <- function(id, x0, y0, r, fill, alpha = 0.35) {
  t <- seq(0, 2 * pi, length.out = 60)
  add("poly", tibble(x = x0 + r * cos(t), y = y0 + r * sin(t), group = id, fill = fill, alpha = alpha))
}

# Two overlapping circles: the shared part are the shared edges
venn <- function(id, x0, y0, r, fills) {
  circle(paste(id, 1), x0 - 0.55 * r, y0, r, fills[[1]])
  circle(paste(id, 2), x0 + 0.55 * r, y0, r, fills[[2]])
}

title <- function(y0, number, text, question) {
  label(2, y0 + 77, paste0(number, " · ", text), size = 8.5, face = "bold")
  label(2, y0 + 72.5, question, size = 6.5, colour = INK_SECONDARY)
}

# Edges of the example networks: CT and AD share 5 of 13 edges;
# halves A and B share 8 of 10 (same group, different people)
EDGES_CT <- rbind(c(8, 1), c(8, 3), c(8, 5), c(1, 2), c(3, 4), c(5, 6), c(6, 7), c(2, 3), c(7, 1))
EDGES_AD <- rbind(c(8, 1), c(8, 3), c(1, 2), c(3, 4), c(4, 5), c(8, 6), c(2, 7), c(5, 7), c(6, 7))
EDGES_A <- EDGES_CT
EDGES_B <- rbind(c(8, 1), c(8, 3), c(8, 5), c(1, 2), c(3, 4), c(5, 6), c(2, 3), c(7, 1), c(4, 5))


# ---- Panel a: analysis 1 -------------------------------------------------------

y0 <- PANEL_Y[["a"]]
title(y0, "1", "AD vs CT con el mismo n  (análisis principal)",
      "¿Cuánto cambia la red de co-expresión entre AD y CT en cada región?")

# People
label(4, y0 + 66, "CT: se usan todas las personas", face = "bold", colour = COL_CT)
label(4, y0 + 62.5, "ej. DLPFC: las 141", colour = INK_SECONDARY)
crowd("ct1", 6, y0 + 46, 8, 2, COL_CT)

label(4, y0 + 37.5, "AD: se sortean tantas personas como en CT", face = "bold", colour = COL_AD)
label(4, y0 + 34, "ej. DLPFC: 141 de 289, con la misma proporción de mujeres que CT", colour = INK_SECONDARY)
drawn <- c(1, 3, 4, 6, 9, 10, 13, 15, 18, 19, 21, 24)
crowd("ad1", 6, y0 + 9, 8, 3, if_else(seq_len(24) %in% drawn, COL_AD, FADED))
label(4, y0 + 3.5, "En gris: personas que no salieron en ese sorteo", size = 5.5, colour = INK_SECONDARY)

# ARACNe-AP
line(42, y0 + 51, 58, y0 + 51)
label(50, y0 + 53.5, "ARACNe-AP", size = 5.5, hjust = 0.5, colour = INK_SECONDARY)
network(67, y0 + 51, 6, EDGES_CT, COL_CT)
label(67, y0 + 42.5, "1 red CT", hjust = 0.5, face = "bold", colour = COL_CT)

line(42, y0 + 17, 58, y0 + 17)
label(50, y0 + 19.5, "ARACNe-AP", size = 5.5, hjust = 0.5, colour = INK_SECONDARY)
label(50, y0 + 14, "× 10 sorteos", size = 6, hjust = 0.5, face = "bold", colour = COL_AD)
walk(2:0, \(k) box(59.5 + 1.3 * k, 74.5 + 1.3 * k, y0 + 9.5 + 1.3 * k, y0 + 24.5 + 1.3 * k,
                   fill = "white", colour = COL_AD))
network(67, y0 + 17, 5.5, EDGES_AD, COL_AD)
label(67, y0 + 5.5, "10 redes AD", hjust = 0.5, face = "bold", colour = COL_AD)

# Comparison
line(76, y0 + 49, 92, y0 + 40)
line(77, y0 + 24, 94, y0 + 31)
label(104, y0 + 48, "Se compara cada red AD con la red CT", face = "bold", hjust = 0.5)
venn("v1", 104, y0 + 36, 7, c(COL_CT, COL_AD))
label(97.5, y0 + 36, "CT", hjust = 0.5, face = "bold", colour = COL_CT)
label(110.5, y0 + 36, "AD", hjust = 0.5, face = "bold", colour = COL_AD)
label(104, y0 + 25.5, "J = enlaces compartidos / enlaces en total", size = 6, hjust = 0.5)
label(104, y0 + 22, "J = 1: redes idénticas · J = 0: nada en común", size = 5.5, hjust = 0.5, colour = INK_SECONDARY)

box(84, 128, y0 + 2, y0 + 17, fill = "grey96", colour = "grey80")
label(86, y0 + 14.5, "¿Por qué el mismo n?", size = 6, face = "bold")
label(86, y0 + 11.5, "Con más personas, ARACNe detecta más enlaces.", size = 5.5)
label(86, y0 + 8.5, "Si AD tuviera 289 personas y CT 141, las redes", size = 5.5)
label(86, y0 + 5.5, "serían distintas solo por eso.", size = 5.5)

# Result: 10 J values per region (illustrative)
line(119, y0 + 36, 126, y0 + 36)
label(128, y0 + 62, "Resultado: 10 valores de J por región", face = "bold")
j_to_y <- \(j) y0 + 18 + j / 0.3 * 38
line(137, j_to_y(0), 137, j_to_y(0.3), colour = INK, head = FALSE, linewidth = PUB_LINE_WIDTH)
walk(c(0, 0.1, 0.2, 0.3), \(j) {
  line(136.2, j_to_y(j), 137, j_to_y(j), colour = INK, head = FALSE, linewidth = PUB_LINE_WIDTH)
  label(135.6, j_to_y(j), format(j, nsmall = 1), size = 5, hjust = 1, colour = INK_SECONDARY)
})
label(129.5, j_to_y(0.15), "J (AD vs CT)", size = 5.5, hjust = 0.5, angle = 90)
example_j <- c(DLPFC = 0.25, ACC = 0.21, PCC = 0.23, TCX = 0.13, CBE = 0.17)
jitter_x <- c(-1.2, -0.6, 0, 0.6, 1.2, -0.9, -0.3, 0.3, 0.9, 0)
jitter_j <- c(-0.008, 0.004, 0.009, -0.003, 0.006, 0.001, -0.006, 0.007, -0.002, -0.009)
iwalk(example_j, \(j, region) {
  x <- 143 + 8 * (match(region, names(example_j)) - 1)
  dot(x + jitter_x, j_to_y(j + jitter_j), colour = COL_REGION[[region]], size = 1)
  label(x, j_to_y(0) - 3, region, size = 5.5, hjust = 0.5)
})
label(158, y0 + 9, "valores ilustrativos", size = 5.5, hjust = 0.5, colour = INK_SECONDARY)


# ---- Panel b: analysis 2 -------------------------------------------------------

y0 <- PANEL_Y[["b"]]
title(y0, "2", "Mismo n en todas las regiones  (robustez)",
      "En el análisis 1 cada región tiene un n distinto. ¿El orden de las regiones se mantiene si todas tienen el mismo n?")

n_to_w <- \(n) n * 0.2 # bar length (mm) per person
bars <- function(top, title_text, n) {
  label(4, top, title_text, face = "bold")
  iwalk(n, \(value, region) {
    y <- top - 6 - 5 * (match(region, names(n)) - 1)
    label(14, y, region, size = 6, hjust = 1)
    box(16, 16 + n_to_w(value), y - 1.7, y + 1.7, fill = COL_REGION[[region]], colour = NA)
    label(17 + n_to_w(value), y, value, size = 5.5, colour = INK_SECONDARY)
  })
}
bars(y0 + 63, "Análisis 1: n distinto por región", c(DLPFC = 141, ACC = 94, PCC = 99))
line(16 + n_to_w(90), y0 + 44.5, 16 + n_to_w(90), y0 + 60, colour = INK_SECONDARY, head = FALSE, linetype = "dashed", linewidth = 0.3)
label(16 + n_to_w(90), y0 + 42.5, "n = 90", size = 5.5, hjust = 0.5, colour = INK_SECONDARY)
bars(y0 + 34, "Análisis 2: n = 90 en todas", c(DLPFC = 90, ACC = 90, PCC = 90))
label(4, y0 + 9, "90 es apenas menos que el grupo\nmás chico (ACC CT: 94 personas)", size = 5.5, colour = INK_SECONDARY, vjust = 0.5)

line(56, y0 + 40, 64, y0 + 40)

# In each region and each draw: 90 CT and 90 AD
label(68, y0 + 63, "En cada región y en cada sorteo:", face = "bold")
label(68, y0 + 57.5, "se sortean 90 CT", colour = COL_CT, face = "bold")
crowd("ct2", 71, y0 + 47, 5, 1, COL_CT)
label(68, y0 + 35.5, "se sortean 90 AD", colour = COL_AD, face = "bold")
crowd("ad2", 71, y0 + 25, 5, 1, COL_AD)
line(90, y0 + 49.5, 94, y0 + 49.5)
line(90, y0 + 27.5, 94, y0 + 27.5)
network(99, y0 + 50, 4.5, EDGES_CT, COL_CT)
network(99, y0 + 28, 4.5, EDGES_AD, COL_AD)
line(104, y0 + 48, 109, y0 + 42)
line(104, y0 + 30, 109, y0 + 36)
label(110, y0 + 39, "J(AD, CT)", face = "bold")
label(68, y0 + 14, "× 10 sorteos, en cada una de las 3 regiones de ROSMAP\n(DLPFC, ACC y PCC): 60 redes", size = 5.5, colour = INK_SECONDARY)

line(123, y0 + 39, 129, y0 + 39)

# Ranking: same order in both analyses? (illustrative)
label(131, y0 + 63, "¿Se mantiene el orden?", face = "bold")
label(140, y0 + 57, "Análisis 1", size = 6, hjust = 0.5, face = "bold")
label(171, y0 + 57, "Análisis 2", size = 6, hjust = 0.5, face = "bold")
order_1 <- c("PCC", "ACC", "DLPFC")
order_2 <- c("PCC", "ACC", "DLPFC")
rank_y <- \(r) y0 + 51 - 8 * (r - 1)
walk(1:3, \(r) {
  label(131, rank_y(r), paste0(r, "."), size = 6)
  box(134, 146, rank_y(r) - 2.5, rank_y(r) + 2.5, fill = COL_REGION[[order_1[r]]], colour = NA)
  label(140, rank_y(r), order_1[r], size = 6, hjust = 0.5, face = "bold", colour = "white")
  box(165, 177, rank_y(r) - 2.5, rank_y(r) + 2.5, fill = COL_REGION[[order_2[r]]], colour = NA)
  label(171, rank_y(r), order_2[r], size = 6, hjust = 0.5, face = "bold", colour = "white")
  line(146.5, rank_y(r), 164.5, rank_y(match(order_1[r], order_2)), colour = COL_REGION[[order_1[r]]], head = FALSE)
})
label(155.5, y0 + 30, "mismo orden", size = 6, hjust = 0.5, face = "bold", colour = NPG[["green"]])
label(131, y0 + 22, "1 = la región con más rewiring", size = 5.5, colour = INK_SECONDARY)
label(131, y0 + 15, "Resultado: % de sorteos con el mismo orden\ny τ de Kendall contra el orden del análisis 1", size = 6)
label(131, y0 + 7, "orden ilustrativo", size = 5.5, colour = INK_SECONDARY)


# ---- Panel c: analysis 3 -------------------------------------------------------

y0 <- PANEL_Y[["c"]]
title(y0, "3", "J de referencia  (¿cuánto se parecen dos redes sin ninguna diferencia biológica?)",
      "Dos redes del mismo grupo nunca salen idénticas, porque se hacen con personas distintas. ¿Cuánto se parecen solo por eso?")

# One group split into two halves with nobody in common
label(4, y0 + 65, "Un solo grupo", face = "bold", colour = COL_AD)
label(4, y0 + 61.5, "ej. DLPFC AD (289 personas)", size = 6, colour = INK_SECONDARY)
crowd("grp", 6, y0 + 27, 4, 4, rep(c(COL_HALF[["B"]], COL_HALF[["A"]]), each = 8))
line(23, y0 + 45, 33, y0 + 52)
line(23, y0 + 31, 33, y0 + 24)

label(43, y0 + 64.5, "Mitad A (n = 100)", size = 6, hjust = 0.5, face = "bold", colour = COL_HALF[["A"]])
crowd("hA", 37, y0 + 48, 4, 2, COL_HALF[["A"]])
label(43, y0 + 37, "nadie está\nen las dos", size = 5.5, hjust = 0.5, colour = INK_SECONDARY)
crowd("hB", 37, y0 + 15, 4, 2, COL_HALF[["B"]])
label(43, y0 + 10.5, "Mitad B (n = 100)", size = 6, hjust = 0.5, face = "bold", colour = COL_HALF[["B"]])

line(53, y0 + 53, 59, y0 + 53)
line(53, y0 + 20, 59, y0 + 20)
network(65, y0 + 53, 5, EDGES_A, COL_HALF[["A"]])
network(65, y0 + 20, 5, EDGES_B, COL_HALF[["B"]])
label(65, y0 + 45, "red A", size = 6, hjust = 0.5, face = "bold", colour = COL_HALF[["A"]])
label(65, y0 + 12, "red B", size = 6, hjust = 0.5, face = "bold", colour = COL_HALF[["B"]])

line(71, y0 + 50, 80, y0 + 41)
line(71, y0 + 23, 80, y0 + 31)
venn("v3", 87, y0 + 36, 6, COL_HALF)
label(87, y0 + 26, "J(A, B)", hjust = 0.5, face = "bold")
label(87, y0 + 22.5, "× 5 sorteos", size = 5.5, hjust = 0.5, colour = INK_SECONDARY)

line(95, y0 + 36, 101, y0 + 36)

# Reference curve J_ref(n) (illustrative)
n_to_x <- \(n) 110 + (n - 50) * 0.4
jref_to_y <- \(j) y0 + 15 + j / 0.4 * 37
label(104, y0 + 65, "Curva de J de referencia", face = "bold")
line(n_to_x(50), jref_to_y(0), n_to_x(150), jref_to_y(0), colour = INK, head = FALSE, linewidth = PUB_LINE_WIDTH)
line(n_to_x(50), jref_to_y(0), n_to_x(50), jref_to_y(0.4), colour = INK, head = FALSE, linewidth = PUB_LINE_WIDTH)
walk(c(60, 100, 140), \(n) {
  line(n_to_x(n), jref_to_y(0), n_to_x(n), jref_to_y(0) - 0.8, colour = INK, head = FALSE, linewidth = PUB_LINE_WIDTH)
  label(n_to_x(n), jref_to_y(0) - 2.8, n, size = 5, hjust = 0.5, colour = INK_SECONDARY)
})
walk(c(0, 0.2, 0.4), \(j) {
  line(n_to_x(50) - 0.8, jref_to_y(j), n_to_x(50), jref_to_y(j), colour = INK, head = FALSE, linewidth = PUB_LINE_WIDTH)
  label(n_to_x(50) - 1.4, jref_to_y(j), format(j, nsmall = 1), size = 5, hjust = 1, colour = INK_SECONDARY)
})
label(n_to_x(100), jref_to_y(0) - 7, "n (personas por red)", size = 5.5, hjust = 0.5)
label(103.5, jref_to_y(0.2), "J_ref", size = 5.5, hjust = 0.5, angle = 90)

calibration <- tibble(n = c(60, 100, 140), j = c(0.19, 0.28, 0.36))
fit <- lm(j ~ poly(n, 2), data = calibration)
curve_n <- seq(55, 145, by = 2)
curve_j <- predict(fit, tibble(n = curve_n))
line(n_to_x(head(curve_n, -1)), jref_to_y(head(curve_j, -1)), n_to_x(curve_n[-1]), jref_to_y(curve_j[-1]),
     colour = COL_AD, head = FALSE, linewidth = 0.5)
dot(n_to_x(calibration$n), jref_to_y(calibration$j), colour = COL_AD, size = 1.8)
validation <- tibble(n = c(70, 85, 75), j = c(0.217, 0.256, 0.228))
dot(n_to_x(validation$n), jref_to_y(validation$j), colour = INK, size = 1.6, shape = 21, fill = "white")

dot(111, y0 + 60.5, colour = COL_AD, size = 1.6)
label(113, y0 + 60.5, "calibración: DLPFC AD, n = 60, 100 y 140", size = 5.5)
dot(111, y0 + 57, colour = INK, size = 1.5, shape = 21, fill = "white")
label(113, y0 + 57, "validación: DLPFC CT, ACC AD y PCC AD", size = 5.5)
label(130, y0 + 3, "valores ilustrativos", size = 5.5, hjust = 0.5, colour = INK_SECONDARY)

# Rewiring index
box(153, 181, y0 + 6, y0 + 68, fill = "grey96", colour = "grey80")
label(155, y0 + 64, "Índice de rewiring", size = 6.5, face = "bold")
label(155, y0 + 59.5, "J(AD, CT) / J_ref(n)", size = 6.5)
label(155, y0 + 52, "≈ 1: sin rewiring", size = 6, face = "bold")
label(155, y0 + 48.5, "AD y CT se parecen tanto\ncomo dos redes del\nmismo grupo", size = 5.5, vjust = 1, colour = INK_SECONDARY)
label(155, y0 + 36, "≪ 1: hay rewiring", size = 6, face = "bold", colour = COL_AD)
label(155, y0 + 32.5, "AD y CT difieren más de\nlo que explica tener\npersonas distintas", size = 5.5, vjust = 1, colour = INK_SECONDARY)
label(155, y0 + 19.5, "Ejemplo hipotético", size = 5.5, face = "bold")
label(155, y0 + 16, "PCC, n = 99:\n0.087 / 0.30 = 0.29\n→ hay rewiring", size = 5.5, vjust = 1, colour = INK_SECONDARY)


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
  geom_polygon(data = L$poly, aes(x = x, y = y, group = group, fill = fill, alpha = alpha), colour = NA) +
  geom_segment(
    data = filter(L$seg, !head),
    aes(x = x, y = y, xend = xend, yend = yend, colour = colour, linetype = linetype, linewidth = linewidth)
  ) +
  geom_segment(
    data = filter(L$seg, head),
    aes(x = x, y = y, xend = xend, yend = yend, colour = colour, linetype = linetype, linewidth = linewidth),
    arrow = arrow(length = unit(1.2, "mm"), type = "closed")
  ) +
  geom_point(data = L$point, aes(x = x, y = y, colour = colour, size = size, shape = shape, fill = fill), stroke = 0.4) +
  geom_text(
    data = L$text,
    aes(x = x, y = y, label = text, size = size, fontface = face, colour = colour, hjust = hjust, vjust = vjust, angle = angle),
    family = PUB_FONT, lineheight = 0.95
  ) +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_alpha_identity() +
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
