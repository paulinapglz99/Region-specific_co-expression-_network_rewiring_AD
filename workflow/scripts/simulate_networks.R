# Synthetic networks with a known answer, to test the network analyses
# (docs/experimental_design.md) before running ARACNe-AP.
#
# Not part of the Snakemake workflow: run it by hand, from the repository root,
# with the environment workflow/envs/network_analysis.yaml:
#   Rscript workflow/scripts/simulate_networks.R
#
# Needs the subsample tables (results/subsamples/) and Sage's DE tables (resources/DE/).
# Writes, under results/simulated/:
#   DE/{study}_Differential_Expression_(diagnosis).tsv
#       Sage's DE rows of the simulated genes only (same columns)
#   {scenario}/aracne/{design}/{table}/{network}/network.txt
#       one file per network, same layout and format as results/aracne/
#   {scenario}/truth_rates.tsv   planted global rewiring per region
#   {scenario}/truth_genes.tsv   planted genes per region and their role
#
# Two scenarios:
#   planted: AD differs from CT by known changes (positive control)
#   null:    AD and CT have the same true network (negative control)
#
# Model:
#   1. One true network per study: a scale-free graph (Barabási-Albert), so it
#      has hubs. Each region rewires 20% of it, so regions differ. This is the
#      CT truth of each region.
#   2. AD truth (planted scenario) = CT truth plus:
#      - global rewiring: a fraction of edges moved at random, different per
#        region (TCX the most, CBE the least: the prediction of analysis 6);
#      - ego rewiring: 10 top DE up and 10 top DE down genes lose 60% of their
#        edges and connect to DE genes of the same direction instead (analysis 4);
#      - 5 gained hubs and 5 lost hubs (analysis 5).
#   3. Each network observes its truth with noise that shrinks with n, like
#      ARACNe: each true edge is missed with probability NOISE / sqrt(n), and
#      false edges are added (half as many as the missed ones). So two networks
#      of the same truth get more similar as n grows (the J_ref curve of
#      analysis 3), and networks with more samples have more edges.

library(dplyr)
library(purrr)
library(vroom)
library(igraph)
source("workflow/scripts/R/subsample_utils.R") # set_seed_for()


# ---- Parameters ---------------------------------------------------------------

SEED <- 2026
OUT_DIR <- "results/simulated"

STUDY_OF <- c(DLPFC = "ROSMAP", ACC = "ROSMAP", PCC = "ROSMAP", TCX = "Mayo", CBE = "Mayo")
DESIGNS <- c("matched", "equal_n", "curve")

# DE gene definition and top DE genes: same values as config/config.yaml (de:)
DE_ADJ_P <- 0.05
DE_FC <- 1.2
TOP_DE <- 25

N_GENES <- 1500 # genes per study (small, so the files stay light)
N_EXTRA_DE <- 300 # extra DE genes in the gene set, as ego rewiring partners
BA_M <- 4 # edges added per gene in the Barabási-Albert graph (mean degree ~ 2 x 4)
REGION_SHIFT <- 0.20 # fraction of the study network rewired in each region

# Planted changes in AD (planted scenario)
GLOBAL_REWIRING <- c(DLPFC = 0.20, ACC = 0.15, PCC = 0.15, TCX = 0.30, CBE = 0.10)
EGO_GENES <- 10 # per direction (up, down), taken from the top DE genes
EGO_REWIRED <- 0.60 # fraction of their edges moved to DE genes of the same direction
N_HUBS <- 5 # gained hubs and lost hubs per region
HUB_NEW_EDGES <- 40 # edges added to each gained hub
HUB_LOST_EDGES <- 0.80 # fraction of edges removed from each lost hub

# Observation noise
NOISE <- 4 # each true edge is missed with probability NOISE / sqrt(n)
FALSE_EDGE_RATIO <- 0.5 # false edges added = 0.5 x missed edges


# ---- Edge helpers ------------------------------------------------------------
# A network is a tibble of undirected edges (a, b) with a < b.

no_edges <- tibble(a = character(), b = character())

canon <- function(edges) {
  edges |>
    transmute(lo = pmin(a, b), hi = pmax(a, b)) |>
    rename(a = lo, b = hi) |>
    filter(a != b) |>
    distinct()
}

# k random new edges between `genes`, none of them already in `exclude`
random_edges <- function(genes, k, exclude) {
  if (k <= 0) {
    return(no_edges)
  }
  tibble(a = sample(genes, 3 * k, replace = TRUE), b = sample(genes, 3 * k, replace = TRUE)) |>
    canon() |>
    anti_join(exclude, by = c("a", "b")) |>
    slice_head(n = k)
}

neighbors_of <- function(edges, gene) {
  c(edges$b[edges$a == gene], edges$a[edges$b == gene])
}

degrees <- function(edges) {
  bind_rows(tibble(gene = edges$a), tibble(gene = edges$b)) |> count(gene, name = "degree")
}

# Move a fraction of the edges to random new pairs (same number of edges)
rewire_global <- function(edges, frac, genes) {
  kept <- slice_sample(edges, prop = 1 - frac)
  bind_rows(kept, random_edges(genes, nrow(edges) - nrow(kept), edges))
}

# A gene loses a fraction of its edges and connects to `partners` instead
rewire_ego <- function(edges, gene, frac, partners) {
  own <- filter(edges, a == gene | b == gene)
  lost <- slice_sample(own, prop = frac)
  candidates <- setdiff(partners, c(gene, neighbors_of(edges, gene)))
  new <- tibble(a = gene, b = sample(candidates, min(nrow(lost), length(candidates))))
  edges |>
    anti_join(lost, by = c("a", "b")) |>
    bind_rows(canon(new))
}

gain_hub <- function(edges, gene, k, genes) {
  candidates <- setdiff(genes, c(gene, neighbors_of(edges, gene)))
  bind_rows(edges, canon(tibble(a = gene, b = sample(candidates, k))))
}

lose_hub <- function(edges, gene, frac) {
  own <- filter(edges, a == gene | b == gene)
  anti_join(edges, slice_sample(own, prop = frac), by = c("a", "b"))
}

# What one network "sees" of its truth with n samples
observe <- function(truth, n_samples, genes) {
  kept <- filter(truth, runif(n()) > NOISE / sqrt(n_samples))
  n_false <- round(FALSE_EDGE_RATIO * (nrow(truth) - nrow(kept)))
  bind_rows(kept, random_edges(genes, n_false, truth))
}

# ARACNe-AP format. With all genes as regulators, ARACNe can report a pair in
# both directions: here half of the edges are written twice (A->B and B->A).
to_aracne <- function(edges) {
  edges <- edges |>
    mutate(
      MI = round(runif(n(), 0.05, 0.8), 4),
      pvalue = signif(runif(n(), 0, 0.05), 3),
      flip = runif(n()) < 0.5,
      both = runif(n()) < 0.5
    )
  bind_rows(
    edges |> transmute(Regulator = if_else(flip, b, a), Target = if_else(flip, a, b), MI, pvalue),
    edges |> filter(both) |> transmute(Regulator = if_else(flip, a, b), Target = if_else(flip, b, a), MI, pvalue)
  )
}

jaccard <- function(x, y) {
  shared <- nrow(semi_join(x, y, by = c("a", "b")))
  shared / (nrow(x) + nrow(y) - shared)
}


# ---- DE tables and gene sets -------------------------------------------------

read_de <- function(study) {
  vroom(sprintf("resources/DE/%s_Differential_Expression_(diagnosis).tsv", study), show_col_types = FALSE) |>
    mutate(region = sub("^AD_(\\w+) - .*$", "\\1", Comparison)) |>
    filter(region %in% names(STUDY_OF)[STUDY_OF == study], Comparison == sprintf("AD_%s - CT_%s", region, region))
}

classify_de <- function(de) {
  de |>
    mutate(
      is_de = adj.P.Val < DE_ADJ_P & abs(logFC) > log2(DE_FC),
      direction = if_else(logFC > 0, "up", "down")
    )
}

# Top DE genes of each region and direction, by |t|
top_de_genes <- function(de) {
  de |>
    filter(is_de) |>
    group_by(region, direction) |>
    slice_max(abs(t), n = TOP_DE, with_ties = FALSE) |>
    ungroup()
}

# Genes of one study: all its top DE genes, extra DE genes, and random genes
gene_set <- function(de, top, study) {
  set_seed_for(SEED, paste("genes", study))
  must <- unique(top$ensembl_gene_id)
  extra_de <- de |>
    filter(is_de, !ensembl_gene_id %in% must) |>
    distinct(ensembl_gene_id) |>
    slice_sample(n = N_EXTRA_DE) |>
    pull()
  rest <- setdiff(unique(de$ensembl_gene_id), c(must, extra_de))
  c(must, extra_de, sample(rest, N_GENES - length(must) - length(extra_de)))
}


# ---- True networks -------------------------------------------------------------

study_truth <- function(genes, study) {
  set_seed_for(SEED, paste("truth", study))
  graph <- sample_pa(length(genes), m = BA_M, directed = FALSE)
  ends <- as_edgelist(graph, names = FALSE)
  shuffled <- sample(genes)
  canon(tibble(a = shuffled[ends[, 1]], b = shuffled[ends[, 2]]))
}

# CT truth and planted AD truth of one region, plus the planted genes
region_truths <- function(region, base, genes, de, top) {
  set_seed_for(SEED, paste("truth", region))
  ct <- rewire_global(base, REGION_SHIFT, genes)

  de_region <- filter(de, region == .env$region, ensembl_gene_id %in% genes)
  ego <- top |>
    filter(region == .env$region) |>
    group_by(direction) |>
    slice_max(abs(t), n = EGO_GENES, with_ties = FALSE) |>
    ungroup() |>
    select(gene = ensembl_gene_id, direction)
  de_partners <- split(de_region$ensembl_gene_id[de_region$is_de], de_region$direction[de_region$is_de])

  deg <- degrees(ct)
  top_hubs <- deg |> slice_max(degree, n = ceiling(0.01 * length(genes)), with_ties = FALSE) |> pull(gene)
  low_degree <- deg |>
    filter(degree <= median(degree), !gene %in% top$ensembl_gene_id) |>
    pull(gene)
  gained <- sample(low_degree, N_HUBS)
  lost <- sample(top_hubs, N_HUBS)

  ad <- rewire_global(ct, GLOBAL_REWIRING[[region]], genes)
  ad <- reduce2(ego$gene, ego$direction, \(e, g, d) rewire_ego(e, g, EGO_REWIRED, de_partners[[d]]), .init = ad)
  ad <- reduce(gained, \(e, g) gain_hub(e, g, HUB_NEW_EDGES, genes), .init = ad)
  ad <- reduce(lost, \(e, g) lose_hub(e, g, HUB_LOST_EDGES), .init = ad)

  planted <- bind_rows(
    tibble(region = region, gene = ego$gene, role = paste0("ego_rewired_", ego$direction)),
    tibble(region = region, gene = gained, role = "hub_gained"),
    tibble(region = region, gene = lost, role = "hub_lost")
  )
  list(CT = ct, AD = ad, planted = planted)
}


# ---- Networks to simulate (same names and n as the real subsamples) ---------

read_networks <- function(file, design) {
  table <- sub("\\.tsv$", "", basename(file))
  vroom(file, col_types = cols(.default = "c")) |>
    count(network, name = "n_samples") |>
    mutate(design = design, table = table, .before = 1)
}

catalog <- DESIGNS |>
  map(\(d) map(list.files(file.path("results/subsamples", d), full.names = TRUE), read_networks, design = d)) |>
  list_flatten() |>
  list_rbind() |>
  mutate(
    region = if_else(design == "curve", sub("_.*$", "", table), table),
    dx = if_else(design == "curve", sub("^[^_]+_([^_]+)_n.*$", "\\1", table), sub("_.*$", "", network)),
    study = STUDY_OF[region]
  )
message("Networks per scenario: ", nrow(catalog))


# ---- Build the truths ----------------------------------------------------------

truths <- list()
planted_genes <- list()

for (study in unique(STUDY_OF)) {
  de <- classify_de(read_de(study))
  top <- top_de_genes(de)
  genes <- gene_set(de, top, study)
  base <- study_truth(genes, study)

  dir.create(file.path(OUT_DIR, "DE"), recursive = TRUE, showWarnings = FALSE)
  de |>
    filter(ensembl_gene_id %in% genes) |>
    select(-region, -is_de, -direction) |>
    vroom_write(file.path(OUT_DIR, "DE", sprintf("%s_Differential_Expression_(diagnosis).tsv", study)))

  for (region in names(STUDY_OF)[STUDY_OF == study]) {
    rt <- region_truths(region, base, genes, de, top)
    truths[[region]] <- list(genes = genes, CT = rt$CT, AD = rt$AD)
    planted_genes[[region]] <- rt$planted
    message(
      region, ": ", length(genes), " genes | true edges CT=", nrow(rt$CT), " AD=", nrow(rt$AD),
      " | J(true AD, true CT)=", round(jaccard(rt$AD, rt$CT), 3)
    )
  }
}


# ---- Write the networks --------------------------------------------------------

simulate_scenario <- function(scenario) {
  networks <- catalog |>
    mutate(edges = pmap(list(design, table, network, region, dx, n_samples), \(design, table, network, region, dx, n_samples) {
      truth <- truths[[region]]
      true_net <- if (scenario == "null") truth$CT else truth[[dx]]
      set_seed_for(SEED, paste(scenario, design, table, network))
      edges <- observe(true_net, n_samples, truth$genes)

      dir <- file.path(OUT_DIR, scenario, "aracne", design, table, network)
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      vroom_write(to_aracne(edges), file.path(dir, "network.txt"))
      edges
    }))

  rates <- tibble(
    region = names(GLOBAL_REWIRING),
    global_rewiring = if (scenario == "null") 0 else unname(GLOBAL_REWIRING)
  )
  genes_out <- if (scenario == "null") {
    tibble(region = character(), gene = character(), role = character())
  } else {
    list_rbind(planted_genes)
  }
  vroom_write(rates, file.path(OUT_DIR, scenario, "truth_rates.tsv"))
  vroom_write(genes_out, file.path(OUT_DIR, scenario, "truth_genes.tsv"))

  networks
}


# ---- Quick check: does the simulation behave as intended? -------------------

# Mean J(AD_repk, CT_all) per region in the matched design, and J(A, B) of the curve
check_scenario <- function(networks, scenario) {
  edges_of <- \(d, t, n) networks$edges[[which(networks$design == d & networks$table == t & networks$network == n)]]

  matched <- networks |>
    filter(design == "matched", dx == "AD") |>
    mutate(J = map2_dbl(table, edges, \(t, e) jaccard(e, edges_of("matched", t, "CT_all")))) |>
    group_by(region) |>
    summarise(n = first(n_samples), mean_J = round(mean(J), 3), sd_J = round(sd(J), 3))

  curve <- networks |>
    filter(design == "curve", grepl("_A$", network)) |>
    mutate(J = pmap_dbl(list(table, network, edges), \(t, n, e) jaccard(e, edges_of("curve", t, sub("_A$", "_B", n))))) |>
    group_by(table) |>
    summarise(n = first(n_samples), J_ref = round(mean(J), 3))

  message("\n== ", scenario, ": J(AD, CT), matched design")
  print(matched)
  message("== ", scenario, ": J(A, B), reference curve")
  print(curve)
}

for (scenario in c("planted", "null")) {
  networks <- simulate_scenario(scenario)
  check_scenario(networks, scenario)
}
message("\nDone: ", OUT_DIR)
