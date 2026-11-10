# Analysis 5: are the hubs that change between AD and CT DE genes? (all regions)
#   1. PageRank of every gene in each network (undirected, unweighted).
#      Hub of a network = top `top_pagerank` of the genes of the study
#      (k = top_pagerank x genes of the DE table, ~190 for 1%).
#   2. Hub type of each gene, from CT_all and the 10 AD replicates:
#        shared   hub in CT and in >= min_replicates AD replicates
#        gained   not hub in CT, hub in >= min_replicates AD replicates
#        lost     hub in CT, hub in <= max_replicates_lost AD replicates
#        unstable hub in CT and in between: not classified
#   3. Overlap test of each hub type x DE (de_utils.R): DE labels shuffled
#      within deciles of AveExpr. Not within degree: hubs are defined by their
#      connections, so matching on degree would remove what is being tested.
#      BH across regions and hub types.
# Universe: genes of the DE table with >= 1 edge in at least one network of the region.

source(snakemake@input[["network_utils"]])
source(snakemake@input[["de_utils"]])
library(tidyr)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

study_of <- unlist(snakemake@params[["study_of"]])
top_pagerank <- snakemake@params[["top_pagerank"]]
min_reps <- snakemake@params[["min_replicates"]]
max_reps_lost <- snakemake@params[["max_replicates_lost"]]
n_perm <- snakemake@params[["n_perm"]]
bins <- snakemake@params[["strata"]]
set.seed(snakemake@params[["seed"]])

de_files <- unlist(snakemake@input[["de"]])

# Network files: .../matched/{region}/{network}/network.txt
networks <- tibble(file = unlist(snakemake@input[["networks"]])) |>
  mutate(
    network = basename(dirname(file)),
    region = basename(dirname(dirname(file)))
  )

pagerank <- function(file) {
  graph <- graph_from_data_frame(read_network(file), directed = FALSE)
  tibble(gene = V(graph)$name, pagerank = page_rank(graph)$vector)
}

hub_table <- function(region) {
  de <- read_de(de_file_of(de_files, study_of[[region]]), region, snakemake@params[["adj_p"]], snakemake@params[["fold_change"]])
  k <- round(top_pagerank * nrow(de))

  scores <- networks |>
    filter(region == .env$region) |>
    mutate(scores = map(file, pagerank)) |>
    select(network, scores) |>
    unnest(scores) |>
    filter(gene %in% de$gene) |>
    mutate(hub = rank(-pagerank, ties.method = "first") <= k, .by = network)

  scores |>
    summarise(
      pagerank_ct = sum(pagerank[network == "CT_all"]),
      ct_hub = any(hub & network == "CT_all"),
      ad_hub_replicates = sum(hub & network != "CT_all"),
      .by = gene
    ) |>
    mutate(hub_type = case_when(
      ct_hub & ad_hub_replicates >= min_reps ~ "shared",
      !ct_hub & ad_hub_replicates >= min_reps ~ "gained",
      ct_hub & ad_hub_replicates <= max_reps_lost ~ "lost",
      ct_hub ~ "unstable",
      TRUE ~ "not hub"
    )) |>
    inner_join(de, by = "gene") |>
    filter(!is.na(AveExpr)) |>
    add_strata("AveExpr", bins) |>
    mutate(region = region, k = k, .before = 1)
}

genes <- map(unique(networks$region), hub_table) |> list_rbind()

tests <- expand_grid(region = unique(genes$region), hub_type = c("shared", "gained", "lost")) |>
  pmap(\(region, hub_type) {
    df <- genes[genes$region == region, ]
    overlap_test(df$hub_type == hub_type, df$is_de, df$stratum, n_perm) |>
      mutate(region = region, hub_type = hub_type, .before = 1)
  }) |>
  list_rbind() |>
  mutate(p_adj = p.adjust(p_perm, method = "BH"))

vroom_write(select(genes, -stratum), snakemake@output[["genes"]])
vroom_write(tests, snakemake@output[["tests"]])

print(count(genes, region, hub_type) |> pivot_wider(names_from = hub_type, values_from = n))
print(mutate(tests, across(where(is.double), \(v) signif(v, 3))), n = Inf, width = Inf)
