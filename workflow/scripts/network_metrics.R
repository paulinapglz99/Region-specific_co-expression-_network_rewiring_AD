# Compare the networks of one subsample table, pair by pair (analyses 1-3):
#   matched/{region}: AD_repk vs CT_all
#   equal_n/{region}: AD_repk vs CT_repk
#   curve/{table}:    repk_A  vs repk_B
# One row per pair: edge Jaccard, Spearman rho of degree, NMI of Leiden
# communities, and the size of both networks.

source(snakemake@input[["utils"]])

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

design <- snakemake@wildcards[["design"]]
files <- unlist(snakemake@input[["networks"]])
resolution <- snakemake@params[["resolution"]]
seed <- snakemake@params[["seed"]]

# Samples per network, from the subsample table
n_samples <- vroom(snakemake@input[["subsamples"]], col_types = cols(.default = "c")) |>
  count(network, name = "n")

# Networks and their communities, by network name (the folder name)
networks <- set_names(map(files, read_network), basename(dirname(files)))
communities <- map(networks, \(edges) leiden_communities(edges, resolution, seed))

# Which network is compared with which
ad_networks <- grep("^AD_", names(networks), value = TRUE)
pairs <- switch(design,
  matched = tibble(net_x = ad_networks, net_y = "CT_all"),
  equal_n = tibble(net_x = ad_networks, net_y = sub("^AD_", "CT_", ad_networks)),
  curve = tibble(net_x = grep("_A$", names(networks), value = TRUE)) |>
    mutate(net_y = sub("_A$", "_B", net_x)),
  stop("Unknown design: ", design)
)

compare_pair <- function(net_x, net_y) {
  x <- networks[[net_x]]
  y <- networks[[net_y]]
  genes_x <- nrow(degrees(x))
  genes_y <- nrow(degrees(y))
  tibble(
    net_x = net_x,
    net_y = net_y,
    jaccard = jaccard(x, y),
    spearman_degree = degree_spearman(x, y),
    nmi = communities_nmi(communities[[net_x]], communities[[net_y]]),
    edges_x = nrow(x),
    edges_y = nrow(y),
    genes_x = genes_x,
    genes_y = genes_y,
    mean_degree_x = 2 * nrow(x) / genes_x, # over the genes with at least one edge
    mean_degree_y = 2 * nrow(y) / genes_y
  )
}

metrics <- pairs |>
  pmap(compare_pair) |>
  list_rbind() |>
  mutate(
    design = design,
    table = snakemake@wildcards[["table"]],
    replicate = sub("^(AD|CT)_|_[AB]$", "", net_x),
    .before = 1
  ) |>
  left_join(rename(n_samples, net_x = network, n_x = n), by = "net_x") |>
  left_join(rename(n_samples, net_y = network, n_y = n), by = "net_y")

vroom_write(metrics, snakemake@output[[1]])

message(design, "/", snakemake@wildcards[["table"]], ": ", nrow(metrics), " pairs")
print(summarise(metrics, across(c(jaccard, spearman_degree, nmi), \(v) round(mean(v), 3))))
