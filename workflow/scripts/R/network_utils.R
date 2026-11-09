# Helpers to read and compare ARACNe-AP networks.
# A network is a tibble of undirected edges (a, b) with a < b: with all genes as
# regulators, ARACNe can report A->B and B->A, and both count as one edge.

library(dplyr)
library(purrr)
library(vroom)
library(igraph)

# Read a network.txt of ARACNe-AP (columns Regulator, Target, MI, pvalue)
read_network <- function(file) {
  vroom(file, col_types = cols(Regulator = "c", Target = "c", .default = "d")) |>
    transmute(lo = pmin(Regulator, Target), hi = pmax(Regulator, Target)) |>
    rename(a = lo, b = hi) |>
    distinct()
}

# Degree of every gene with at least one edge
degrees <- function(edges) {
  bind_rows(tibble(gene = edges$a), tibble(gene = edges$b)) |>
    count(gene, name = "degree")
}

# Edges in both / edges in at least one
jaccard <- function(x, y) {
  shared <- nrow(semi_join(x, y, by = c("a", "b")))
  shared / (nrow(x) + nrow(y) - shared)
}

# Spearman correlation of the degree, over the genes with edges in both networks
degree_spearman <- function(x, y) {
  both <- inner_join(degrees(x), degrees(y), by = "gene")
  cor(both$degree.x, both$degree.y, method = "spearman")
}

# Leiden communities (modularity objective). Returns one row per gene with edges
leiden_communities <- function(edges, resolution, seed) {
  graph <- graph_from_data_frame(edges, directed = FALSE)
  set.seed(seed)
  partition <- cluster_leiden(
    graph,
    objective_function = "modularity",
    resolution = resolution,
    n_iterations = -1 # repeat until the partition stops changing
  )
  tibble(gene = V(graph)$name, community = as.integer(membership(partition)))
}

# NMI of two partitions, over the genes with edges in both networks
# (isolated genes would form one community each and inflate the NMI)
communities_nmi <- function(cx, cy) {
  both <- inner_join(cx, cy, by = "gene")
  compare(both$community.x, both$community.y, method = "nmi")
}
