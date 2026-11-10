# Analysis 4: GO over-representation of the "rewired only" genes (rewired, not
# DE) of each region. Universe = the genes of that region that pass the degree
# filter of rewiring_vs_de.R, so a term is enriched compared with genes that
# could have been called rewired. BH correction within each region (enrichGO).

library(dplyr)
library(purrr)
library(vroom)
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

genes <- vroom(snakemake@input[["genes"]], col_types = cols(region = "c", gene = "c", group = "c"))

enrich_region <- function(df, region) {
  result <- enrichGO(
    gene = df$gene[df$group == "rewired only"],
    universe = df$gene,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = snakemake@params[["ontology"]],
    minGSSize = snakemake@params[["min_size"]],
    maxGSSize = snakemake@params[["max_size"]],
    pvalueCutoff = 1, # keep every tested term; the significance filter is p_adj
    qvalueCutoff = 1
  )
  if (is.null(result)) {
    return(NULL)
  }
  as_tibble(result) |>
    transmute(region = region, term = ID, description = Description, gene_ratio = GeneRatio,
              bg_ratio = BgRatio, count = Count, p = pvalue, p_adj = p.adjust, genes = geneID)
}

go <- genes |>
  group_by(region) |>
  group_map(\(df, key) enrich_region(df, key$region)) |>
  list_rbind()

if (nrow(go) == 0) {
  go <- tibble(region = character(), term = character(), description = character(), gene_ratio = character(),
               bg_ratio = character(), count = integer(), p = double(), p_adj = double(), genes = character())
}

vroom_write(go, snakemake@output[[1]])

message("Terms with p_adj < ", snakemake@params[["max_adj_p"]], " per region:")
print(count(filter(go, p_adj < snakemake@params[["max_adj_p"]]), region))
