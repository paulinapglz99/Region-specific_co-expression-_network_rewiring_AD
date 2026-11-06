# Samples of one region and one diagnosis, keeping one sample per donor.
#
# Some donors were sequenced twice in the same region (different batches).
# We keep the sample with the highest RIN (ties: first specimenID alphabetically).

library(dplyr)
library(vroom)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

region <- snakemake@wildcards[["region"]]
dx <- snakemake@wildcards[["dx"]]

# IDs as text: Mayo donor IDs look like numbers
meta <- vroom(
  snakemake@input[["metadata"]],
  col_types = cols(specimenID = "c", individualID = "c", age_death = "c", apoe4_allele = "c")
)

group <- meta |>
  filter(tissue == region, diagnosis == dx)

if (nrow(group) == 0) {
  stop("No samples found for tissue=", region, ", diagnosis=", dx)
}

selected <- group |>
  arrange(individualID, desc(RIN), specimenID) |>
  distinct(individualID, .keep_all = TRUE) |>
  select(specimenID, individualID, sex, age_death, apoe4_allele, RIN, pmi)

vroom_write(selected, snakemake@output[[1]])

message("region=", region, " diagnosis=", dx)
message("samples=", nrow(group), " donors=", nrow(selected))
dropped <- setdiff(group$specimenID, selected$specimenID)
if (length(dropped) > 0) message("dropped duplicates: ", paste(dropped, collapse = ", "))
