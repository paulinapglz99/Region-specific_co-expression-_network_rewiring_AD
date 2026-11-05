# Samples of one region and one diagnosis, keeping one sample per donor.
#
# Some donors were sequenced twice in the same region (different batches).
# We keep the sample with the highest RIN (ties: first specimenID alphabetically).

library(data.table)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

region <- snakemake@wildcards[["region"]]
dx <- snakemake@wildcards[["dx"]]

meta <- fread(
  snakemake@input[["metadata"]],
  colClasses = c(age_death = "character", apoe4_allele = "character")
)
group <- meta[tissue == region & diagnosis == dx]

if (nrow(group) == 0) {
  stop("No samples found for tissue=", region, ", diagnosis=", dx)
}

setorder(group, individualID, -RIN, specimenID)
selected <- group[!duplicated(individualID)]

fwrite(
  selected[, .(specimenID, individualID, sex, age_death, apoe4_allele, RIN, pmi)],
  snakemake@output[[1]],
  sep = "\t"
)

message("region=", region, " diagnosis=", dx)
message("samples=", nrow(group), " donors=", nrow(selected))
dropped <- group[!specimenID %in% selected$specimenID, specimenID]
if (length(dropped) > 0) message("dropped duplicates: ", paste(dropped, collapse = ", "))
