# Design 3 (reference J): two disjoint subsamples of size n from the same group
# ("rep01_A", "rep01_B", ...). No donor is in both A and B.
# The female proportion is the one of the group itself.

source(snakemake@input[["utils"]])

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

pool <- read_samples(snakemake@input[["samples"]])
n <- as.integer(snakemake@wildcards[["n"]])
dx <- snakemake@wildcards[["dx"]]
reps <- rep_label(seq_len(snakemake@params[["replicates"]]))

if (2 * n > nrow(pool)) {
  stop("Two disjoint subsamples of n = ", n, " need ", 2 * n, " donors; the group has ", nrow(pool))
}

set_seed_for(
  snakemake@params[["seed"]],
  paste("curve", snakemake@wildcards[["region"]], dx, n)
)
female_prop <- mean(pool$sex == "female")

# Two disjoint halves (A, B) per replicate
subsamples <- reps |>
  map(\(rep) sample_two_disjoint(pool, n, female_prop) |> bind_rows(.id = "half")) |>
  set_names(reps) |>
  bind_rows(.id = "replicate") |>
  mutate(network = paste0(replicate, "_", half), diagnosis = dx) |>
  select(network, diagnosis, specimenID, individualID, sex, RIN)

vroom_write(subsamples, snakemake@output[[1]])

message("n=", n, " x 2 | group=", nrow(pool), " (", round(100 * female_prop, 1), "% female)")
print(summarise_networks(subsamples), n = Inf)
