#!/usr/bin/env Rscript

library(tidyverse)

# helper to find repetitions on a vector
find_reps <- function(x) {
  tib <- tibble(X = x)
  tib |>
    group_by(X) |>
    summarize(n = n()) |>
    filter(n != 1)
}


IN_neigh <- "data/01_neigh_bsub.tsv"
IN_genomes <- "data/02_genomes_bsub.tsv"
IN_iscan <- "data/03_iscan.tsv"
IN_seqs <- "data/04_endo_neigh.tsv"

neigh <- read_tsv(IN_neigh)

# info to append to neighbors
endo <- neigh |>
  filter(q_alias == "endo", Nseq == 0) |>
  mutate(endo_neigh = str_c(genome, "_", contig)) |>
  relocate(endo_neigh) |>
  rename(pid_endo = pid,
         strand_endo = strand,
         order_endo = order,
         start_endo = start,
         end_endo = end) |>
  select(endo_neigh, pid_endo, strand_endo,
         order_endo, start_endo, end_endo)


# generate key to cross
endo_neighs <- neigh |>
  filter(q_alias == "endo", Nseq != 0) |>
  mutate(endo_neigh = str_c(genome, "_", contig))
  

# annotate neighbors
endo_neighs <- left_join(endo_neighs, endo, join_by(endo_neigh))