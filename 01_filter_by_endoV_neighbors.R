#!/usr/bin/env Rscript

library(tidyverse)

IN <- "data/01_neigh_bsub.tsv"
neigh <- read_tsv(IN)

neigh_endo <- neigh |>
  filter(q_alias == "endo") |>
  filter(Nseq != 0)

neigh_endo |>
  pull(pid) |>
  unique() |>
  sort() |>
  write_lines("pids_endo_neigh.txt")

