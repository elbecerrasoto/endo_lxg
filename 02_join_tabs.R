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
  rename(
    pid_endo = pid,
    locus_endo = locus_tag,
    strand_endo = strand,
    order_endo = order,
    start_endo = start,
    end_endo = end
  ) |>
  select(
    endo_neigh, pid_endo,
    locus_endo, strand_endo,
    order_endo, start_endo, end_endo
  )


# generate key to cross
endo_neighs <- neigh |>
  filter(q_alias == "endo", Nseq != 0) |>
  mutate(endo_neigh = str_c(genome, "_", contig))


# annotate neighbors
endo_neighs <- left_join(endo_neighs, endo, join_by(endo_neigh))

endo_neighs <- endo_neighs |>
  mutate(
    same_strand = strand == strand_endo,
    distance2lxg = start - start_endo
  )


#---- iscan


LXG <- "IPR006829"

ISCAN_NAMES <- c(
  "pid",
  "md5",
  "length",
  "analysis",
  "memberDB",
  "memberDB_txt",
  "start",
  "end",
  "score",
  "recommended",
  "date",
  "interpro",
  "interpro_txt",
  "GO",
  "residue"
)

iscan <- read_tsv(IN_iscan,
  col_names = ISCAN_NAMES,
  na = c("-", "NA", "")
)

domains <- iscan |>
  filter(recommended, !is.na(interpro)) |>
  distinct(pid, interpro, .keep_all = T) |>
  select(pid, interpro, interpro_txt) |>
  mutate(lxg = interpro == LXG)


domains_wide <- domains |>
  group_by(pid) |>
  reframe(domains = str_flatten(sort(unique(interpro)), collapse = ";"))


domains2cross <- domains |>
  left_join(domains_wide, join_by(pid)) |>
  filter(lxg)

#---- endo neighs annotate
endo_neighs <- semi_join(
  endo_neighs, domains2cross,
  join_by(pid)
)

endo_neighs <- left_join(
  endo_neighs, domains2cross,
  join_by(pid)
)



#---- add genome, seqs, annots

genomes <- read_tsv(IN_genomes)
seqs <- read_tsv(IN_seqs, col_names = c("pid", "seq"))


endo_neighs <- left_join(endo_neighs, genomes, join_by(genome))
endo_neighs <- left_join(endo_neighs, seqs, join_by(pid))


COUTS <- c("Nseq", "distance2lxg", "strand", "strand_endo", "gene", "product", "pid", "locus_tag", "genome", "contig", "pid_endo", "locus_endo", "domains", "seq")

endo_neighs <- endo_neighs |>
  select(all_of(COUTS)) |>
  rename(genes2lxg = Nseq)

endo_neighs <- endo_neighs |>
  arrange(genome)

write_tsv(endo_neighs, "endo_lxg_neighs_sig.tsv")

endo_neighs_abs <- endo_neighs |>
  mutate(
    genes2lxg = abs(genes2lxg),
    distance2lxg = abs(distance2lxg)
  ) |>
  arrange(desc(genes2lxg), desc(distance2lxg))

write_tsv(endo_neighs_abs, "endo_lxg_neighs.tsv")
