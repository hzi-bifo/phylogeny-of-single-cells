log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)

raxml_ng_genotype_lookup_table <- tribble(
   ~a,  ~b, ~iupac,
  "A", "A",    "A",
  "C", "C",    "C",
  "G", "G",    "G",
  "T", "T",    "T",
  "A", "C",    "M",
  "A", "G",    "R",
  "A", "T",    "W",
  "C", "A",    "M",
  "C", "G",    "S",
  "C", "T",    "Y",
  "G", "A",    "R",
  "G", "C",    "S",
  "G", "T",    "K",
  "T", "A",    "W",
  "T", "C",    "Y",
  "T", "G",    "K"
)

iupac_of <- function(a, b, table) {
  return table %>%
    filter( a == a, b == b) %>%
    pull(iupac)
}

to_raxml_ng_genotype <- function(genotype, ref, alt, lookup_table) {
  if ( genotype == "hom_ref" ) {
    iupac <- iupac_of(ref, ref, lookup_table)
  } else if ( genotype == "het" ) {
    iupac <- iupac_of(ref, alt, lookup_table)
  } else if ( genotype == "hom_alt" ) {
    iupac <- iupac_of(alt, alt, lookup_table)
  } else {
    cli_abort("Unknown genotype name. Please only use 'hom_ref', 'het', or 'hom_alt'.")
  }
  iupac
}

raxml_ng_genotype_order <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K")

to_raxml_ng_genotype_likelihoods <- function(ref, alt, lookup_table, order, hom_ref_likelihood, het_likelihood, hom_alt_likelihood) {
  likelihoods = rep(0.0, 10)

  hom_ref_index <- which(iupac_of(ref, ref, lookup_table), order)
  likelihoods[hom_ref_index] <- hom_ref_likelihood

  het_index <- which(iupac_of(ref, alt, lookup_table), order)
  likelihoods[het_index] <- het_likelihood

  hom_alt_index <- which(iupac_of(alt, alt, lookup_table), order)
  likelihoods[hom_alt_index] <- hom_alt_likelihood

  likelihoods
}

read_genotype <- function(genotype) {
  return read_tsv(snakemake@input[[genotype]]) %>%
    mutate(
      ml_genotype = to_raxml_ng_genotype(genotype, REF, ALT, raxml_ng_genotype_lookup_table)
      likelihoods = to_raxml_ng_genotype_likelihoods(REF, ALT, raxml_ng_genotype_lookup_table, HOM_REF, HET, HOM_ALT)
    ) %>%
    select(
      !c(HOM_REF, HET, HOM_ALT)
    )
}

all_genotypes <- read_genotype("hom_ref") %>%
  bind_rows(
    read_genotype("het")
  ) %>%
  bind_rows(
    read_genotype("hom_alt")
  ) %>%
  add_column(cell = snakemake@wildcards.sc)

write_tsv(all_genotypes, snakemake@output[[1]])