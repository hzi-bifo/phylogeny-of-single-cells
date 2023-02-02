log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)

to_scelestial_genotype <- function(genotype, ref, alt) {
  if ( genotype == "hom_ref" ) {
    gt <- str_c(ref, "/", ref)
  } else if ( genotype == "het" ) {
    gt <- str_c(ref, "/", alt)
  } else if ( genotype == "hom_alt" ) {
    gt <- str_c(alt, "/", alt)
  } else {
    cli_abort("Unknown genotype name. Please only use 'hom_ref', 'het', or 'hom_alt'.")
  }
  gt
}

read_genotype <- function(genotype) {
  return read_tsv(snakemake@input[[genotype]]) %>%
    mutate(
      genotype = to_scelestial_genotype(genotype, REF, ALT)
    ) %>%
    select(
      c(position, genotype)
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