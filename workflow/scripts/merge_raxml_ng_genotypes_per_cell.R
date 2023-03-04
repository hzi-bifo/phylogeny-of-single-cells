log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)

lookup_table <- c(
  "AA" = "A",
  "CC" = "C",
  "GG" = "G",
  "TT" = "T",
  "AC" = "M",
  "AG" = "R",
  "AT" = "W",
  "CA" = "M",
  "CG" = "S",
  "CT" = "Y",
  "GA" = "R",
  "GC" = "S",
  "GT" = "K",
  "TA" = "W",
  "TC" = "Y",
  "TG" = "K"
)

to_raxml_ng_genotype <- function(ref, alt, gt, iupac_of) {
  if ( gt == "hom_ref" ) {
    iupac <- iupac_of[[str_c(ref, ref)]]
  } else if ( gt == "het" ) {
    iupac <- iupac_of[[str_c(ref, alt)]]
  } else if ( gt == "hom_alt" ) {
    iupac <- iupac_of[[str_c(alt, alt)]]
  } else {
    cli_abort("Unknown genotype name. Please only use 'hom_ref', 'het', or 'hom_alt'.")
  }
  iupac
}

raxml_ng_genotype_order <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K")

to_raxml_ng_genotype_likelihoods <- function(ref, alt, hom_ref_likelihood, het_likelihood, hom_alt_likelihood, iupac_of, order) {
  likelihoods = rep(0.0, 10)

  hom_ref_index <- which(iupac_of[[str_c(ref, ref)]] == order)
  likelihoods[hom_ref_index] <- hom_ref_likelihood

  het_index <- which(iupac_of[[str_c(ref, alt)]] == order)
  likelihoods[het_index] <- het_likelihood

  hom_alt_index <- which(iupac_of[[str_c(alt, alt)]] == order)
  likelihoods[hom_alt_index] <- hom_alt_likelihood

  str_c(likelihoods, collapse=",")
}

read_genotype <- function(genotype) {
  read_tsv(snakemake@input[[genotype]]) %>%
    mutate(
      ml_genotype = pmap_chr(
        .,
        function(REF, ALT, ... )
        to_raxml_ng_genotype(
          ref= REF,
          alt = ALT,
          gt = genotype,
          iupac_of = lookup_table
        )
      ),
      likelihoods = pmap_chr(
        .,
        function(REF, ALT, HOM_REF, HET, HOM_ALT, ...)
        to_raxml_ng_genotype_likelihoods(
          ref = REF,
          alt = ALT, 
          hom_ref_likelihood = HOM_REF,
          het_likelihood = HET,
          hom_alt_likelihood = HOM_ALT,
          iupac_of = lookup_table,
          order = raxml_ng_genotype_order
        )
      )
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
  add_column(cell = snakemake@wildcards[["sc"]])

write_tsv(all_genotypes, snakemake@output[[1]])