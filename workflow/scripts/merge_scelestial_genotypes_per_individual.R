log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)

individual_genotype_likelihoods <- snakemake@input %>% 
  map(~ read_tsv(.)) %>%
  reduce(bind_rows)

all_vars_all_cells <- individual_genotype_likelihoods %>%
  expand(position, cell)

flat_likelihood <- str_c(rep(0.1, 10), collapse=",")

all_individual_genotype_likelihoods <- individual_genotype_likelihoods %>%
  right_join(
    all_vars_all_cells
  ) %>%
  replace_na(
    list(
      genotype = "./."
    )
  ) %>%
  pivot_wider(
    names_from = cell,
    values_from = genotype
  )

write_tsv(
  all_individual_genotype_likelihoods,
  snakemake@output[[1]],
  col_names = FALSE
)