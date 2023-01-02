log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)

individual_genotype_likelihoods <- snakemake@input %>% 
  map(~ read_tsv(.)) %>%
  reduce(bind_rows)

number_of_taxa <- individual_genotype_likelihoods %>%
  select(cell) %>%
  n_distinct()

all_vars_all_cells <- individual_genotype_likelihoods %>%
  expand(CHROM, POS, REF, ALT, cell)

flat_likelihood <- str_c(rep(0.1, 10), collapse=",")

all_individual_genotype_likelihoods <- individual_genotype_likelihoods %>%
  right_join(
    all_vars_all_cells
  ) %>%
  replace_na(
    list(
      ml_genotype = "N",
      likelihoods = flat_likelihood
    )
  ) %>%
  pivot_wider(
    names_from = cell,
    values_from = c(ml_genotype, likelihoods)
  ) %>%
  unite(
    starts_with("ml_genotype_"),
    sep = "",
    col = "consensus_genotypes"
  ) %>%
  rename_with(
    ~ str_replace(., "likelihoods_", ""),
    starts_with("likelihoods_"),
  ) %>%
  select(
    !c(CHROM, POS, REF, ALT)
  )

alignment_sites <- all_individual_genotype_likelihoods %>% nrow()

row_one <- tibble_row(n_taxa = number_of_taxa, alignment_sites = alignment_sites)
write_tsv(
  row_one,
  snakemake@output[[1]],
  col_names = FALSE
)

row_two <- all_individual_genotype_likelihoods %>%
  select(!consensus_genotypes) %>%
  colnames() %>%
  as_tibble() %>%
  pivot_wider(names_from = value)
write_tsv(
  row_two, snakemake@output[[1]],
  append = TRUE,
  col_names = FALSE
)

write_tsv(
  all_individual_genotype_likelihoods,
  snakemake@output[[1]],
  append = TRUE,
  col_names = FALSE
)