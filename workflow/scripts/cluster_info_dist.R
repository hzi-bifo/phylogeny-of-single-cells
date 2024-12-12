log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library("TreeDist")
library("tidyverse")

trees <- ape::read.tree(snakemake@input[["trees"]])

pairwise_cid <- ClusteringInfoDistance(trees, trees, normalize=TRUE)

dedup_cids <- pairwise_cid[lower.tri(pairwise_cid)]

annotated_tibble <- dedup_cids |>
  as_tibble_col(column_name = "cluster_information_distance") |>
  add_column(
    model=snakemake@wildcards[["model"]],
    max_missing=snakemake@wildcards[["n_missing_cells"]],
    tree_type=snakemake@wildcards[["type"]],
    software=snakemake@wildcards[["software"]]
  )

write_tsv(
  x = annotated_tibble,
  file = snakemake@output[["cluster_info_dist"]]
)

