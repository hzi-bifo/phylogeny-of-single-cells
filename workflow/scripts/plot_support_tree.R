log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library("treeio")
library("ggtree")
library("tidyverse")

support_tree <- read.newick(snakemake@input[["support"]])

rooting_group <- read_tsv(snakemake@input[["samples"]]) |>
  filter(
    cell_type == "Endothel",
    individual == snakemake@wildcards[["individual"]],
    sample_type == "single_cell",
    sample_name %in% (support_tree |> as_tibble() |> pull(label))
  ) |>
  pull(sample_name)

print("Rooting group: ")
print(str_flatten(rooting_group, ", "))

safe_root <- safely(root)

rooted_support_tree <- safe_root(
  support_tree,
  outgroup = rooting_group,
  resolve.root = TRUE
)

if ( is_null(rooted_support_tree$result) ) {
  print("Could not root tree: ")
  print(rooted_support_tree$error$message)
  print("Will plot unrooted tree.")
  tree <- support_tree
} else {
  tree <- rooted_support_tree$result
}

tree_plot <- ggtree(
  tree
) + 
  geom_tiplab() +
  geom_nodelab(
    hjust = 0
  )

ggsave(
  snakemake@output[["support"]],
  plot = tree_plot,
  width = 10,
  height = 12
)