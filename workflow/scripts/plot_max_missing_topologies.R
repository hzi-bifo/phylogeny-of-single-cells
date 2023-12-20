log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")

rfs_and_topos <- read_tsv(
  snakemake@input[["tsv"]]
) |>
  pivot_longer(
    c(
      total_trees,
      unique_topologies,
      relative_rf_distance,
      absolute_rf_distance
    ),
    names_to = "label",
    values_to = "value"
  ) |>
  mutate(
    value_type = case_match(
      label,
      c("total_trees", "unique_topologies") ~ "count",
      "relative_rf_distance" ~ "relative distance",
      "absolute_rf_distance" ~ "absolute distance"
    ),
    count_type = case_match(
      label,
      "total_trees" ~ "number of trees",
      "unique_topologies" ~ "number of topologies"
    ),
    tree_type = case_match(
      tree_type,
      "startTree" ~ "starting trees",
      "mlTrees" ~ "maximum likelihood trees",
      "bootstraps" ~ "bootstrap trees"
    )
  )

plot <- ggplot(
  rfs_and_topos,
  aes(
    x = max_missing,
    y = value,
    color = tree_type,
    shape = count_type
  )
) +
  scale_color_brewer(
    palette = "Dark2"
  ) +
  geom_line() +
  geom_point() +
  facet_grid(
    rows = vars(value_type),
    cols = vars(model),
    scales = "free_y"
  )

ggsave(
  snakemake@output[["plot"]],
  plot = plot,
)