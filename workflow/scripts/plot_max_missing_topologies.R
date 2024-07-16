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
      convergence,
      relative_rf_distance,
      absolute_rf_distance
    ),
    names_to = "label",
    values_to = "value"
  ) |>
  mutate(
    value_type = case_match(
      label,
      c("total_trees", "unique_topologies", "convergence") ~ "count",
      "relative_rf_distance" ~ "Robinson-Foulds distance, relative",
      "absolute_rf_distance" ~ "Robinson-Foulds distance, absolute"
    ),
    `count of` = case_match(
      label,
      "total_trees" ~ "total trees",
      "unique_topologies" ~ "distinct topologies",
      "convergence" ~ "bootstrap trees for convergence"
    ),
    `source of tree set` = case_match(
      tree_type,
      "startTree" ~ "start of maximum likelihood search",
      "mlTrees" ~ "after maximum likelihood search",
      "bootstraps" ~ "bootstrapping"
    )
  )

plot <- ggplot(
  rfs_and_topos,
  aes(
    x = max_missing,
    y = value,
    color = `source of tree set`,
    shape = `count of`
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