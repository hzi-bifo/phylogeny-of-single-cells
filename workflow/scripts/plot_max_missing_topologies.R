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
      absolute_rf_distance,
      usable_sites,
      aic_score,
      aicc_score,
      bic_score
    ),
    names_to = "label",
    values_to = "value"
  ) |>
  filter(
    label != "relative_rf_distance"
  ) |>
  mutate(
    value_type = case_match(
      label,
      "usable_sites" ~ "usable sites",
      c("aic_score", "aicc_score", "bic_score") ~ "information criteria scores",
      c("total_trees", "unique_topologies", "convergence") ~ "number of trees",
#      "relative_rf_distance" ~ "Robinson-Foulds distance, relative",
      "absolute_rf_distance" ~ "Robinson-Foulds distance"
    ),
    value_type = factor(
      value_type,
      levels = c(
        "usable sites",
        "information criteria scores",
        "Robinson-Foulds distance",
        "number of trees"
      )
    ),
    `number of trees` = case_match(
      label,
      c(
        "total_trees",
        "absolute_rf_distance",
        "usable_sites",
        "aic_score",
        "aicc_score",
        "bic_score"
      ) ~ "total",
      "unique_topologies" ~ "distinct topologies",
      "convergence" ~ "bootstrapping converges"
    ),
    `trees from` = case_match(
      tree_type,
      "startTree" ~ "start of ML search",
      "mlTrees" ~ "result of ML search",
      "bootstraps" ~ "bootstrapping"
    ),
    `information criteria scores` = case_match(
      label,
      "aic_score" ~ "Akaike",
      "aicc_score" ~ "Akaike corrected",
      "bic_score" ~ "Bayesian",
      .default = "non-score lines"
    ),
    `information criteria scores` = factor(
      `information criteria scores`,
      levels = c(
      "non-score lines",
      "bic_score" ~ "Bayesian",
      "aicc_score" ~ "Akaike corrected",
      "aic_score" ~ "Akaike"
      )
    )
  )

plot <- ggplot(
  rfs_and_topos,
  aes(
    x = max_missing,
    y = value,
    color = `number of trees`,
    linetype = `information criteria scores`,
    shape = `trees from`
  )
) +
  scale_color_brewer(
    palette = "Dark2"
  ) +
  geom_line(
  ) +
  geom_point(
  ) +
  facet_grid(
    rows = vars(value_type),
    cols = vars(factor(model)),
    scales = "free_y"
  ) +
  ylim(0,NA) +
  theme_bw()

number_of_models <- rfs_and_topos |> distinct(model) |> count() |> pull(n)
plot_width = 3 + number_of_models * 3.5


ggsave(
  snakemake@output[["plot"]],
  plot = plot,
  width = plot_width,
  limitsize = FALSE
)