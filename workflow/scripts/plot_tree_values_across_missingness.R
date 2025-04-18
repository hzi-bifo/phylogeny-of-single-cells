log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("treeio")
library("tidyverse")

rlang::global_entrace()

best_trees_support_values <- snakemake@input[["support_trees"]] |>
  map(
    \(filename)
    read.newick(
      filename,
      node.label = "support"
    ) |>
    as_tibble() |>
    drop_na(support) |>
    add_column(
      max_missing = str_extract(
        filename,
        "max_(\\d+)_missing\\.support",
        group = 1
      ),
      model = str_extract(
        filename,
        "results/([^/]+)/max_",
        group = 1
      ),
      tree_type = str_extract(
        filename,
        "support\\.([^.]+)\\.raxml\\.collapsed.support$",
        group = 1
      )
    ) |>
  mutate(max_missing = as.integer(max_missing)) |>
  select(support, tree_type, max_missing, model)) |>
  list_rbind() |>
  mutate(
    tree_type = str_c(tree_type, " bootstrap support")
  ) |>
  rename(
    value = support,
    value_type = tree_type
  )
  
cluster_info_dist <- read_tsv(snakemake@input[["cluster_info_dist"]]) |>
  mutate(
    tree_type = str_c(tree_type, " cluster info dist")
  ) |>
  rename(
    value = cluster_information_distance,
    value_type = tree_type
  )

all_values <- bind_rows(
    best_trees_support_values,
    cluster_info_dist
  )

highest_mean_besttree_bootstrap_support = all_values |>
  filter(
    value_type == "bestTree bootstrap support"
  ) |>
  mutate(
    full_model = str_c(
      value_type,
      max_missing,
      model,
      sep = "_"
    )
  ) |>
  group_by(full_model) |>
  summarise(mean = mean(value)) |>
  slice_max(mean, n = 1) |>
  pull(full_model)

highest_mean_consensustree_bootstrap_support = all_values |>
  filter(
    value_type == "consensusTree bootstrap support"
  ) |>
  mutate(
    full_model = str_c(
      value_type,
      max_missing,
      model,
      sep = "_"
    )
  ) |>
  group_by(full_model) |>
  summarise(mean = mean(value)) |>
  slice_max(mean, n = 1) |>
  pull(full_model)

lowest_mean_mltrees_cluster_info_dist = all_values |>
  filter(
    value_type == "mlTrees cluster info dist"
  ) |>
  mutate(
    full_model = str_c(
      value_type,
      max_missing,
      model,
      sep = "_"
    )
  ) |>
  group_by(full_model) |>
  summarise(mean = mean(value)) |>
  slice_min(mean, n = 1) |>
  pull(full_model)


support_across_missingness <- ggplot(
    all_values,
    aes(
      x=model, 
      y=value
    )
  ) + 
  geom_violin(
    draw_quantiles = c(0.25, 0.5, 0.75),
    aes(
      color=factor(
        str_c(
          value_type,
          max_missing,
          model,
          sep="_"
        )
      )
    )
  ) +
  scale_color_manual(
    str_wrap(
      "best model for respective value",
      width = 12
    ),
    values = setNames(
      c("blue", "blue", "blue"),
      c(
        highest_mean_besttree_bootstrap_support,
        highest_mean_consensustree_bootstrap_support,
        lowest_mean_mltrees_cluster_info_dist
      )
    )
  ) +
  stat_summary(
    fun.data="mean_se",
    color="red"
  ) +
  facet_grid(
    cols = vars(factor(max_missing)),
    rows = vars(factor(value_type))
  ) +
theme_bw() +
  theme(
    # text = element_text(size=rel(5.5)),
    # x-axis facet labels do not seem to inherit from text above
    # strip.text.x = element_text(size=rel(5.5)),
    axis.text.x = element_text(
      angle=45,
      vjust=0.9,
      hjust=0.9
    )
  )


number_of_models <- best_trees_support_values |> distinct(model) |> count() |> pull(n)
plot_width = 2 + number_of_models * 5

ggsave(
  filename = snakemake@output[["support_plot"]],
  plot = support_across_missingness,
  width = plot_width,
  height = 4 * 5 + 1,
  limitsize = FALSE
)
