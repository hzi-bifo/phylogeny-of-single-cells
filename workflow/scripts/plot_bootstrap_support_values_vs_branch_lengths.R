log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library("treeio")
library("tidyverse")
library("hexbin")

best_trees_support_vs_branch_length <- snakemake@input[["support_trees"]] |>
  map(
    \(filename)
    read.newick(
      filename,
      node.label = "support"
    ) |>
    as_tibble() |>
    drop_na(
      branch.length,
      support
    ) |>
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
      )
    ) |>
    mutate(max_missing = as.integer(max_missing)) |>
    select(
      branch.length,
      support,
      max_missing,
      model
    )
  ) |>
  list_rbind()

number_of_models <- best_trees_support_vs_branch_length |> distinct(model) |> count() |> pull(n)
plot_height = number_of_models * 3
plot_width = 2 + plot_height

data_plot <- ggplot(
    data = best_trees_support_vs_branch_length,
    mapping = aes(
      x=branch.length,
      y=support
    )
  ) +
  geom_hex(
    bins = 10
  ) + 
  scale_fill_viridis_c() +
  geom_point(
    color = "orange",
    alpha = .6
  ) +
  facet_grid(
    cols=vars(model),
    rows=vars(max_missing)
  ) +
  theme_bw() +
  theme(
    text = element_text(size=rel(5.5)),
    # x-axis facet labels do not seem to inherit from text above
    strip.text.x = element_text(size=rel(5.5)),
    axis.text.x = element_text(
      angle=45,
      vjust=0.9,
      hjust=0.9
    )
  )


ggsave(
  filename = snakemake@output[["data_plot"]],
  plot = data_plot,
  width = plot_width,
  height = plot_height
)


best_trees_2d_summary <- best_trees_support_vs_branch_length |>
  group_by(model, max_missing) |>
  summarise(
    branch_length_summary = mean_se(branch.length),
    support_summary = mean_se(support),
    .groups = "drop"
  ) |>
  unnest(
    cols = c(
      branch_length_summary,
      support_summary
    ),
    names_sep="_"
  )

summary_plot <- ggplot(
   data = best_trees_2d_summary,
    mapping = aes(
      x=branch_length_summary_y,
      y=support_summary_y
    )
  ) +
  geom_point(
     color="red"
  ) +
  geom_errorbarh(
    mapping = aes(
      xmin=branch_length_summary_ymin,
      xmax=branch_length_summary_ymax
    ),
    color="red"
  ) +
  geom_errorbar(
    mapping = aes(
      ymin=support_summary_ymin,
      ymax=support_summary_ymax
    ),
    color="red"
  ) +
  facet_grid(
    cols=vars(model),
    rows=vars(max_missing)
  ) +
  theme_bw() +
  theme(
    text = element_text(size=rel(5.5)),
    # x-axis facet labels do not seem to inherit from text above
    strip.text.x = element_text(size=rel(5.5)),
    axis.text.x = element_text(
      angle=45,
      vjust=0.9,
      hjust=0.9
    )
  )

ggsave(
  filename = snakemake@output[["summary_plot"]],
  plot = summary_plot,
  width = plot_width,
  height = plot_height
)
