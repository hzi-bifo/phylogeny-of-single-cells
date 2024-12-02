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

usable_sites <- read_tsv(
    snakemake@input[["tsv"]]
  ) |>
  select(model, max_missing, usable_sites) |>
  distinct()


number_of_models <- best_trees_support_vs_branch_length |> distinct(model) |> count() |> pull(n)
plot_height = number_of_models * 3
plot_width = 2 + plot_height

non_na_branch_lengths <- best_trees_support_vs_branch_length |>
  drop_na(branch.length)

branch_length_hist <- ggplot(
    data = non_na_branch_lengths
  ) +
  geom_histogram(
    aes(x = branch.length)
  ) +
  facet_grid(
    cols=vars(model),
    rows=vars(max_missing)
  )

ggsave(
  filename = snakemake@output[["branch_length_hist"]],
  plot = branch_length_hist,
  width = plot_width,
  height = plot_height * 0.75,
  limitsize = FALSE
)

normalized_branch_lengths <- non_na_branch_lengths |>
  left_join(
    usable_sites,
    by = join_by(
      model,
      max_missing
    )
  )

branch_length_ecdf <- ggplot(
    data = normalized_branch_lengths,
    aes(
      x = branch.length
    )
  ) +
  stat_ecdf(geom = "point") +
  facet_grid(
    cols=vars(model),
    rows=vars(max_missing)
  )

ggsave(
  filename = snakemake@output[["branch_length_ecdf"]],
  plot = branch_length_ecdf,
  width = plot_width,
  height = plot_height * 0.75,
  limitsize = FALSE
)



support_hist <- ggplot(
    data = best_trees_support_vs_branch_length |>
      drop_na(support)
  ) +
  geom_histogram(
    aes(x = support)
  ) +
  facet_grid(
    cols=vars(model),
    rows=vars(max_missing)
  )

ggsave(
  filename = snakemake@output[["support_hist"]],
  plot = support_hist,
  width = plot_width,
  height = plot_height * 0.75,
  limitsize = FALSE
)


data_plot <- ggplot(
    data = best_trees_support_vs_branch_length |>
      drop_na(
        branch.length,
        support
      ),
    mapping = aes(
      x=support,
      y=branch.length
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
  height = plot_height,
  limitsize = FALSE
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
      x=support_summary_y,
      y=branch_length_summary_y
    )
  ) +
  geom_point(
     color="red"
  ) +
  geom_errorbar(
    mapping = aes(
      ymin=branch_length_summary_ymin,
      ymax=branch_length_summary_ymax
    ),
    color="red"
  ) +
  geom_errorbarh(
    mapping = aes(
      xmin=support_summary_ymin,
      xmax=support_summary_ymax
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
  height = plot_height,
  limitsize = FALSE
)
