log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library("treeio")
library("tidyverse")

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
      )
    ) |>
    mutate(max_missing = as.integer(max_missing)) |>
    select(support, max_missing, model)) |>
    list_rbind()

support_across_missingness <- ggplot(
    best_trees_support_values,
    aes(
      x=model, 
      y=support
    )
  ) + 
  geom_violin() +
  geom_jitter(width=0.25) +
  stat_summary(color="red") +
  facet_grid(
    cols = vars(factor(max_missing))
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
  limitsize = FALSE
)
