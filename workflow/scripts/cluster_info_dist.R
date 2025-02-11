log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

rlang::global_entrace()

library("TreeDist")
library("tidyverse")
library("cluster")
library("protoclust")

# cluster information distance of all input trees

trees <- ape::read.tree(snakemake@input[["trees"]])

all_trees_cid <- ClusteringInfoDistance(trees, trees, normalize=TRUE)

all_trees_cid_dedup <- all_trees_cid[lower.tri(all_trees_cid)]

all_cid_annotated <- all_trees_cid_dedup |>
  as_tibble_col(column_name = "cluster_information_distance") |>
  add_column(
    model=snakemake@wildcards[["model"]],
    max_missing=snakemake@wildcards[["n_missing_cells"]],
    tree_type=snakemake@wildcards[["type"]],
    software=snakemake@wildcards[["software"]]
  )

write_tsv(
  x = all_cid_annotated,
  file = snakemake@output[["all_cid"]]
)

# select best clustering

number_of_clusters <- as_tibble_col(2:20, column_name = "number_of_clusters")

hierarchical_clustering_tree <- protoclust(distances)

clusterings <- number_of_clusters |> 
  mutate(
    clustering_partitioning_around_medoids = map(
      number_of_clusters,
      \(x) cluster::pam(!!distances, x)
    ),
    clustering_hierarchical = map(
      number_of_clusters,
      \(x) protocut(!!hierarchical_clustering_tree, k = x)
    ),
    clustering_k_means = map(
      number_of_clusters, 
      \(x) KMeansPP(!!distances, k = x)
    )
  ) |>
  pivot_longer(
    starts_with("clustering_"),
    names_prefix = "clustering_",
    names_to = "clustering_method",
    values_to = "clustering_object"
  ) |>
  mutate(
    clustering = map2(
      clustering_method,
      clustering_object,
      \(m, o)
        case_match(
          m,
          c("partitioning_around_medoids", "k_means")  ~ list(o$cluster),
          "hierarchical" ~ list(o$cl)
        )
      )
  ) |>
  unnest(clustering) |>
  mutate(
    mean_silhouette = map(
      clustering,
      \(x)
        mean(cluster::silhouette(x, !!distances)[, "sil_width"])
    ) 
  ) |>
  unnest(mean_silhouette)

clustering_silhouette_scores_plot <- ggplot(clusterings) +
  geom_point(
    aes(
      x=number_of_clusters,
      y = mean_silhouette,
      shape=clustering_method,
      color=clustering_method
    )
  )

ggsave(
  snakemake@output[["silhouette_scores"]],
  clustering_silhouette_scores_plot
)

best_clustering <- clusterings |>
  slice_max(mean_silhouette) |>
  pull(clustering) |>
  first()

best_clustering_tidy <- best_clustering |>
  enframe(
    name = "tree",
    value = "cluster"
  ) |>
  mutate(
    cluster = factor(cluster)
  )

# select best cluster from best clustering

tree_log_likelihoods <- left_join(
    best_clustering_tidy,
    read_tsv(snakemake@input[["likelihoods"]], col_names = c("tree", "log_likelihood")),
    by = "tree"
  )

cluster_log_likelihoods_plot <- ggplot(
    tree_log_likelihoods,
    aes(
      x=cluster,
      y=log_likelihood
    )
  ) +
  geom_violin() +
  geom_jitter() +
  theme_bw()

ggsave(
  snakemake@output[["best_clustering_log_likelihoods"]],
  cluster_log_likelihoods_plot,
  width = max(best_clustering) * 1.2,
  height = 8
)

best_clustering_silhouette_scores <- cluster::silhouette(best_clustering, distances)[, "sil_width"] |>
  enframe(name = "tree", value = "sil_width") |>
  left_join(
    tree_log_likelihoods,
    by = "tree"
  )

cluster_silhouette_scores_plot <- ggplot(
    best_clustering_silhouette_scores,
    aes(
      x=cluster,
      y=sil_width
    )
  ) +
  geom_violin() +
  geom_jitter() +
  theme_bw()

ggsave(
  snakemake@output[["best_clustering_silhouette_scores"]],
  cluster_log_likelihoods_plot,
  width = max(best_clustering) * 1.2,
  height = 8
)

pcoa_tidy <- cmdscale(distances, k = 8) |>
  as_tibble() |>
  rownames_to_column("tree") |>
  pivot_longer(
    cols = starts_with("V"),
    names_to = "dimension"
  ) |>
  mutate(
    dimension = str_replace(dimension, "V", ""),
    tree = as.integer(tree)
  )

all_combinations <- full_join(
    pcoa_tidy,
    pcoa_tidy,
    by = "tree",
    relationship="many-to-many"
  )

cluster_stats <- best_clustering_silhouette_scores |>
  group_by(cluster) |>
  summarize(
    n = n(),
    mean_log_likelihood = mean(log_likelihood),
    median_log_likelihood = median(log_likelihood),
    max_log_likelihood = max(log_likelihood),
    mean_silhouette_score = mean(sil_width)
  ) |>
  mutate(
    cluster = factor(cluster)
  )

all_annotated <- all_combinations |>
  left_join(
    tree_log_likelihoods,
    by="tree"
  ) |>
  left_join(
    cluster_stats,
    by="cluster"
  )

total_trees <- cluster_stats |>
  summarise(total = sum(n)) |>
  pull()

mapping_plot_cluster_and_log_likelihood <- ggplot(
    all_annotated |>
      filter(n > .02 * total_trees)
  ) +
  geom_point(
    aes(
      x=value.x,
      y=value.y,
      fill = log_likelihood,
      colour = cluster
    ),
    shape = 21,
    size = 1.5,
    stroke = 1
  ) +
  facet_grid(
    cols=vars(dimension.x),
    rows=vars(dimension.y)
  ) +
  scale_color_hue(aesthetics = "colour") +
  scale_fill_viridis_c(aesthetics = "fill") +
  theme_bw()

ggsave(
  snakemake@output[["mapping_cluster_likelihood"]],
  mapping_plot_cluster_and_log_likelihood,
  width = 8 * 3 + 2,
  height = 8 * 2.5 + 1
)

cluster_likelihoods <- tree_log_likelihoods |>
  group_by(cluster) |>
  summarise(
    n = n(),
    mean = mean(log_likelihood),
    max = max(log_likelihood),
    median = num(median(as.numeric(log_likelihood)), digits=6)
  )

best_cluster <- cluster_likelihoods |>
  # ignore very small clusters, assuming they are not well supported
  filter(n > .02 * total_trees) |>
  slice_max(mean) |>
  pull(cluster)

best_cluster_trees <- trees[
    best_clustering |>
    filter(cluster == best_cluster) |>
    pull(tree)
  ]

ape::write.tree(
  best_cluster_trees,
  snakemake@output[["best_cluster_trees"]]
)

