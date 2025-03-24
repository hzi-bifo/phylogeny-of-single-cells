# Changelog

## 1.0.0 (2025-03-24)


### Features

* add collapsed-tree support calculations via gotree rules ([6391a79](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/6391a79b7e7fedbd5ad3757ede0d3bb4a225a56d))
* add further plots to compare across models and missingness ([e8cebab](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/e8cebab5c1c36bb06daeb11f72e4681892dc2b6e))
* also look at consensus tree of all maximum likelihood trees ([33c2dd7](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/33c2dd70ce01a7a4bd0ee07fe2a65220c39ceaea))
* clean up and improve max_missing topologies summary plot ([70cd27d](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/70cd27dacd4444f238fccafac02a3809f48245c9))
* complete cluster_info_dist.R script ([cce6e1d](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/cce6e1da2a9116710d1e8cf504c861b1a03e279d))
* compute and plot cluster information distance for all relevant tree sets ([9917ebb](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/9917ebbb81195cd3a8b5fa9d98c1da44bcefd7c4))
* generate and collect bootstrap convergence and rfdist metrics ([ff634fb](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/ff634fb892285962236405f7ebbd150ced066cdf))
* implement cluster selection and qc for more focused consensus tree ([2bf3be7](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/2bf3be710febf00fb98e0e747028bbdc595a1d6b))
* implement consensus tree building from CID cluster among mlTrees with highest average likelihood ([c3eff03](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/c3eff03c509ad54be07a0bb3fe838656c09ef851))
* improve error messages backtraces in cluster_info_dist.R ([5497d31](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/5497d31dbbb0aeb97c40ce2474e58e69c716a086))
* overhaul summary plots (add information criteria, add ecdf and histogram plots) ([7f111f6](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/7f111f6c773470c8b6a6421b9d7c333168a5e2f5))
* switch axes and facet labels in large overview plots ([89183f7](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/89183f72f6a4b3a651c610bc5992aaf79c8d4f49))


### Bug Fixes

* adjust file name downstream of previous change ([4ef281d](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/4ef281dd04973e7d45f24c24047db8ade91c0626))
* adjust file name extraction pattern ([7cf439e](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/7cf439e77cd8343dedac86b853548e54b49b2362))
* allow big plots ([e79a8ce](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/e79a8ce27994ad93ddd4cc7e744a1a2847869839))
* always compute support based on raxml-ng bootstraps with gotree ([8b08105](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/8b081059bac1192c23f86a3573298e4ed49bdcb8))
* amend previous commit ([441d5c5](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/441d5c5522f1d54e31b89369bda4eb4b7d243da6))
* better error backtrace ([7ebc0bb](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/7ebc0bb3b0149959fce8a7f89fbf6ecaa7b8c0c6))
* file path and resources adjustments ([04e63ef](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/04e63efa196452607b939ea593c4056bfd18aa54))
* hard-coded output file names will give you missing output exceptions ([c200c06](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/c200c065c0e8dc2dec8bb065740dd798c42e76c6))
* highlight model setup with best respective score(s) ([a31e600](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/a31e600ffafa25f053aab9b7510f7a8d65287273))
* include new (and required) r-protoclust in treedist script dependencies ([6ea1cc6](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/6ea1cc6aeaced5ba491c6d82ae01f9a96f6e3961))
* load protoclust in cluster_info_dist.R ([c606784](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/c6067847af0acc2dc2bba8d2cdb95cfcfbf422b4))
* revert to working plot layout ([09aaf43](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/09aaf4378a8e11ad00eea8fb5e665880240aded1))
* update names as needed ([27e4609](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/27e46094294c86612bd48cbc318e999b371579e4))
* use extended regular expressions for sed extraction ([3a72bf7](https://www.github.com/hzi-bifo/phylogeny-of-single-cells/commit/3a72bf7316bdf93dd05c40618f52a13b53e652e2))
