#!/usr/bin/env Rscript
# Clustering module (Seurat-backed) for omnibenchmark: Seurat::FindClusters on the
# connectivities graph of a *_neighbors.h5 file (UMAP-style affinities from the scanpy
# and rapids modules, the SNN graph from this module's knn.R).
#
# Phases (obkit-events.jsonl):
#   load     h5 -> sparse graph -> Seurat object
#   cluster  FindClusters (--flavor: louvain_original, Seurat's default;
#            louvain_multilevel_refinement; slm; leiden)
#   write    TSV serialization
#
# Output: {name}_clusters.tsv (cell_id, cluster).

suppressPackageStartupMessages({
  library(HDF5Array)
  library(rhdf5)
  library(Seurat)
  library(Matrix)
  library(jsonlite)
})

source("src/read_neighbors.R")  # read_neighbors(): a CSR graph from *_neighbors.h5
source("src/flavor.R")          # flavor_algorithm_id(): --flavor -> FindClusters algorithm id
source("src/common/cli.R")
source("src/obkit_logger.R")
source("src/phases.R")

parse_cluster_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  p <- arg_parser("CLUST module")
  p <- add_base_args(p)            # --output_dir, --name
  p <- add_stage_args(p, "CLUST")  # --neighbors_h5
  p <- add_argument(p, "--flavor", type = "character",
                    help = "Clustering algorithm (louvain_original|louvain_multilevel_refinement|slm|leiden)")
  p <- add_argument(p, "--resolution", type = "numeric", help = "Clustering resolution")
  p <- add_argument(p, "--random_seed", type = "integer", help = "Random seed")
  parse_args(p, argv)
}

main <- function() {
  args <- parse_cluster_args()
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in c("output_dir", "name", "neighbors_h5", "flavor", "resolution", "random_seed"))
    cat(sprintf("  %s: %s\n", k, format(args[[k]])))
  algorithm <- flavor_algorithm_id(args$flavor)  # errors on a missing or unknown flavour

  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
  logger_init(args$output_dir)
  set.seed(args$random_seed)

  so <- phase("load", function(attrs) {
    g <- as.Graph(read_neighbors(args$neighbors_h5, group = "connectivities"))
    so <- CreateSeuratObject(counts = g, assay = "RNA")
    so@graphs$neighbors <- g
    attrs$n_cells <- ncol(g)
    attrs$connectivities_nnz <- length(g@x)
    so
  })

  so <- phase("cluster", function(attrs) {
    attrs$flavor <- args$flavor
    attrs$algorithm <- algorithm
    attrs$resolution <- args$resolution
    so <- FindClusters(so, algorithm = algorithm, resolution = args$resolution, graph.name = "neighbors",
                       random.seed = args$random_seed, verbose = FALSE)
    attrs$n_clusters <- length(unique(so$seurat_clusters))
    so
  })

  phase("write", function(attrs) {
    out <- file.path(args$output_dir, paste0(args$name, "_clusters.tsv"))
    write.table(data.frame(cell_id = colnames(so), cluster = as.character(so$seurat_clusters)),
                file = out, sep = "\t", quote = FALSE, row.names = FALSE)
    attrs$path <- out
    cat(sprintf("  wrote: %s\n", out))
  })
}

if (sys.nframe() == 0L) {
  main()
}
