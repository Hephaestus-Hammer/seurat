#!/usr/bin/env Rscript
# PCA module (Seurat-backed) for omnibenchmark.
#
# Solvers:
#   approximate — Seurat::RunPCA(approx = TRUE)  [IRLBA]
#   exact       — Seurat::RunPCA(approx = FALSE) [full SVD]
#
# --scale false (default): ScaleData(do.scale = FALSE), mean-centre only, to match
#   scanpy/scrapper/rapids, which don't divide by per-gene std.
# --scale true: Seurat's own default ScaleData (centre, divide by the standard
#   deviation, clip at --scale_max, default 10), i.e. what `RunPCA` sees in a default
#   Seurat workflow.
#
# Phases (obkit-events.jsonl):
#   load      h5 -> dgCMatrix -> Seurat object (everything that gets the input
#             into the framework's expected in-memory form).
#   densify   ScaleData — allocates a dense (n_genes × n_cells) double matrix in the
#             `scale.data` layer because that's what RunPCA reads from. This is the
#             memory elephant for Seurat: ~n_genes * n_cells * 8 bytes.
#   pca       RunPCA (IRLBA or full SVD).
#   write     TSV serialization.
#
# Outputs: {name}_embedding.tsv (cell_id, PC1..PCn) and {name}_loadings.tsv
# (gene_id, PC1..PCn).

suppressPackageStartupMessages({
  library(Matrix)
  library(HDF5Array)
  library(Seurat)
  library(data.table)
  library(jsonlite)  # logger emits via toJSON; satisfied transitively by Seurat
})

source("src/common/cli.R")
source("src/obkit_logger.R")
source("src/phases.R")
source("src/args.R")

parse_pca_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  p <- arg_parser("PCA module")
  p <- add_base_args(p)          # --output_dir, --name
  p <- add_stage_args(p, "PCA")  # --normalized_selected_h5
  p <- add_argument(p, "--solver", type = "character", default = "approximate",
                    help = "approximate (IRLBA) or exact (full SVD)")
  p <- add_argument(p, "--n_components", type = "integer", default = 50L, help = "number of PCs")
  p <- add_argument(p, "--random_seed", type = "integer", default = 42L, help = "random seed")
  p <- add_argument(p, "--scale", type = "character", default = "false",
                    help = "true: Seurat's default ScaleData (unit variance, clipped); false: centre only")
  p <- add_argument(p, "--scale_max", type = "numeric", default = 10, help = "clip for scaled values (--scale true)")
  args <- parse_args_outputs(p, argv)  # argparser, minus --output id=PATH (api 0.8)
  if (!(args$solver %in% c("exact", "approximate")))
    stop("Invalid --solver: ", args$solver, " (valid: exact, approximate)")
  args$scale <- as_flag(args$scale, "--scale")
  args
}

run_pca <- function(so, X, args) {
  # so: Seurat object with `counts` + `data` layers populated from X.
  # X:  the same gene-by-cell sparse dgCMatrix (for gene/cell ids + PCA features).
  set.seed(args$random_seed)
  approx <- args$solver == "approximate"

  so <- phase("densify", function(attrs) {
    attrs$do_scale        <- args$scale
    attrs$dense_bytes_est <- as.numeric(nrow(X)) * as.numeric(ncol(X)) * 8
    ScaleData(so, features = rownames(X), do.scale = args$scale, do.center = TRUE,
              scale.max = args$scale_max, verbose = FALSE)
  })

  so <- phase("pca", function(attrs) {
    attrs$solver       <- args$solver
    attrs$approx       <- approx
    attrs$n_components <- args$n_components
    RunPCA(so,
           features = rownames(X),
           npcs = args$n_components,
           approx = approx,
           seed.use = args$random_seed,
           verbose = FALSE)
  })

  dimred    <- so[["pca"]]
  embedding <- Embeddings(dimred)           # (n_cells, n_components)
  loadings  <- Loadings(dimred)             # (n_genes, n_components)
  rownames(embedding) <- colnames(X)
  colnames(embedding) <- paste0("PC", seq_len(ncol(embedding)))
  rownames(loadings)  <- rownames(X)
  colnames(loadings)  <- paste0("PC", seq_len(ncol(loadings)))
  list(embedding = embedding, loadings = loadings)
}

main <- function() {
  args <- parse_pca_args()
  cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
  for (k in c("output_dir", "name", "normalized_selected_h5", "solver", "n_components", "random_seed", "scale"))
    cat(sprintf("  %s: %s\n", k, format(args[[k]])))

  out_emb  <- output_path(args, "embedding_tsv", "_embedding.tsv")  # --output id=PATH (api 0.8)
  out_load <- output_path(args, "loadings_tsv", "_loadings.tsv")
  dir.create(dirname(out_emb), showWarnings = FALSE, recursive = TRUE)
  logger_init(dirname(out_emb))

  loaded <- phase("load", function(attrs) {
    # TENxMatrix is lazy/disk-backed; as(., dgCMatrix) materializes it in memory.
    # CreateSeuratObject + SetAssayData wrap the sparse matrix in Seurat's container.
    m_lazy <- TENxMatrix(args$normalized_selected_h5, group = "matrix")
    m_mem  <- as(m_lazy, "dgCMatrix")
    so     <- CreateSeuratObject(counts = m_mem)
    so     <- SetAssayData(so, layer = "data", new.data = m_mem)
    attrs$n_genes <- nrow(m_mem)
    attrs$n_cells <- ncol(m_mem)
    attrs$nnz     <- as.integer(length(m_mem@x))
    list(X = m_mem, so = so)
  })
  cat(sprintf("  matrix (genes x cells): %d x %d\n", nrow(loaded$X), ncol(loaded$X)))

  res <- run_pca(loaded$so, loaded$X, args)
  cat(sprintf("  embedding: %d x %d, loadings: %d x %d\n",
              nrow(res$embedding), ncol(res$embedding), nrow(res$loadings), ncol(res$loadings)))

  phase("write", function(attrs) {
    fwrite(data.table(cell_id = rownames(res$embedding), res$embedding), out_emb, sep = "\t", quote = FALSE)
    fwrite(data.table(gene_id = rownames(res$loadings), res$loadings), out_load, sep = "\t", quote = FALSE)
    attrs$embedding <- out_emb
    attrs$loadings  <- out_load
    cat(sprintf("  wrote: %s\n  wrote: %s\n", out_emb, out_load))
  })
}

if (sys.nframe() == 0L) {
  main()
}
