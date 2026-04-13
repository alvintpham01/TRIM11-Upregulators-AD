# Author: Alvin Pham
# Project: TRIM11 Upregulator Analysis
# Description: Bayesian network inference from post-Pearson gene expression data

suppressPackageStartupMessages({
  library(bnlearn)
  library(pheatmap)
})

# =========================================================
# 1. Helper function
# =========================================================

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

# =========================================================
# 2. Read and format post-Pearson gene count matrix
# =========================================================

prepare_bn_input <- function(input_csv) {
  genecounts <- read.csv(input_csv, check.names = FALSE)
  
  data_mat <- as.matrix(genecounts)
  data_mat <- t(data_mat)
  data_df <- as.data.frame(data_mat, stringsAsFactors = FALSE)
  
  colnames(data_df) <- data_df[1, ]
  data_df <- data_df[-1, , drop = FALSE]
  
  data_df[] <- lapply(data_df, as.numeric)
  
  return(data_df)
}

# =========================================================
# 3. Build Bayesian network and weighted adjacency matrix
# =========================================================

run_bayesian_network <- function(
    data_df,
    score_type = "bic-g"
) {
  bn_structure <- hc(data_df, score = score_type)
  bn_fitted <- bn.fit(bn_structure, data_df)
  
  weighted_adj_matrix <- matrix(
    0,
    nrow = ncol(data_df),
    ncol = ncol(data_df),
    dimnames = list(colnames(data_df), colnames(data_df))
  )
  
  for (node in names(bn_fitted)) {
    if (!is.null(bn_fitted[[node]]$coefficients)) {
      parents <- bn_fitted[[node]]$parents
      coefficients <- bn_fitted[[node]]$coefficients
      
      for (parent in parents) {
        if (parent %in% names(coefficients)) {
          weighted_adj_matrix[parent, node] <- coefficients[parent]
        }
      }
    }
  }
  
  list(
    structure = bn_structure,
    fitted = bn_fitted,
    weighted_adj_matrix = weighted_adj_matrix
  )
}

# =========================================================
# 4. Save outputs
# =========================================================

save_bn_outputs <- function(
    weighted_adj_matrix,
    target_gene,
    cluster_name,
    output_dir
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  safe_cluster_name <- sanitize_name(cluster_name)
  safe_target_gene <- sanitize_name(target_gene)
  
  result_csv <- file.path(
    output_dir,
    paste0("bayesian_weighted_adjacency_", safe_cluster_name, ".csv")
  )
  
  result_heatmap <- file.path(
    output_dir,
    paste0("bayesian_heatmap_", safe_cluster_name, ".pdf")
  )
  
  result_upreg <- file.path(
    output_dir,
    paste0("bayesian_upregulators_", safe_target_gene, "_", safe_cluster_name, ".csv")
  )
  
  write.csv(weighted_adj_matrix, file = result_csv, row.names = TRUE)
  
  pheatmap(
    weighted_adj_matrix,
    main = "Weighted Adjacency Matrix of Bayesian Network",
    display_numbers = FALSE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 8,
    fontsize_col = 8,
    color = colorRampPalette(c("white", "blue"))(50),
    filename = result_heatmap,
    width = 15,
    height = 15
  )
  
  if (!target_gene %in% colnames(weighted_adj_matrix)) {
    stop(target_gene, " not found in weighted adjacency matrix.")
  }
  
  regulators_with_weights <- weighted_adj_matrix[
    weighted_adj_matrix[, target_gene] != 0,
    target_gene,
    drop = FALSE
  ]
  
  write.csv(regulators_with_weights, file = result_upreg, row.names = TRUE)
  
  invisible(list(
    adjacency_csv = result_csv,
    heatmap_pdf = result_heatmap,
    upregulator_csv = result_upreg
  ))
}

# =========================================================
# 5. Wrapper function
# =========================================================

run_bn_pipeline <- function(
    input_csv,
    cluster_name,
    target_gene = "TRIM11",
    output_dir = "results/bayesian",
    score_type = "bic-g"
) {
  data_df <- prepare_bn_input(input_csv)
  
  if (ncol(data_df) < 2) {
    stop("Bayesian network input must contain at least two genes.")
  }
  
  bn_results <- run_bayesian_network(
    data_df = data_df,
    score_type = score_type
  )
  
  save_bn_outputs(
    weighted_adj_matrix = bn_results$weighted_adj_matrix,
    target_gene = target_gene,
    cluster_name = cluster_name,
    output_dir = output_dir
  )
  
  return(bn_results)
}

# =========================================================
# 6. Example usage
# =========================================================

# cluster_name <- "Immune-like / myeloid receptor-high (CD200R1L+; unclear subtype)"
#
# input_csv <- "results/GSE173731/postpearson_gene_counts/post_pearson_TRIM11_Immune_like_myeloid_receptor_high_CD200R1L_unclear_subtype.csv"
#
# bn_results <- run_bn_pipeline(
#   input_csv = input_csv,
#   cluster_name = cluster_name,
#   target_gene = "TRIM11",
#   output_dir = "results/GSE173731/bayesian_results",
#   score_type = "bic-g"
# )
