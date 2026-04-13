# Author: Alvin Pham
# Project: TRIM11 Upregulator Analysis
# Description: Pearson correlation analysis within a specified Seurat cluster

suppressPackageStartupMessages({
  library(Seurat)
})

# =========================================================
# 1. Helper function: sanitize names for files
# =========================================================

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

# =========================================================
# 2. Correlate target gene with all genes in one cluster
# =========================================================

correlate_gene_within_cluster <- function(
    seurat_object,
    cluster_name,
    target_gene,
    output_dir = NULL,
    assay_name = "SCT",
    data_slot = "data",
    save_cluster_rds = FALSE
) {
  # Subset the selected cluster
  cluster_subset <- subset(seurat_object, idents = cluster_name)
  
  # Optionally save cluster object
  safe_cluster_name <- sanitize_name(cluster_name)
  safe_target_gene <- sanitize_name(target_gene)
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (save_cluster_rds && !is.null(output_dir)) {
    rds_path <- file.path(output_dir, paste0("cluster_", safe_cluster_name, ".rds"))
    saveRDS(cluster_subset, file = rds_path)
  }
  
  # Extract normalized expression matrix
  expr_mat <- GetAssayData(cluster_subset, assay = assay_name, slot = data_slot)
  
  # Check that target gene is present
  if (!target_gene %in% rownames(expr_mat)) {
    stop(target_gene, " was not found in assay ", assay_name, " slot ", data_slot, ".")
  }
  
  target_expression <- as.numeric(expr_mat[target_gene, ])
  genes <- rownames(expr_mat)
  
  correlations <- rep(NA_real_, length(genes))
  p_values <- rep(NA_real_, length(genes))
  
  for (i in seq_along(genes)) {
    gene <- genes[i]
    
    if (gene == target_gene) {
      next
    }
    
    gene_expression <- as.numeric(expr_mat[gene, ])
    
    test_result <- suppressWarnings(
      cor.test(target_expression, gene_expression, method = "pearson")
    )
    
    correlations[i] <- unname(test_result$estimate)
    p_values[i] <- test_result$p.value
  }
  
  correlation_results <- data.frame(
    gene = genes,
    correlation = correlations,
    p_value = p_values,
    stringsAsFactors = FALSE
  )
  
  correlation_results <- subset(correlation_results, gene != target_gene)
  correlation_results <- correlation_results[order(correlation_results$p_value), ]
  
  if (!is.null(output_dir)) {
    output_csv <- file.path(
      output_dir,
      paste0(safe_target_gene, "_pearson_correlation_", safe_cluster_name, ".csv")
    )
    write.csv(correlation_results, output_csv, row.names = FALSE)
  }
  
  return(correlation_results)
}

# =========================================================
# 3. Extract significant genes after multiple-testing correction
# =========================================================

extract_significant_genes <- function(
    correlation_df,
    seurat_object,
    target_gene,
    cluster_name,
    output_dir = NULL,
    assay_name = "SCT",
    data_slot = "data",
    adj_pval_threshold = 0.1,
    adjustment_method = "BH"
) {
  valid <- correlation_df[!is.na(correlation_df$p_value), , drop = FALSE]
  
  valid$p_adj <- p.adjust(valid$p_value, method = adjustment_method)
  
  sig_genes <- valid$gene[valid$p_adj < adj_pval_threshold]
  
  expr_mat <- GetAssayData(seurat_object, assay = assay_name, slot = data_slot)
  
  if (!target_gene %in% rownames(expr_mat)) {
    stop(target_gene, " was not found in assay ", assay_name, " slot ", data_slot, ".")
  }
  
  target_expr <- expr_mat[target_gene, , drop = FALSE]
  
  if (length(sig_genes) > 0) {
    sig_expr <- expr_mat[sig_genes, , drop = FALSE]
    combined_expr <- rbind(target_expr, sig_expr)
  } else {
    combined_expr <- target_expr
  }
  
  rownames(combined_expr) <- make.unique(rownames(combined_expr))
  
  safe_cluster_name <- sanitize_name(cluster_name)
  safe_target_gene <- sanitize_name(target_gene)
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    output_csv <- file.path(
      output_dir,
      paste0("post_pearson_", safe_target_gene, "_", safe_cluster_name, ".csv")
    )
    
    write.csv(as.data.frame(combined_expr), output_csv)
  }
  
  return(valid[valid$p_adj < adj_pval_threshold, , drop = FALSE])
}

# =========================================================
# 4. Example usage
# =========================================================

# seurat_object <- readRDS("results/GSE173731/GSE173731_seurat_obj.rds")
# target_gene <- "TRIM11"
# cluster_name <- "Immune-like / myeloid receptor-high (CD200R1L+; unclear subtype)"
#
# correlation_output_dir <- "results/GSE173731/pearson_results"
# postpearson_output_dir <- "results/GSE173731/postpearson_gene_counts"
#
# correlation_results <- correlate_gene_within_cluster(
#   seurat_object = seurat_object,
#   cluster_name = cluster_name,
#   target_gene = target_gene,
#   output_dir = correlation_output_dir,
#   assay_name = "SCT",
#   data_slot = "data",
#   save_cluster_rds = FALSE
# )
#
# sig_results <- extract_significant_genes(
#   correlation_df = correlation_results,
#   seurat_object = subset(seurat_object, idents = cluster_name),
#   target_gene = target_gene,
#   cluster_name = cluster_name,
#   output_dir = postpearson_output_dir,
#   assay_name = "SCT",
#   data_slot = "data",
#   adj_pval_threshold = 0.1,
#   adjustment_method = "BH"
# )
