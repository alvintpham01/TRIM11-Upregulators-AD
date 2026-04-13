# Author: Alvin Pham
# Project: TRIM11 Upregulator Analysis
# Description: General Seurat preprocessing and clustering workflow for GEO sc/snRNA-seq data

suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(dplyr)
})

# =========================================================
# 1. User-defined inputs
# =========================================================

counts_path <- "data/GSE173731_counts.rds"
metadata_path <- "data/GSE173731_metadata.rds"

output_dir <- "results/GSE173731"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_name <- "GSE173731"

min_features <- 300
max_features <- 10000
min_counts <- 500
max_counts <- 1e5
max_percent_mt <- 10

use_sctransform <- TRUE
n_pcs <- 10
cluster_resolution <- 0.2

target_gene <- "TRIM11"

# Optional: set to NULL if you do not want to rename clusters
new_cluster_names <- NULL
# Example:
# new_cluster_names <- c(
#   "0" = "Cluster A",
#   "1" = "Cluster B"
# )

# =========================================================
# 2. Load data
# =========================================================

count_data <- readRDS(counts_path)
meta_data <- readRDS(metadata_path)

seurat_obj <- CreateSeuratObject(
  counts = count_data,
  meta.data = meta_data,
  project = sample_name
)

# =========================================================
# 3. Quality control
# =========================================================

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

qc_violin <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

qc_scatter_1 <- FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

qc_scatter_2 <- FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

qc_scatter_combined <- qc_scatter_1 + qc_scatter_2

ggsave(
  filename = file.path(output_dir, paste0(sample_name, "_qc_violin.png")),
  plot = qc_violin,
  width = 10,
  height = 4,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, paste0(sample_name, "_qc_scatter.png")),
  plot = qc_scatter_combined,
  width = 10,
  height = 4,
  dpi = 300
)

seurat_obj <- subset(
  seurat_obj,
  subset =
    nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    nCount_RNA > min_counts &
    nCount_RNA < max_counts &
    percent.mt < max_percent_mt
)

# =========================================================
# 4. Normalization and scaling
# =========================================================

options(future.globals.maxSize = 10 * 1024^3)

if (use_sctransform) {
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
} else {
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
}

# Extra scaling after SCTransform is usually not necessary,
# so only do it if you specifically need it.
if (!use_sctransform) {
  all_genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all_genes, verbose = FALSE)
}

# =========================================================
# 5. Dimensional reduction and clustering
# =========================================================

if (use_sctransform) {
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
} else {
  seurat_obj <- RunPCA(
    seurat_obj,
    features = VariableFeatures(seurat_obj),
    verbose = FALSE
  )
}

pca_plot <- DimPlot(seurat_obj, reduction = "pca") + NoLegend()
elbow_plot <- ElbowPlot(seurat_obj)

ggsave(
  filename = file.path(output_dir, paste0(sample_name, "_pca.png")),
  plot = pca_plot,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, paste0(sample_name, "_elbow.png")),
  plot = elbow_plot,
  width = 6,
  height = 5,
  dpi = 300
)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = cluster_resolution, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs, verbose = FALSE)

umap_clusters <- DimPlot(seurat_obj, label = TRUE, repel = TRUE)

ggsave(
  filename = file.path(output_dir, paste0(sample_name, "_umap_clusters.png")),
  plot = umap_clusters,
  width = 8,
  height = 6,
  dpi = 300
)

# =========================================================
# 6. Marker detection
# =========================================================

seurat_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)

marker_hits <- seurat_markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1)

top10_markers <- seurat_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(
  seurat_markers,
  file = file.path(output_dir, paste0(sample_name, "_all_markers.csv")),
  row.names = FALSE
)

write.csv(
  marker_hits,
  file = file.path(output_dir, paste0(sample_name, "_marker_hits_log2fc_gt1.csv")),
  row.names = FALSE
)

write.csv(
  top10_markers,
  file = file.path(output_dir, paste0(sample_name, "_top10_markers.csv")),
  row.names = FALSE
)

# =========================================================
# 7. Optional cluster renaming
# =========================================================

if (!is.null(new_cluster_names)) {
  seurat_obj <- RenameIdents(seurat_obj, new_cluster_names)
  seurat_obj$celltype <- Idents(seurat_obj)

  umap_annotated <- DimPlot(seurat_obj, label = TRUE, repel = TRUE)

  ggsave(
    filename = file.path(output_dir, paste0(sample_name, "_umap_annotated.png")),
    plot = umap_annotated,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# =========================================================
# 8. Target gene visualization
# =========================================================

if (target_gene %in% rownames(seurat_obj)) {
  trim11_plot <- FeaturePlot(seurat_obj, features = target_gene)

  ggsave(
    filename = file.path(output_dir, paste0(sample_name, "_", target_gene, "_featureplot.png")),
    plot = trim11_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
} else {
  message(target_gene, " not found in the Seurat object.")
}

# =========================================================
# 9. Save object
# =========================================================

saveRDS(
  seurat_obj,
  file = file.path(output_dir, paste0(sample_name, "_seurat_obj.rds"))
)

writeLines(capture.output(sessionInfo()),
           con = file.path(output_dir, paste0(sample_name, "_sessionInfo.txt")))
