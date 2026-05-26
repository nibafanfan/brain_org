#!/usr/bin/env Rscript
# Test cell-type segmentation on the 100k subset, clustering on the reused scVI latent.
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
set.seed(0)

dir <- "data/seurat_subset_100k"
outdir <- file.path(dir, "seurat_out")
dir.create(outdir, showWarnings = FALSE)

cat("[1] Read10X + metadata...\n")
mat  <- Read10X(dir)
obj  <- CreateSeuratObject(mat)
meta <- read.csv(file.path(dir, "metadata.csv"), row.names = 1)
obj  <- AddMetaData(obj, meta[colnames(obj), ])
cat("    cells:", ncol(obj), " genes:", nrow(obj), "\n")

cat("[2] load scVI latent as DimReduc...\n")
lat <- read.csv(file.path(dir, "scvi_latent.csv"), row.names = 1)
lat <- as.matrix(lat[colnames(obj), ])
colnames(lat) <- paste0("scvi_", seq_len(ncol(lat)))
obj[["scvi"]] <- CreateDimReducObject(embeddings = lat, key = "scvi_", assay = "RNA")

cat("[3] normalize (marker viz), cluster on scVI, UMAP...\n")
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "scvi", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 1.0, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "scvi", dims = 1:30, verbose = FALSE)

cat("[4] save plots...\n")
ggsave(file.path(outdir, "umap_clusters.png"),
       DimPlot(obj, label = TRUE) + NoLegend(), width = 7, height = 6, dpi = 120)
ggsave(file.path(outdir, "umap_celltype_origin.png"),
       DimPlot(obj, group.by = "cell_type_origin"), width = 9, height = 6, dpi = 120)
ggsave(file.path(outdir, "umap_gsm.png"),
       DimPlot(obj, group.by = "gsm") + NoLegend(), width = 7, height = 6, dpi = 120)
ggsave(file.path(outdir, "umap_protocol.png"),
       DimPlot(obj, group.by = "protocol") + NoLegend(), width = 7, height = 6, dpi = 120)

cat("[5] cluster x cell_type_origin table...\n")
tab <- table(cluster = obj$seurat_clusters, origin = obj$cell_type_origin)
write.csv(as.data.frame.matrix(tab), file.path(outdir, "cluster_vs_origin.csv"))

# batch-mixing sanity: how many GSMs contribute to each cluster (higher = better mixing)
gsm_per_cluster <- tapply(obj$gsm, obj$seurat_clusters, function(x) length(unique(x)))
cat("    n clusters:", nlevels(obj$seurat_clusters), "\n")
cat("    GSMs per cluster (min/median/max):",
    min(gsm_per_cluster), median(gsm_per_cluster), max(gsm_per_cluster), "\n")

saveRDS(obj, file.path(outdir, "obj.rds"))
cat("DONE ->", outdir, "\n")
print(table(obj$seurat_clusters))
