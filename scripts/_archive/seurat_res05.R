#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
dir <- "data/seurat_subset_100k"; outdir <- file.path(dir,"seurat_out")
obj <- readRDS(file.path(outdir,"obj.rds"))

obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
cat("n clusters @ res 0.5:", nlevels(obj$seurat_clusters), "\n")

ggsave(file.path(outdir,"umap_clusters_res05.png"),
       DimPlot(obj, label=TRUE)+NoLegend(), width=7,height=6,dpi=120)

tab <- table(cluster=obj$seurat_clusters, class=obj$CellClass_cal)
write.csv(as.data.frame.matrix(tab), file.path(outdir,"cluster_vs_CellClass_res05.csv"))
purity <- apply(tab,1,function(r) max(r)/sum(r))
domclass <- apply(tab,1,function(r) names(which.max(r)))
gsm_pc <- tapply(obj$gsm, obj$seurat_clusters, function(x) length(unique(x)))
cat("Per-cluster dominant CellClass + purity + #GSMs:\n")
for(c in rownames(tab)) cat(sprintf("  cl%2s n=%5d  %-14s %4.0f%%   GSMs=%d\n",
    c, sum(tab[c,]), domclass[c], 100*purity[c], gsm_pc[c]))
cat(sprintf("\nWeighted mean cluster purity: %.1f%%\n",
    100*sum(purity*rowSums(tab))/sum(tab)))
cat(sprintf("GSMs/cluster min/median/max: %d / %g / %d\n",
    min(gsm_pc), median(gsm_pc), max(gsm_pc)))
saveRDS(obj, file.path(outdir,"obj.rds"))
