#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
dir <- "data/seurat_subset_100k"; outdir <- file.path(dir,"seurat_out")
obj <- readRDS(file.path(outdir,"obj.rds"))
lab <- read.csv(file.path(dir,"braun_labels.csv"), row.names=1)
obj <- AddMetaData(obj, lab[colnames(obj),])

ggsave(file.path(outdir,"umap_CellClass.png"),
       DimPlot(obj, group.by="CellClass_cal", label=TRUE, repel=TRUE)+NoLegend(),
       width=7,height=6,dpi=120)
ggsave(file.path(outdir,"umap_Region.png"),
       DimPlot(obj, group.by="Region"), width=9,height=6,dpi=120)

# purity of each Seurat cluster wrt calibrated CellClass
tab <- table(cluster=obj$seurat_clusters, class=obj$CellClass_cal)
write.csv(as.data.frame.matrix(tab), file.path(outdir,"cluster_vs_CellClass.csv"))
purity <- apply(tab,1,function(r) max(r)/sum(r))
domclass <- apply(tab,1,function(r) names(which.max(r)))
cat("Per-cluster dominant CellClass + purity:\n")
for(c in rownames(tab)) cat(sprintf("  cl%2s n=%5d  %-14s %4.0f%%\n",
    c, sum(tab[c,]), domclass[c], 100*purity[c]))
cat(sprintf("\nMean cluster purity (weighted): %.1f%%\n",
    100*sum(purity*rowSums(tab))/sum(tab)))

# overall agreement: assign each cluster its majority label, count matches
maj <- domclass[as.character(obj$seurat_clusters)]
cat(sprintf("Cells whose cluster-majority == their own label: %.1f%%\n",
    100*mean(maj==obj$CellClass_cal)))
saveRDS(obj, file.path(outdir,"obj.rds"))
