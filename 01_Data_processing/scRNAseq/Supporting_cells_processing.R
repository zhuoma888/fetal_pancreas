library(Seurat)
library(ggplot2)
library(ggsci)
library(tidyverse)

sup <- readRDS("/home/marzon/Rproject/fetal_panc/W7_11.rds")
sup <- subset(sup, subset = cell_class %in% c("Mesenchymal","Endothelial","Neural","Immune") )
sup <- NormalizeData(sup)
sup <- FindVariableFeatures(sup,nfeatures = 2000)
sup <- ScaleData(sup, vars.to.regress = c("S.Score","G2M.Score","orig.ident"))
sup <- RunPCA(sup)
sup <- RunUMAP(sup, dims = 1:20, reduction = "pca")
sup <- FindNeighbors(sup, reduction = "pca", dims = 1:20)
sup <- FindClusters(sup, resolution = 0.2)


sup@meta.data$cell_type <- "Fibroblast"
sup@meta.data$cell_type[which(sup@meta.data$seurat_clusters %in% c(5,6))] <- "Mesothelial"
sup@meta.data$cell_type[which(sup@meta.data$seurat_clusters == 3)] <- "Pericyte"
sup@meta.data$cell_type[which(sup@meta.data$seurat_clusters == 10)] <- "Immune"
sup@meta.data$cell_type[which(sup@meta.data$seurat_clusters %in% c(8,11))] <- "Neural"
sup@meta.data$cell_type[which(sup@meta.data$seurat_clusters %in% c(2,12))] <- "Endothelial"

Idents(sup) <- sup@meta.data$cell_type

DotPlot(sup, features = c("UPK3B","KRT19","PDGFRA","PDGFRB","CSPG4","ACTA2","PECAM1","ASCL1","PTPRC"), group.by = "cell_type") + 
  theme_classic() + RotatedAxis() + 
  scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)

saveRDS(sup, "./sup.rds")