library(Seurat)
library(ggplot2)
library(ggsci)
library(tidyverse)

set.seed(20220512)
#Load data
W456 <- readRDS("./W456.rds")
W7_11 <- readRDS("./W7_11.rds")

#Extract pancreatic epithelial cells
epi_W7_11 <- subset(W7_11, subset = cell_class %in% c("Non-endocrine","Endocrine"))
epi_W456 <- subset(W456, subset = organ == "Pancreas")
epi <- merge(epi_W456, epi_W7_11)

#Clustering
epi <- NormalizeData(epi)
epi <- FindVariableFeatures(epi, selection.method = "vst", nfeatures = 3000)
epi <- CellCycleScoring(epi, s.features = cc.genes.updated.2019$s.genes, 
                        g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
epi <- ScaleData(epi, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","nFeature_RNA","orig.ident"))
epi <- RunPCA(epi,npcs = 50)
epi <- RunUMAP(epi,reduction = "pca",dims = 1:40)
epi <- FindNeighbors(epi, reduction = "pca", dims = 1:40)
epi <- FindClusters(epi, resolution = 2.5)

#Annotation
DotPlot(epi, features = c("TOP2A","FOXA2","PDX1","GATA4","SOX9","NKX6-1","NR2F1","CPA2","RBPJL","CLPS","CTRB1",
                          "TTYH1","HES4","DCDC2","CFTR","ASCL2","SAT1","INS","GCG","SST","PPY","GHRL","NEUROG3"), cluster.idents = T) + RotatedAxis()

epi@meta.data$cell_type <- "Early tip"
epi@meta.data$cell_type[which(epi@meta.data$seurat_clusters == 27)] <- "Ventral MP"
epi@meta.data$cell_type[which(epi@meta.data$seurat_clusters == 25)] <- "Dorsal MP"
epi@meta.data$cell_type[which(epi@meta.data$seurat_clusters %in% c(9,13))] <- "Trunk"
epi@meta.data$cell_type[which(epi@meta.data$seurat_clusters %in% c(3,4,5,17,19,22))] <- "Early trunk"
epi@meta.data$cell_type[which(epi@meta.data$seurat_clusters %in% c(8,11,15,18,20,21))] <- "Tip"
epi@meta.data$cell_type[which(epi@meta.data$seurat_clusters %in% c(7,26,32))] <- "Duct"
epi@meta.data$cell_type[which(epi@meta.data$seurat_clusters %in% c(10,34))] <- "Acinar"
epi@meta.data$cell_type[which(epi@meta.data$seurat_clusters == 24)] <- "EP"
epi@meta.data$cell_type[which(epi@meta.data$seurat_clusters == 29)] <- "Beta"
epi@meta.data$cell_type[which(epi@meta.data$seurat_clusters == 28)] <- "Alpha"
epi@meta.data$cell_type[which(epi@meta.data$seurat_clusters == 33)] <- "Delta"
epi@meta.data$cell_type[which(epi@meta.data$seurat_clusters == 35)] <- "Epsilon"

epi@meta.data$cell_type <- factor(epi@meta.data$cell_type, levels = c("Dorsal MP","Ventral MP","Early tip","Tip","Acinar",
                                                                      "Early trunk","Trunk","Duct","EP","Beta","Alpha","Delta","Epsilon"))
epi@meta.data$PCW <- factor(epi@meta.data$PCW, levels = c("W4","W5","W6","W7","W8","W9","W10","W11"))

Idents(epi) <- epi@meta.data$cell_type

saveRDS(epi, file = "./epi.rds")