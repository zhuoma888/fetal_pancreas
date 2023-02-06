library(Seurat)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(harmony)

set.seed(123)
#Load data
epi <- readRDS("./epi.rds")

#Extract pancreatic endocrine cells
EP <- subset(epi, subset = cell_type %in% c("EP","Beta","Alpha","Delta","Epsilon"))

#Clustering
EP <- NormalizeData(EP)
EP <- FindVariableFeatures(EP,nfeatures = 2000)
EP <- CellCycleScoring(EP, s.features = cc.genes.updated.2019$s.genes, 
                       g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
EP <- ScaleData(EP, vars.to.regress = c("G2M.Score","S.Score","nCount_RNA","nFeature_RNA","orig.ident") )
EP <- RunPCA(EP, npcs = 50)
EP <- RunHarmony(
  EP,
  group.by.vars = "orig.ident",
  reduction = 'pca',
  dims.use = 1:50,
  project.dim = F,
  max.iter.harmony = 20,
  max.iter.cluster = 50,
  sigma = 0.2
)
EP <- RunUMAP(EP, reduction = "harmony", dims = 1:35)
EP <- FindNeighbors(EP, reduction = "harmony", dims = 1:35)
EP <- FindClusters(EP, resolution = 2)

#Annotation
DotPlot(EP, features = c("SOX9","ID3","RBPJ","TEAD2","NEUROG3","HES6","FEV","PAX4","ARX",
                         "ISL1","ETV1","GCG","PPY","IRX2","GHRL","ONECUT2","INS","NKX6-1","MAFA","MAFB","MEIS2","PLAGL1",
                         "SST","HHEX","CHGA"), group.by = "cell_type") + theme_classic() + RotatedAxis()
EP@meta.data$cell_type <- "Beta"
EP@meta.data$cell_type[which(EP@meta.data$seurat_clusters == 9)] <- "EP early"
EP@meta.data$cell_type[which(EP@meta.data$seurat_clusters == 6)] <- "EP mid"
EP@meta.data$cell_type[which(EP@meta.data$seurat_clusters == 5)] <- "EP late"
EP@meta.data$cell_type[which(EP@meta.data$seurat_clusters == 3)] <- "EP beta"
EP@meta.data$cell_type[which(EP@meta.data$seurat_clusters == 1)] <- "Delta"
EP@meta.data$cell_type[which(EP@meta.data$seurat_clusters == 4)] <- "Epsilon"
EP@meta.data$cell_type[which(EP@meta.data$seurat_clusters == 8)] <- "EP alpha"
EP@meta.data$cell_type[which(EP@meta.data$seurat_clusters == 0)] <- "Alpha/PP"

Idents(EP) <- EP@meta.data$cell_type

saveRDS(EP, file = "./EP.rds")

#Calculate DEG
marker <- FindAllMarkers(EP, logfc.threshold = 0.25, only.pos = T)