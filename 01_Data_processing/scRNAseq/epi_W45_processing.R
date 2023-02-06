library(Seurat)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(harmony)
library(EnhancedVolcano)
library(clusterProfiler)

set.seed(123)
#Load data
W456 <- readRDS("./W456.Rds")

epi_W45 <- subset(W456, subset = cell_class == "Epithelial" &
                PCW %in% c("W4","W5") )

#Clustering
epi_W45 <- NormalizeData(epi_W45)
epi_W45 <- FindVariableFeatures(epi_W45,nfeatures = 2000)

epi_W45 <- CellCycleScoring(epi_W45,s.features = cc.genes.updated.2019$s.genes,
                        g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)

epi_W45 <- ScaleData(epi_W45, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA"))
epi_W45 <- RunPCA(epi_W45, npcs = 50)

epi_W45 <- RunHarmony(
  epi_W45,
  group.by.vars = "orig.ident",
  reduction = 'pca',
  dims.use = 1:50,
  project.dim = F,
  max.iter.harmony = 20,
  max.iter.cluster = 50,
  sigma = 0.2
)
epi_W45 <- RunUMAP(epi_W45, reduction = "harmony", dims = 1:30)
epi_W45 <- FindNeighbors(epi_W45, reduction = "harmony", dims = 1:30)
epi_W45 <- FindClusters(epi_W45,resolution = 0.7)

#Annotation
DotPlot(epi_W45, features = c("PDX1","SOX9","PTF1A","AFP","ALB","HNF4A","HNF1B","ONECUT1",
                          "KRT19","NR2F1","NKX6-2","TOP2A","TBX3","SPP1","ONECUT3","FOXA1","CPA2",
                          "GATA4","CDX2","LGALS3","KCTD12","SHH","NTS","HES6","SOX2","ONECUT2",
                          "NKX6-3"), cluster.idents = T) + RotatedAxis()

epi_W45@meta.data$cell_type <- "Enterocyte"
epi_W45@meta.data$cell_type[which(epi_W45@meta.data$seurat_clusters == 4)] <- "Hepatoblast"
epi_W45@meta.data$cell_type[which(epi_W45@meta.data$seurat_clusters %in% c(2,8))] <- "Dorsal MP"
epi_W45@meta.data$cell_type[which(epi_W45@meta.data$seurat_clusters == 3)] <- "Ventral MP"
epi_W45@meta.data$cell_type[which(epi_W45@meta.data$seurat_clusters == 5)] <- "PB"
epi_W45@meta.data$cell_type[which(epi_W45@meta.data$seurat_clusters == 6)] <- "EHBD"

Idents(epi_W45) <- epi_W45@meta.data$cell_type
epi_W45@meta.data$cell_type <- factor(epi_W45@meta.data$cell_type, levels = c("Dorsal MP","Ventral MP","PB","EHBD","Enterocyte","Hepatoblast"))

saveRDS(epi_W45, "./epi_W45.rds")
