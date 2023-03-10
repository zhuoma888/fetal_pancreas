library(Seurat)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(SoupX)

set.seed(123)

####Load data
S9 <- load10X("./data/S9/outs/")
S13 <- load10X("./data/S13/outs/")
S14 <- load10X("./data/S14/outs/")
S27 <- load10X("./data/S27/outs/")
S35 <- load10X("./data/S35/outs/")
S37 <- load10X("./data/S37/outs/")
S57 <- load10X("./data/S57/outs/")
S66 <- load10X("./data/S66/outs/")
S76 <- load10X("./data/S76/outs/")
S82 <- load10X("./data/S82/outs/")
S85 <- load10X("./data/S85/outs/")
S86 <- load10X("./data/S86/outs/")
S88 <- load10X("./data/S88/outs/")
S91 <- load10X("./data/S91/outs/")
S113 <- load10X("./data/S113/outs/")
S120 <- load10X("./data/S120/outs/")
S129 <- load10X("./data/S129/outs/")

panc.list <- list(S9,S13,S14,S27,S35,S37,S57,S66,S76,
                  S82,S85,S86,S88,S91,S113,S120,S129)

#remove ambient RNA contamination
panc.soupx <- lapply(panc.list, function(x) autoEstCont(x))
panc.soupx <- lapply(panc.soupx, function(x) adjustCounts(x, roundToInt = T))

#creat seurat object
colnames(panc.soupx[[1]]) <- paste(colnames(panc.soupx[[1]]),"S9",sep="_")
colnames(panc.soupx[[2]]) <- paste(colnames(panc.soupx[[2]]),"S13",sep="_")
colnames(panc.soupx[[3]]) <- paste(colnames(panc.soupx[[3]]),"S14",sep="_")
colnames(panc.soupx[[4]]) <- paste(colnames(panc.soupx[[4]]),"S27",sep="_")
colnames(panc.soupx[[5]]) <- paste(colnames(panc.soupx[[5]]),"S35",sep="_")
colnames(panc.soupx[[6]]) <- paste(colnames(panc.soupx[[6]]),"S37",sep="_")
colnames(panc.soupx[[7]]) <- paste(colnames(panc.soupx[[7]]),"S57",sep="_")
colnames(panc.soupx[[8]]) <- paste(colnames(panc.soupx[[8]]),"S66",sep="_")
colnames(panc.soupx[[9]]) <- paste(colnames(panc.soupx[[9]]),"S76",sep="_")
colnames(panc.soupx[[10]]) <- paste(colnames(panc.soupx[[10]]),"S82",sep="_")
colnames(panc.soupx[[11]]) <- paste(colnames(panc.soupx[[11]]),"S85",sep="_")
colnames(panc.soupx[[12]]) <- paste(colnames(panc.soupx[[12]]),"S86",sep="_")
colnames(panc.soupx[[13]]) <- paste(colnames(panc.soupx[[13]]),"S88",sep="_")
colnames(panc.soupx[[14]]) <- paste(colnames(panc.soupx[[14]]),"S91",sep="_")
colnames(panc.soupx[[15]]) <- paste(colnames(panc.soupx[[15]]),"S113",sep="_")
colnames(panc.soupx[[16]]) <- paste(colnames(panc.soupx[[16]]),"S120",sep="_")
colnames(panc.soupx[[17]]) <- paste(colnames(panc.soupx[[17]]),"S129",sep="_")

panc <- cbind(panc.soupx[[1]],panc.soupx[[2]],panc.soupx[[3]],panc.soupx[[4]],
              panc.soupx[[5]],panc.soupx[[6]],panc.soupx[[7]],panc.soupx[[8]],
              panc.soupx[[9]],panc.soupx[[10]],panc.soupx[[11]],panc.soupx[[12]],
              panc.soupx[[13]],panc.soupx[[14]],panc.soupx[[15]],panc.soupx[[16]],
              panc.soupx[[17]])
panc <- CreateSeuratObject(panc, names.field = 2, names.delim = "_")

#Add metadata
panc@meta.data[["PCW"]] <- "W5"
panc@meta.data$PCW[which(panc@meta.data$orig.ident %in% c("S57","S113"))] <- "W4"
panc@meta.data$PCW[which(panc@meta.data$orig.ident %in% c("S88","S91","S37"))] <- "W6"
panc@meta.data$PCW[which(panc@meta.data$orig.ident %in% c("S5","S120","S129"))] <- "W7"
panc@meta.data$PCW[which(panc@meta.data$orig.ident %in% c("S9","S27","S66","S86"))] <- "W8"
panc@meta.data$PCW[which(panc@meta.data$orig.ident == "S85")] <- "W9"
panc@meta.data$PCW[which(panc@meta.data$orig.ident == "S35")] <- "W10"
panc@meta.data$PCW[which(panc@meta.data$orig.ident %in% c("S13","S14"))] <- "W11"
panc@meta.data[["carnegie_stage"]] <- "fetal_stage"
panc@meta.data$carnegie_stage[which(panc@meta.data$orig.ident %in% c("S57","S113"))] <- "CS13"
panc@meta.data$carnegie_stage[which(panc@meta.data$orig.ident == "S76")] <- "CS14"
panc@meta.data$carnegie_stage[which(panc@meta.data$orig.ident %in% c("S82"))] <- "CS15"
panc@meta.data$carnegie_stage[which(panc@meta.data$orig.ident %in% c("S88","S91"))] <- "CS16"
panc@meta.data$carnegie_stage[which(panc@meta.data$orig.ident == "S37")] <- "CS17"
panc@meta.data$carnegie_stage[which(panc@meta.data$orig.ident == "S120")] <- "CS19"
panc@meta.data$carnegie_stage[which(panc@meta.data$orig.ident == "S129")] <- "CS18"
panc@meta.data$carnegie_stage[which(panc@meta.data$orig.ident == "S9")] <- "CS21"
panc@meta.data$carnegie_stage[which(panc@meta.data$orig.ident %in% c("S27","S66","S86"))] <- "CS23"

#Quality control
panc[["percent.mito"]] <- PercentageFeatureSet(panc,pattern = "^MT-")
panc <- subset(panc, subset = 
                 nFeature_RNA > 500 & 
                 nCount_RNA > 4000 &
                 nFeature_RNA < 8000 & 
                 nCount_RNA < 50000 & 
                 percent.mito < 15)

saveRDS(panc, file = "./panc_all.rds")

####W456 sample clustering
W456 <- subset(panc, subset =  PCW %in% c("W4","W5","W6"))
W456 <- NormalizeData(W456)
W456 <- FindVariableFeatures(W456, selection.method = "vst", nfeatures = 3000)
W456 <- CellCycleScoring(W456, s.features = cc.genes.updated.2019$s.genes, 
                         g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
W456 <- ScaleData(W456, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","nFeature_RNA","orig.ident"))
W456 <- RunPCA(W456,npcs = 50)
W456 <- RunUMAP(W456,reduction = "pca",dims = 1:50)
W456 <- FindNeighbors(W456, reduction = "pca", dims = 1:50)
W456 <- FindClusters(W456, resolution = 1)

#First annotation
DotPlot(W456, features = c("EPCAM","PDX1","CPA2","CHGA","NKX2-2","NKX6-1","PROX1","SPP1","CDX2","TOP2A",
                           "COL3A1","PECAM1","ASCL1","HBA1","PTPRC","SOX9","SOX17","GATA4","AFP","ALB"), cluster.idents = T)

W456@meta.data$cell_class <- "Mesenchymal"
W456@meta.data$cell_class[which(W456@meta.data$seurat_clusters %in% c(5,8,16,20,27))] <- "Epithelial"
W456@meta.data$cell_class[which(W456@meta.data$seurat_clusters == 24)] <- "Neural"
W456@meta.data$cell_class[which(W456@meta.data$seurat_clusters == 28)] <- "Erythroid"
W456@meta.data$cell_class[which(W456@meta.data$seurat_clusters == 25)] <- "Immune"
W456@meta.data$cell_class[which(W456@meta.data$seurat_clusters %in% c(14,33))] <- "Endothelial"

#remove doublets
W456@meta.data["mix"] <- "False"
W456@meta.data["mes"] <- W456@assays$RNA@data["COL3A1",]
W456@meta.data["hb"] <- W456@assays$RNA@data["HBA1",]
W456@meta.data["endo"] <- W456@assays$RNA@data["PECAM1",]
W456@meta.data["immune"] <- W456@assays$RNA@data["PTPRC",]
W456@meta.data["epi"] <- W456@assays$RNA@data["EPCAM",]
W456@meta.data["neuro"] <- W456@assays$RNA@data["ASCL1",]

W456@meta.data$mix[which(W456@meta.data$cell_class == "Epithelial" &
                           W456@meta.data$mes > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Epithelial" &
                           W456@meta.data$hb > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Epithelial" &
                           W456@meta.data$endo > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Epithelial" &
                           W456@meta.data$immune > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Epithelial" &
                           W456@meta.data$neuro > 0.2)] <- "True"

W456@meta.data$mix[which(W456@meta.data$cell_class == "Mesenchymal" &
                           W456@meta.data$epi > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Mesenchymal" &
                           W456@meta.data$hb > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Mesenchymal" &
                           W456@meta.data$endo > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Mesenchymal" &
                           W456@meta.data$immune > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Mesenchymal" &
                           W456@meta.data$neuro > 0.2)] <- "True"


W456@meta.data$mix[which(W456@meta.data$cell_class == "Endothelial" &
                           W456@meta.data$epi > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Endothelial" &
                           W456@meta.data$hb > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Endothelial" &
                           W456@meta.data$mes > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Endothelial" &
                           W456@meta.data$immune > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Endothelial" &
                           W456@meta.data$neuro > 0.2)] <- "True"


W456@meta.data$mix[which(W456@meta.data$cell_class == "Neural" &
                           W456@meta.data$epi > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Neural" &
                           W456@meta.data$hb > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Neural" &
                           W456@meta.data$endo > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Neural" &
                           W456@meta.data$immune > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Neural" &
                           W456@meta.data$mes > 0.2)] <- "True"


W456@meta.data$mix[which(W456@meta.data$cell_class == "Immune" &
                           W456@meta.data$epi > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Immune" &
                           W456@meta.data$hb > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Immune" &
                           W456@meta.data$endo > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Immune" &
                           W456@meta.data$mes > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Immune" &
                           W456@meta.data$neuro > 0.2)] <- "True"

W456@meta.data$mix[which(W456@meta.data$cell_class == "Erythroid" &
                           W456@meta.data$epi > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Erythroid" &
                           W456@meta.data$mes > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Erythroid" &
                           W456@meta.data$endo > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Erythroid" &
                           W456@meta.data$immune > 0.2)] <- "True"
W456@meta.data$mix[which(W456@meta.data$cell_class == "Erythroid" &
                           W456@meta.data$neuro > 0.2)] <- "True"

W456 <- subset(W456, subset = mix == "False")

#Further annote epithelial cells
W456@meta.data$cell_class[which(W456@meta.data$seurat_clusters == 27)] <- "Liver"
W456@meta.data$cell_class[which(W456@meta.data$seurat_clusters %in% c(8,20))] <- "Pancreas"
W456@meta.data$cell_class[which(W456@meta.data$seurat_clusters == 16)] <- "EHBD"
W456@meta.data$cell_class[which(W456@meta.data$seurat_clusters == 5)] <- "Duodenum"

saveRDS(W456, file = "./W456.rds")

####W7_11 sample clustering
W7_11 <- subset(panc, subset = PCW %in% c("W7","W8","W9","W10","W11") )
W7_11 <- NormalizeData(W7_11)
W7_11 <- FindVariableFeatures(W7_11, selection.method = "vst", nfeatures = 2000)
W7_11 <- CellCycleScoring(W7_11, s.features = cc.genes.updated.2019$s.genes, 
                          g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
W7_11 <- ScaleData(W7_11, vars.to.regress = c("S.Score","G2M.Score","orig.ident","nCount_RNA","nFeature_RNA"))
W7_11 <- RunPCA(W7_11,npcs = 50)
W7_11 <- RunUMAP(W7_11,reduction = "pca",dims = 1:50)
W7_11 <- FindNeighbors(W7_11, reduction = "pca", dims = 1:50)
W7_11 <- FindClusters(W7_11, resolution = 0.3)

#First annotation
DotPlot(W7_11, features = c("EPCAM","PDX1","CPA2","CHGA","NKX2-2","NEUROD1","COL3A1","COL1A1","PECAM1","PTPRC","ASCL1","UPK3B",
                            "HBA1","POU5F1"), cluster.idents = T) + RotatedAxis()

W7_11@meta.data$cell_class <- "Mesenchymal"
W7_11@meta.data$cell_class[which(W7_11@meta.data$seurat_clusters %in% c(0,3,6))] <- "Non-endocrine"
W7_11@meta.data$cell_class[which(W7_11@meta.data$seurat_clusters == 11)] <- "Endocrine"
W7_11@meta.data$cell_class[which(W7_11@meta.data$seurat_clusters %in% c(4,16))] <- "Endothelial"
W7_11@meta.data$cell_class[which(W7_11@meta.data$seurat_clusters %in% c(14,15))] <- "Erythroid"
W7_11@meta.data$cell_class[which(W7_11@meta.data$seurat_clusters == 12)] <- "Immune"
W7_11@meta.data$cell_class[which(W7_11@meta.data$seurat_clusters == 8)] <- "Neural"

#remove doublets
W7_11@meta.data["mix"] <- "False"
W7_11@meta.data["mes"] <- W7_11@assays$RNA@data["COL3A1",]
W7_11@meta.data["hb"] <- W7_11@assays$RNA@data["HBA1",]
W7_11@meta.data["endo"] <- W7_11@assays$RNA@data["PECAM1",]
W7_11@meta.data["immune"] <- W7_11@assays$RNA@data["PTPRC",]
W7_11@meta.data["epi"] <- W7_11@assays$RNA@data["EPCAM",]
W7_11@meta.data["neuro"] <- W7_11@assays$RNA@data["ASCL1",]
W7_11@meta.data["CPA2"] <- W7_11@assays$RNA@data["CPA2",]
W7_11@meta.data["NEUROD1"] <- W7_11@assays$RNA@data["NEUROD1",]
table(W7_11@meta.data$mix)

W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Non-endocrine" &
                            W7_11@meta.data$mes > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Non-endocrine" &
                            W7_11@meta.data$hb > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Non-endocrine" &
                            W7_11@meta.data$endo > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Non-endocrine" &
                            W7_11@meta.data$immune > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Non-endocrine" &
                            W7_11@meta.data$neuro > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Non-endocrine" &
                            W7_11@meta.data$NEUROD1 > 0.2)] <- "True"

W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Endocrine" &
                            W7_11@meta.data$mes > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Endocrine" &
                            W7_11@meta.data$hb > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Endocrine" &
                            W7_11@meta.data$endo > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Endocrine" &
                            W7_11@meta.data$immune > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Endocrine" &
                            W7_11@meta.data$neuro > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Endocrine" &
                            W7_11@meta.data$CPA2 > 0.2)] <- "True"

W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Mesenchymal" &
                            W7_11@meta.data$NEUROD1 > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Mesenchymal" &
                            W7_11@meta.data$hb > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Mesenchymal" &
                            W7_11@meta.data$endo > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Mesenchymal" &
                            W7_11@meta.data$immune > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Mesenchymal" &
                            W7_11@meta.data$neuro > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Mesenchymal" &
                            W7_11@meta.data$epi > 0.2)] <- "True"

W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Endothelial" &
                            W7_11@meta.data$NEUROD1 > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Endothelial" &
                            W7_11@meta.data$hb > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Endothelial" &
                            W7_11@meta.data$mes > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Endothelial" &
                            W7_11@meta.data$immune > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Endothelial" &
                            W7_11@meta.data$neuro > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Endothelial" &
                            W7_11@meta.data$epi > 0.2)] <- "True"

W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Erythroid" &
                            W7_11@meta.data$NEUROD1 > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Erythroid" &
                            W7_11@meta.data$endo > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Erythroid" &
                            W7_11@meta.data$mes > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Erythroid" &
                            W7_11@meta.data$immune > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Erythroid" &
                            W7_11@meta.data$neuro > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Erythroid" &
                            W7_11@meta.data$epi > 0.2)] <- "True"

W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Immune" &
                            W7_11@meta.data$NEUROD1 > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Immune" &
                            W7_11@meta.data$endo > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Immune" &
                            W7_11@meta.data$mes > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Immune" &
                            W7_11@meta.data$hb > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Immune" &
                            W7_11@meta.data$neuro > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Immune" &
                            W7_11@meta.data$epi > 0.2)] <- "True"

W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Neural" &
                            W7_11@meta.data$NEUROD1 > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Neural" &
                            W7_11@meta.data$endo > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Neural" &
                            W7_11@meta.data$mes > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Neural" &
                            W7_11@meta.data$hb > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Neural" &
                            W7_11@meta.data$immune > 0.2)] <- "True"
W7_11@meta.data$mix[which(W7_11@meta.data$cell_class == "Neural" &
                            W7_11@meta.data$epi > 0.2)] <- "True"

W7_11 <- subset(W7_11, subset = mix == "False")

saveRDS(W7_11, file = "./W7_11.rds")