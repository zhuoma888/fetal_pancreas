library(Seurat)
library(data.table)
library(tibble)
library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(biomaRt)
library(pheatmap)

set.seed(123)

####Comare with exocrine cells
#Load data from GSE115931
Smartseq2 <- read.table("./GSE115931_SmartSeq2.TPM.txt",header = T,sep = "\t",row.names = 1)
Smartseq2[c("ENSMUSG00000063897","ENSMUSG00000095041","ENSMUSG00000099353",
            "ENSMUSG00000099689","ENSMUSG00000100030","ENSMUSG00000100855"),1] <- c("Dhrsx","Pisd","1700031M16Rik","Zfp383","","Gm6976")
Smartseq2 <- Smartseq2[-which(Smartseq2$Symbol == ""),] 
Smartseq2 <- aggregate(.~Symbol,max,data=Smartseq2)
rownames(Smartseq2) <- Smartseq2$Symbol
Smartseq2$Symbol <- NULL

panc_smart <- CreateSeuratObject(Smartseq2, project = "pancreas", names.field = 3)
panc_smart@meta.data$Sample <- rownames(panc_smart@meta.data)

meta <- read.csv("./Smart_meta.csv", header = T)
meta <- meta[which(meta$Note == "/"),]
rownames(meta) <- meta$Sample

panc_smart <- subset(panc_smart, subset = Sample %in% meta$Sample)
panc_smart <- AddMetaData(panc_smart, metadata = meta)

names(panc_smart@meta.data)[which(names(panc_smart@meta.data)=="Embryonic_Day")] <- "time"
panc_smart@meta.data <- panc_smart@meta.data[,c(1,2,3,6,7,16)]

panc_smart <- subset(panc_smart, subset = Batch %in% c("1","2","3","4"))
panc_smart <- subset(panc_smart, subset = Putative_Cell_Type %in% c(
  "Acinar","alpha-1st-early cell","alpha-1st-late cell",
  "alpha-2nd cell","beta cell","Duct","EP1","EP2","EP3",
  "EP4","MP-early","MP-late","Pre-alpha-1st cell","Tip","Trunk"))

exo_embo <- subset(panc_smart, subset = Putative_Cell_Type %in% c("Acinar","Duct","MP-early","MP-late","Tip","Trunk"))

names(exo_embo@meta.data)[which(names(exo_embo@meta.data)=="Putative_Cell_Type")] <- "cell_type"
names(exo_embo@meta.data)[which(names(exo_embo@meta.data)=="speices")] <- "species"

#Load our data
epi <- readRDS("./epi.rds")

exo <- subset(epi, subset = cell_type %in% c("Dorsal MP","Ventral MP","Early tip","Tip",
                                             "Early trunk","Trunk","Acinar","Duct"))

names(exo@meta.data)[which(names(exo@meta.data)=="PCW")] <- "time"

#Sample our dataset
exo.list <- SplitObject(exo, split.by = "cell_type")
subcells <- vector("list", 8)
subcells <- lapply(exo.list, function(x) sample(Cells(x), size=200, replace=F) )

sub <- subcells[[1]]
for (i in 2:length(subcells)) {
  sub <- c(sub, subcells[[i]])
}

exo <- subset(exo, cells= sub)
exo@meta.data$species <- "human"
exo@meta.data$dataset <- "Our data"
exo@meta.data$Batch <- exo@meta.data$orig.ident

#Convert to human gene names
gene <- read.table("./mart_export.txt", header = T, sep = ",") #from BioMart
gene <- gene[which(gene$Human.homology.type == "ortholog_one2one"),]
gene <- gene %>% distinct(Gene.name, .keep_all = T)
gene <- gene[which(gene$Human.gene.name != ""), ]

mouse_exo <- exo_embo@assays$RNA@counts
gene_mouse <- rownames(mouse_exo)

for (i in 1:length(gene_mouse)) { 
  if (gene_mouse[i] %in% gene[,5]) {
    gene_mouse[i] <- gene[which(gene[,5] == gene_mouse[i]),4]
  }
}

rownames(mouse_exo) <- gene_mouse
meta_mouse <- exo_embo@meta.data
meta_mouse <- meta_mouse[,-c(1,2,3)]

mouse_exo <- CreateSeuratObject(mouse_exo, project = "mouse_exo", names.field = 3)
mouse_exo <- AddMetaData(mouse_exo, metadata = meta_mouse)

#Merge two datasets
panc.list <- list(exo, mouse_exo)

panc.list <- lapply(panc.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

panc.anchors<- FindIntegrationAnchors(object.list = panc.list, dims = 1:30, anchor.features = 2000)
exo_merge <- IntegrateData(anchorset = panc.anchors, dims = 1:30 )

DefaultAssay(exo_merge) <- "integrated"
exo_merge <- CellCycleScoring(exo_merge, s.features = cc.genes.updated.2019$s.genes, 
                              g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)

exo_merge <- ScaleData(exo_merge, vars.to.regress = c("nCount_RNA","nFeature_RNA","S.Score","G2M.Score"))
exo_merge <- RunPCA(exo_merge)
exo_merge <- RunUMAP(exo_merge, reduction = "pca", dims = 1:50)
exo_merge <- FindNeighbors(exo_merge, reduction = "pca", dims = 1:50)
exo_merge <- FindClusters(exo_merge, resolution = 1)

exo_merge@meta.data$type <- paste0(exo_merge@meta.data$cell_type,"_",exo_merge@meta.data$species)
exo_merge@meta.data$type <- factor(exo_merge@meta.data$type, levels = c(
  "MP-early_mouse","MP-late_mouse","Tip_mouse","Acinar_mouse","Trunk_mouse",
  "Duct_mouse",
  "Dorsal MP_human","Ventral MP_human","Early tip_human","Tip_human","Acinar_human",
  "Early trunk_human","Trunk_human","Duct_human"))
Idents(exo_merge) <- exo_merge@meta.data$type

saveRDS(exo_merge, file = "./exo_merge.rds")
#Fig. 7a
exo_cols <- setNames(c("#1f77b4","#aec7e8","#fdd0a2","#ff7f0e","#e6550d",
                       "#c7e9c0","#74c476","#31a354",
                       "#c6dbef","#3182bd"),
                     c("Dorsal MP","Ventral MP","Early tip","Tip","Acinar",
                       "Early trunk","Trunk","Duct",
                       "MP-early","MP-late") )

DimPlot(exo_merge, reduction = "umap", label = T,shuffle = T, group.by = "cell_type", split.by = "species") +
  scale_color_manual(values = exo_cols)

###Correlation
panc.list <- SplitObject(exo_merge, split.by = "type")

exp.matrix <- lapply(panc.list, function(x)  GetAssayData(x,slot = "data", assay = "integrated" ))
exp.mean <- lapply(exp.matrix, function(x) apply(x, 1, mean))
for(i in 1:length(exp.mean)) {
  exp.mean[[i]] <- as.data.frame(exp.mean[[i]])
  colnames(exp.mean[[i]]) <- names(exp.mean[i])
}

mean <- do.call(cbind, exp.mean)
mean <- mean[,c("Dorsal MP_human","Ventral MP_human","Early tip_human","Tip_human","Acinar_human",
                "Early trunk_human","Trunk_human","Duct_human",
                "MP-early_mouse","MP-late_mouse","Tip_mouse","Acinar_mouse","Trunk_mouse",
                "Duct_mouse")]
cor <- cor(mean)
cor <- cor[9:14,1:8]

#Fig. 7c (left)
pheatmap(cor,cellwidth = 8,cellheight = 8,border_color = NA, cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#Fig. 7d
DotPlot(exo_merge, features = c("GATA5","ONECUT2","RFX6","NEUROG3","NKX6-2","NR2F2",
                                "NR2F1","HOXB2",
                                "SIM1","TBX3","SOX6","NKX6-1","ID4",
                                "STAT1","PTF1A"),
        idents = c("Dorsal MP_human","Ventral MP_human","MP-early_mouse","MP-late_mouse")) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)

DotPlot(exo_merge, features = c("EGR1","NR4A1","MAFF","EPAS1",
                                "GATA4","NR3C1","PTF1A","NFATC1","NFIL3",
                                "BHLHA15"),
        idents = c("Early tip_human","Tip_human","Acinar_human","Tip_mouse","Acinar_mouse")) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)

DotPlot(exo_merge, features = c("HEY1","HES4","ASCL2","ID1","ID4",
                                "HNF1B","NKX6-1","NEUROG3","SOX5","SOX9"
),
idents = c("Early trunk_human","Trunk_human","Duct_human","Trunk_mouse","Duct_mouse")) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)

####Comare with endocrine cells
#Load data from GSE139627

Smartseq2 <- read.table("./GSE139627_TPM.txt",header = T,sep = "\t",row.names = 1)
Smartseq2[c("ENSMUSG00000063897","ENSMUSG00000095041","ENSMUSG00000099353",
            "ENSMUSG00000099689","ENSMUSG00000100030","ENSMUSG00000100855"),1] <- c("Dhrsx","Pisd","1700031M16Rik","Zfp383","","Gm6976")
Smartseq2 <- Smartseq2[-which(Smartseq2$Symbol == ""),] 
Smartseq2 <- aggregate(.~Symbol,max,data=Smartseq2)
rownames(Smartseq2) <- Smartseq2$Symbol
Smartseq2$Symbol <- NULL

panc_smart2 <- CreateSeuratObject(Smartseq2, project = "pancreas", names.field = 3)
panc_smart2@meta.data$Sample <- rownames(panc_smart2@meta.data)

meta2 <- read.csv("./meta_smart2.csv", header = T, row.names = 1)
meta2 <- meta2[which(meta2$Note == ""),]
panc_smart2 <- subset(panc_smart2,  subset = Sample %in% rownames(meta2))

panc_smart2 <- AddMetaData(panc_smart2, metadata = meta2)

names(panc_smart2@meta.data)[which(names(panc_smart2@meta.data)=="Developmental.Time")] <- "time"
names(panc_smart2@meta.data)[which(names(panc_smart2@meta.data)=="Putative.Cell.Type")] <- "cell_type"
panc_smart2@meta.data <- panc_smart2@meta.data[,c(1,2,3,4,6,8,15)]

panc_smart2 <- subset(panc_smart2, subset = cell_type %in% c(
  "alpha","alpha/PP-Pro-I (epsilon)",
  "alpha/PP-Pro-II","betaearly","betalate","delta","EP1","EP2","EP3",
  "EP4early","EP4late","PP"))
panc_smart2@meta.data$species <- 'mouse'
panc_smart2@meta.data$dataset <- "Yu et al., 2021"

gene <- read.table("/home/marzon/Rproject/integration/mart_export.txt", header = T, sep = ",") #from BioMart
gene <- gene[which(gene$Human.homology.type == "ortholog_one2one"),]
gene <- gene %>% distinct(Gene.name, .keep_all = T)
gene <- gene[which(gene$Human.gene.name != ""), ]

mouse_endo <- panc_smart2@assays$RNA@counts
gene_mouse <- rownames(mouse_endo)

for (i in 1:length(gene_mouse)) { 
  if (gene_mouse[i] %in% gene[,5]) {
    gene_mouse[i] <- gene[which(gene[,5] == gene_mouse[i]),4]
  }
}

rownames(mouse_endo) <- gene_mouse
meta_mouse <- panc_smart2@meta.data
meta_mouse <- meta_mouse[,-c(1,2,3)]

mouse_endo <- CreateSeuratObject(mouse_endo, project = "mouse_endo", names.field = 3)
mouse_endo <- AddMetaData(mouse_endo, metadata = meta_mouse)

mouse_endo <- subset(mouse_endo, subset = time %in% c("E13.5","E14.5","E15.5","E16.5","E17.5"))

#Load our endocrine dataset
EP <- readRDS("./EP.rds")
names(EP@meta.data)[which(names(EP@meta.data)=="PCW")] <- "time"
EP@meta.data$species <- "human"
EP@meta.data$dataset <- "Our data"

mouse_endo@meta.data$cell_type[which(mouse_endo@meta.data$cell_type == "betaearly")] <- "Beta"
mouse_endo@meta.data$cell_type[which(mouse_endo@meta.data$cell_type == "betalate")] <- "Beta"
mouse_endo@meta.data$cell_type[which(mouse_endo@meta.data$cell_type == "EP4early")] <- "EP4"
mouse_endo@meta.data$cell_type[which(mouse_endo@meta.data$cell_type == "EP4late")] <- "EP4"

panc.list <- list(EP, mouse_endo)

panc.list <- lapply(panc.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

panc.anchors<- FindIntegrationAnchors(object.list = panc.list, dims = 1:30, anchor.features = 2000)
endo_merge <- IntegrateData(anchorset = panc.anchors, dims = 1:30 )

DefaultAssay(endo_merge) <- "integrated"
endo_merge <- CellCycleScoring(endo_merge, s.features = cc.genes.updated.2019$s.genes, 
                               g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
endo_merge <- ScaleData(endo_merge, vars.to.regress = c("nCount_RNA","nFeature_RNA","S.Score","G2M.Score"))
endo_merge <- RunPCA(endo_merge)
endo_merge <- RunUMAP(endo_merge, reduction = "pca", dims = 1:50)
endo_merge <- FindNeighbors(endo_merge, reduction = "pca", dims = 1:50)
endo_merge <- FindClusters(endo_merge, resolution = 0.8)

endo_cols <- setNames(c("#c6dbef","#6baed6","#3182bd",
                        "#fdae6b","#e6550d","#636363",
                        "#a1d99b","#31a354","#756bb1",
                        "#c6dbef","#9ecae1","#6baed6","#3182bd",
                        "#e6550d","#8c564b",
                        "#fdae6b","#756bb1"),
                      c("EP early","EP mid","EP late",
                        "EP alpha","Alpha/PP","Epsilon",
                        "EP beta","Beta","Delta",
                        "EP1","EP2","EP3","EP4",
                        "alpha","PP",
                        "Alpha/PP−Pro","delta") )
endo_merge@meta.data$cell_type[which(endo_merge@meta.data$cell_type == "alpha/PP-Pro-I (epsilon)")] <- "Epsilon"
endo_merge@meta.data$cell_type[which(endo_merge@meta.data$cell_type == "alpha/PP-Pro-II")] <- "Alpha/PP−Pro"

#Fig. 7b
DimPlot(endo_merge, reduction = "umap", label = T,shuffle = T, group.by = "cell_type", split.by = "species") +
  scale_color_manual(values = endo_cols)

endo_merge@meta.data$type <- paste0(endo_merge@meta.data$cell_type,"_",endo_merge@meta.data$species)
endo_merge@meta.data$type <- factor(endo_merge@meta.data$type, levels = c(
  "EP1_mouse","EP2_mouse","EP3_mouse","EP4_mouse","Epsilon_mouse",
  "Alpha/PP−Pro_mouse","alpha_mouse","PP_mouse","Beta_mouse","delta_mouse",
  "EP early_human","EP mid_human","EP late_human","Epsilon_human","EP alpha_human",
  "Alpha/PP_human","EP beta_human","Beta_human","Delta_human"))
Idents(endo_merge) <- endo_merge@meta.data$type
saveRDS(endo_merge, file = "./endo_merge.rds")

#Fig. 7c (right)
panc.list <- SplitObject(endo_merge, split.by = "type")

exp.matrix <- lapply(panc.list, function(x)  GetAssayData(x,slot = "data", assay = "integrated" ))
exp.mean <- lapply(exp.matrix, function(x) apply(x, 1, mean))
for(i in 1:length(exp.mean)) {
  exp.mean[[i]] <- as.data.frame(exp.mean[[i]])
  colnames(exp.mean[[i]]) <- names(exp.mean[i])
}

mean <- do.call(cbind, exp.mean)
mean <- mean[,c("EP early_human","EP mid_human","EP late_human","Epsilon_human","EP alpha_human",
                "Alpha/PP_human","EP beta_human","Beta_human","Delta_human",
                "EP1_mouse","EP2_mouse","EP3_mouse","EP4_mouse","Epsilon_mouse",
                "Alpha/PP−Pro_mouse","alpha_mouse","PP_mouse","Beta_mouse","delta_mouse")]
cor <- cor(mean)
cor <- cor[10:19,1:9]

pheatmap(cor,cellwidth = 8,cellheight = 8,border_color = NA, cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#Fig. 7e
DotPlot(endo_merge, features = c("NEUROG3","FOXA3","ESR1","LHX1","AKNA",
                                 "ZBTB18","SETBP1","SOX11","INSM1","ZNF503",
                                 "DACH1","ARID5B","POU2F2"),
        idents = c("EP early_human","EP mid_human","EP late_human","EP1_mouse","EP2_mouse",
                   "EP3_mouse","EP4_mouse")) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)

DotPlot(endo_merge, features = c("SOX4","SOX6",
                                 "LRRFIP1","ZFHX3","ETV1","SLC38A5","GAST",
                                 "VTN","FGF14","VSTM2L","CHL1",
                                 "POU6F2","IRX2","MAFB","PAX6"),
        idents = c("EP alpha_human","Alpha/PP_human","Epsilon_human",
                   "Alpha/PP−Pro_mouse","alpha_mouse","PP_mouse","Epsilon_mouse")) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)

DotPlot(endo_merge, features = c("PLAGL1","ASCL2","ID4","MNX1","SAMD11","MAFA","UCN3","SLC2A2","SLC2A1","IAPP","G6PC2","SIX2",
                                 "NFIL3","NPY","GIP",
                                 "SOX4","PBX3","MAFB","NKX6-1","POU2F2",
                                 "EHF"),
        idents = c("EP beta_human","Beta_human","Delta_human",
                   "Beta_mouse","delta_mouse")) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)

