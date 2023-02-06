library(Seurat)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(harmony)
library(EnhancedVolcano)
library(clusterProfiler)

#Load data
epi_W45 <- readRDS("./epi_W45.rds")

#Calculate DEG
marker <- FindAllMarkers(epi_W45, logfc.threshold = 0.25, only.pos = T)

##Visulation
#Fig. 2a
DV <- FindMarkers(epi_W45, ident.1 = "Dorsal MP", ident.2 = "Ventral MP")
DV$gene <- rownames(DV)
DV <- DV %>% filter(p_val_adj < 0.01)

group<-ifelse(
  DV$avg_log2FC<(-0.5)& DV$p_val_adj<1e-2,'#1f77b4',
  ifelse(DV$avg_log2FC>(0.5) & DV$p_val_adj<1e-2, '#aec7e8',
         'grey'))

names(group)[group=='#1f77b4']<-'Dorsal MP'
names(group)[group=='grey']<-'Nodiff'
names(group)[group=='#aec7e8']<-'Ventral MP'

EnhancedVolcano(VP_EHBD, lab = VP_EHBD$gene,
                x = "avg_log2FC",
                y = "p_val_adj", pCutoff = 1E-2, FCcutoff = 0.5,
                colAlpha= 1, xlim = c(-2, 2),
                labSize = 2, pointSize = 0.5,colCustom = group) + theme_classic()

#Fig. 2b
dorsal <- DV$gene[which(DV$avg_log2FC > 0)]
ventral <- DV$gene[which(DV$avg_log2FC < 0)]
dorsal_EntrezID <- bitr(dorsal, fromType="SYMBOL", toType="ENTREZID",
                        OrgDb="org.Hs.eg.db")
ventral_EntrezID <- bitr(ventral, fromType="SYMBOL", toType="ENTREZID",
                         OrgDb="org.Hs.eg.db")

dorsal_GO <- enrichGO(dorsal_EntrezID$ENTREZID,,OrgDb=org.Hs.eg.db, ont='BP',
                      pvalueCutoff=0.1,readable = T,
                      keyType='ENTREZID')
ventral_GO <- enrichGO(ventral_EntrezID$ENTREZID,,OrgDb=org.Hs.eg.db, ont='BP',
                       pvalueCutoff=0.1,readable = T,
                       keyType='ENTREZID')
dorsal_GO <- as.data.frame(dorsal_GO)
ventral_GO <- as.data.frame(ventral_GO)

#Fig.2c
epi_W45@meta.data["PDX1"] <- epi_W45@assays$RNA@data["PDX1",]
epi_W45@meta.data["PTF1A"] <- epi_W45@assays$RNA@data["PTF1A",]
epi_W45@meta.data["NKX6-1"] <- epi_W45@assays$RNA@data["NKX6-1",]
epi_W45@meta.data["SOX9"] <- epi_W45@assays$RNA@data["SOX9",]
epi_W45@meta.data["NR2F1"] <- epi_W45@assays$RNA@data["NR2F1",]
epi_W45@meta.data["SIM1"] <- epi_W45@assays$RNA@data["SIM1",]
epi_W45@meta.data["TBX3"] <- epi_W45@assays$RNA@data["TBX3",]
epi_W45@meta.data["SOX6"] <- epi_W45@assays$RNA@data["SOX6",]
epi_W45@meta.data["GPC3"] <- epi_W45@assays$RNA@data["GPC3",]
epi_W45@meta.data["FRZB"] <- epi_W45@assays$RNA@data["FRZB",]
epi_W45@meta.data["FZD5"] <- epi_W45@assays$RNA@data["FZD5",]
epi_W45@meta.data["ID1"] <- epi_W45@assays$RNA@data["ID1",]
epi_W45@meta.data["ID2"] <- epi_W45@assays$RNA@data["ID2",]
epi_W45@meta.data["ID3"] <- epi_W45@assays$RNA@data["ID3",]

exp <- epi_W45@meta.data[,c("barcode","PDX1","PTF1A","NKX6-1","SOX9","NR2F1","SIM1","TBX3","SOX6",
                            "GPC3","FRZB","FZD5",
                            "ID1","ID2","ID3","cell_type")]
exp <- exp[which(exp$cell_type %in% c("Dorsal MP","Ventral MP")),]

exp <- melt(exp, id.vars=c("barcode","cell_type"),variable.name="gene", value.name="Expression")
exp$gene <- factor(exp$gene, levels = c("PDX1","PTF1A","NKX6-1","SOX9",
                                        "NR2F1","SIM1","TBX3","SOX6",
                                        "GPC3","FRZB","FZD5","ID1","ID2","ID3"))

DV_col <- setNames(c("#1f77b4","#aec7e8"),
                   c("Dorsal MP","Ventral MP"))

ggplot(exp, aes(x=cell_type, y=Expression, fill = cell_type)) + geom_boxplot(outlier.size = 0.1) + facet_wrap(~ gene, ncol = 8) +
  geom_signif(comparisons = list(c("Dorsal MP","Ventral MP")),test = "wilcox.test", textsize = 2) +
  theme_classic() + scale_fill_manual(values = DV_col) + theme(axis.text.x = element_text(angle = 45, hjust =0.5, vjust =0.5))

#Fig.2d
W45_col <- setNames(c("#1f77b4","#aec7e8","#ff7f0e","#9467bd","#2ca02c","#d62728"),
                    c("Dorsal MP","Ventral MP","PB","EHBD","Enterocyte","Hepatoblast"))

DimPlot(epi_W45, reduction = "umap",pt.size = 0.5,shuffle = T, group.by = "cell_type") + scale_color_manual(values = W45_col)
DimPlot(epi_W45, reduction = "umap",pt.size = 0.5,shuffle = T, group.by = "carnegie_stage") + scale_color_viridis_d()

#Fig.2e
gene <- c("PDX1","PTF1A","NKX6-1","NR2F1","SOX6","TBX3","ISL1","NKX6-2","HHEX","SULT1E1","SPP1",
          "LGALS3","CDX2","RFX6","ALB","APOA2")

DotPlot(epi_W45, features = gene, group.by = "cell_type") + theme_classic() + RotatedAxis() +
  scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0) 

#Fig.2f
FeaturePlot(epi_W45, features = c("ISL1","SULT1E1","HHEX","NKX6-2"), reduction = "umap", ncol = 2, cols = c("grey","red") )
