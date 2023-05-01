library(Seurat)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(harmony)
library(EnhancedVolcano)
library(clusterProfiler)
library(ggsignif)
library(reshape2)

#Load data
epi <- readRDS("./epi.rds")
epi_W45 <- readRDS("./epi_W45.rds")


##Visulation
#Fig. 2a
DV <- FindMarkers(epi, ident.1 = "Dorsal MP", ident.2 = "Ventral MP")
DV$gene <- rownames(DV)
DV <- DV %>% filter(p_val_adj < 0.01)

group<-ifelse(
  DV$avg_log2FC<(-0.5)& DV$p_val_adj<1e-2,'#1f77b4',
  ifelse(DV$avg_log2FC>(0.5) & DV$p_val_adj<1e-2, '#aec7e8',
         'grey'))

names(group)[group=='#1f77b4']<-'Dorsal MP'
names(group)[group=='grey']<-'Nodiff'
names(group)[group=='#aec7e8']<-'Ventral MP'

EnhancedVolcano(DV, lab = DV$gene,
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
dorsal_GO <- dorsal_GO[which(dorsal_GO$ID %in% c("GO:0016055","GO:0007389","GO:0034329","GO:0050807")),]
dorsal_GO <- dorsal_GO %>% arrange(desc(pvalue))
dorsal_GO$pvalue_log <- -log10(dorsal_GO$pvalue)

ggbarplot(dorsal_GO, y = "pvalue_log", x = "Description",
             color = "#1f77b4", fill = "#1f77b4", rotate = T,
             xlab = "Dorsal MP",
             ylab = "-log10(p-value)") 

ventral_GO <- as.data.frame(ventral_GO)
ventral_GO <- ventral_GO[which(ventral_GO$ID %in% c("GO:0042255","GO:0045661","GO:0034249","GO:0060537")),]
ventral_GO <- ventral_GO %>% arrange(desc(pvalue))
ventral_GO$pvalue_log <- -log10(ventral_GO$pvalue)

ggbarplot(ventral_GO, y = "pvalue_log", x = "Description",
             color = "#aec7e8", fill = "#aec7e8", rotate = T,
             xlab = "Ventral MP",
             ylab = "-log10(p-value)") 
#Fig.2c
epi@meta.data["PDX1"] <- epi@assays$RNA@data["PDX1",]
epi@meta.data["PTF1A"] <- epi@assays$RNA@data["PTF1A",]
epi@meta.data["NKX6-1"] <- epi@assays$RNA@data["NKX6-1",]
epi@meta.data["SOX9"] <- epi@assays$RNA@data["SOX9",]
epi@meta.data["NR2F1"] <- epi@assays$RNA@data["NR2F1",]
epi@meta.data["SIM1"] <- epi@assays$RNA@data["SIM1",]
epi@meta.data["TBX3"] <- epi@assays$RNA@data["TBX3",]
epi@meta.data["SOX6"] <- epi@assays$RNA@data["SOX6",]
epi@meta.data["GPC3"] <- epi@assays$RNA@data["GPC3",]
epi@meta.data["FRZB"] <- epi@assays$RNA@data["FRZB",]
epi@meta.data["FZD5"] <- epi@assays$RNA@data["FZD5",]
epi@meta.data["ID1"] <- epi@assays$RNA@data["ID1",]
epi@meta.data["ID2"] <- epi@assays$RNA@data["ID2",]
epi@meta.data["ID3"] <- epi@assays$RNA@data["ID3",]

exp <- epi@meta.data[,c("barcode","PDX1","PTF1A","NKX6-1","SOX9","NR2F1","SIM1","TBX3","SOX6",
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

#Fig.2g
marker <- FindAllMarkers(epi_W45, logfc.threshold = 0.25, only.pos = T)
marker <- marker %>% filter(p_val_adj < 0.05)

PB_marker <- marker$gene[which(marker$cluster == "PB")]

PB_EntrezID <- bitr(PB_marker, fromType="SYMBOL", toType="ENTREZID",
                    OrgDb="org.Hs.eg.db")

PB_GO <- enrichGO(PB_EntrezID$ENTREZID,,OrgDb=org.Hs.eg.db, ont='BP',
                  pvalueCutoff=0.1,readable = T,
                  keyType='ENTREZID')
PB_GO <- as.data.frame(PB_GO)
PB_GO <- PB_GO[which(PB_GO$ID %in% c("GO:0010975","GO:0048732","GO:0016055","GO:0090287","GO:0007409")),]
PB_GO <- PB_GO %>% arrange(desc(pvalue))
PB_GO$pvalue_log <- -log10(PB_GO$pvalue)

ggbarplot(PB_GO, y = "pvalue_log", x = "Description",
          color = "#ff7f0e", fill = "#ff7f0e", rotate = T,
          xlab = "PB progenitors",
          ylab = "-log10(p-value)") 