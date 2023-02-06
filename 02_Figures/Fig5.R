library(Seurat)
library(ggplot2)
library(ggsci)
library(Seurat)
library(tidyverse)
library(monocle3)
library(circlize)
library(ComplexHeatmap)
library(clusterProfiler)
library(ggsignif)
library(Signac)
library(SeuratWrappers)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(chromVARmotifs)
library(Matrix)
library(Matrix.utils)

#Load data
epi <- readRDS("./epi.rds")
cds <- readRDS("./cds.rds")
pseudotime <- as.data.frame(pseudotime(cds))
names(pseudotime) <- "pseudotime"
epi <- AddMetaData(epi, metadata = pseudotime)

#Fig.5a
ep <- subset(epi, subset = cell_type %in%  c("EP","Trunk") &
               pseudotime > 8 )
duct <- subset(epi, subset = cell_type %in%  c("Duct","Trunk") &
                 pseudotime > 8)
nrow(ep@meta.data)
nrow(duct@meta.data)

marker <- read.table("./epi_marker_genes.txt", header = T, sep = "\t")
genes <- unique(marker$gene[which(marker$cluster %in% c("Trunk","Duct","EP") &
                                    marker$avg_log2FC > 0.4)])

ep@meta.data <- ep@meta.data %>% arrange(-pseudotime)
duct@meta.data <- duct@meta.data %>% arrange(pseudotime)

mat1 <- ep@assays$RNA@data[genes,rownames(ep@meta.data)]
mat2 <- duct@assays$RNA@data[genes,rownames(duct@meta.data)]
mat <- cbind(mat1,mat2)
mat <- t(apply(mat,1,function(x){smooth.spline(x,df=3)$y}))
mat <- t(apply(mat,1,function(x){(x-mean(x))/sd(x)}))

meta1 <- ep@meta.data[,c("cell_type","PCW")]
meta2 <- duct@meta.data[,c("cell_type","PCW")]
meta <- rbind(meta1,meta2)
meta$cell_type <- factor(meta$cell_type, levels = c("Trunk","Duct","EP"))

type_col <- setNames(c("#74c476","#31a354","#17becf"),
                     c("Trunk","Duct","EP"))
pcw_col <- setNames(c("#440154","#46337E","#365C8D","#277F8E","#1FA187",
                      "#4AC16D","#9FDA3A","#FDE725"),
                    c("W4","W5","W6","W7","W8","W9","W10","W11"))

col_anno <- HeatmapAnnotation(Cell_type=meta$cell_type,PCW=meta$PCW,
                              col=list(Cell_type = type_col,
                                       PCW = pcw_col))
set.seed(234)
htmp <- Heatmap(
  mat,
  name                         = "z-score",
  km = 4,
  col                          = colorRamp2(seq(from=-2,to=2,length=11),colorRampPalette(c("navy", "white", "firebrick3"))(11)),
  show_row_names               = F,
  show_column_names            = FALSE,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation = col_anno,
  column_split = c(rep("EP",1418),rep("Duct",2179)) )
dht <- draw(htmp)
clusters <- row_order(dht)

#Fig.5b
gene <- vector("list",4)
gene[[1]] <- genes[clusters$`1`]
gene[[2]] <- genes[clusters$`2`]
gene[[3]] <- genes[clusters$`3`]
gene[[4]] <- genes[clusters$`4`]

gene_id <- lapply(gene, function(x) bitr(x, fromType="SYMBOL", toType="ENTREZID",
                                         OrgDb="org.Hs.eg.db") )

GO <- vector("list",4)
for (i in 1:length(gene)) {
  GO[[i]] <- enrichGO(gene_id[[i]]$ENTREZID,OrgDb=org.Hs.eg.db, ont='BP',pvalueCutoff=0.05,
                      keyType='ENTREZID', readable = T)
}

GO_df <- lapply(GO, function(x) as.data.frame(x))

GO_df[[1]]$group <- "1"
GO_df[[2]]$group <- "2"
GO_df[[3]]$group <- "3"
GO_df[[4]]$group <- "4"

GO_df <- rbind(GO_df[[1]],GO_df[[2]],GO_df[[3]],GO_df[[4]])

#Fig.5c
epi <- subset(epi, subset = cell_type %in% c("Trunk","Duct","EP"))
epi@meta.data["ASCL2"] <- epi@assays$RNA@data["ASCL2",]
exp <- epi@meta.data[,c("barcode","ASCL2","cell_type")]
exp_ascl2 <- melt(exp, id.vars=c("barcode","cell_type"), variable.name="gene", value.name="Expression")

ggplot(exp_ascl2, aes(x=cell_type, y=Expression, fill = cell_type)) + geom_boxplot(outlier.size = 0.1) + facet_wrap(~ gene, ncol = 1) +
  theme_classic() + scale_fill_manual(values = type_col) + theme(axis.text.x = element_text(angle = 45, hjust =0.5, vjust =0.5)) +
  geom_signif(comparisons = list(c("Trunk","EP"),
                                 c("Trunk","Duct")),
              step_increase = 0.1,
              test = "wilcox.test")

#Fig.5e
epi@meta.data["HDAC2"] <- epi@assays$RNA@data["HDAC2",]
epi@meta.data["KDM1A"] <- epi@assays$RNA@data["KDM1A",]
epi@meta.data["RCOR2"] <- epi@assays$RNA@data["RCOR2",]

exp <- epi@meta.data[,c("barcode","HDAC2","KDM1A","RCOR2","cell_type")]
exp <- melt(exp, id.vars=c("barcode","cell_type"),variable.name="gene", value.name="Expression")
ggplot(exp_other, aes(x=cell_type, y=Expression, fill = cell_type)) + geom_boxplot(outlier.size = 0.1) + facet_wrap(~ gene, ncol = 1) +
  theme_classic() + scale_fill_manual(values = type_col) + theme(axis.text.x = element_text(angle = 45, hjust =0.5, vjust =0.5)) +
  geom_signif(comparisons = list(c("Trunk","EP"),
                                 c("EP","Duct")),
              step_increase = 0.1,
              test = "wilcox.test")

####scATAC
#Fig.5f
w1 <- c("INSM1","ONECUT2","FOXA2")
split1 <- c("YBX1","YBX1","HES4")
w2 <- c("GLIS3","TCF12","KLF5","TCF7L1","ID1","TCF7L2","HES4","MAFA","ELF3")
split2 <- c("HES4","HES4","HES4","HES1","HES1","HMGA2","ID4","ASCL2","ASCL2")
avg <- AverageExpression(rna,assays = "RNA",group.by = "celltype",features = c(w1,w2))
avg <- avg[[1]]
avg <- avg[,3:5]
avg_scaled = t(scale(t(avg)))
splite <- c(split1,split2)
row_ha = rowAnnotation(type=c("p",'p','n','p','p','p',"p",'p','p',"p",'p','n'),
                       col = list(type = c("p" = "red", "n" = "black")))
tfscoret <- data.frame(tf=c("YBX1","ASCL2","HES4","HES1","HMGA2","ID4"))
rownames(tfscoret) <- tfscoret$tf
tfscoret$score <- 0
for (i in 1:nrow(tfscoret)) {
  tfscoret$score[i] <- mean(trunk_footprints_cl$center[which(trunk_footprints_cl$tf==tfscoret$tf[i])])
}
tfscoret2 <- tfscoret[splite,]
row_anno <- rowAnnotation(foo = anno_barplot(
  tfscoret2$score
))
ComplexHeatmap::Heatmap(avg_scaled,
                        heatmap_legend_param = list(title='scRNA-seq',
                                                    title_position='leftcenter-rot'),
                        cluster_rows = F,
                        cluster_columns = F,
                        show_column_names = T,
                        show_row_names = T,
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 6),
                        col = colorRamp2(seq(from=(-1.5),to=1.5,length=11),colorRampPalette(c("navy", "white", "firebrick3"))(11)),
                        column_names_rot = 90,
                        heatmap_width = unit(7, "cm"), 
                        heatmap_height = unit(5, "cm"),
                        column_names_centered = T,
                        row_split = factor(splite,levels = c("YBX1","ASCL2","HES4","HES1","HMGA2","ID4")),
                        row_title_gp = gpar(fontsize = 8),
                        row_title_rot = 0,
                        right_annotation = row_ha,
                        left_annotation = row_anno,
                        row_gap = unit(2, "mm")
)

#Fig.5g
endogenes <- marker_rna$gene[which(marker_rna$cluster%in%c("EP","Endocrine")&
                                     marker_rna$avg_log2FC>0&
                                     marker_rna$p_val_adj<0.05)]
ductgenes <- marker_rna$gene[which(marker_rna$cluster%in%c("Duct")&
                                     marker_rna$avg_log2FC>0&
                                     marker_rna$p_val_adj<0.05)]
data <- matrix(ncol = 5)
for(i in c(1:3)){
  y1 <- endo_Reg_motif_cl[which(endo_Reg_motif_cl$TFs %in%w1[i]),"Target"]
  y1 <- y1[y1%in%endogenes]
  print(w1[i])
  print(length(y1))
  
  avg <- AverageExpression(rna,assays = "RNA",group.by = "celltype",features = y1)
  avg <- avg[[1]]
  
  data <- rbind(data,avg)
  
  
}
for(i in c(1:9)){
  y1 <- duct_Reg_motif_cl[which(duct_Reg_motif_cl$TFs %in%w2[i]),"Target"]
  y1 <- y1[y1%in%ductgenes]
  print(w2[i])
  print(length(y1))
  avg <- AverageExpression(rna,assays = "RNA",group.by = "celltype",features = y1)
  avg <- avg[[1]]
  data <- rbind(data,avg)
}
data <- data[,c("Trunk","Duct","Endocrine")]
avg_scaled = t(scale(t(data)))
mark_gene <- c("PAX6","FEV","NEUROD1","CRYBA2","SLC37A1","HES4","MAFK",
               "ANXA2","YAP1",w1,w2)
gene_pos <- which(rownames(avg_scaled) %in% mark_gene)
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                                 labels = rownames(avg_scaled)[gene_pos],
                                                 labels_gp = gpar(fontsize = 6)))

ComplexHeatmap::Heatmap(avg_scaled,
                        heatmap_legend_param = list(title='scRNA-seq',
                                                    title_position='leftcenter-rot'),
                        cluster_rows = F,
                        cluster_columns = F,
                        show_column_names = T,
                        show_row_names = F,
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8),
                        col = colorRamp2(seq(from=-2,to=2,length=11),colorRampPalette(c("navy", "white", "firebrick3"))(11)),
                        column_names_rot = 0,
                        heatmap_width = unit(7, "cm"), 
                        heatmap_height = unit(15, "cm"),
                        column_names_centered = T,
                        row_title_gp = gpar(fontsize = 8),
                        row_title_rot = 0,
                        right_annotation = row_anno,
                        row_gap = unit(2, "mm"))

#Fig.5h
FeaturePlot(rna_all,features = c("YBX1","INSM1"),blend = T,blend.threshold = 1)
FeaturePlot(rna_all,features = c("INSM1","FEV"),blend = T,blend.threshold = 1)
