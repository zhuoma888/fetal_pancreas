library(Signac)
library(Seurat)
library(SeuratWrappers)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(harmony)
library(ArchR)
library(monocle3)
library(chromVARmotifs)
library(tidyverse)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggalluvial)

rna <- readRDS("/epi_RNA.Rds")
rna <- subset(rna,subset=PCW%in%c("W8","W9","W10","W11")&carnegie_stage!="CS21")
rna$celltype[which(rna$celltype%in%c("Alpha/PP","Beta","Delta","Epsilon"))] <- "Endocrine"
rna$celltype[which(rna$celltype=="Early tip")] <- "Tip"
rna$celltype[which(rna$celltype=="Early trunk")] <- "Trunk"
rna$lineage <- "endocrine"
rna$lineage[which(rna$celltype%in%c("Tip","Acinar"))] <- "acinar"
rna$lineage[which(rna$celltype%in%c("Trunk","Duct"))] <- "duct"
atac <- readRDS("/epi_ATAC.Rds")

#Figure 3B
DimPlot(epi,group.by = "celltype")+
  scale_color_manual(values = c("#ff7f0e","#e6550d","#74c476","#31a354","#17becf","#756bb1"))+
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank())
DimPlot(epi,group.by = "pcw")+
  scale_color_manual(values = hcl.colors(4))+
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank())

#Figure 3C
ncells_atac <- table(atac$celltype,atac$PCW)
freqs_atac <- prop.table(ncells_atac,margin = 2)
freqs_atac <- as.data.frame(freqs_atac)
colnames(freqs_atac) <- c("celltype","PCW","percentage")

ncells_rna <- table(rna$celltype,rna$PCW)
freqs_rna <- prop.table(ncells_rna,margin = 2)
freqs_rna <- as.data.frame(freqs_rna)
colnames(freqs_rna) <- c("celltype","PCW","percentage")
p1 <- ggplot(freqs_atac,aes(x=PCW,y=percentage,fill=celltype))+ 
  scale_fill_manual(values = c("#ff7f0e","#e6550d","#74c476","#31a354","#17becf","#756bb1"))+
  geom_bar(stat = "identity",position = "fill") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor =element_blank())+ 
  ylab("Percentage")+
  ggtitle("scATAC-seq")
p2 <- ggplot(freqs_rna,aes(x=PCW,y=percentage,fill=celltype))+ 
  scale_fill_manual(values = c("#ff7f0e","#e6550d","#74c476","#31a354","#17becf","#756bb1"))+
  geom_bar(stat = "identity",position = "fill") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor =element_blank())+ 
  ylab("Percentage")+
  ggtitle("scRNA-seq")

#Figure 3D
plot_cells(monoclecds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = F,label_roots = F)+
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank())+ 
  scale_color_viridis_c()

plot_cells(monoclecds, color_cells_by = "celltype", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = F,label_roots = F)+
  scale_color_manual(values = c("#ff7f0e","#e6550d","#74c476","#31a354","#17becf","#756bb1"))+
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank())

#Figure 3E
sigpeaks <- marker_peaks$gene[which(marker_peaks$gene%in%p2g$peaks & 
                                      marker_peaks$p_val_adj<0.05 &
                                      marker_peaks$pct.1>marker_peaks$pct.2)]
avg <- AverageExpression(epi,group.by = "celltype",features = sigpeaks,assays = "peaks")
avg <- avg[[1]][sigpeaks, ]
avg_scaled = t(scale(t(avg)))
color = c("#ff7f0e","#e6550d","#74c476","#31a354","#17becf","#756bb1")
names(color) <- c("Tip","Acinar","Trunk","Duct","EP","Endocrine")
top_anno <- HeatmapAnnotation(cluster=anno_block(gp=gpar(fill=color),
                                                 labels = c("Tip","Acinar","Trunk","Duct","EP","Endocrine"),
                                                 labels_gp = gpar(cex=0.5,color='white',fontsize = 18)))
ComplexHeatmap::Heatmap(avg_scaled,
                        heatmap_legend_param = list(title='peak accessibility',
                                                    title_position='leftcenter-rot'),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        show_column_names = T,
                        show_row_names = F,
                        col = paletteContinuous("solarExtra"),
                        column_title = "Cell-type-specific peaks",
                        column_title_gp = gpar(fontsize = 8, fontface = "bold"))

closestfeatures <- ClosestFeature(epi,regions = marker_peaks$gene)
closest <- intersect(unique(closestfeatures$gene_name),rownames(rna@assays$RNA@counts))
siggenes <- marker_rna$gene[which(marker_rna$gene%in%closest & 
                                    marker_rna$p_val_adj<0.05 &
                                    marker_rna$pct.1>marker_rna$pct.2)]
siggenes <- unique(siggenes)

avg <- AverageExpression(rna,group.by = "celltype",features = siggenes)
avg <- avg[[1]][siggenes, ]
avg_scaled = t(scale(t(avg)))
mark_gene <- c("CPA2","RBPJL","PTF1A","LGR4","XBP1","GP2","ASCL2",
               "CFTR","CTGF","ATP1A1","CHGA","NEUROD1","FEV","ZBTB7A")
gene_pos <- which(rownames(avg_scaled) %in% mark_gene)
anno <-  HeatmapAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                                 labels = rownames(avg_scaled)[gene_pos],
                                                 labels_gp = gpar(fontsize = 6)))

ComplexHeatmap::Heatmap(t(avg_scaled),
                        heatmap_legend_param = list(title='scRNA-seq',
                                                    #legend_direction = "horizontal",
                                                    title_position='leftcenter-rot'),
                        cluster_rows = F,
                        cluster_columns = F,
                        show_column_names = F,
                        show_row_names = T,
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8),
                        col = viridis(100),
                        column_names_rot = 0,
                        heatmap_width = unit(15, "cm"), 
                        heatmap_height = unit(7, "cm"),
                        column_names_centered = T,
                        row_title_gp = gpar(fontsize = 8),
                        row_title_rot = 0,
                        top_annotation = anno,
                        row_gap = unit(2, "mm")
                        
)