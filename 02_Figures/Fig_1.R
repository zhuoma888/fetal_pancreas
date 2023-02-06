library(Seurat)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(ComplexHeatmap)
library(monocle3)

#Load data
epi <- readRDS("./epi.rds")

###Fig.1
#Calculate DEG
marker <- FindAllMarkers(epi, logfc.threshold = 0.25, only.pos = T)
marker <- marker %>% filter(p_val_adj < 0.05)
marker$cluster <- factor(marker$cluster, levels = c("Epsilon","EP","Beta","Delta","Alpha",
                                                    "Duct","Trunk","Ventral MP","Dorsal MP","Early trunk","Early tip",
                                                    "Tip","Acinar" ) )

#Pseudotime and trajectory analysis
data <- GetAssayData(epi, assay = 'RNA', slot = 'counts')
cell_metadata <- epi@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
colnames(colData(cds))
cds <- reduce_dimension(cds, preprocess_method = "PCA")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(epi, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds,k = 20)
cds <- learn_graph(cds,use_partition = F, learn_graph_control = list(minimal_branch_len = 4))
cds <- order_cells(cds) #Choose root
#Set colors
type_col <- setNames(c("#1f77b4","#aec7e8","#fdd0a2","#ff7f0e","#e6550d",
                       "#c7e9c0","#74c476","#31a354","#17becf",
                       "#756bb1","#9e9ac8","#bcbddc","#dadaeb"),
                     c("Ventral MP","Dorsal MP","Early tip","Tip","Acinar",
                       "Early trunk","Trunk","Duct","EP",
                       "Beta","Alpha","Delta","Epsilon"))

#Fig. 1b
DimPlot(epi, reduction = "umap", pt.size = 0.5, group.by = "cell_type",shuffle = T) + scale_color_manual(values = type_col)
DimPlot(epi, reduction = "umap",group.by = "PCW") + scale_color_viridis_d()

#Fig. 1c
ncells <- table(epi$PCW, epi$cell_type)
per <- as.data.frame(ncells)
colnames(per) <- c("PCW","cell_type","freq")

ggplot(per,aes(x=PCW,y=freq,fill=cell_type))+
  geom_bar(stat = "identity",position = "fill") + scale_fill_manual(values = type_col) +
  theme_classic() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor =element_blank()) + ylab("Percentage")

#Fig. 1d
marker <- marker %>% arrange(cluster,-avg_log2FC)
top100 <- marker %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)

epi.list <- SplitObject(epi, split.by = "cell_type")
exp.matrix <- lapply(epi.list, function(x)  GetAssayData(x,slot = "data", assay = "RNA" ))
exp.mean <- lapply(exp.matrix, function(x) apply(x, 1, mean))
for(i in 1:length(exp.mean)) {
  exp.mean[[i]] <- as.data.frame(exp.mean[[i]])
  colnames(exp.mean[[i]]) <- names(exp.mean[i])
}
mean <- do.call(cbind, exp.mean)
mean <- mean[,c("Dorsal MP","Ventral MP","Early tip","Tip","Acinar",
                "Early trunk","Trunk","Duct","EP","Beta","Alpha","Delta","Epsilon")]

rowanno <- subset(marker, select = c(cluster,gene))
rowanno <- rowanno%>%arrange(cluster)
rownames(rowanno) <- paste0(rowanno$gene,"_",rowanno$cluster)

colors=list(cell_type=type_col)

marker_adj$cluster <- factor(marker_adj$cluster, levels = c("Dorsal MP","Ventral MP","Early tip","Tip","Acinar",
                                                            "Early trunk","Trunk","Duct","EP","Beta","Alpha","Delta","Epsilon"))
marker_adj <- marker_adj%>%arrange(cluster) 

mat <- as.matrix(mean[rowanno$gene,])
rownames(mat) <- rownames(rowanno)
mat <- t(mat)
mat <- scale(mat)
mat <- t(mat)

cell_type <- as.data.frame(unique(epi@meta.data$cell_type))
names(cell_type) <- "celltype"
rownames(cell_type) <- cell_type$celltype

cell_num <- as.data.frame(table(epi@meta.data$cell_type))
rownames(cell_num) <- cell_num$Var1
names(cell_num) <- c("celltype","cellnum")


colors <- list(celltype = type_col)
col_ha <- HeatmapAnnotation(celltype=cell_num$celltype, cellnum=anno_barplot(cell_num$cellnum,
                                                                             gp = gpar(fill = type_col)),
                            col = colors)

Heatmap(mat, 
        name = "z-score",
        col = colorRamp2(seq(from=-2,to=2,length=11),colorRampPalette(c("navy", "white", "firebrick3"))(11)),
        show_row_names = F,
        show_column_names = T,
        cluster_rows = F,
        cluster_columns = T,
        top_annotation = col_ha
)

#Fig.1e
FeaturePlot(epi, reduction = "umap", features = c("GATA4","RBPJL","DCDC2","NEUROD1"), 
            ncol=2, cols = c("grey","red")) 

#Fig.1f
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_roots = F,cell_size =0.5,
           label_leaves = F,  label_branch_points = F,trajectory_graph_color = "grey") +  scale_color_viridis_c()

