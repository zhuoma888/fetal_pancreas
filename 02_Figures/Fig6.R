library(Seurat)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(SCENIC)
library(SCopeLoomR)
library(Seurat)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(clusterProfiler)

#Load data
EP <- readRDS(file = "./EP.rds")

#Set colors
EP_type_col <- setNames(c("#c6dbef","#6baed6","#3182bd",
                           "#fdae6b","#e6550d","#636363",
                           "#a1d99b","#31a354","#756bb1"),
                         c("EP early","EP mid","EP late",
                           "EP alpha","Alpha/PP","Epsilon",
                           "EP beta","Beta","Delta") )

#Fig.6a
DimPlot(EP, reduction = "umap", pt.size = 0.5, group.by = "cell_type",shuffle = T) + scale_color_manual(values = EP_type_col)
DimPlot(EP, reduction = "umap",group.by = "PCW") + scale_color_viridis_d()

#Fig.6b
ncells <- table(EP$PCW, EP$cell_type)
per <- as.data.frame(ncells)
colnames(per) <- c("PCW","cell_type","freq")

ggplot(per,aes(x=PCW,y=freq,fill=cell_type))+
  geom_bar(stat = "identity",position = "fill") + scale_fill_manual(values = EP_type_col) +
  theme_classic() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor =element_blank()) + ylab("Percentage")

#Fig.6c
DotPlot(EP, features = c("SOX9","ID3","RBPJ","TEAD2","NEUROG3","HES6","FEV","PAX4","ARX",
                         "ISL1","ETV1","GCG","PPY","IRX2","GHRL","ONECUT2","INS","NKX6-1","MAFA","MAFB","MEIS2","PLAGL1",
                         "SST","HHEX","CHGA"), group.by = "cell_type") + theme_classic() + RotatedAxis() +
  scale_color_gradient2(high="red",mid = "lightgrey",low ="darkblue", midpoint = 0)

#Fig.6d
loom <- open_loom('./out_SCENIC_EP.loom')#output file from pySCENIC
 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 

cellInfo <- EP@meta.data[,c("barcode","cell_type")]
cellInfo$barcode <- NULL
n <- t(scale(t(getAUC(regulonAUC[,] ))))
AUC <- t(getAUC(regulonAUC[,]))

EP@meta.data$barcode <- rownames(EP@meta.data)
dim(AUC)
EP[["AUC"]] <- CreateAssayObject(t(AUC))
DefaultAssay(EP) <- "AUC"

da_regulon <- FindAllMarkers(EP, assay = "AUC", slot = "count", logfc.threshold = 0.05, only.pos = T)
da_regulon$cluster <- factor(da_regulon$cluster, levels = c("EP early","EP mid","EP late","Epsilon","Beta",
                                                            "Delta","EP alpha","EP beta","Alpha/PP"))
da_regulon <- da_regulon %>% arrange(cluster)
write.table(da_regulon, file = "./da_regulon_EP.txt", sep = "\t") # Differentially expressed regulons


#Fig.6f
data <- GetAssayData(EP, assay = 'RNA', slot = 'counts')
cell_metadata <- EP@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
colnames(colData(cds))
cds <- reduce_dimension(cds, preprocess_method = "PCA")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(EP, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds,k = 20)
cds <- learn_graph(cds,use_partition = F, learn_graph_control = list(minimal_branch_len = 4, ncenter =100))
cds <- order_cells(cds) #Set root

plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_roots = F,cell_size =0.5,
           label_leaves = F,  label_branch_points = F, scale_to_range = F) 

#Fig.6g
pseudotime <- as.data.frame(pseudotime(cds))
names(pseudotime) <- "pseudotime"
EP <- AddMetaData(EP, metadata = pseudotime)
cells <- EP@meta.data$barcode[which(EP@meta.data$cell_type %in% c("EP early","EP mid","EP late"))]

ep <- subset(EP, subset = cell_type %in% c("EP early","EP mid","EP late") )
genes <- unique(marker_adj$gene[which(marker_adj$cluster %in% c("EP early","EP mid","EP late"))])

ep@meta.data <- ep@meta.data %>% arrange(pseudotime)
View(ep@meta.data)
mat <- ep@assays$RNA@data[genes,rownames(ep@meta.data)]

mat <- t(apply(mat,1,function(x){smooth.spline(x,df=3)$y}))
mat <- t(apply(mat,1,function(x){(x-mean(x))/sd(x)}))

meta <- ep@meta.data[,c("cell_type","PCW")]

pcw <- meta
pcw$cell_type <- NULL
cell_type <- meta
cell_type$PCW <- NULL

col_anno <- HeatmapAnnotation(Cell_type=cell_type$cell_type,
                              col=list(Cell_type = EP_type_col))
set.seed(123)
htmp <- Heatmap(
  mat,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),colorRampPalette(c("navy", "white", "firebrick3"))(11)),
  show_row_names               = F,
  show_column_names            = F,
  km = 4,
  cluster_rows                 = T,
  cluster_row_slices           = F,
  cluster_columns              = F,
  top_annotation = col_anno )
dht <- draw(htmp)
clusters <- row_order(dht)

#Fig.6h
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

GO_df[[1]]$region <- "1"
GO_df[[2]]$region <- "2"
GO_df[[3]]$region <- "3"
GO_df[[4]]$region <- "4"
GO <- rbind(GO_df[[1]],GO_df[[2]],GO_df[[3]],GO_df[[4]])