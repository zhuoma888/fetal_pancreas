library(Seurat)
library(ggplot2)
library(ggsci)
library(Seurat)
library(tidyverse)
library(monocle3)
library(circlize)
library(ComplexHeatmap)
library(clusterProfiler)
library(Cellchat)
library(Signac)
library(SeuratWrappers)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(chromVARmotifs)
library(Matrix)
library(Matrix.utils)
library(ggpubr)
library(ggalluvial)

#Load data
epi <- readRDS("./epi.rds")
cds <- readRDS("./cds.rds")
pseudotime <- as.data.frame(pseudotime(cds))
names(pseudotime) <- "pseudotime"
epi <- AddMetaData(epi, metadata = pseudotime)

#Fig.4a
tip_trunk <- FindMarkers(epi, ident.1 = "Tip", ident.2 = "Trunk")
tip_trunk <- subset(epi, subset = cell_type %in% c("Tip","Trunk"))

epi.list <- SplitObject(tip_trunk, split.by = "cell_type")
exp.matrix <- lapply(epi.list, function(x)  GetAssayData(x,slot = "data", assay = "RNA" ))
exp.mean <- lapply(exp.matrix, function(x) apply(x, 1, mean))
for(i in 1:length(exp.mean)) {
  exp.mean[[i]] <- as.data.frame(exp.mean[[i]])
  colnames(exp.mean[[i]]) <- names(exp.mean[i])
}
mean <- do.call(cbind, exp.mean)

tip_trunk_gene <- FindAllMarkers(tip_trunk,logfc.threshold = 0.25,only.pos = T)
mean$gene <- rownames(mean)
tip_gene <- tip_trunk_gene$gene[which(tip_trunk_gene$cluster == "Tip")]
trunk_gene <- tip_trunk_gene$gene[which(tip_trunk_gene$cluster == "Trunk")]
mean$group <- "NS"
mean$group[which(mean$gene %in% tip_gene)] <- "Tip"
mean$group[which(mean$gene %in% trunk_gene)] <- "Trunk"
mean$label <- ""
mean$label[which(mean$group != "NS")] <- mean$gene
col <- setNames(c("#ff7f0e","#74c476","grey"),
                c("Tip","Trunk","NS"))

ggplot(mean,aes(x=Tip,y=Trunk,col=group)) + geom_point(size = 0.5) + 
  scale_color_manual(values = col) + theme_classic()

#Fig.4b
tip <- subset(epi, subset = cell_type %in%  c("Early tip","Tip","Acinar"))
trunk <- subset(epi, subset = cell_type %in%  c("Early trunk","Trunk","Duct") )
nrow(tip@meta.data)
nrow(trunk@meta.data)

marker <- FindAllMarkers(epi, logfc.threshold = 0.25, only.pos = T)
genes <- unique(marker$gene[which(marker$cluster %in% c("Early tip","Tip","Acinar",
                                                        "Early trunk","Trunk","Duct") &
                                    marker$avg_log2FC > 0.4)])

tip@meta.data <- tip@meta.data %>% arrange(-pseudotime)
trunk@meta.data <- trunk@meta.data %>% arrange(pseudotime)

mat1 <- tip@assays$RNA@data[genes,rownames(tip@meta.data)]
mat2 <- trunk@assays$RNA@data[genes,rownames(trunk@meta.data)]
mat <- cbind(mat1,mat2)
mat <- t(apply(mat,1,function(x){smooth.spline(x,df=3)$y}))
mat <- t(apply(mat,1,function(x){(x-mean(x))/sd(x)}))

meta1 <- tip@meta.data[,c("cell_type","PCW")]
meta2 <- trunk@meta.data[,c("cell_type","PCW")]
meta <- rbind(meta1,meta2)
meta$cell_type <- factor(meta$cell_type, levels = c("Early tip","Tip","Acinar",
                                                    "Early trunk","Trunk","Duct"))

type_col <- setNames(c("#fdd0a2","#ff7f0e","#e6550d",
                       "#c7e9c0","#74c476","#31a354"),
                     c("Early tip","Tip","Acinar",
                       "Early trunk","Trunk","Duct"))

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
  column_split = c(rep("tip",9793),rep("trunk",6036)) )
dht <- draw(htmp)
clusters <- row_order(dht)

#Fig.4c
gene <- vector("list",4)
gene[[1]] <- genes[clusters$`1`]
gene[[2]] <- genes[clusters$`2`]
gene[[3]] <- genes[clusters$`3`]
gene[[4]] <- genes[clusters$`4`]

gene_id <- lapply(gene, function(x) bitr(x, fromType="SYMBOL", toType="ENTREZID",
                                         OrgDb="org.Hs.eg.db") )

GO <- vector("list",4)
for (i in 1:length(gene)) {
  GO[[i]] <- enrichGO(gene_id[[i]]$ENTREZID,OrgDb=org.Hs.eg.db, ont='BP',pvalueCutoff=0.1,
                      keyType='ENTREZID', readable = T)
}

GO_df <- lapply(GO, function(x) as.data.frame(x))
GO_df[[1]]$group <- "1"
GO_df[[2]]$group <- "2"
GO_df[[3]]$group <- "3"
GO_df[[4]]$group <- "4"

GO_df <- rbind(GO_df[[1]],GO_df[[2]],GO_df[[3]],GO_df[[4]])

#Fig.4d
epi@meta.data["HES4"] <- epi@assays$RNA@data["HES4",]
epi@meta.data["HES1"] <- epi@assays$RNA@data["HES1",]
epi@meta.data["HEY1"] <- epi@assays$RNA@data["HEY1",]

exp <- epi@meta.data[,c("barcode","HES4","HES1","HEY1","cell_type")]
exp <- melt(exp, id.vars=c("barcode","cell_type"),variable.name="gene", value.name="Expression")
ggplot(exp_other, aes(x=cell_type, y=Expression, fill = cell_type)) + geom_boxplot(outlier.size = 0.1) + facet_wrap(~ gene, ncol = 1) +
  theme_classic() + scale_fill_manual(values = type_col) + theme(axis.text.x = element_text(angle = 45, hjust =0.5, vjust =0.5)) 

#Fig.4f
sup <- readRDS("./sup.rds") #Load supporting cells
new <- merge(epi, sup)
new <- subset(new, subset = cell_type %in% c("Acinar","Duct","Tip","Trunk","Mesothelial",
                                             "Pericyte","Fibroblast","Endothelial","Immune","Neural"))

new@meta.data$cell_type[which(new@meta.data$cell_type %in% c("Tip","Acinar"))] <- "Tip & Acinar"
new@meta.data$cell_type[which(new@meta.data$cell_type %in% c("Trunk","Duct"))] <- "Trunk & Duct"

#Cellchat
cellchat <- createCellChat(object = new, group.by = "cell_type")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat,raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)

all_col <- setNames(c("#1f77b4","#7f7f7f","#e377c2","#d62728","#9467bd","#8c564b",
                      "#e6550d","#31a354"),
                    c("Fibroblast","Mesothelial","Pericyte","Immune","Neural","Endothelial",
                      "Tip & Acinar","Trunk & Duct"))

netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4,5,6,7,8), targets.use = c(7,8), signaling = "NOTCH", legend.pos.x = 0, color.use = all_col)

netVisual_aggregate(cellchat, signaling = "NOTCH", layout = "circle", color.use = c("#8c564b","#1f77b4","#d62728","#7f7f7f","#9467bd","#e377c2","#e6550d","#31a354"),
                    sources.use = c(1,2,3,4,5,6), targets.use = c(7,8), edge.curved = 0.1)

#Fig.4g
pairLR <- extractEnrichedLR(cellchat, signaling = c("FGF","NT","HGF"), geneLR.return = FALSE)
pairLR <- pairLR[c(2,4,7,8,9),]
pairLR <- as.data.frame(pairLR)
names(pairLR) <- "interaction_name"

netVisual_bubble(cellchat, sources.use = c(1,2,3,4,5,6), targets.use = c(7,8), remove.isolate = FALSE, pairLR.use = pairLR) 


####scATAC
#Fig.4h
tfscoredl <- data.frame(tf=dltfs)
rownames(tfscoredl) <- tfscoredl$tf
#set footprints score of acinar lineage TFs to 0 for none footprints were enriched in duct lineages
tfscoredl$al_score <- 0
tfscoredl$dl_score <- 0
for (i in 1:nrow(tfscoredl)) {
  tfscoredl$dl_score[i] <- mean(dl_footprints_cl$center[which(dl_footprints_cl$tf==tfscoredl$tf[i])])
}
tfscoreal <- data.frame(tf=altfs)
rownames(tfscoreal) <- tfscoreal$tf
tfscoreal$al_score <- 0
tfscoreal$dl_score <- 0
for (i in 1:nrow(tfscoreal)) {
  tfscoreal$al_score[i] <- mean(al_footprints_cl$center[which(al_footprints_cl$tf==tfscoreal$tf[i])])
}
tfscore <- data.frame(tf=common)
rownames(tfscore) <- tfscore$tf
tfscore$al_score <- 0
tfscore$dl_score <- 0
for (i in 1:nrow(tfscore)) {
  tfscore$al_score[i] <- mean(al_footprints_cl$center[which(al_footprints_cl$tf==tfscore$tf[i])])
  tfscore$dl_score[i] <- mean(dl_footprints_cl$center[which(dl_footprints_cl$tf==tfscore$tf[i])])
}
tfscore2 <- rbind(tfscoreal,tfscoredl,tfscore)
tfscore2$acinar <- abs(tfscore2$al_score)
tfscore2$duct <- abs(tfscore2$dl_score)
tfscore2$specificity <- "common"
tfscore2$specificity[which(tfscore2$tf%in%dltfs)] <- "duct_specific"
tfscore2$specificity[which(tfscore2$tf%in%altfs)] <- "acinar_specific"
p1=ggscatter(data=tfscore2,x="acinar",y='duct',size=1,color = "specificity",
             palette = c("#ff7f0e","black","#74c476"),
             xlab="acinar lineage footprint score",ylab="ductal lineage footprint score",
             label =rownames(tfscore2),label.select = c(x,x2,common,"SMARCC2","ELK4","REST","ETS1","KLF10"),
             repel = T,font.label = c(10),label.rectangle = F)
p1

#Fig.4i
dl.peak <- da_peaks$gene[which(da_peaks$p_val_adj < 0.005 & da_peaks$cluster%in%c("Trunk","Duct"))]
al.peak <- da_peaks$gene[which(da_peaks$p_val_adj < 0.005 & da_peaks$cluster%in%c("Tip","Acinar"))]
dmotif <- FindMotifs(
  object = epi,
  features = dl.peak,
  background = al.peak
)
amotif <- FindMotifs(
  object = epi,
  features = al.peak,
  background = dl.peak
)
amotif <- amotif[order(amotif$fold.enrichment,decreasing = T),]
amotif <- amotif[!duplicated(amotif$motif.name),]
rownames(amotif) <- amotif$motif.name
amotif <- amotif[altfs,]
dmotif <- dmotif[order(dmotif$fold.enrichment,decreasing = T),]
dmotif <- dmotif[!duplicated(dmotif$motif.name),]
rownames(dmotif) <- dmotif$motif.name
dmotif <- dmotif[dltfs,]
aldata <- amotif[,c("motif.name","fold.enrichment")]
dldata <- dmotif[,c("motif.name","fold.enrichment")]

alrna <- subset(marker_rna,subset=cluster%in%c("Tip","Acinar"))
alrna <- subset(alrna,subset=gene%in%altfs)
alrna <- alrna[order(alrna$avg_log2FC,decreasing = T),]
alrna <- alrna[!duplicated(alrna$gene),]
rownames(alrna) <- alrna$gene
alrna <- alrna[rownames(aldata),]
aldata$rna <- alrna$avg_log2FC


dlrna <- subset(marker_rna,subset=cluster%in%c("Trunk","Duct"))
dlrna <- subset(dlrna,subset=gene%in%dltfs)
dlrna <- dlrna[order(dlrna$avg_log2FC,decreasing = T),]
dlrna <- dlrna[!duplicated(dlrna$gene),]
rownames(dlrna) <- dlrna$gene
dlrna <- dlrna[rownames(dldata),]
dldata$rna <- -dlrna$avg_log2FC

aldata$lineage <- "acinar lineage"
dldata$lineage <- "duct lineage"

data <- rbind(aldata,dldata)
colnames(data) <- c("gene","atac","rna")
data$threshold <- "no"
data$threshold[which(abs(data$rna)>=0.5)] <- "yes"
p1=ggscatter(data=data,x="rna",y='atac',color='celltype',size=1,
             palette=c("#ff7f0e","#74c476"),xlab="log2FC of gene expression",ylab="motif enrichment",
             label =rownames(data),label.select = data$gene[which(data$threshold=="yes")],
             repel = T,font.label = c(10),label.rectangle = F)
p1 <- p1 + theme(axis.line = element_line(linetype = "solid"))+xlim(-2,2)
p1=p1+geom_vline(xintercept=c((-0.5),0.5),lty=3,col="azure4",lwd=1)
p1=p1+geom_hline(yintercept =0.58,lty=3,col="azure4",lwd=1)  
p1 <- p1 + theme(legend.position = "right")
p1

#Fig.4j
rna <- AddModuleScore(rna,features = list(common),name = "common")
rna <- AddModuleScore(rna,features = list(altfs),name = "altf")
rna <- AddModuleScore(rna,features = list(dltfs),name = "dltf")
t <- rna@meta.data[,c("celltype","common1","altf1","dltf1")]
colnames(t) <- c("celltype","common","acinar_TF","ductal_TF")
t <- melt(t)
t
colnames(t) <- c("celltype","TF_specificity","value")
p <- ggplot(t)+
  geom_boxplot(aes(x=celltype,y=value,fill=TF_specificity),outlier.size = 0.1)+
  scale_fill_manual(values = c("grey","#ff7f0e","#74c476"))
p

#Fig.4k
sigal <- c("PTF1A","JUNB","JUND","XBP1","STAT3","MYC","EPAS1","ONECUT1","CEBPD")
sigdl <- c("GLIS3","HES1","EHF","NFIB","HMGA2","HEY1","ID4","ASCL2","HES4","TCF4")

split <- character()
datad <- matrix(ncol = 3)
for(i in 1:length(sigdl)){
  y <- dl_Reg_motif_cl[which(dl_Reg_motif_cl$TFs %in%sigdl[i]),"Target"]
  y <- y[y%in%dl_genes]
  print(x[i])
  print(length(y))
  split <- c(split,rep(as.character(sigdl[i]),length(y)))
  avg <- AverageExpression(rna,assays = "RNA",group.by = "lineage",features = y)
  avg <- avg[[1]]
  datad <- rbind(datad,avg)
}
datad <- as.data.frame(datad)
datad$gene <- rownames(datad)
datad$tf <- split
datad$tfclass <- "duct_specific"
datad$targetclass <- "duct"

split <- character()
datad <- matrix(ncol = 3)
for(i in 1:length(sigal)){
  y <- al_Reg_motif_cl[which(al_Reg_motif_cl$TFs %in%sigal[i]),"Target"]
  y <- y[y%in%al_genes]
  split <- c(split,rep(as.character(sigal[i]),length(y)))
  avg <- AverageExpression(rna,assays = "RNA",group.by = "lineage",features = y)
  avg <- avg[[1]]
  dataa <- rbind(dataa,avg)
}
dataa <- as.data.frame(dataa)
dataa$gene <- rownames(dataa)
dataa$tf <- split
dataa$tfclass <- "acinar_specific"
dataa$targetclass <- "acinar"

split <- character()
datac <- matrix(ncol = 3)
for(i in 1:length(common)){
  y1 <- al_Reg_motif_cl[which(al_Reg_motif_cl$TFs %in%common[i]),"Target"]
  y2 <- dl_Reg_motif_cl[which(dl_Reg_motif_cl$TFs %in%common[i]),"Target"]
  y1 <- y1[y1%in%al_genes]
  y2 <- y2[y2%in%dl_genes]
  y <- c(y1,y2)
  y <- unique(y)
  split <- c(split,rep(as.character(common[i]),length(y)))
  avg <- AverageExpression(rna,assays = "RNA",group.by = "lineage",features = y)
  avg <- avg[[1]]
  datac <- rbind(datac,avg)
}
datac <- as.data.frame(datac)
datac$gene <- rownames(datac)
datac$tf <- split
datac$tfclass <- "common"
datac$targetclass <- "duct"
datac$targetclass[which(datac$gene%in%al_genes)] <- "acinar"
data <- rbind(dataa,datad)
data <- rbind(data,datac)
p <- ggplot(data = data, aes(axis1 = tfclass, axis2 = tf, axis3 = targetclass)) +
  scale_x_discrete(limits = c("tfclass","tf" ,"targetclass"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = targetclass)) +
  scale_fill_manual(values = c("#ff7f0e","#74c476"))+
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            color = "black",check_overlap = T,size=2)+
  theme(element_blank())

