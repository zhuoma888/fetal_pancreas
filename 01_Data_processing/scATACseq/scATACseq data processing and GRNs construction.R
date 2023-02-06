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

#quality control
counts <- Read10X_h5(filename = "/cellranger-atac/fetal_panc_atac/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "/cellranger-atac/fetal_panc_atac/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '/cellranger-atac/fetal_panc_atac/outs/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
samples <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(samples) <- annotations
samples <- NucleosomeSignal(samples)
samples <- TSSEnrichment(samples, fast = FALSE)
samples$blacklist_fraction <- FractionCountsInRegion(
  object = samples, 
  assay = 'peaks',
  regions = blacklist_hg38
)
#choose high-quality samples and then perform quality control
samples <- subset(
  x = samples,
  subset = orig.ident %in% c("S14","S35","S96","S98") &
    blacklist_fraction < 0.02 &
    nucleosome_signal < 4 &
    TSS.enrichment > 4 &
    TSS.enrichment <20 &
    total < 50000 
)
samples <- SplitObject(samples,split.by = "orig.ident")


#build windows assay to identify cell classes and extract epithelial cells
#References
#Granja,J.M.et al.Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia. Nat Biotechnol 37, 1458-1465 (2019).
#https://github.com/GreenleafLab/MPAL-Single-Cell-2019
genome <- BSgenome.Hsapiens.UCSC.hg38
chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
windows <- unlist(tile(chromSizes, width = 2500))
windows <- subsetByOverlaps(windows,blacklist_hg38,invert = T)
samples[["S14"]]$barcode <- gsub("2","1",rownames(samples[["S14"]]@meta.data))
samples[["S35"]]$barcode <- gsub("3","1",rownames(samples[["S35"]]@meta.data))
samples[["S96"]]$barcode <- gsub("4","1",rownames(samples[["S96"]]@meta.data))
samples[["S98"]]$barcode <- gsub("5","1",rownames(samples[["S98"]]@meta.data))

frag14 <- CreateFragmentObject("/cellranger-atac/S14/outs/fragments.tsv.gz",
                               cells = samples[["S14"]]$barcode)
frag35 <- CreateFragmentObject("/cellranger-atac/S35/outs/fragments.tsv.gz",
                               cells = samples[["S35"]]$barcode$barcode)
frag96 <- CreateFragmentObject("/cellranger-atac/S96/outs/fragments.tsv.gz",
                               cells = samples[["S96"]]$barcode$barcode)
frag98 <- CreateFragmentObject("/cellranger-atac/S98/outs/fragments.tsv.gz",
                               cells = samples[["S98"]]$barcode$barcode)
windows14 <- FeatureMatrix(
  fragments=frag14,
  features = windows,
  cells = Cells(frag14),
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
)
windows35 <- FeatureMatrix(
  fragments=frag35,
  features = windows,
  cells = Cells(frag35),
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
)
windows96 <- FeatureMatrix(
  fragments=frag96,
  features = windows,
  cells = Cells(frag96),
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
)
windows98 <- FeatureMatrix(
  fragments=frag98,
  features = windows,
  cells = Cells(frag98),
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
)
colnames(windows14) <- gsub("1","2",colnames(windows14))
colnames(windows35) <- gsub("1","3",colnames(windows35))
colnames(windows96) <- gsub("1","4",colnames(windows96))
colnames(windows98) <- gsub("1","5",colnames(windows98))

whole <- merge.Matrix(windows14,windows35,
                      by.x = rownames(windows14),
                      by.y = rownames(windows35))
whole <- merge.Matrix(whole,windows96,
                      by.x = rownames(whole),
                      by.y = rownames(windows96))
whole <- merge.Matrix(whole,windows98,
                      by.x = rownames(whole),
                      by.y = rownames(windows98))

chrom_assay <- CreateChromatinAssay(
  counts = whole,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = list(frag14,frag35,frag96,frag98),
  min.cells = 10,
  min.features = 200
)
meta <- rbind(samples[["S14"]]@meta.data,samples[["S35"]]@meta.data)
meta <- rbind(meta,samples[["S96"]]@meta.data)
meta <- rbind(meta,samples[["S98"]]@meta.data)
panc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "windows",
  meta.data = meta
)
panc <- BinarizeCounts(panc)
panc <- RunTFIDF(panc,method = 3)
panc <- RunSVD(panc)
panc <- RunHarmony(
  object = panc,
  group.by.vars = 'orig.ident',
  reduction = 'lsi',
  assay.use = 'windows',
  project.dim = F
)
panc <- RunUMAP(panc, dims = 2:30, reduction = 'harmony',
                n.neighbors = 50L,min.dist = 0.5)
gene.activities <- GeneActivity(panc)
panc[['RNA']] <- CreateAssayObject(counts = gene.activities)

#remove clusters coexpressing markers of more than one cell classes or expressing none of the markers of interested cell classes
panc <- subset(panc,subset=seurat_clusters!="0")
panc <- subset(panc,subset=seurat_clusters!="15")
panc$cellclass <- ""
panc$cellclass[which(panc$seurat_clusters %in% c('1','2','13'))] <- "Non-endocrine"
panc$cellclass[which(panc$seurat_clusters %in% c('3','4','7'))] <- "Mesenchymal"
panc$cellclass[which(panc$seurat_clusters %in% c('9','11'))] <- "Neural"
panc$cellclass[which(panc$seurat_clusters == '10')] <- "Immune"
panc$cellclass[which(panc$seurat_clusters == '8')] <- "Endocrine"
panc$cellclass[which(panc$seurat_clusters == '14')] <- "Erythroid"
panc$cellclass[which(panc$seurat_clusters == '6')] <- "Endothelial"

epi <- subset(panc,subset=cellclass %in% c("Non-endocrine","Endocrine"))

#call uniform peaks and build peaks to identify cell types
epi <- FindNeighbors(epi, reduction = "harmony", dims = 2:30)
epi <- FindClusters(epi, verbose = T, resolution = 0.8)
peaks <- CallPeaks(epi,assay = "windows",group.by = "seurat_clusters",
                   macs2.path = "/home/hupk/.local/bin/macs2",
                   broad = FALSE,format = "BED",
                   outdir = tempdir(),fragment.tempdir = tempdir(),
                   combine.peaks = F,effective.genome.size = 2.7e+09,
                   extsize = as.character(150),shift = as.character(-150/2) ,
                   additional.args = "--call-summits --nolambda --qval 5e-2 --keep-dup all" ,
                   name = Project(object),
                   cleanup = F,verbose = TRUE)
extendedPeakSet <- function(df, BSgenome = NULL, extend = 250, blacklist = NULL, nSummits = 100000){
  #Helper Functions
  readSummits <- function(file){
    df <- suppressMessages(data.frame(readr::read_tsv(file, col_names = c("chr","start","end","name","score"))))
    df <- df[,c(1,2,3,5)] #do not keep name column it can make the size really large
    return(GenomicRanges::makeGRangesFromDataFrame(df=df,keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE))
  }
  nonOverlappingGRanges <- function(gr, by = "score", decreasing = TRUE, verbose = FALSE){
    stopifnot(by %in% colnames(mcols(gr)))
    clusterGRanges <- function(gr, filter = TRUE, by = "score", decreasing = TRUE){
      gr <- sort(sortSeqlevels(gr))
      r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
      o <- findOverlaps(gr,r)
      mcols(gr)$cluster <- subjectHits(o)
      gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
      gr <- gr[!duplicated(mcols(gr)$cluster),]
      gr <- sort(sortSeqlevels(gr))
      mcols(gr)$cluster <- NULL
      return(gr)
    }
    if(verbose){
      message("Converging", appendLF = FALSE)
    }
    i <-  0
    gr_converge <- gr
    while(length(gr_converge) > 0){
      if(verbose){
        message(".", appendLF = FALSE)
      }
      i <-  i + 1
      gr_selected <- clusterGRanges(gr = gr_converge, filter = TRUE, by = by, decreasing = decreasing)
      gr_converge <- subsetByOverlaps(gr_converge ,gr_selected, invert=TRUE) #blacklist selected gr
      if(i == 1){ #if i=1 then set gr_all to clustered
        gr_all <- gr_selected
      }else{
        gr_all <- c(gr_all, gr_selected)
      } 
    }
    if(verbose){
      message("\nSelected ", length(gr_all), " from ", length(gr))
    }
    gr_all <- sort(sortSeqlevels(gr_all))
    return(gr_all)
  }
  #Check-------
  stopifnot(extend > 0)
  stopifnot("samples" %in% colnames(df))
  stopifnot("groups" %in% colnames(df))
  stopifnot("summits" %in% colnames(df))
  stopifnot(!is.null(BSgenome))
  stopifnot(all(apply(df,1,function(x){file.exists(paste0(x[3]))})))
  #------------
  #Deal with blacklist
  if(is.null(blacklist)){
    blacklist <- GRanges()
  }else if(is.character(blacklist)){
    blacklist <- rtracklayer::import.bed(blacklist)
  }
  stopifnot(inherits(blacklist,"GenomicRanges"))
  #------------
  #Time to do stuff
  chromSizes <- GRanges(names(seqlengths(BSgenome)), IRanges(1, seqlengths(BSgenome)))
  chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
  groups <- unique(df$groups)
  groupGRList <- GenomicRanges::GRangesList(lapply(seq_along(groups), function(i){
    message("groupno.1")
    df_group = df[which(df$groups==groups[i]),]
    grList <- GenomicRanges::GRangesList(lapply(paste0(df_group$summits), function(x){
      message("no.2")
      extended_summits <- readSummits(x) %>%
        resize(., width = 2 * extend + 1, fix = "center") %>%     
        subsetByOverlaps(.,chromSizes,type="within") %>%
        subsetByOverlaps(.,blacklist,invert=TRUE) %>%
        nonOverlappingGRanges(., by="score", decreasing=TRUE)
      extended_summits <- extended_summits[order(extended_summits$score,decreasing=TRUE)]
      if(!is.null(nSummits)){
        extended_summits <- head(extended_summits, nSummits)
      }
      mcols(extended_summits)$scoreQuantile <- trunc(rank(mcols(extended_summits)$score))/length(mcols(extended_summits)$score)
      extended_summits
    }))
    #Non Overlapping
    grNonOverlapping <- nonOverlappingGRanges(unlist(grList), by = "scoreQuantile", decreasing = TRUE)
    #Free Up Memory
    remove(grList)
    gc()
    grNonOverlapping
  }))
  grFinal <- nonOverlappingGRanges(unlist(groupGRList), by = "scoreQuantile", decreasing = TRUE)
  grFinal <- sort(sortSeqlevels(grFinal))
  return(grFinal)
}
df <- data.frame(
  samples = gsub("\\_summits.bed","",list.files(tempdir(), pattern = "\\_summits.bed", full.names = FALSE)),
  groups = "scATAC",
  summits = list.files(tempdir(), pattern = "\\_summits.bed", full.names = TRUE)
)
genome <- BSgenome.Hsapiens.UCSC.hg38
unionPeaks <- extendedPeakSet(
  df = df,
  BSgenome = genome, 
  extend = 250,
  blacklist = blacklist_hg38,
  nSummits = 200000
)
unionPeaks <- unionPeaks[seqnames(unionPeaks) %in% paste0("chr",c(1:22,"X"))]
unionPeaks <- keepSeqlevels(unionPeaks, paste0("chr",c(1:22,"X")))
epi$barcode <- rownames(epi@meta.data)
samples[["S14"]]$barcode_sample <- rownames(samples[["S14"]]@meta.data)
samples[["S35"]]$barcode_sample <- rownames(samples[["S35"]]@meta.data)
samples[["S96"]]$barcode_sample <- rownames(samples[["S96"]]@meta.data)
samples[["S98"]]$barcode_sample <- rownames(samples[["S98"]]@meta.data)

samples[["S14"]] <- subset(samples[["S14"]],subset=barcode_sample%in%epi$barcode)
samples[["S35"]] <- subset(samples[["S35"]],subset=barcode_sample%in%epi$barcode)
samples[["S96"]] <- subset(samples[["S96"]],subset=barcode_sample%in%epi$barcode)
samples[["S98"]] <- subset(samples[["S98"]],subset=barcode_sample%in%epi$barcode)
frag14 <- CreateFragmentObject("/data/cellranger-atac/S14/outs/fragments.tsv.gz",
                               cells = samples[["S14"]]$barcode)
frag35 <- CreateFragmentObject("/data/cellranger-atac/S35/outs/fragments.tsv.gz",
                               cells = samples[["S35"]]$barcode)
frag96 <- CreateFragmentObject("/data/cellranger-atac/S96/outs/fragments.tsv.gz",
                               cells = samples[["S96"]]$barcode)
frag98 <- CreateFragmentObject("/data/cellranger-atac/S98/outs/fragments.tsv.gz",
                               cells = samples[["S98"]]$barcode)
peak14 <- FeatureMatrix(
  fragments=frag14,
  features = unionPeaks,
  cells = Cells(frag14),
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
)
peak35 <- FeatureMatrix(
  fragments=frag35,
  features = unionPeaks,
  cells = Cells(frag35),
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
)
peak96 <- FeatureMatrix(
  fragments=frag96,
  features = unionPeaks,
  cells = Cells(frag96),
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
)
peak98 <- FeatureMatrix(
  fragments=frag98,
  features = unionPeaks,
  cells = Cells(frag98),
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
)
colnames(peak14) <- gsub("1","2",colnames(peak14))
colnames(peak35) <- gsub("1","3",colnames(peak35))
colnames(peak96) <- gsub("1","4",colnames(peak96))
colnames(peak98) <- gsub("1","5",colnames(peak98))
whole <- merge.Matrix(peak14,peak35,
                      by.x = rownames(peak14),
                      by.y = rownames(peak35))
whole <- merge.Matrix(whole,peak96,
                      by.x=rownames(whole),
                      by.y=rownames(peak96))
whole <- merge.Matrix(whole,peak98,
                      by.x=rownames(whole),
                      by.y=rownames(peak98))
chrom_assay <- CreateChromatinAssay(
  counts = whole,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = list(frag14,frag35,frag96,frag98),
  min.cells = 10,
  min.features = 200
)
epi[["peaks"]] <- chrom_assay
epi <- BinarizeCounts(epi,assay = "peaks500")
epi <- RunTFIDF(epi,method = 3,scale.factor = 1e4)
epi <- FindTopFeatures(epi)
epi <- RunSVD(epi)
epi <- RunHarmony(
  object = epi,
  group.by.vars = 'batch',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)
epi <- RunUMAP(epi,reduction = "harmony",dims = 2:25,
               n.neighbors = 50L,min.dist = 0.5)
epi <- FindNeighbors(epi, reduction = 'harmony', dims = 2:30)
epi <- FindClusters(epi, verbose = T,resolution = 1.5)
epi$celltype <- "no"
epi$celltype[which(epi$seurat_clusters=="13")] <- "EP"
epi$celltype[which(epi$seurat_clusters%in%c("12","14"))] <- "Endocrine"
epi$celltype[which(epi$seurat_clusters%in%c("0","6","15"))] <- "Trunk"
epi$celltype[which(epi$seurat_clusters%in%c("1","7"))] <- "Duct"
epi$celltype[which(epi$seurat_clusters%in%c("3","4","5","8","9","10","11"))] <- "Tip"
epi$celltype[which(epi$seurat_clusters%in%c("2"))] <- "Acinar"

#motif enrichment
epi <- AddMotifs(
  object = epi,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = human_pwms_v2
)

#pseudotime analysis
monoclecds <- as.cell_data_set(epi)
monoclecds <- cluster_cells(monoclecds, reduction_method = "UMAP")
monoclecds <- learn_graph(monoclecds, use_partition = TRUE, close_loop = F,learn_graph_control = list(ncenter=330,minimal_branch_len=15,geodesic_distance_ratio=1/5),verbose = FALSE)
monoclecds <- order_cells(monoclecds, reduction_method = "UMAP")

#Build a ArchR project and put above Signac results into ArchR
epiarch <- ArchRProject(ArrowFiles = c("./S14.arrow",
                                       "./S35.arrow",
                                       "./S96.arrow",
                                       "./S98.arrow"),
                        copyArrows = F,
                        outputDirectory = "./")
epiarch <- addIterativeLSI(
  ArchRProj = epiarch,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

epiarch <- addClusters(
  input = epiarch,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)
epiarch <- addHarmony(
  ArchRProj = epiarch,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = T
)

epiarch <- addUMAP(
  ArchRProj = epiarch, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = T
)

row <- rownames(epiarch@reducedDims$Harmony$matDR)
cellnames <- data.frame(row.names = row,
                        archr=row,
                        sample=substr(row,1,3))
cellnames$signac <- substr(cellnames$archr,5,22)
cellnames <- subset(cellnames,subset=signac%in%epi$barcode)
epiarch2 <- subsetArchRProject(epiarch,cells = cellnames$archr,
                               outputDirectory = "./new/",force = T)

epiarch2 <- addUMAP(
  ArchRProj = epiarch2, 
  reducedDims = "Harmony", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = T
)

harmony <- epi@reductions$harmony@cell.embeddings
harmony <- as.data.frame(harmony)
harmony <- harmony[cellnames$signac,]
colnames(epiarch2@reducedDims$Harmony$matDR)
colnames(harmony)
harmony <- harmony[,1:30]
rownames(harmony) <- cellnames$archr
colnames(harmony) <- colnames(epiarch2@reducedDims$Harmony$matDR)
harmony.m <- as.matrix(harmony)
epiarch2@reducedDims$Harmony$matDR <- harmony.m
umap <- epi@reductions$umap@cell.embeddings
umap <- as.data.frame(umap)
umap <- umap[cellnames$signac,]
rownames(umap) <- cellnames$archr
colnames(umap)<- colnames(epiarch2@embeddings$UMAP$df)
epiarch2@embeddings$UMAP$df <- umap
meta <- epi@meta.data[,c("seurat_clusters","peak_region_fragments","celltype")]
meta <- meta[cellnames$signac,]
epiarch2$celltype <- meta$celltype
epiarch2$peak_region_fragments <- meta$peak_region_fragments
epiarch2$Clusters <- meta$seurat_clusters
p1 <- plotEmbedding(epiarch2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(epiarch2, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
plotPDF(p1,p2, name = "umap.pdf", ArchRProj = epiarch2, addDOC = FALSE, width = 5, height = 5)

#transport RNA assay in Signac to ArchR project
rna_signac <- as.SingleCellExperiment(epi,assay = "RNA")
rna_signac@assays@data[["logcounts"]] <- NULL
rna_signac <- as(rna_signac,"RangedSummarizedExperiment")
epiarch2 <- addGeneExpressionMatrix(
  input = epiarch2,
  seRNA = rna_signac,
  chromSizes = getChromSizes(epiarch2),
  excludeChr = c("chrM", "chrY"),
  scaleTo = 10000,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addGeneExpressionMatrix")
)

seRNA <- readRDS("/epi_RNA.Rds")
seRNA <- subset(seRNA,subset=PCW%in%c("W8","W9","W10","W11")&carnegie_stage!="CS21")
seRNA$celltype[which(seRNA$celltype%in%c("Alpha/PP","Beta","Delta","Epsilon"))] <- "Endocrine"
seRNA$celltype[which(seRNA$celltype=="Early tip")] <- "Tip"
seRNA$celltype[which(seRNA$celltype=="Early trunk")] <- "Trunk"

epiarch2 <- addGeneIntegrationMatrix(
  ArchRProj = epiarch2, 
  useMatrix = "GeneExpressionMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = seRNA,
  addToArrow = T,
  groupRNA = "celltype",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  force = T,
  threads = 1
)
epiarch2 <- addPeak2GeneLinks(
  ArchRProj = epiarch2,
  reducedDims = "Harmony",
  dimsToUse = 2:30
)

#get gene regulatory networks using IReNA2
#References
#Lyu, P. et al. Gene regulatory networks controlling temporal patterning, neurogenesis, and cell-fate specification in mammalian retina. Cell Rep 37, 109994 (2021)
#following https://github.com/Pinlyu3/IReNA-v2
Get_p2g_fun <- function(x){
  corCutOff = 0.20
  FDRCutOff = 1e-6
  varCutOffATAC = 0.7
  varCutOffRNA = 0.3
  p2g <- metadata(x@peakSet)$Peak2GeneLinks
  p2g <- p2g[which(abs(p2g$Correlation) >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]
  if(!is.null(varCutOffATAC)){
    p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
  }
  if(!is.null(varCutOffRNA)){
    p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
  }
  mATAC <- readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ]
  mRNA <- readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ]
  p2g$peak <- paste0(rowRanges(mATAC))
  p2g$gene <- rowData(mRNA)$name
  return(p2g)
}
p2g = Get_p2g_fun(epiarch2)
allgenes <- unique(markers$gene)
p2g2 <- subset(p2g,subset=gene%in%allgenes)
Selection_peaks_for_one <- function(All_peaks_list,All_genes,p2g,distance_F=100000,mm10_TSS_GR_all){
  All_genes = All_genes
  #### first to find the TSS Peaks ##########
  Genes1 = All_genes
  TSS_peak_list = All_peaks_list$TSS
  TSS_peak_list$peaks = as.character(TSS_peak_list)
  #### filter out TSS peaks which are not in the enriched genes' TSS ####
  Genes1_TSS1 = TSS_peak_list[which(TSS_peak_list$gene_name %in% Genes1 == T)]
  m = match(Genes1_TSS1$gene_name,Genes1)
  Genes1_TSS1 = Genes1_TSS1[order(m)]
  #### Second to find the Body peaks ######
  ##### Gene Body list ##############
  Body_peak_list = All_peaks_list$GeneBody
  Body_peak_list$peaks = as.character(Body_peak_list)####
  ### need p2g correlation ###############
  p2g$G_P = paste(p2g$gene,p2g$peak,sep='@')
  Body_peak_list$G_P = paste(Body_peak_list$gene_name,Body_peak_list$peaks,sep='@')
  #### filter the body peaks, add the p2g correlations #####
  k = which(Body_peak_list$G_P %in% p2g$G_P == T)
  Body_peak_list_cl = Body_peak_list[k]
  m = match(Body_peak_list_cl$G_P,p2g$G_P)
  Body_peak_list_cl$Correlation = p2g$Correlation[m]
  ##### filter out Body peaks which are not in the enriched genes' TSS#####
  Genes1_Body1 = Body_peak_list_cl[which(Body_peak_list_cl$gene_name %in% Genes1 == T)]
  m = match(Genes1_Body1$gene_name,Genes1)
  Genes1_Body1 = Genes1_Body1[order(m)]
  ##### print(Genes1_Body1$gene_name[!duplicated(Genes1_Body1$gene_name)]) #####
  ###### third to find the integenetic peaks ################
  Inter_peak_list = All_peaks_list$Intergenic
  Inter_peak_list$peaks = as.character(Inter_peak_list)
  #### first filter p2g data.frame, keep the p2g with enriched genes and Intergenic peaks ######
  p2g_cl = p2g[which(p2g$gene %in% c(Genes1) == T),]
  p2g_clcl = p2g_cl[which(p2g_cl$peak %in% Inter_peak_list$peaks == T),]
  ### p2g to GRanges #######
  p2g_clcl_GR = GRanges(p2g_clcl$peak)
  p2g_clcl_GR$peaks = p2g_clcl$peak
  p2g_clcl_GR$gene_name = p2g_clcl$gene
  p2g_clcl_GR$Correlation = p2g_clcl$Correlation
  #### find the TSS region for each gene_name ###########
  mm10_TSS_GR_all$peaks = as.character(mm10_TSS_GR_all)
  ####
  m = match(p2g_clcl_GR$gene_name,mm10_TSS_GR_all$gene_name)
  p2g_clcl_GR$TSS_region = mm10_TSS_GR_all$peaks[m]
  p2g_clcl_GR$distance = distance(GRanges(p2g_clcl_GR$peaks),GRanges(p2g_clcl_GR$TSS_region))
  #### filter intergenic peaks with distance to TSS ####
  p2g_clcl_GR = p2g_clcl_GR[which(p2g_clcl_GR$distance < distance_F)]
  ########  #####
  Gene1_inter1 = p2g_clcl_GR[which(p2g_clcl_GR$gene_name %in% Genes1 == T)]
  m = match(Gene1_inter1$gene_name,Genes1)
  Gene1_inter1 = Gene1_inter1[order(m)] #######
  ###### print(Gene1_inter1$gene_name[!duplicated(Gene1_inter1$gene_name)])
  #### Output the list !!!! ########
  Out_list = list(Genes1_TSS1,Genes1_Body1,Gene1_inter1)
  names(Out_list) <- c('Gene1_TSS','Gene1_Body','Gene1_Inter')
  return(Out_list)
}
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
pos <- subset(annotations, strand == "+")
pos <- pos[order(start(pos)),] 
pos <- pos[!duplicated(pos$gene_name),] # remove all but the first exons per transcript
pos$st <- start(pos)
pos$en <- pos$st + 1
neg <- subset(annotations, strand == "-")
neg <- neg[order(start(neg), decreasing = TRUE),] 
neg <- neg[!duplicated(neg$gene_name),] # remove all but the first exons per transcript
neg$en <- end(neg)
neg$st <- neg$en-1
gene_annotation_sub <- c(pos, neg)
gene_annotation_sub$gene <- gene_annotation_sub$gene_name
gene_anno_df <- DataFrame(row.names = rownames(gene_annotation_sub),
                          seqnames=gene_annotation_sub@seqnames,
                          start=gene_annotation_sub$st,
                          end=gene_annotation_sub$en,
                          gene=gene_annotation_sub$gene_name)
gene_anno_df <- as.data.frame(gene_anno_df)
TSS_GR_all <- makeGRangesFromDataFrame(gene_anno_df,keep.extra.columns = T,ignore.strand = T,
                                    seqnames.field = "seqnames",start.field = "start",end.field = "end")
peak_gene_list = Selection_peaks_for_one(peaklist,allgenes,p2g2,distance_F=100000,TSS_GR_all)
cellnames_cl <- cellnames[which(cellnames$signac %in% epi$barcode[which(epi$celltype%in%c("Tip","Acinar","Trunk","Duct","EP","endocrine"))]),]
epi_cl <- subset(epi,subset=barcode%in%cellnames_cl$signac)
new <- CreateFragmentObject("/fetal_panc_atac/outs/fragments.tsv.gz",
                            cells = cellnames_cl$signac)
chrom_assay <- CreateChromatinAssay(
  counts = epi_cl@assays$peaks@counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments =new,
  min.cells = 10,
  min.features = 200
)
forfrag <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = epi_cl@meta.data
)
SplitFragments(
  forfrag,
  assay = "peaks",
  group.by = "lineage",
  outdir = getwd(),
  file.suffix = "",
  append = TRUE,
  buffer_length = 256L,
  verbose = TRUE
)
Read_fragment_to_GR <- function(x){
  tab <- read_tsv(x,col_names = F)
  GR <- GRanges(seqnames=tab$X1,ranges=IRanges(start=tab$X2,end=tab$X3),index=tab$X4,dup=tab$X5)
  return(GR)
}
al_tab_GR <- Read_fragment_to_GR("./acinar_lineage.bed")
Fragment_list <- list()
Fragment_list[[1]] <- al_tab_GR
names(Fragment_list) <- c("al")
Output_to_bed_files = function(Fragment_list_cl,folder){
  tags = names(Fragment_list_cl)
  ### to bam files ###
  fragments_cl_bamGR_list = list()
  for(i in 1:length(Fragment_list_cl)){
    print(i)
    temp_Fragment = Fragment_list_cl[[i]]
    temp_Fragment$Reads = paste(tags[i],1:length(temp_Fragment),sep='_')
    temp_Fragment_bamGR = fragments_to_bam_GR(temp_Fragment)
    fragments_cl_bamGR_list = c(fragments_cl_bamGR_list,list(temp_Fragment_bamGR))
  }
  ###
  for(i in 1:length(fragments_cl_bamGR_list)){
    print(i)
    library('readr')
    library(tibble)
    FN = paste(tags[i],"fragments_cl_bamGR_pe.bed",sep='_')
    print(FN)
    fragments_cl_bamGR_tab = Convert_to_bedpe(fragments_cl_bamGR_list[[i]])
    print(head(fragments_cl_bamGR_tab))
    write_tsv(fragments_cl_bamGR_tab, path=FN, col_names= FALSE)
  }
}
Output_to_bed_files(Fragment_list,'./')
for(i in names(Fragment_list)){
  print(i)
  shell= paste("bedtools bedpetobam -i ",i,"_fragments_cl_bamGR_pe.bed -g /usr/share/bedtools/genomes/human.hg38.genome > ",i,"_fragments_cl_bamGR_pe.bam" ,sep="")
  system(shell,wait=TRUE)
}

for(i in names(Fragment_list)){
  print(i)
  shell= (paste("samtools sort -o ",i,"_fragments_cl_bamGR_pe_s.bam ",i,"_fragments_cl_bamGR_pe.bam",sep=""))
  system(shell,wait=TRUE)
}
for(i in names(Fragment_list)){
  print(i)
  shell= paste("samtools index ",i,"_fragments_cl_bamGR_pe_s.bam",sep="")
  system(shell,wait=TRUE)
}
#run TOBIAS 
#nohup TOBIAS ATACorrect --read_shift 0 0 --bam /data/al_tab_GR_fragments_cl_bamGR_pe_s.bam --genome /data/hg38.fa --peaks /data/peak.bed --blacklist /data/hg38-blacklist.v2.bed --outdir /data/ --cores 4 &
Check_normalized_Signal <- function(file,savefile){
  temp_bw = import.bw(file)
  print(summary(width(temp_bw)))
  print(summary(temp_bw$score))
  saveRDS(temp_bw,file=savefile)
}
file = './al_tab_GR_fragments_cl_bamGR_pe_s_corrected.bw'
savefile = './al_signal'
Check_normalized_Signal(file,savefile)
Total_footprint_Motif = matchMotifs(human_pwms_v2,peak,genome = BSgenome.Hsapiens.UCSC.hg38,out='positions',p.cutoff = 5e-05)
Total_footprint_Motif_GR = Must_to_GR(Total_footprint_Motif)
al_signal  = readRDS('./al_signal')
motifs <- names(human_pwms_v2)
motifs <- substr(motifs,17,100)
tfs <- str_match(motfs, "_\\s*(.*?)\\s*_")
out_all_ext <- data.frame(Motif=motifs,TFs=tfs[,2])
Calculate_footprint <- function(footprint_GR,Signal){
  width = width(footprint_GR)
  #####
  left_s = start(footprint_GR)-width*3
  left_e = start(footprint_GR)-1
  #####
  left_GR = footprint_GR
  start(left_GR) = left_s
  end(left_GR) = left_e
  ######
  right_s = end(footprint_GR)+1
  right_e = end(footprint_GR)+width*3
  ######
  right_GR = footprint_GR
  start(right_GR) = right_s
  end(right_GR) = right_e
  ######
  center_GR = footprint_GR
  ######
  ######
  print('Calculate_footprint_center')
  footprint_GR$center = Calculate_signal_bw(center_GR,Signal)
  print('Calculate_footprint_left')
  footprint_GR$left = Calculate_signal_bw(left_GR,Signal)
  print('Calculate_footprint_right')
  footprint_GR$right = Calculate_signal_bw(right_GR,Signal)
  ######
  return(footprint_GR)
}
Calculate_footprint_celltypes <- function(footprint_GR,Signal,TFs,out_all_ext){
  motifs_need = out_all_ext$Motif[which(out_all_ext$TFs %in% TFs == T)]
  ####
  print(length(motifs_need))
  ####
  footprint_GR_cl = footprint_GR[which(footprint_GR$motifs %in% motifs_need == T)]
  ####
  footprint_GR_cl_out = Calculate_footprint(footprint_GR_cl,Signal)
  ####
  return(footprint_GR_cl_out)
}
Calculate_signal_bw <- function(GR,Signal){
  ####
  res = findOverlaps(GR,Signal)
  res_score = Signal$score[subjectHits(res)]
  ####
  res_score_merge = tapply(res_score,queryHits(res),sum)
  ####
  GR$score = 0
  ####
  GR$score[as.numeric(names(res_score_merge))] = as.numeric(res_score_merge)
  ####
  return(GR$score)
}
Filter_footprints <- function(footprint_GR_cl_out,delta=0.1){
  ####
  left_delta = (footprint_GR_cl_out$left)/3 - footprint_GR_cl_out$center
  right_delta = (footprint_GR_cl_out$right)/3 - footprint_GR_cl_out$center
  ####
  k = which(footprint_GR_cl_out$left > delta*3 & footprint_GR_cl_out$right > delta*3 & footprint_GR_cl_out$center < -delta)
  ####
  footprint_GR_cl_out_filtered = footprint_GR_cl_out[k]
  ####
  return(footprint_GR_cl_out_filtered)
}
Early_Diff_Genes = markers[which(markers$avg_log2FC > 0 & markers$p_val_adj < 0.01),]
Early_Diff_Genes_tab = Process_DEGs_to_Celltypes(Early_Diff_Genes,index=index)
al_genes = Early_Diff_Genes_tab$genes[which(Early_Diff_Genes_tab$Tip >0 |Early_Diff_Genes_tab$Acinar>0)]

al_footprints = Calculate_footprint_celltypes(Total_footprint_Motif_GR,al_signal,al_genes,out_all_ext)
al_footprints_cl = Filter_footprints(al_footprints,delta=0.1)

random_cells_by_celltypes = function(x,celltypes){
  x$cell_id = colnames(x)
  ####
  x_cl = subset(x,subset= celltype %in% celltypes == T)
  print(dim(x_cl@meta.data))
  ####
  x_cl$celltype <- droplevels(x_cl$celltype)
  min_numbers = min(table(x_cl$celltype))
  print(min_numbers)
  #### sample for each cell types #####
  cell_list = c()
  ####
  for(i in 1:length(celltypes)){
    tmp_celltypes = celltypes[i]
    print(tmp_celltypes)
    tmp_cells = colnames(x_cl)[which(x_cl$celltype == celltypes[i])]
    print(length(tmp_cells))
    index_sample = sample(1:length(tmp_cells),min_numbers,replace=F)
    ###
    tmp_cells_choose = tmp_cells[index_sample]
    ##
    cell_list = c(cell_list,tmp_cells_choose)
    ###
    print(length(cell_list))
  }
  print(head(cell_list))
  #k = which(x_cl$cell_id %in% cell_list == T)
  x_cl_choose = subset(x_cl,subset = cell_id %in% cell_list)
  print(table(x_cl_choose$celltype))
  ####
  return(x_cl_choose)
  ####
}
sparse.cor3 <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  cSums <- colSums(x)
  # Calculate the population covariance matrix.
  # There's no need to divide by (n-1) as the std. dev is also calculated the same way.
  # The code is optimized to minize use of memory and expensive operations
  covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
  crossp <- as.matrix(crossprod(x))
  covmat <- covmat+crossp
  sdvec <- sqrt(diag(covmat)) # standard deviations of columns
  covmat/crossprod(t(sdvec)) # correlation matrix
}
Reg_one_cells_RPC_MG <- function(footprint_res,peak_gene_list,out_all_ext,diff_list){
  ###### process the list ##########
  ###### generate peak-gene tables  ######
  peak_gene_list$Gene1_TSS$Correlation = 'TSS'
  peaks = c()
  genes = c()
  PtoG = c()
  for(i in 1:length(peak_gene_list)){
    peaks = c(peaks,peak_gene_list[[i]]$peaks)
    print(length(peak_gene_list[[i]]$peaks))
    genes = c(genes,peak_gene_list[[i]]$gene_name)
    print(length(peak_gene_list[[i]]$gene_name))
    PtoG = c(PtoG,peak_gene_list[[i]]$Correlation)
  }
  #### ###### ###### ###### ###### #####
  ### k = which(genes %in% rownames(diff_list) == T)
  #print(length(k)/length(genes))
  #print(length(genes[!duplicated(genes)]))
  ####
  Out_Dat = data.frame(Peaks =peaks ,Target=genes,PtoG=PtoG)
  dim(Out_Dat)
  ####
  #### Next filter the motifs according to gene expression ############
  print(length(diff_list))
  out_all_ext_cl = out_all_ext[which(out_all_ext$TFs %in% diff_list == T),]
  Motif_need = out_all_ext_cl$Motif
  #### filter the footprint ###################
  k = which(footprint_res$motifs %in% Motif_need == T)
  footprint_res_cl = footprint_res[k]
  #### merge Out_Dat with footprint_res_cl ######
  Out_Dat_GR = GRanges(Out_Dat$Peaks,Target=Out_Dat$Target,PtoG=PtoG)
  #### remove their not overlapped regions ####################
  k1 = which(countOverlaps(Out_Dat_GR,footprint_res_cl) >0)
  k2 = which(countOverlaps(footprint_res_cl,Out_Dat_GR) >0)
  #### keep overlapped ######
  Out_Dat_GR_Overlap = Out_Dat_GR[k1]
  footprint_res_cl_Overlap = footprint_res_cl[k2]
  #### merge Out_Dat_GR_Overlap and footprint_res_cl_Overlap ####
  Res_find = data.frame(findOverlaps(footprint_res_cl_Overlap,Out_Dat_GR_Overlap))
  ####
  Res_find$Motifs = footprint_res_cl_Overlap$motifs[Res_find$queryHits]
  Res_find$footprint = as.character(footprint_res_cl_Overlap)[Res_find$queryHits]
  ####
  Res_find$Target = Out_Dat_GR_Overlap$Target[Res_find$subjectHits]
  Res_find$peaks = as.character(Out_Dat_GR_Overlap)[Res_find$subjectHits]
  Res_find$PtoG = Out_Dat_GR_Overlap$PtoG[Res_find$subjectHits]
  ##### add TF names to it #######
  ##### one motif may corresponding to multiple TF genes #####
  colnames(out_all_ext) = c('Motifs','TFs')
  out_all_ext_cl = out_all_ext[which(out_all_ext$TFs %in% diff_list == T),]
  #####
  Res_find_merge = merge(Res_find,out_all_ext_cl)
  ####
  Res_find_merge = Res_find_merge[,c(8,1,4,6,7,5)]
  ####
  return(Res_find_merge)
}
RNA_seurat_choose = random_cells_by_celltypes(rna,c("Tip","Acinar","Trunk","Duct","EP","endocrine"))
Early_Corr = RNA_Corr_Add_cutoff(RNA_seurat_choose[['RNA']]@data)
al_Reg_motif = Reg_one_cells_RPC_MG(al_footprints_cl,peak_gene_list,out_all_ext,All_genes_test)
al_Reg_motif = Add_Cor_to_GRN_network_and_Filter(Early_Corr,al_Reg_motif,All_genes_test)
f1 = which(al_Reg_motif$TFs %in% al_genes == T)
al_Reg_motif_cl = al_Reg_motif[f1,]
al_Reg_motif_cl <- al_Reg_motif_cl[!duplicated(al_Reg_motif_cl$index),]

#repeat the above analysis for duct lineage cells,trunk cells, duct cells and endocrine lineage cells


#common and specific TFs in acinar and duct lineages
common <- unique(intersect(al$TFs,dl$TFs))

altfs <- unique(al_Reg_motif_cl$TFs)
altfs <- altfs[altfs%ni%common]
dltfs <- unique(dl_Reg_motif_cl$TFs)
dltfs <- dltfs[dltfs%ni%common]

#screen common TFs using scRNA-seq data and choose TFs that had targets in both lineages
dl_genes <- marker_rna$gene[which(marker_rna$cluster%in%c("Trunk","Duct")&
                                    marker_rna$avg_log2FC>0&
                                    marker_rna$p_val_adj<0.05)]

al_genes <- marker_rna$gene[which(marker_rna$cluster%in%c("Tip","Acinar")&
                                     marker_rna$avg_log2FC>0&
                                     marker_rna$p_val_adj<0.05)]
common_df <- DataFrame(row.names = common,
                     genes=common)
common_df$dlTg <- 0
for (i in 1:31) {
  y <- dl_Reg_motif_cl[which(dl_Reg_motif_cl$TFs %in%common[i]),"Target"]
  y <- y[y%in%dl_genes]
  common_df[i,"dlTg"] <- length(y)
}
common_df$alTg <- 0
for (i in 1:31) {
  y <- al_Reg_motif_cl[which(al_Reg_motif_cl$TFs %in%common[i]),"Target"]
  y <- y[y%in%al_genes]
  common_df[i,"alTg"] <- length(y)
}
common <- commons$genes[which(common_df$alTg>0&common_df$dlTg>0)]


#choose GRNs regulating bi-potent trunk cells differentiation and GRNs in duct and endocrine cells
whole <- list()
whole[[1]] <- trunk_Reg_motif_cl
whole[[2]] <- duct_Reg_motif_cl
whole[[3]] <- endo_Reg_motif_cl

class1_endo <- whole[[1]][which(whole[[1]]$Target%in%whole[[3]]$TFs),]
class1_endo <- subset(class1_endo,subset=class1_endo$index%ni%whole[[2]]$index)
class2_endo <- whole[[1]][which(whole[[1]]$index%in%whole[[3]]$index),]
bptoendo <- rbind(class1_endo,class2_endo)
bptoendo$direction <- "endo"
bptoendo <- bptoendo[!duplicated(bptoendo$index),]
class1_duct <- whole[[1]][which(whole[[1]]$Target%in%whole[[2]]$TFs),]
class1_duct<- subset(class1_duct,subset=class1_duct$index%ni%whole[[3]]$index)
class2_duct <- whole[[1]][which(whole[[1]]$index%in%whole[[2]]$index),]
class2_duct <- subset(class2_duct,subset=class2_duct$index%ni%whole[[3]]$index)
bptoduct <- rbind(class1_duct,class2_duct)
bptoduct$direction <- "duct"
bptoendo$origin <- "bp"
bptoduct$origin <- "bp"
bptobp <- whole[[1]][which(whole[[1]]$Target%in%whole[[1]]$TFs),]
bptobp$direction <- "bp"
bptobp$origin <- "bp"
duplicated(bptobp$index)
dtod <- whole[[2]][which(whole[[2]]$Target%in%whole[[2]]$TFs),]
dtod$direction <- "duct"
dtod$origin <- "duct"
dtod <- dtod[!duplicated(dtod$index),]
etoe <- whole[[3]][which(whole[[3]]$Target%in%whole[[3]]$TFs),]
etoe$direction <- "endo"
etoe$origin <- "endo"
etoe <- etoe[!duplicated(etoe$index),]

bptoendo$from <- "bptoendo"
bptoduct$from <- "bptoduct"
bptobp$from <- "bptobp"
dtod$from <- "dtod"
etoe$from <- "etoe"
etoe <- etoe[!duplicated(etoe$index),]
dtod <- dtod[!duplicated(dtod$index),]

inter <- rbind(bptoendo,bptoduct,bptobp,dtod,etoe)
inter <- subset(inter,subset=Target%in%out_all_ext$TFs)

inters <- inter[,c("TFs","Target","index","PtoG","Cor","Tag","direction","origin")]
inters <- inters[!duplicated(inters$index),]

name<-data.frame(c(inters$TFs,inters$Target))
nodes<-name%>%
  distinct()
nodes$groups <- "no"
colnames(nodes) <- c("genes","groups")

nodes$groups[which(nodes$genes%in%inters$Target)] <- "Target"
nodes$groups[which(nodes$genes%in%inters$TFs)] <- "TF"
nodes$origins <- "no"
nodes$origins[which(nodes$genes%in%c(dtod$TFs,dtod$Target))] <- "duct"
nodes$origins[which(nodes$genes%in%c(etoe$TFs,etoe$Target))] <- "ep"
nodes$origins[which(nodes$genes%in%c(bptobp$TFs,bptoendo$TFs,bptoduct$TFs))] <- "bp"
hl <- names(degree(net)[order(degree(net),decreasing = T)][1:10])
edge <- data.frame(inters$TFs,inters$Target,inters$index,inters$origin,inters$direction,inters$Tag)
edge$type <- "no"
edge$type[which(edge$inters.index%in%c(bptoduct$index,bptoendo$index,bptobp$index))] <- "yes"
edge$type[which(edge$inters.TFs%in%hl|edge$inters.Target%in%hl)] <- "type2"
edge$type[which(edge$inters.Target%in%trunk_genes|edge$inters.TFs%in%trunk_genes)] <- "yes"

edges<-edge%>%
  rename(from=inters.TFs,to=inters.Target,origin=inters.origin,tag=inters.Tag,type=type)

net<-graph_from_data_frame(
  d=edges,vertices=nodes,
  directed=T)

write.graph(net,"./0511/0516/graph.gml",format = "gml")

#the network were visualized with Cytoscape and the interested TF-Target pairs were chosen manually and checked with scRNA-seq expression