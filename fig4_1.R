library(Seurat)
library(ggplot2)
library(sctransform)
library(harmony)
library(Rcpp)
library(dplyr)
library(BuenColors)
set.seed(1234)

rnaAggr <- readRDS("PKDContAggr.rds")

#The target cell type subset using the annotations on an integrated dataset (Fig. 1b) was then further subset with the annotations on each dataset to extract cell type with high confidence. 
#(Supplementary Fig. 1d, Supplementary Fig. 2d)

Idents(rnaAggr) <- "celltype"
FIB <- subset(rnaAggr,idents = c("PKD-FIB1","PKD-FIB2","PKD-FIB3","Cont-FIB"))
Idents(FIB) <- "celltype_all"
FIB <- subset(FIB,idents = "FIB")

DefaultAssay(FIB) <- "RNA"
FIB <- NormalizeData(FIB)
FIB <- FindVariableFeatures(FIB, selection.method = "vst", nfeatures = 3000)
FIB <- ScaleData(FIB, verbose = FALSE)
FIB <- RunPCA(FIB, verbose = FALSE)
FIB <- RunHarmony(FIB,group.by.vars = c("orig.ident"))
FIB <- RunUMAP(FIB, reduction = "harmony", dims = 1:20)
FIB <- FindNeighbors(FIB, reduction = "harmony", dims = 1:20)
FIB <- FindClusters(FIB, resolution = 0.2)

fig4a_1 <- DimPlot(FIB,label = T)+NoLegend() #600x500
Idents(FIB) <- "disease"
levels(FIB) <- c("PKD","control")
fig4a_2 <- DimPlot(FIB)+NoLegend() #600x500

Idents(FIB) <- "seurat_clusters"
new.cluster.ids <- c("MyoFib","PKD-FIB","Unknown1","FIB1","FIB2","Unknown2","Fat")
names(new.cluster.ids) <- levels(FIB)
FIB <- RenameIdents(FIB, new.cluster.ids)
levels(FIB) <- c("FIB1","FIB2","PKD-FIB","MyoFib","Fat","Unknown1","Unknown2")
FIB@meta.data[["subtype"]] <- FIB@active.ident
saveRDS(FIB,"FIB.rds")

marker <- FindAllMarkers(FIB,only.pos = T,logfc.threshold = 0.5)

FIB2 <- subset(FIB,idents = c("Unknown1","Unknown2"),invert=T)

features <- c("NTRK3","RGS6","CSPG4","PDGFRB","EDNRA","HS3ST3A1","AGTR1","REN","PDGFRA","FGF14","SEMA4A","IL6","ACTA2","FN1","VCAN","PLIN1","LPL","ADIPOQ")
features2 <- c("GLI1","GLI2","GLI3","GAS1","CDON","BOC","PTCH1","PTCH2")

levels(FIB2) <- rev(levels(FIB2))
fig4b <- DotPlot(FIB2, features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #750x400
fig4g <- DotPlot(FIB2, features = features2, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #600x400
levels(FIB2) <- rev(levels(FIB2))

fig4d <- FeaturePlot(rnaAggr,"GLI1",raster = F,pt.size = 0.2,
            order=T,split.by = "disease") #900x440

figS5_1 <- FeaturePlot(rnaAggr,"IL6",order=T,raster=F) #600x530
figS5_2 <- FeaturePlot(FIB,"IL6",order=T,raster=F) #600x530
figS5_3 <- FeaturePlot(rnaAggr,"TNF",order=T,raster = F,pt.size = 0.5) #600x530
figS5_4 <- FeaturePlot(FIB,"GLI1",raster = F,pt.size = 0.2,
            order=T)

##################################### ATAC-seq data ########################################

library(Seurat)
library(Signac) 
library(monocle3)
library(cicero)
library(BuenColors)
library(EnsDb.Hsapiens.v86)
set.seed(1234)
sub_atac <- readRDS("pkd_contATAC_Aggr_sub80.rds")

FIB_ATAC <- subset(sub_atac,idents = "FIB")

DefaultAssay(FIB_ATAC) <- "peaks"
Idents(FIB_ATAC) <- "disease"
levels(FIB_ATAC) <- c("PKD","Cont")

fig4e <- CoveragePlot(
  object = FIB_ATAC,
  region = "chr12-57458717-57463141",
  annotation = T,
  peaks = T,
  extend.upstream = 3000,
  extend.downstream = 3000,
) #800x800 GLI1 promoter

DefaultAssay(FIB_ATAC) <- "chromvar"
fig4f <- VlnPlot(FIB_ATAC,"MA0734.1",pt.size = 0) #GLI2, 500x500
