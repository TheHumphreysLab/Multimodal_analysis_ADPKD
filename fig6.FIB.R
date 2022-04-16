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

Idents(FIB) <- "seurat_clusters"
new.cluster.ids <- c("MyoFib","PKD-FIB","Unknown1","FIB1","FIB2","Unknown2","Fat")
names(new.cluster.ids) <- levels(FIB)
FIB <- RenameIdents(FIB, new.cluster.ids)
levels(FIB) <- c("FIB1","FIB2","PKD-FIB","MyoFib","Fat","Unknown1","Unknown2")
FIB@meta.data[["subtype"]] <- FIB@active.ident
saveRDS(FIB,"FIB.rds")