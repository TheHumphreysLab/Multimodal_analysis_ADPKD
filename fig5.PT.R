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
PTC <- subset(rnaAggr,idents = c("PKD-PT","Cont-PT1","Cont-PT2"))
Idents(PTC) <- "celltype_all"
PTC <- subset(PTC,idents = c("PT","FR-PTC"))

DefaultAssay(PTC) <- "RNA"
PTC <- NormalizeData(PTC)
PTC <- FindVariableFeatures(PTC, selection.method = "vst", nfeatures = 3000)
PTC <- ScaleData(PTC, verbose = FALSE)
PTC <- RunPCA(PTC, verbose = FALSE)
PTC <- RunHarmony(PTC,group.by.vars = c("orig.ident"))
PTC <- RunUMAP(PTC, reduction = "harmony", dims = 1:10)
PTC <- FindNeighbors(PTC, reduction = "harmony", dims = 1:10)
PTC <- FindClusters(PTC, resolution = 0.1)

Idents(PTC) <- "seurat_clusters"
new.cluster.ids <- c("FR-PTC","N-PTC")
names(new.cluster.ids) <- levels(PTC)
PTC <- RenameIdents(PTC, new.cluster.ids)
levels(PTC) <- c("N-PTC","FR-PTC")
PTC@meta.data[["subtype"]] <- PTC@active.ident
saveRDS(PTC,"PTC.rds")