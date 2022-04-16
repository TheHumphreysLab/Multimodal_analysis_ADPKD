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
CNT_PC <- subset(rnaAggr,idents = c("PKD-CNT_PC","Cont-CNT_PC","PKD-CTC"))
Idents(CNT_PC) <- "celltype_all"
CNT_PC <- subset(CNT_PC,idents = "CNT_PC")

CNT_PC <- NormalizeData(CNT_PC)
CNT_PC <- FindVariableFeatures(CNT_PC, selection.method = "vst", nfeatures = 3000)
CNT_PC <- ScaleData(CNT_PC, verbose = FALSE)
CNT_PC <- RunPCA(CNT_PC, verbose = FALSE)
CNT_PC <- RunHarmony(CNT_PC,group.by.vars = c("orig.ident"))
CNT_PC <- RunUMAP(CNT_PC, reduction = "harmony", dims = 1:20)
CNT_PC <- FindNeighbors(CNT_PC, reduction = "harmony", dims = 1:20)
CNT_PC <- FindClusters(CNT_PC, resolution = 0.3)
fig7a <- DimPlot(CNT_PC,label = T,repel=T)+NoLegend() #550x440

Idents(CNT_PC) <- "seurat_clusters"
new.cluster.ids <- c("PKD-CTC1","N-CNT","N-PC","LowQC1","PKD-CTC2","PKD-CNT","LowQC2","PKD-CTC3")
names(new.cluster.ids) <- levels(CNT_PC)
CNT_PC <- RenameIdents(CNT_PC, new.cluster.ids)
levels(CNT_PC) <- c("N-CNT","PKD-CNT","N-PC","PKD-CTC1","PKD-CTC2","PKD-CTC3","LowQC1","LowQC2")
CNT_PC@meta.data[["subtype"]] <- CNT_PC@active.ident
saveRDS(CNT_PC,"CNT_PC.rds")

Idents(CNT_PC) <- "subtype"
fig7a_2 <- DimPlot(CNT_PC,label = T)+NoLegend()
Idents(CNT_PC) <- "disease"
fig7a_3 <- DimPlot(CNT_PC) #450x330
