library(Seurat)
library(Signac) 
library(monocle3)
library(cicero)
library(BuenColors)
library(EnsDb.Hsapiens.v86)
library(here)
set.seed(1234)

sub_atac <- readRDS("pkd_contATAC_Aggr_sub80.rds")
CNT_PCatac <- subset(sub_atac,idents = "CNT_PC")

DefaultAssay(CNT_PCatac) <- "peaks"
CNT_PCatac <- RunUMAP(CNT_PCatac, dims = 2:20, reduction = 'harmony')
CNT_PCatac <- FindNeighbors(object = CNT_PCatac, reduction = "harmony", dims = 2:20)
CNT_PCatac <- FindClusters(object = CNT_PCatac, verbose = FALSE, algorithm = 3,resolution = 0.2)

#Sub_clustering of CNT_PC in snATAC-seq / Annotation of subtypes
Idents(CNT_PCatac) <- "seurat_clusters"
new.cluster.ids <- c("N-CNT","N-PC","PKD-PC")
names(new.cluster.ids) <- levels(CNT_PCatac)
CNT_PCatac <- RenameIdents(CNT_PCatac, new.cluster.ids)
CNT_PCatac@meta.data[["subtype"]] <- CNT_PCatac@active.ident

Idents(CNT_PCatac) <- "subtype"
fig8a_1 <- DimPlot(CNT_PCatac)+NoLegend() #550x440
Idents(CNT_PCatac) <- "disease"
levels(CNT_PCatac) <- c("PKD","Cont")
fig8a_2 <- DimPlot(CNT_PCatac)+NoLegend()

#Dot plot of gene activities
Idents(CNT_PCatac) <- "subtype"
DefaultAssay(CNT_PCatac) <- "RNA"
features <- c("SLC8A1","FXYD2","CALB1","AQP2","SCNN1G",
              "HRH1","CD44","GPRC5A","ROR1")

levels(CNT_PCatac) <- rev(levels(CNT_PCatac))
fig8b <- DotPlot(CNT_PCatac, features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #600x450
levels(CNT_PCatac) <- rev(levels(CNT_PCatac))