library(Seurat)
library(ggplot2)
library(sctransform)
library(harmony)
library(Rcpp)
library(dplyr)
set.seed(1234)

sub_atac <- readRDS("pkd_contATAC_Aggr_sub80.rds")

Idents(sub_atac) <- "seurat_clusters"
fig2a <- DimPlot(sub_atac) #750x630
Idents(sub_atac) <- "disease"
levels(sub_atac) <- c("PKD","Cont")
fig2b <- DimPlot(sub_atac)+NoLegend() #700x630

Idents(sub_atac) <- "disease"
cont <- subset(sub_atac,idents = "Cont")
pkd <- subset(sub_atac,idents = "PKD")
Idents(cont) <- "celltype"
Idents(pkd) <- "celltype"

features <- c("SLC34A1","LRP2","SLC5A2","SLC5A1","HAVCR1","CFH","NPHS2",
              "SLC12A1","SLC12A3","AQP2","SLC26A7",
              "SLC26A4","EMCN","ACTA2","PTPRC")
levels(cont) <- rev(levels(cont))
fig2c <- DotPlot(cont, features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #800x600
levels(cont) <- rev(levels(cont)) #800x600

features <- c("ACSM2A","CDH6","CFH","NPHS1","SLC12A1",
              "SLC12A3","SCNN1G","SLC26A7","SLC26A4",
              "FLT1","EMCN","ACTA2","MYH11","PTPRC","SLC11A1")

levels(pkd) <- rev(levels(pkd))
fig2c <- DotPlot(pkd, features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #800x600
levels(pkd) <- rev(levels(pkd))