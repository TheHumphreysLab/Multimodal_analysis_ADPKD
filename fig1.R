library(Seurat)
library(ggplot2)
library(sctransform)
library(harmony)
library(Rcpp)
library(dplyr)
set.seed(1234)

rnaAggr <- readRDS("PKDContAggr.rds")

Idents(rnaAggr) <- "seurat_clusters"
fig1a <- DimPlot(rnaAggr,raster = F) #750x630
Idents(rnaAggr) <- "disease"
levels(rnaAggr) <- c("PKD","control")
fig1b <- DimPlot(rnaAggr,raster = F)+NoLegend() #700x630

Idents(rnaAggr) <- "disease"
pkd <- subset(rnaAggr,idents = "PKD")
cont <- subset(rnaAggr,idents = "control")
Idents(cont) <- "celltype_all"
Idents(pkd) <- "celltype_all"

features <- c("SLC34A1","LRP2","RHEX","CFH","SLC12A1","CD44","SLC12A3",
              "SLC8A1","AQP2","SLC26A7","SLC26A4","NPHS2",
              "EMCN","PDGFRB","ACTA2")

cont <- subset(cont,idents = "LEUK",invert=T)
levels(cont) <- rev(levels(cont))
fig1c <- DotPlot(cont, features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #800x600
levels(cont) <- rev(levels(cont)) #800x600

features <- c("CDH6","TGFB2","REHX","CFH","KIRREL3","COL7A1","SLC12A1",
              "SLC12A3","SLC14A2","SCNN1G","SLC26A7","SLC26A4",
              "NPHS1","FLT1","GPM6A","ACTA2","MYH11","MEG3","PTPRC","SLC11A1","UPK3A","PSCA")

levels(pkd) <- rev(levels(pkd))
fig1d <- DotPlot(pkd, features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #850x600
levels(pkd) <- rev(levels(pkd))