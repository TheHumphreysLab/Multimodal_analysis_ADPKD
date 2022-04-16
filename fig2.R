library(Seurat) #3.2.0
library(Signac)
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

sub_atac@meta.data$celltype_disease <- paste0(sub_atac@meta.data$celltype,"_",sub_atac@meta.data$disease)
Idents(sub_atac) <- "celltype_disease"
levels(sub_atac) <- c("PCT_Cont","PCT_PKD","PST_Cont","PST_PKD","FR-PTC_Cont","FR-PTC_PKD",
                     "PEC_PODO_Cont","PEC_PODO_PKD","TAL_Cont","TAL_PKD",
                     "DCT_Cont","DCT_PKD","CNT_PC_Cont","CNT_PC_PKD",
                     "IC_Cont","IC_PKD","ENDO_Cont","ENDO_PKD",
                     "FIB_Cont","FIB_PKD","LEUK_Cont","LEUK_PKD")

features <- c("SLC34A1","SLC36A2","SLC5A2","SLC5A1","SLC6A18","VCAM1","CRYAB","CFH","NPHS1","NPHS2",
              "SLC12A1","UMOD","PAPPA2","TCF24","SLC12A3","SALL3","AQP2","SCNN1G","SCNN1B","FAM20A","SLC26A4","SLC26A7",
              "KIT","FLT1","EMCN","ACTA2","MYH11","PDGFRA","GLI2","GLI1","PTPRC","SLC11A1","CSF2RA")

DefaultAssay(sub_atac) <- "RNA"
levels(sub_atac) <- rev(levels(sub_atac))
fig2d <- DotPlot(sub_atac, features = features, cols = c("lightyellow","#992FA0")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
levels(sub_atac) <- rev(levels(sub_atac)) #1400x750
