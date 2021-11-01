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
ENDO <- subset(rnaAggr,idents = c("PKD-ENDO1","PKD-ENDO2","PKD-ENDO3","Cont-ENDO1","Cont-ENDO2","Cont-ENDO3"))
Idents(ENDO) <- "celltype_all"
ENDO <- subset(ENDO,idents = "ENDO")

DefaultAssay(ENDO) <- "RNA"
ENDO <- NormalizeData(ENDO)
ENDO <- FindVariableFeatures(ENDO, selection.method = "vst", nfeatures = 3000)
ENDO <- ScaleData(ENDO, verbose = FALSE)
ENDO <- RunPCA(ENDO, verbose = FALSE)
ENDO <- RunHarmony(ENDO,group.by.vars = c("orig.ident"))
ENDO <- RunUMAP(ENDO, reduction = "harmony", dims = 1:10)
ENDO <- FindNeighbors(ENDO, reduction = "harmony", dims = 1:10)
ENDO <- FindClusters(ENDO, resolution = 0.2)
fig4h_1 <- DimPlot(ENDO,label = T,repel=T)+NoLegend() #600x500
Idents(ENDO) <- "disease"
levels(ENDO) <- c("PKD","control")
fig4h_2 <- DimPlot(ENDO)+NoLegend() #600x500

Idents(ENDO) <- "seurat_clusters"
new.cluster.ids <- c("EC_SELE","AEC","LowQC","CEC","VEC","LEC")
names(new.cluster.ids) <- levels(ENDO)
ENDO <- RenameIdents(ENDO, new.cluster.ids)
levels(ENDO) <- c("CEC","VEC","AEC","EC_SELE","LEC","LowQC")
ENDO@meta.data[["subtype"]] <- ENDO@active.ident
saveRDS(ENDO,"ENDO.rds")

#Vision
library(VISION)
signatures <- "/vision/h.all.v6.2.symbols.gmt"
vision.obj <- Vision(ENDO, signatures = signatures)
vision.obj <- analyze(vision.obj)

sigScores <- getSignatureScores(vision.obj)
sigScores <- as.data.frame(sigScores)

ENDO <- AddMetaData(ENDO,sigScores)
figS7_b <- FeaturePlot(ENDO,"HALLMARK_INFLAMMATORY_RESPONSE",cols =jdb_palette("brewer_yes"),order = T)

#Remove low QC cluster
ENDO2 <- subset(ENDO,idents = c("LowQC"),invert=T)

features <- c("EMCN","HECW2","ITGA8","NRG3","DNASE1L3","CEACAM1","ADGRB3",
               "ADAMTS6","VEGFC","PDE3A","NKAIN2","ANO2","VWF","HIF1A","SELE","VCAM1",
              "MMRN1","PROX1","RELN","CCL21")

levels(ENDO2) <- rev(levels(ENDO2))
fig4i <- DotPlot(ENDO2, features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #750x450

levels(ENDO2) <- rev(levels(ENDO2))

figS7a_1 <- FeaturePlot(ENDO,"SELE",pt.size = 0.3) #600x500
figS7a_2 <- FeaturePlot(ENDO,"VCAM1",pt.size = 0.3) #600x500
figS7a_3 <- FeaturePlot(ENDO,"ICAM1",pt.size = 0.3) #600x500

fig4j_1 <- FeaturePlot(rnaAggr,"DHH",pt.size = 0.5,order=T,raster=F) #600x500
fig4j_2 <- FeaturePlot(ENDO,"DHH",pt.size = 1,order=T) #600x500

##################################### ATAC-seq data ########################################

library(Seurat)
library(Signac) 
library(monocle3)
library(cicero)
library(BuenColors)
library(EnsDb.Hsapiens.v86)
set.seed(1234)
sub_atac <- readRDS("pkd_contATAC_Aggr_sub80.rds")
ENDO_ATAC <- subset(sub_atac,idents = "ENDO")

Idents(ENDO_ATAC) <- "disease"
DefaultAssay(ENDO_ATAC) <- "peaks"
ENDO_ATAC_peakMarker <- FindAllMarkers(ENDO_ATAC,only.pos = T)
test <- ClosestFeature(ENDO_ATAC, regions = ENDO_ATAC_peakMarker$gene)
ENDO_ATAC_peakMarker <- cbind(ENDO_ATAC_peakMarker,test[,c(2,8)])

DefaultAssay(ENDO_ATAC) <- "chromvar"
ENDO_ATAC_chromvarMarker <- FindAllMarkers(ENDO_ATAC,only.pos = T, logfc.threshold = 0.25)

library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(TFBSTools)

pfm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = T)
)

x <- NULL                        
for (i in ENDO_ATAC_chromvarMarker$gene) { 
  x <- c(x,pfm@listData[[i]]@name)
}
ENDO_ATAC_chromvarMarker$gene_name <- x

Idents(ENDO_ATAC) <- "disease"
levels(ENDO_ATAC) <- c("PKD","Cont")
figS7c_rela <- VlnPlot(ENDO_ATAC,"MA0107.1",pt.size=0) #RELA 500x500 
figS7c_fosjun <- VlnPlot(ENDO_ATAC,"MA0099.3",pt.size=0) #STAT3 500x500
