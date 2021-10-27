#Individual snRNA-seq data preprocessing
#For control 1-5 and ADPKD 1-8

library(SoupX)
library(DoubletFinder)
library(Seurat)
library(ggplot2)
library(dplyr)

#Correction of ambient RNA with SoupX -------------------------------------------------------------------------------------------------------

toc = Read10X_h5("/control_1/outs/filtered_feature_bc_matrix.h5")
tod = Read10X_h5("/control_1/outs/raw_feature_bc_matrix.h5")

control_1 <- CreateSeuratObject(counts = toc, min.cells = 10, project = "RNA")
control_1 <- NormalizeData(control_1)
control_1 <- ScaleData(control_1)
control_1 <- FindVariableFeatures(control_1, selection.method = "vst", nfeatures = 2000)
control_1 <- RunPCA(control_1)
control_1 <- FindNeighbors(control_1, dims = 1:20)
control_1 <- FindClusters(control_1)
control_1 <- RunUMAP(control_1, dims = 1:20)

sc = SoupChannel(tod,toc)
sc = setClusters(sc,setNames(control_1@meta.data[["seurat_clusters"]], rownames(control_1@meta.data)))
sc = autoEstCont(sc,forceAccept=TRUE)
out = adjustCounts(sc)


#Quality control step ------------------------------------------------------------------------------------------------------------------------
#remove low-quality nuclei (nuclei with top 5% and bottom 1% in the distribution of feature counts or RNA count, 
#or those with %Mitochondrial genes > 0.25). 

control_1 <- CreateSeuratObject(counts = out, min.cells = 10, min.features = 500, project = "RNA")
control_1 <- PercentageFeatureSet(control_1, pattern = "^MT-", col.name = "percent.mt")
control_1 <- PercentageFeatureSet(control_1, pattern = "^RPL", col.name = "percent.rpl")
control_1 <- PercentageFeatureSet(control_1, pattern = "^RPS", col.name = "percent.rps")

LowQC1 <- rownames(top_frac(control_1@meta.data, 0.05, nFeature_RNA))
LowQC2 <- rownames(top_frac(control_1@meta.data, 0.05, nCount_RNA))
LowQC3 <- rownames(top_frac(control_1@meta.data, -0.01, nFeature_RNA))
LowQC4 <- rownames(top_frac(control_1@meta.data, -0.01, nCount_RNA))
LowQC <- unique(c(LowQC1,LowQC2,LowQC3,LowQC4))

control_1 <- subset(control_1, cells =  LowQC,invert = T)
control_1 <- subset(control_1, percent.mt < 0.25)

#Doublet detection  ---------------------------------------------------------------------------------------------------------------------------

control_1 <- NormalizeData(control_1)
control_1 <- ScaleData(control_1)
control_1 <- FindVariableFeatures(control_1, selection.method = "vst", nfeatures = 2000)
control_1 <- RunPCA(control_1)
control_1 <- FindNeighbors(control_1, dims = 1:10)
control_1 <- FindClusters(control_1)
control_1 <- RunUMAP(control_1, dims = 1:10)
DimPlot(control_1,label = T)+NoLegend()

## pK Identification  
sweep.res.list_kidney <- paramSweep_v3(control_1, PCs = 1:10, sct = F)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)
#determine pK for individual dataset as instructed in tutorial

#Heterotypic doublet rate is set to 8%
nExp_poi <- round(0.08*length(control_1@active.ident))

## Run DoubletFinder
control_1 <- doubletFinder_v3(control_1, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#Use pK determined above, pK=0.005 for control_1

## Add detected doublet to metadata of Seurat object
doubletdata <- control_1@meta.data[["DF.classifications_0.25_0.005_561"]]
names(doubletdata) <- rownames(control_1@meta.data)
doubletdata <- as.data.frame(doubletdata)
control_1 <- AddMetaData(control_1,doubletdata)
saveRDS(control_1,"control_1.rds")
