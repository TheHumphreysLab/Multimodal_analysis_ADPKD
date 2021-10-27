#Data integration for control, ADPKD and all datasets)

library(Seurat)
library(ggplot2)
library(sctransform)
library(harmony)
library(Rcpp)
library(dplyr)
set.seed(1234)

######### control dataset ##########

control_1 <- readRDS("control_1.rds")
control_2 <- readRDS("control_2.rds")
control_3 <- readRDS("control_3.rds")
control_4 <- readRDS("control_4.rds")
control_5 <- readRDS("control_5.rds")

control_1@meta.data[["orig.ident"]] <- "control1"
control_2@meta.data[["orig.ident"]] <- "control2"
control_3@meta.data[["orig.ident"]] <- "control3"
control_4@meta.data[["orig.ident"]] <- "control4"
control_5@meta.data[["orig.ident"]] <- "control5"

# Integration (control dataset)
control.list <- list(control_1,control_2,control_3,control_4,control_5)
control.list <- lapply(X = control.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = control.list, nfeatures = 3000)
control.list <- PrepSCTIntegration(object.list = control.list, anchor.features = features)
control.anchors <- FindIntegrationAnchors(object.list = control.list, normalization.method = "SCT", 
                                          anchor.features = features)
cont <- IntegrateData(anchorset = control.anchors, normalization.method = "SCT")
DefaultAssay(cont) <- "integrated"

# visualization and clustering
cont <- ScaleData(cont, verbose = FALSE)
cont <- RunPCA(cont, npcs = 30, verbose = FALSE)
cont <- RunUMAP(cont, reduction = "pca", dims = 1:20)
cont <- FindNeighbors(cont, reduction = "pca", dims = 1:20)
cont <- FindClusters(cont, resolution = 0.3)
DimPlot(cont, reduction = "umap", label = TRUE)

# annotation
Idents(cont) <- "seurat_clusters"
new.cluster.ids <- c("PT1","DCT","TAL1","TAL2","TAL3","CNT_PC","ENDO1","PT2","Doublet",
                     "PEC","ICA","ICB","PT3","PODO","FIB","ATL","ENDO2")
#Cluster 8 expressed both PT/TAL markers, and gene expression pattern was not specific.
#Cluster 8 enriched doublets predicted by doubletfinder.
#Cluster 8 was thought to be composed of heterotypic doublets, annotated as "Doublets" 
#See Supplementary Fig.2 for detail

names(new.cluster.ids) <- levels(cont)
cont <- RenameIdents(cont, new.cluster.ids)

levels(cont) <- c("PT1","PT2","PT3","PEC","ATL","TAL1","TAL2","TAL3",
                  "DCT","CNT_PC","ICA","ICB","PODO","ENDO1","ENDO2","FIB","Doublet")

#Remove doublets from control dataset
cont <- subset(cont,idents = "Doublet",invert=T)
Idents(cont) <- "doubletdata"
cont <- subset(cont,idents = "Singlet")

# visualization and clustering after clean up
DefaultAssay(cont) <- "integrated"
cont <- ScaleData(cont, verbose = FALSE)
cont <- RunPCA(cont, npcs = 30, verbose = FALSE)
cont <- RunUMAP(cont, reduction = "pca", dims = 1:20)
cont <- FindNeighbors(cont, reduction = "pca", dims = 1:20)
cont <- FindClusters(cont, resolution = 0.3)
DimPlot(cont, reduction = "umap", label = TRUE)

#annotation
Idents(cont) <- "seurat_clusters"
new.cluster.ids <- c("DCT","TAL2","TAL1","PT1","PT1","TAL1","CNT_PC","PEC","PT1",
                     "ICA","ENDO2","ICB","PT2","PODO","FIB","ENDO1","ATL","ENDO3")
names(new.cluster.ids) <- levels(cont)
cont <- RenameIdents(cont, new.cluster.ids)

levels(cont) <- c("PT1","PT2","PEC","ATL","TAL1","TAL2","DCT","CNT_PC","ICA","ICB","PODO",
                  "ENDO1","ENDO2","ENDO3","FIB")

cont@meta.data[["celltype"]] <- cont@active.ident

saveRDS(cont,"controlAggr_filtered.rds")

######### PKD dataset ##########

PKD_1 <- readRDS("PKD_1.rds")
PKD_2 <- readRDS("PKD_2.rds")
PKD_3 <- readRDS("PKD_3.rds")
PKD_4 <- readRDS("PKD_4.rds")
PKD_5 <- readRDS("PKD_5.rds")
PKD_6 <- readRDS("PKD_6.rds")
PKD_7 <- readRDS("PKD_7.rds")
PKD_8 <- readRDS("PKD_8.rds")

PKD_1@meta.data[["orig.ident"]] <- "PKD1"
PKD_2@meta.data[["orig.ident"]] <- "PKD2"
PKD_3@meta.data[["orig.ident"]] <- "PKD3"
PKD_4@meta.data[["orig.ident"]] <- "PKD4"
PKD_5@meta.data[["orig.ident"]] <- "PKD5"
PKD_6@meta.data[["orig.ident"]] <- "PKD6"
PKD_7@meta.data[["orig.ident"]] <- "PKD7"
PKD_8@meta.data[["orig.ident"]] <- "PKD8"

# Integration (ADPKD dataset)
PKD.list <- list(PKD_1,PKD_2,PKD_3,PKD_4,PKD_5,PKD_6,PKD_7,PKD_8)
PKD.list <- lapply(X = PKD.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = PKD.list, nfeatures = 3000)
PKD.list <- PrepSCTIntegration(object.list = PKD.list, anchor.features = features)
PKD.anchors <- FindIntegrationAnchors(object.list = PKD.list, normalization.method = "SCT", 
                                      anchor.features = features)
pkd <- IntegrateData(anchorset = PKD.anchors, normalization.method = "SCT")

DefaultAssay(pkd) <- "integrated"

# visualization and clustering
pkd <- ScaleData(pkd, verbose = FALSE)
pkd <- RunPCA(pkd, npcs = 30, verbose = FALSE)
pkd <- RunUMAP(pkd, reduction = "pca", dims = 1:20)
pkd <- FindNeighbors(pkd, reduction = "pca", dims = 1:20)
pkd <- FindClusters(pkd, resolution = 0.4)
DimPlot(pkd, reduction = "umap", label = TRUE)

#annotation
Idents(pkd) <- "seurat_clusters"
new.cluster.ids <- c("TAL","LowQC","CNT_PC","PT","FIB","DCT",
                     "Transient","FIB","LEUK","ENDO","PT","PC","IC",
                     "PODO","LEUK","ENDO","PEC","FIB")
names(new.cluster.ids) <- levels(pkd)
pkd <- RenameIdents(pkd, new.cluster.ids)
#Cluster 1 did not have specific marker gene expressions, and they enriched mitochondrial genes.
#Cluster 1 was thought to be composed of low-quality nuclei, annotated as "LowQC".
#See Supplementary Fig.3 for detail

levels(pkd) <- c("PT","PEC","Transient","TAL","DCT","CNT_PC","PC","IC","PODO",
                     "ENDO","FIB","LEUK","LowQC")

#Remove LowQC & doublets from ADPKD dataset
pkd <- subset(pkd,idents = "LowQC",invert=T) 
Idents(pkd) <- "doubletdata"
pkd <- subset(pkd,idents = "Singlet")

DefaultAssay(pkd) <- "integrated"

# visualization and clustering after clean up
pkd <- ScaleData(pkd, verbose = FALSE)
pkd <- RunPCA(pkd, npcs = 30, verbose = FALSE)
pkd <- RunUMAP(pkd, reduction = "pca", dims = 1:20)
pkd <- FindNeighbors(pkd, reduction = "pca", dims = 1:20)
pkd <- FindClusters(pkd, resolution = 0.4)
DimPlot(pkd, reduction = "umap", label = TRUE)

#annotation
Idents(pkd) <- "seurat_clusters"
new.cluster.ids <- c("PT","TAL","CNT_PC","FIB1","PC","DCT","LEUK1",
                     "Transitional","FIB2","ENDO1","IC","PODO","ENDO2",
                     "PEC","LEUK2","FIB3","ENDO3")
names(new.cluster.ids) <- levels(pkd)
pkd <- RenameIdents(pkd, new.cluster.ids)

levels(pkd) <- c("PT","PEC","Transitional","TAL","DCT","CNT_PC","PC","IC","PODO",
                 "ENDO1","ENDO2","ENDO3","FIB1","FIB2","FIB3","LEUK1","LEUK2")

pkd@meta.data[["celltype"]] <- pkd@active.ident

saveRDS(pkd,"PKDAggr_filtered.rds")


######### All dataset (control + ADPKD) ##########

Idents(cont) <- "celltype"
Idents(pkd) <- "celltype"
cont@meta.data[["celltype"]] <- paste("Cont",cont@active.ident,sep = "-")
pkd@meta.data[["celltype"]] <- paste("PKD",pkd@active.ident,sep = "-")

#merging ADPKD and control data
rnaAggr <- merge(pkd,cont,add.cell.ids = c("PKD","Cont"))

Idents(rnaAggr) <- "orig.ident"
current.ids<-c("PKD1","PKD2","PKD3","PKD4","PKD5","PKD6","PKD7","PKD8","control1","control2","control3","control4","control5")
new.ids<-c('female','female','female','male','female','male','male','male','male','male','female','male','female')
rnaAggr@meta.data$gender <-plyr::mapvalues(rnaAggr@meta.data$orig.ident,from = current.ids,to=new.ids)

current.ids<-c("PKD1","PKD2","PKD3","PKD4","PKD5","PKD6","PKD7","PKD8","control1","control2","control3","control4","control5")
new.ids<-c('PKD','PKD','PKD','PKD','PKD','PKD','PKD','PKD','control','control','control','control','control')
rnaAggr@meta.data$disease <-plyr::mapvalues(rnaAggr@meta.data$orig.ident,from = current.ids,to=new.ids)

#Integration with Harmony on assay "RNA"
DefaultAssay(rnaAggr) <- "RNA"
rnaAggr <- NormalizeData(rnaAggr)
rnaAggr <- FindVariableFeatures(rnaAggr, selection.method = "vst", nfeatures = 3000)
rnaAggr <- ScaleData(rnaAggr, verbose = FALSE)
rnaAggr <- RunPCA(rnaAggr, verbose = FALSE)
rnaAggr <- RunHarmony(rnaAggr,group.by.vars = c("orig.ident"))
rnaAggr <- RunUMAP(rnaAggr, reduction = "harmony", dims = 1:20)
rnaAggr <- FindNeighbors(rnaAggr, reduction = "harmony", dims = 1:20)
rnaAggr <- FindClusters(rnaAggr, resolution = 0.3)

Idents(rnaAggr) <- "seurat_clusters"
new.cluster.ids <- c("TAL","PT","CNT_PC","DCT","FIB","ENDO","LEUK","FR-PTC","ICA","PEC","PODO","ATL","ICB","URO1","URO2")
names(new.cluster.ids) <- levels(rnaAggr)
rnaAggr <- RenameIdents(rnaAggr, new.cluster.ids)
levels(rnaAggr) <- c("PT","FR-PTC","PEC","ATL","TAL","DCT","CNT_PC","ICA","ICB","PODO","ENDO","FIB","LEUK","URO1","URO2")
#Almost all of "URO1"/"URO2" clusters are from PKD8

rnaAggr@meta.data[["celltype_all"]] <- rnaAggr@active.ident

saveRDS(rnaAggr,"PKDContAggr.rds")

#> sessionInfo()

#R version 4.0.4 (2021-02-15)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS High Sierra 10.13.6

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
#LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] dplyr_1.0.5        harmony_1.0        Rcpp_1.0.6         sctransform_0.3.2  ggplot2_3.3.3     
#[6] Seurat_4.0.0       SeuratObject_4.0.0

#loaded via a namespace (and not attached):
#[1] nlme_3.1-152         matrixStats_0.58.0   RcppAnnoy_0.0.18     RColorBrewer_1.1-2  
#[5] httr_1.4.2           tools_4.0.4          utf8_1.1.4           R6_2.5.0            
#[9] irlba_2.3.3          rpart_4.1-15         KernSmooth_2.23-18   uwot_0.1.10         
#[13] mgcv_1.8-34          DBI_1.1.1            lazyeval_0.2.2       colorspace_2.0-0    
#[17] withr_2.4.1          tidyselect_1.1.0     gridExtra_2.3        compiler_4.0.4      
#[21] plotly_4.9.3         labeling_0.4.2       scales_1.1.1         lmtest_0.9-38       
#[25] spatstat.data_2.0-0  ggridges_0.5.3       pbapply_1.4-3        spatstat_1.64-1     
#[29] goftest_1.2-2        stringr_1.4.0        digest_0.6.27        spatstat.utils_2.0-0
#[33] pkgconfig_2.0.3      htmltools_0.5.1.1    parallelly_1.23.0    limma_3.44.3        
#[37] fastmap_1.1.0        htmlwidgets_1.5.3    rlang_0.4.10         shiny_1.6.0         
#[41] farver_2.1.0         generics_0.1.0       zoo_1.8-9            jsonlite_1.7.2      
#[45] ica_1.0-2            magrittr_2.0.1       patchwork_1.1.1      Matrix_1.3-2        
#[49] munsell_0.5.0        fansi_0.4.2          abind_1.4-5          reticulate_1.18     
#[53] lifecycle_1.0.0      stringi_1.5.3        MASS_7.3-53.1        Rtsne_0.15          
#[57] plyr_1.8.6           grid_4.0.4           parallel_4.0.4       listenv_0.8.0       
#[61] promises_1.2.0.1     ggrepel_0.9.1        crayon_1.4.1         miniUI_0.1.1.1      
#[65] deldir_0.2-10        lattice_0.20-41      cowplot_1.1.1        splines_4.0.4       
#[69] tensor_1.5           pillar_1.5.1         igraph_1.2.6         future.apply_1.7.0  
#[73] reshape2_1.4.4       codetools_0.2-18     leiden_0.3.7         glue_1.4.2          
#[77] data.table_1.14.0    vctrs_0.3.6          png_0.1-7            httpuv_1.5.5        
#[81] gtable_0.3.0         RANN_2.6.1           purrr_0.3.4          polyclip_1.10-0     
#[85] tidyr_1.1.3          scattermore_0.7      future_1.21.0        assertthat_0.2.1    
#[89] mime_0.10            xtable_1.8-4         RSpectra_0.16-0      later_1.1.0.1       
#[93] survival_3.2-7       viridisLite_0.3.0    tibble_3.1.0         cluster_2.1.1       
#[97] globals_0.14.0       fitdistrplus_1.1-3   ellipsis_0.3.1       ROCR_1.0-11  
