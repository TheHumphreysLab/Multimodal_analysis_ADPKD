library(Signac)
library(GenomeInfoDb)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(Seurat)

#scRNA-seq data (after preprocessing)
atacAggr <- readRDS("pkd_contATACaggr.rds")

#scRNA-seq data of control / ADPKD kidneys (after clean-up)
contRNA <- readRDS("controlAggr_filtered.rds")
pkdRNA <- readRDS("PKDAggr_filtered.rds")
DefaultAssay(contRNA) <- "RNA"
DefaultAssay(pkdRNA) <- "RNA"
contRNA <- FindVariableFeatures(contRNA, nfeatures = 4000)
pkdRNA <- FindVariableFeatures(pkdRNA, nfeatures = 4000)

# Define low-resolution celltype identities for thresholding snATAC-seq data after label-transfer ------------------------------------------------------
#(e.g., group PT and FR-PTC, or distal nephron together)

## for control snRNA-seq data
Idents(contRNA) <- "celltype"
new.cluster.ids <- c("PT","PT","PT","TAL","TAL","TAL","DCT_CNT_PC","DCT_CNT_PC","IC","IC","PODO","ENDO","ENDO","ENDO","FIB")
names(new.cluster.ids) <- levels(contRNA)
contRNA <- RenameIdents(contRNA, new.cluster.ids)
contRNA@meta.data[["low_res_celltype"]] <- contRNA@active.ident

## for ADPKD snRNA-seq data
Idents(pkdRNA) <- "celltype"
new.cluster.ids <- c("PT","PT","Transitional","TAL","DCT_CNT_PC","DCT_CNT_PC","DCT_CNT_PC",
                     "IC","PODO","ENDO","ENDO","ENDO","FIB","FIB","FIB","LEUK","LEUK")
names(new.cluster.ids) <- levels(pkdRNA)
pkdRNA <- RenameIdents(pkdRNA, new.cluster.ids)
pkdRNA@meta.data[["low_res_celltype"]] <- pkdRNA@active.ident

#Divide snATAC data for label-transferring separately control and ADPKD data  ---------------------------------------------------------------------------
Idents(atacAggr) <- "orig.ident"
contATAC <- subset(atacAggr,idents = c("Cont1","Cont2","Cont3","Cont4","Cont5"))
pkdATAC <- subset(atacAggr,idents = c("PKD1","PKD2","PKD3","PKD4",
                                      "PKD5","PKD6","PKD7","PKD8"))
contATAC@meta.data[["disease"]] <- "control"
pkdATAC@meta.data[["disease"]] <- "PKD"

# Label transfer ----------------------------------------------------------------------------------------------------------------------------------------

##Control data label transfer
DefaultAssay(contATAC) <- "RNA"
contATAC<- NormalizeData(contATAC)
contATAC <- ScaleData(contATAC)
transfer.anchors_cont <- FindTransferAnchors(reference = contRNA, query = contATAC, features = VariableFeatures(object = contRNA), 
                                             reference.assay = "RNA", query.assay = "RNA", reduction = "cca")

##ADPKD data label transfer
DefaultAssay(pkdATAC) <- "RNA"
pkdATAC<- NormalizeData(pkdATAC)
pkdATAC <- ScaleData(pkdATAC)
transfer.anchors_pkd <- FindTransferAnchors(reference = pkdRNA, query = pkdATAC, features = VariableFeatures(object = pkdRNA), 
                                            reference.assay = "RNA", query.assay = "RNA", reduction = "cca")


#Celltype prediction using low resolution and high resolution celltypes of snRNA-seq data -------------------------------------------------------------- 

celltype.predictions_cont_highres <- TransferData(anchorset = transfer.anchors_cont, refdata = contRNA$celltype, 
                                          weight.reduction = contATAC[["lsi"]], dims = 2:30)

celltype.predictions_cont_lowres <- TransferData(anchorset = transfer.anchors_cont, refdata = contRNA$celltype_lowres, 
                                                  weight.reduction = contATAC[["lsi"]], dims = 2:30)

celltype.predictions_pkd_highres <- TransferData(anchorset = transfer.anchors_pkd, refdata = pkdRNA$celltype, 
                                                  weight.reduction = pkdATAC[["lsi"]], dims = 2:30)

celltype.predictions_pkd_lowres <- TransferData(anchorset = transfer.anchors_pkd, refdata = pkdRNA$celltype_lowres, 
                                                 weight.reduction = pkdATAC[["lsi"]], dims = 2:30)


##Dataframe for high/low-resolution predicted celltypes for both control and ADPKD samples
celltype.predictions_all_highres <- rbind(celltype.predictions_cont_highres[,c(1,17)]
                                          ,celltype.predictions_pkd_highres[,c(1,19)])
colnames(celltype.predictions_all_highres) <- c("highres_predicted.id","highres_prediction.score.max")
celltype.predictions_all_lowres <- rbind(celltype.predictions_cont_lowres[,c(1,9)]
                                          ,celltype.predictions_pkd_lowres[,c(1,11)])
colnames(celltype.predictions_all_lowres) <- c("lowres_predicted.id","lowres_prediction.score.max")
prediction <-cbind(celltype.predictions_all_lowres,celltype.predictions_all_highres)

##Add dataframe of predicted celltypes to integrated object as metadata
atacAggr <- AddMetaData(atacAggr,prediction)

#Obtain high-confidence cognate cell type nuclei in snATAC-seq by thresholding low-resolution prediction score ----------------------------------------------
#use max prediction score for low resolution celltypes from snRNA-seq data)
#See Supplementary Fig.4 for detail

sub_atac <- subset(atacAggr,lowres_prediction.score.max > 0.8) #62802 => 50986

# Finalize integration step ---------------------------------------------------------------------------------------------------------------------------------

## Add disease info (control vs ADPKD) to predicted celltypes
current.ids<-levels(sub_atac@meta.data[["orig.ident"]])
new.ids<-c("Cont","Cont","Cont","Cont","Cont","PKD","PKD","PKD","PKD","PKD","PKD","PKD","PKD")
sub_atac@meta.data$disease<-plyr::mapvalues(sub_atac@meta.data$orig.ident,from = current.ids,to=new.ids)

sub_atac@meta.data[["highres_predicted.id"]] <- paste(sub_atac@meta.data[["disease"]],sub_atac@meta.data[["highres_predicted.id"]],sep = "-")
sub_atac@meta.data[["lowres_predicted.id"]] <- paste(sub_atac@meta.data[["disease"]],sub_atac@meta.data[["lowres_predicted.id"]],sep = "-")

## batch correction
library(harmony)

sub_atac <- RunHarmony(
  object = sub_atac,
  group.by.vars = 'orig.ident',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)

## Clustering and visualization
sub_atac <- RunUMAP(sub_atac, dims = 2:20, reduction = 'harmony')
sub_atac <- FindNeighbors(object = sub_atac, reduction = "harmony", dims = 2:20)
sub_atac <- FindClusters(object = sub_atac, verbose = FALSE, algorithm = 3,resolution = 0.3)
DimPlot(sub_atac,  pt.size = 0.1,label = T)+NoLegend() 

## Celltype annotation on unsupervised clustering
Idents(sub_atac) <- "seurat_clusters"
new.cluster.ids <- c("TAL","PCT","CNT_PC","FR-PTC","FIB","PST","DCT","ENDO","IC","LEUK","PEC_PODO")
names(new.cluster.ids) <- levels(sub_atac)
sub_atac <- RenameIdents(sub_atac, new.cluster.ids)
levels(sub_atac) <- c("PCT","PST","FR-PTC","PEC_PODO","TAL","DCT","CNT_PC","IC","ENDO","FIB","LEUK")
sub_atac@meta.data[["celltype"]] <- sub_atac@active.ident

saveRDS(sub_atac,"pkd_contATAC_sub80.rds")
