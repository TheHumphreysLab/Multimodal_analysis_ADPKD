library(Seurat) 
library(Signac) 
library(EnsDb.Hsapiens.v86) 
library(openxlsx) 

Idents(sub_atac) <- "celltype"
idents <- levels(sub_atac)

######## DAR ##########
DefaultAssay(sub_atac) <- "peaks"

DiseaseMarkers <- function(cluster, seurat_obj, test, ref, meta_group) {
  print(paste0("Comparing DAR for: ",cluster))
  atac <- seurat_obj
  group <- dplyr::select(atac@meta.data, all_of(meta_group))[,1]
  atac@meta.data$celltype.stim <- paste0(atac@meta.data$celltype,"_", group)
  Idents(atac) <- "celltype.stim"
  dar <- FindMarkers(atac, 
                     ident.1 = paste0(cluster,"_", test),  
                     ident.2 = paste0(cluster, "_", ref),
                     test.use = 'LR', 
                     latent.vars = "peak_region_fragments")
  cf <- ClosestFeature(atac, rownames(dar))
  return(cbind(dar, gene=cf$gene_name, distance=cf$distance))
}

# compare celltype DAR for ADPKD vs. control
list.disease.dar <- lapply(idents, function(x) {DiseaseMarkers(x, seurat_obj = sub_atac, test="PKD", ref="Cont", meta_group="disease")})
write.xlsx(list.disease.dar, file = "disease.dar.xlsx", sheetName = idents, rowNames = T)


######## Motif ##########
DefaultAssay(sub_atac) <- "chromvar"

library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(TFBSTools)
pfm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = T)
)

DiseaseMarkers <- function(cluster, seurat_obj, test, ref, meta_group) {
  print(paste0("Comparing ADPKD for: ",cluster))
  atac <- seurat_obj
  group <- dplyr::select(atac@meta.data, all_of(meta_group))[,1]
  atac@meta.data$celltype.stim <- paste0(atac@meta.data$celltype,"_", group)
  Idents(atac) <- "celltype.stim"
  dma <- FindMarkers(atac, 
                              ident.1 = paste0(cluster,"_", test),  
                              ident.2 = paste0(cluster, "_", ref),
                              logfc.threshold = 0.1,
                              mean.fxn = rowMeans,
                              fc.name = "avg_diff")
  motifLookup <- rownames(dma)
  motifNames <- sapply(motifLookup, function(x) pfm@listData[[x]]@name)
  return(cbind(dma, gene = motifNames))
}

list.disease.dma <- lapply(idents, function(x) {DiseaseMarkers(x, seurat_obj = sub_atac, test="PKD", ref="Cont", meta_group="disease")})                          
write.xlsx(list.disease.dma, file = "disease.dma.xlsx", sheetName = idents, rowNames = T)


######## Gene activities ##########

DefaultAssay(sub_atac) <- "RNA"
sub_atac <- NormalizeData(sub_atac)

DiseaseMarkers <- function(cluster, seurat_obj, test, ref, meta_group) {
  print(paste0("Comparing DEG for: ",cluster))
  atac <- seurat_obj
  group <- dplyr::select(atac@meta.data, all_of(meta_group))[,1]
  atac@meta.data$celltype.stim <- paste0(atac@meta.data$celltype,"_", group)
  Idents(atac) <- "celltype.stim"
  deg <- FindMarkers(atac, 
                              ident.1 = paste0(cluster,"_", test),  
                              ident.2 = paste0(cluster, "_", ref),
                              logfc.threshold = 0.25)
  return(deg)
}


list.disease.dga <- lapply(idents, function(x) {DiseaseMarkers(x, seurat_obj = sub_atac, test="PKD", ref="Cont", meta_group="disease")})                          
write.xlsx(list.disease.dga, file = "disease_ga.xlsx", sheetName = idents, rowNames = T)
