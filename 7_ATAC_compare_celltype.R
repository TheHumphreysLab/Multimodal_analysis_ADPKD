library(Seurat) 
library(Signac) 
library(EnsDb.Hsapiens.v86) 
library(openxlsx) 

Idents(sub_atac) <- "disease"
cont <- subset(sub_atac,idents = "Cont")
pkd <- subset(sub_atac,idents = "PKD")

Idents(cont) <- "celltype"
Idents(pkd) <- "celltype"
idents_cont <- levels(cont)
idents_pkd <- levels(pkd)

############ Peaks ##############

DefaultAssay(cont) <- "peaks"
DefaultAssay(pkd) <- "peaks"

# wrapper functions for FindMarkers 
CelltypeMarkers <- function(cluster, seurat_obj) {
  print(paste0("Finding DAR for: ",cluster))
  sub_atac <- seurat_obj
  Idents(sub_atac) <- "celltype"
  dar <- FindMarkers(sub_atac, 
                     ident.1 = cluster,    
                     test.use = 'LR', 
                     latent.vars = "peak_region_fragments",
                     min.pct = 0.2) 
  cf <- ClosestFeature(sub_atac, rownames(dar))
  return(cbind(dar, gene=cf$gene_name, distance=cf$distance))
}

# identify celltype DAR

list.celltype.dar_cont <- lapply(idents_cont, function(x) {CelltypeMarkers(x, seurat_obj = cont)})
list.celltype.dar_pkd <- lapply(idents_pkd, function(x) {CelltypeMarkers(x, seurat_obj = pkd)})

write.xlsx(list.celltype.dar_cont, file = "cont_celltype.dar.xlsx", sheetName = idents, rowNames = T)
write.xlsx(list.celltype.dar_pkd, file = "pkd_celltype.dar.xlsx", sheetName = idents, rowNames = T)

############ Motifs ##############

DefaultAssay(cont) <- "chromvar"
DefaultAssay(pkd) <- "chromvar"

library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(TFBSTools)
pfm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = T)
)

# wrapper functions for FindMarkers 
CelltypeMarkers <- function(cluster, seurat_obj) {
  print(paste0("Finding motif for: ",cluster))
  atac <- seurat_obj
  dma <- FindMarkers(atac, 
                     ident.1 = cluster,    
                     logfc.threshold = 0.25,
                     mean.fxn = rowMeans,
                     fc.name = "avg_diff")
  motifLookup <- rownames(dma)
  motifNames <- sapply(motifLookup, function(x) pfm@listData[[x]]@name)
  return(cbind(dma, gene = motifNames))
}

list.cluster.dma_cont <- lapply(idents_cont, function(x) {CelltypeMarkers(x, seurat_obj = cont)})
list.cluster.dma_pkd <- lapply(idents_pkd, function(x) {CelltypeMarkers(x, seurat_obj = pkd)})


write.xlsx(list.cluster.dma_cont, file = "cont_celltype.dma.xlsx", sheetName = idents_cont, rowNames = T)
write.xlsx(list.cluster.dma_pkd, file = "pkd_celltype.dma.xlsx", sheetName = idents_pkd, rowNames = T)


############ Gene activity ##############

DefaultAssay(cont) <- "RNA"
DefaultAssay(pkd) <- "RNA"
cont <- NormalizeData(cont)
pkd <- NormalizeData(pkd)

# wrapper functions for FindMarkers 
CelltypeMarkers <- function(cluster, seurat_obj) {
  print(paste0("Finding DGA for: ",cluster))
  atac <- seurat_obj
  deg <- FindMarkers(atac, 
                     ident.1 = cluster,    
                     min.pct = 0.2)
  return(deg)
}

list.cluster.dga_cont <- lapply(idents_cont, function(x) {CelltypeMarkers(x, seurat_obj = cont)})
list.cluster.dga_pkd <- lapply(idents_pkd, function(x) {CelltypeMarkers(x, seurat_obj = pkd)})

write.xlsx(list.cluster.dga_cont, file = "cont_celltype.dga.xlsx", sheetName = idents_cont, rowNames = T)
write.xlsx(list.cluster.dga_pkd, file = "pkd_celltype.dga.xlsx", sheetName = idents_pkd, rowNames = T)
