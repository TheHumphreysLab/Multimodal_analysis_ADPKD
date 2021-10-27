library(Seurat) 
library(openxlsx)
set.seed(1234)

DiseaseMarkers <- function(cluster, seurat_obj, test, ref, meta_group) {
  print(paste0("Comparing DEG for: ",cluster))
  rnaAggr <- seurat_obj
  group <- dplyr::select(rnaAggr@meta.data, all_of(meta_group))[,1]
  rnaAggr@meta.data$celltype.stim <- paste0(rnaAggr@meta.data$celltype_all,"_", group)
  Idents(rnaAggr) <- "celltype.stim"
  deg <- tryCatch(FindMarkers(rnaAggr, 
                     ident.1 = paste0(cluster,"_", test),  
                     ident.2 = paste0(cluster, "_", ref),
                     logfc.threshold = 0.25),
                  error=function(e) NULL)
  return(deg)
}

idents <- levels(rnaAggr@meta.data$celltype_all)

list.disease.deg <- lapply(idents, function(x) {DiseaseMarkers(x, seurat_obj = rnaAggr, test="PKD", ref="control", meta_group="disease")})                          
write.xlsx(list.disease.deg, file = "disease_deg.xlsx", sheetName = idents, rowNames = T)


Idents(rnaAggr) <- "disease"
cont <- subset(rnaAggr,idents = "control")
pkd <- subset(rnaAggr,idents = "PKD")


Idents(cont) <- "celltype_all"
Idents(pkd) <- "celltype_all"
idents_cont <- levels(cont)
idents_pkd <- levels(pkd)

CelltypeMarkers <- function(cluster, seurat_obj) {
  print(paste0("Finding DGE for: ",cluster))
  rnaAggr <- seurat_obj
  deg <- FindMarkers(rnaAggr, 
                     ident.1 = cluster,    
                     min.pct = 0.2)
  return(deg)
}

list.cluster.deg_cont <- lapply(idents_cont, function(x) {CelltypeMarkers(x, seurat_obj = cont)})
list.cluster.deg_pkd <- lapply(idents_pkd, function(x) {CelltypeMarkers(x, seurat_obj = pkd)})

write.xlsx(list.cluster.deg_cont, file = "cont_celltype.deg.xlsx", sheetName = idents_cont, rowNames = T)
write.xlsx(list.cluster.deg_pkd, file = "pkd_celltype.deg.xlsx", sheetName = idents_pkd, rowNames = T)
