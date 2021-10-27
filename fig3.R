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
PTC <- subset(rnaAggr,idents = c("PKD-PT","Cont-PT1","Cont-PT2"))
Idents(PTC) <- "celltype_all"
PTC <- subset(PTC,idents = c("PT","FR-PTC"))

Idents(rnaAggr) <- "disease"
levels(rnaAggr) <- c("PKD","control")
fig3a_1 <- DimPlot(rnaAggr,raster=FALSE) #450x330

DefaultAssay(PTC) <- "RNA"
PTC <- NormalizeData(PTC)
PTC <- FindVariableFeatures(PTC, selection.method = "vst", nfeatures = 3000)
PTC <- ScaleData(PTC, verbose = FALSE)
PTC <- RunPCA(PTC, verbose = FALSE)
PTC <- RunHarmony(PTC,group.by.vars = c("orig.ident"))
PTC <- RunUMAP(PTC, reduction = "harmony", dims = 1:10)
PTC <- FindNeighbors(PTC, reduction = "harmony", dims = 1:10)
PTC <- FindClusters(PTC, resolution = 0.1)
Idents(PTC) <- "disease"
levels(PTC) <- c("PKD","control")
fig3a_2 <- DimPlot(PTC)+NoLegend() #550x500
fig3a_3 <- FeaturePlot(PTC,"VCAM1",pt.size = 0.2)
fig3h_1 <- FeaturePlot(PTC,"TGFB2",pt.size = 0.2)
fig3e_1 <- FeaturePlot(PTC,"RELB",pt.size = 0.2)
sfig5a <- VlnPlot(PTC,"TGFBR1",pt.size = 0)+NoLegend()
sfig5b <- VlnPlot(PTC,"TGFBR2",pt.size = 0)+NoLegend()


Idents(PTC) <- "seurat_clusters"
new.cluster.ids <- c("FR-PTC","N-PTC")
names(new.cluster.ids) <- levels(PTC)
PTC <- RenameIdents(PTC, new.cluster.ids)
levels(PTC) <- c("N-PTC","FR-PTC")
PTC@meta.data[["subtype"]] <- PTC@active.ident
saveRDS(PTC,"PTC.rds")

features <- c("SLC5A12","SLC4A4","LRP2","TPM1","CDH6","PROM1","VCAM1","HAVCR1","CD24","CCL2")
levels(PTC) <- rev(levels(PTC))
fig3b <- DotPlot(PTC, features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #540x500
levels(PTC) <- rev(levels(PTC))

#pseudotime ordering
library(monocle3)
library(dplyr)
library(Matrix)
library(BuenColors)
set.seed(1234)

count_data <- GetAssayData(PTC, assay = "RNA",slot = "counts")

gene_metadata <- as.data.frame(rownames(count_data))
colnames(gene_metadata) <- "gene_short_name"
rownames(gene_metadata) <- gene_metadata$gene_short_name

cds <- new_cell_data_set(as(count_data, "sparseMatrix"),
                         cell_metadata = PTC@meta.data,
                         gene_metadata = gene_metadata)
cds <- preprocess_cds(cds, num_dim =10)
plot_pc_variance_explained(cds)
cds = align_cds(cds, num_dim = 10, alignment_group = "orig.ident")
cds = reduce_dimension(cds,preprocess_method = "Aligned")
cds <- cluster_cells(cds)

cds <- learn_graph(cds,use_partition = T)
cds <- order_cells(cds) 
fig3g_1 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph=T) #png 600x500

genes <- c("PRICKLE1","CDH6")
lineage_cds <- cds[rowData(cds)$gene_short_name %in% genes,]
fig3g_3 <- plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="disease",
                         min_expr=1,
                         cell_size=0.1,
                         trend_formula = "~ splines::ns(pseudotime, df=4)") #520x520

genes <- c("SLC5A12","SLC4A4")
lineage_cds <- cds[rowData(cds)$gene_short_name %in% genes,]
fig3g_2 <- plot_genes_in_pseudotime(lineage_cds,
                                    color_cells_by="disease",
                                    min_expr=1,
                                    cell_size=0.1,
                                    trend_formula = "~ splines::ns(pseudotime, df=4)") #520x520

#Vision
library(VISION)
signatures <- "vision/h.all.v6.2.symbols.gmt"
vision.obj <- Vision(PTC, signatures = signatures)
vision.obj <- analyze(vision.obj)
sigScores <- getSignatureScores(vision.obj)
sigScores <- as.data.frame(sigScores)
PTC <- AddMetaData(PTC,sigScores)

fig3c_1 <- FeaturePlot(PTC,"HALLMARK_OXIDATIVE_PHOSPHORYLATION",cols =jdb_palette("brewer_yes")) #550x500
fig3c_2 <- FeaturePlot(PTC,"HALLMARK_INFLAMMATORY_RESPONSE",cols =jdb_palette("brewer_yes"))
fig3c_3 <- FeaturePlot(PTC,"HALLMARK_MITOTIC_SPINDLE" ,cols =jdb_palette("brewer_yes"))
fig3c_4 <- FeaturePlot(PTC,"HALLMARK_FATTY_ACID_METABOLISM" ,cols =jdb_palette("brewer_yes"))

##################################### ATAC-seq analysis ########################################

library(Seurat)
library(Signac) 
library(monocle3)
library(cicero)
library(BuenColors)
library(EnsDb.Hsapiens.v86)
library(here)
set.seed(1234)

sub_atac <- readRDS("pkd_contATAC_Aggr_sub80.rds")
PCT_ATAC <- subset(sub_atac,idents = c("PCT","PST","FR-PTC"))

DefaultAssay(PCT_ATAC) <- "peaks"
PCT_ATAC <- RunTFIDF(PCT_ATAC)
PCT_ATAC <- FindTopFeatures(PCT_ATAC, min.cutoff = 'q0')
PCT_ATAC <- RunSVD(
  object = PCT_ATAC
)
PCT_ATAC <- RunHarmony(
  object = PCT_ATAC,
  group.by.vars = 'orig.ident',
 reduction = 'lsi',
  assay.use = 'peaks',
 project.dim = FALSE
)

PCT_ATAC <- RunUMAP(PCT_ATAC, dims = 2:20, reduction = 'harmony')
PCT_ATAC <- FindNeighbors(object = PCT_ATAC, reduction = "harmony", dims = 2:20)
PCT_ATAC <- FindClusters(object = PCT_ATAC, verbose = FALSE, algorithm = 3,resolution = 0.1)

DefaultAssay(PCT_ATAC) <- "peaks"
Idents(PCT_ATAC) <- "disease"
levels(PCT_ATAC) <- c("PKD","Cont")

fig3h_2 <- CoveragePlot(
  object = PCT_ATAC,
  region = "chr1-218344787-218348461",
  annotation = T,
  peaks = T,
  extend.upstream = 1000,
  extend.downstream = 1000,
) #500x600

DefaultAssay(PCT_ATAC) <- "chromvar"
fig3e_2 <- VlnPlot(PCT_ATAC,"MA1117.1",pt.size = 0) #RELB 500x500
