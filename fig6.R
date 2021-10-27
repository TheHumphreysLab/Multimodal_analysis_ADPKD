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
CNT_PC <- subset(rnaAggr,idents = c("PKD-CNT_PC","Cont-CNT_PC","PKD-PC"))
Idents(CNT_PC) <- "celltype_all"
CNT_PC <- subset(CNT_PC,idents = "CNT_PC")

CNT_PC <- NormalizeData(CNT_PC)
CNT_PC <- FindVariableFeatures(CNT_PC, selection.method = "vst", nfeatures = 3000)
CNT_PC <- ScaleData(CNT_PC, verbose = FALSE)
CNT_PC <- RunPCA(CNT_PC, verbose = FALSE)
CNT_PC <- RunHarmony(CNT_PC,group.by.vars = c("orig.ident"))
CNT_PC <- RunUMAP(CNT_PC, reduction = "harmony", dims = 1:20)
CNT_PC <- FindNeighbors(CNT_PC, reduction = "harmony", dims = 1:20)
CNT_PC <- FindClusters(CNT_PC, resolution = 0.3)
fig6a_2 <- DimPlot(CNT_PC,label = T,repel=T)+NoLegend() #550x440

Idents(CNT_PC) <- "seurat_clusters"
new.cluster.ids <- c("PKD-PC1","N-CNT","N-PC","LowQC1","PKD-PC2","PKD-CNT","LowQC2","PKD-PC3")
names(new.cluster.ids) <- levels(CNT_PC)
CNT_PC <- RenameIdents(CNT_PC, new.cluster.ids)
levels(CNT_PC) <- c("N-CNT","PKD-CNT","N-PC","PKD-PC1","PKD-PC2","PKD-PC3","LowQC1","LowQC2")
CNT_PC@meta.data[["subtype"]] <- CNT_PC@active.ident
saveRDS(CNT_PC,"CNT_PC.rds")

Idents(CNT_PC) <- "subtype"
DimPlot(CNT_PC,label = T)+NoLegend()
Idents(CNT_PC) <- "disease"
fig6a_3 <- DimPlot(CNT_PC) #450x330

marker <- FindAllMarkers(CNT_PC,only.pos = T,logfc.threshold = 0.25)

Idents(CNT_PC) <- "subtype"
CNT_PC2 <- subset(CNT_PC,idents = c("LowQC1","LowQC2"),invert=T)

features <- c("SLC8A1","CALB1","MET","AQP2","SCNN1G","CFTR",
              "HRH1","CD44","ALOX5","GPRC5A","MIR31HG",
              "LCN2","MFSD2A","CXCL2",
              "PTGFR","ROBO2","DOCK10","LIN7A")

features2 <- c("CDKN1A","CDKN1B","CDKN1C","CDKN2A","CDKN2B","CDKN2C","CDKN2D","CDKN3")

levels(CNT_PC2) <- rev(levels(CNT_PC2))
fig6b <- DotPlot(CNT_PC2, features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #700x450
figS7a <- DotPlot(CNT_PC2, features = features2, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #640x500
levels(CNT_PC2) <- rev(levels(CNT_PC2))

#Vision
library(VISION)
signatures <- "/Users/mutouyoshiharu/Desktop/PKD_seq/Fig_v3/fig.3/vision/h.all.v6.2.symbols.gmt"
vision.obj <- Vision(CNT_PC, signatures = signatures)
vision.obj <- analyze(vision.obj)

sigScores <- getSignatureScores(vision.obj)
sigScores <- as.data.frame(sigScores)
CNT_PC <- AddMetaData(CNT_PC,sigScores2)

fig6c_1 <- FeaturePlot(CNT_PC,"HALLMARK_GLYCOLYSIS",cols =jdb_palette("brewer_yes")) #550x500
fig6c_2 <- FeaturePlot(CNT_PC,"HALLMARK_OXIDATIVE_PHOSPHORYLATION",cols =jdb_palette("brewer_yes"))
fig6d_3 <- FeaturePlot(CNT_PC,"HALLMARK_INFLAMMATORY_RESPONSE",cols =jdb_palette("brewer_yes"))
fig6d_4 <- FeaturePlot(CNT_PC,"HALLMARK_IL6_JAK_STAT3_SIGNALING" ,cols =jdb_palette("brewer_yes"))

#Pseudotime trajectory
library(monocle3)
library(dplyr)
library(Matrix)
library(BuenColors)
set.seed(1234)

Idents(CNT_PC) <- "disease"
CNT_PC2 <- subset(CNT_PC,idents = "PKD")
Idents(CNT_PC2) <- "subtype"
CNT_PC2 <- subset(CNT_PC2,idents = c("LowQC1","LowQC2"),invert=T)

count_data <- GetAssayData(CNT_PC2, assay = "RNA",slot = "counts")

gene_metadata <- as.data.frame(rownames(count_data))
colnames(gene_metadata) <- "gene_short_name"
rownames(gene_metadata) <- gene_metadata$gene_short_name

cds <- new_cell_data_set(as(count_data, "sparseMatrix"),
                         cell_metadata = CNT_PC2@meta.data,
                         gene_metadata = gene_metadata)
cds <- preprocess_cds(cds, num_dim =20)
plot_pc_variance_explained(cds)
cds = align_cds(cds, num_dim = 20, alignment_group = "orig.ident")
cds = reduce_dimension(cds,preprocess_method = "Aligned")
fig6e_1 <- plot_cells(cds, color_cells_by="subtype", 
                      group_label_size = 4,show_trajectory_graph = F,
                      label_cell_groups =F) #png: 600x500
cds <- cluster_cells(cds)
cds <- learn_graph(cds,use_partition = T)
cds <- order_cells(cds) 
fig6e_2 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph=T) #png 600x500

#subset branch toward PKD-PC1
cds_sub1 <- choose_graph_segments(cds,clear_cds = F)

genes <- c("SLC8A1","GPRC5A")
lineage_cds <- cds_sub1[rowData(cds_sub1)$gene_short_name %in% genes,]
fig6f_1 <- plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="subtype",
                         min_expr=1,
                         cell_size=0.1,
                         trend_formula = "~ splines::ns(pseudotime, df=4)") #480x560

#subset branch toward PKD-PC2
cds_sub2 <- choose_graph_segments(cds,clear_cds = F)

genes <- c("SLC8A1","MIR31HG")
lineage_cds <- cds_sub2[rowData(cds_sub2)$gene_short_name %in% genes,]
fig6f_2 <- plot_genes_in_pseudotime(lineage_cds,
                                    color_cells_by="subtype",
                                    min_expr=1,
                                    cell_size=0.1,
                                    trend_formula = "~ splines::ns(pseudotime, df=4)") #480x560


fig6g <- FeaturePlot(CNT_PC,"GPRC5A") #550x440

figS7b <- FeaturePlot(CNT_PC,"CDKN2A",pt.size = 0.5,order=T) #550x440
figS7c <- FeaturePlot(CNT_PC,"MIR31HG",pt.size = 0.5,order=T) #550x440

FigS6b <- FeaturePlot(rnaAggr,"ROR1",split.by = "disease") #1000x500
Idents(rnaAggr) <- "disease"
levels(rnaAggr) <- c("PKD","control")
FigS6c <- VlnPlot(rnaAggr,"ROR1",pt.size = 0.1) #700x550

FigS8a <- FeaturePlot(rnaAggr,"RDH10") #600x550
FigS8b <- FeaturePlot(CNT_PC,"RDH10")

