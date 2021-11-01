library(Seurat)
library(Signac) 
library(monocle3)
library(cicero)
library(BuenColors)
library(EnsDb.Hsapiens.v86)
library(here)
set.seed(1234)

sub_atac <- readRDS("pkd_contATAC_Aggr_sub80.rds")
CNT_PCatac <- subset(sub_atac,idents = "CNT_PC")

Idents(sub_atac) <- "disease"
levels(sub_atac) <- c("PKD","Cont")
fig6a_1 <- DimPlot(sub_atac,raster=FALSE) #450x330

DefaultAssay(CNT_PCatac) <- "peaks"
CNT_PCatac <- RunUMAP(CNT_PCatac, dims = 2:20, reduction = 'harmony')
CNT_PCatac <- FindNeighbors(object = CNT_PCatac, reduction = "harmony", dims = 2:20)
CNT_PCatac <- FindClusters(object = CNT_PCatac, verbose = FALSE, algorithm = 3,resolution = 0.2)

#Annotation of subtypes
Idents(CNT_PCatac) <- "seurat_clusters"
new.cluster.ids <- c("N-CNT","N-PC","PKD-PC")
names(new.cluster.ids) <- levels(CNT_PCatac)
CNT_PCatac <- RenameIdents(CNT_PCatac, new.cluster.ids)
CNT_PCatac@meta.data[["subtype"]] <- CNT_PCatac@active.ident
saveRDS(CNT_PCatac,"CNT_PCatac.rds")

Idents(CNT_PCatac) <- "subtype"
fig6a_2 <- DimPlot(CNT_PCatac)+NoLegend() #550x440
Idents(CNT_PCatac) <- "disease"
levels(CNT_PCatac) <- c("PKD","Cont")
fig6a_3 <- DimPlot(CNT_PCatac)+NoLegend()

#Dot plot of gene activities
Idents(CNT_PCatac) <- "subtype"
DefaultAssay(CNT_PCatac) <- "RNA"
features <- c("SLC8A1","FXYD2","CALB1","AQP2","SCNN1G",
              "HRH1","CD44","GPRC5A","ROR1")

levels(CNT_PCatac) <- rev(levels(CNT_PCatac))
fig6b <- DotPlot(CNT_PCatac, features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #600x450
levels(CNT_PCatac) <- rev(levels(CNT_PCatac))

DefaultAssay(CNT_PCatac) <- "chromvar"
fig6c_1 <- FeaturePlot(CNT_PCatac,"MA0105.3",cols =jdb_palette("brewer_yes")) #550x500
fig6c_2 <- FeaturePlot(CNT_PCatac,"MA0107.1",cols =jdb_palette("brewer_yes")) #550x500
fig6c_3 <- FeaturePlot(CNT_PCatac,"MA0808.1",cols =jdb_palette("brewer_yes")) #550x500 tead3

fig6f_1 <- FeaturePlot(CNT_PCatac,"MA0018.3",cols =jdb_palette("brewer_yes")) #550x500 CREB1
fig6f_2 <- FeaturePlot(CNT_PCatac,"MA1149.1",cols =jdb_palette("brewer_yes")) #550x500 RARA::RXRG
levels(CNT_PCatac) <- c("PKD","Cont")
fig6f_3 <- VlnPlot(CNT_PCatac,"MA0018.3",pt.size = 0) #550x500 CREB1
fig6f_4 <- VlnPlot(CNT_PCatac,"MA1149.1",pt.size = 0) #550x500 RARA::RXRG

DefaultAssay(CNT_PCatac) <- "peaks"
Idents(CNT_PCatac) <- "subtype"
CNT_PCatac_peakMarker <- FindAllMarkers(CNT_PCatac,only.pos = T,min.pct = 0.05) 
test <- ClosestFeature(CNT_PCatac, regions = CNT_PCatac_peakMarker$gene)
CNT_PCatac_peakMarker <- cbind(CNT_PCatac_peakMarker,test[,c(2,8)])

DefaultAssay(CNT_PCatac) <- "chromvar"
CNT_PCatac_chromvarMarker <- FindAllMarkers(CNT_PCatac,only.pos = T,logfc.threshold = 0.25)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(TFBSTools)
pfm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = F)
)

x <- NULL                        
for (i in CNT_PCatac_chromvarMarker$gene) { 
  x <- c(x,pfm@listData[[i]]@name)
}
CNT_PCatac_chromvarMarker$gene_name <- x

#Coverage plots

DefaultAssay(CNT_PCatac) <- "peaks"
Idents(CNT_PCatac) <- "subtype"

fig6e_2 <- CoveragePlot(
  object = CNT_PCatac,
  region = "chr12-12890233-12893174",
  annotation = T,
  peaks = T,
  group.by = "subtype",
  extend.upstream = 1000,
  extend.downstream = 1000,
) #330x530

fig6e_3 <- CoveragePlot(
  object = CNT_PCatac,
  region = "chr12-12871973-12873059",
  annotation = T,
  peaks = T,
  group.by = "subtype",
  extend.upstream = 1000,
  extend.downstream = 1000,
) #330x530

figS7d_2 <- CoveragePlot(
  object = CNT_PCatac,
  region = "chr9-21994180-21996144",
  annotation = T,
  peaks = T,
  extend.upstream = 3000,
  extend.downstream = 3000,
) #600x600

figS7d_3 <- CoveragePlot(
  object = CNT_PCatac,
  region = "chr9-21684839-21685331",
  annotation = T,
  peaks = T,
  extend.upstream = 3000,
  extend.downstream = 3000,
) #600x600

figS7d_4 <- CoveragePlot(
  object = CNT_PCatac,
  region = "chr9-21559035-21560436",
  annotation = T,
  peaks = T,
  extend.upstream = 3000,
  extend.downstream = 3000,
) #600x600

#CCAN
Idents(sub_atac) <- "disease"
pkd <- subset(sub_atac,idents = "PKD")
DefaultAssay(pkd) <- "peaks"
count_data <- GetAssayData(pkd, slot = "counts")
summ <- summary(count_data)
summ_frame <- data.frame(peak = rownames(count_data)[summ$i],
                         cell.id = colnames(count_data)[summ$j],
                         count = summ$x)

# create cell data set object with cicero constructor
input_cds <- make_atac_cds(summ_frame, binarize = F)
meta <- pkd@meta.data
meta$cells <- rownames(meta)
metanames <- rownames(meta)
meta <- merge(input_cds@colData,meta,by.x="cells", by.y="cells")
rownames(meta) <- meta@listData[["cells"]]
input_cds@colData <- meta
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

contigs <- read.table("/contigLengths.txt")
contigs$V1 <- paste0("chr",contigs$V1)
contigs <- contigs[,c(1,5)]
contigs <- subset(contigs, V1 %in% "chr12")
conns <- run_cicero(cicero_cds, contigs)

gene_annotation <- read.table("/human_hg38_annotation.txt")
gene_annotation$feature <- as.character(gene_annotation$feature)

fig6e_1 <- plot_connections(conns, "chr12", 12870000, 12895000,
                 gene_model = gene_annotation,
                 coaccess_cutoff = .2, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" ) #600x300

#CCAN for p16/MIR31HG (chr9)
contigs <- read.table("/contigLengths.txt")
contigs$V1 <- paste0("chr",contigs$V1)
contigs <- contigs[,c(1,5)]
contigs <- subset(contigs, V1 %in% "chr9")
conns <- run_cicero(cicero_cds, contigs)

figS9d_1 <- plot_connections(conns, "chr9", 21400000,22050000,
                 gene_model = gene_annotation,
                 coaccess_cutoff = .2, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" ) #600x300

