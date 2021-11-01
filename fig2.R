library(Seurat) #3.2.0
library(Signac)
library(ggplot2)
library(sctransform)
library(harmony)
library(Rcpp)
library(dplyr)
set.seed(1234)

sub_atac <- readRDS("pkd_contATAC_Aggr_sub80.rds")

Idents(sub_atac) <- "seurat_clusters"
fig2a <- DimPlot(sub_atac) #750x630
Idents(sub_atac) <- "disease"
levels(sub_atac) <- c("PKD","Cont")
fig2b <- DimPlot(sub_atac)+NoLegend() #700x630

sub_atac@meta.data$celltype_disease <- paste0(sub_atac@meta.data$celltype,"_",sub_atac@meta.data$disease)
Idents(sub_atac) <- "celltype_disease"
levels(sub_atac) <- c("PCT_Cont","PCT_PKD","PST_Cont","PST_PKD","FR-PTC_Cont","FR-PTC_PKD",
                     "PEC_PODO_Cont","PEC_PODO_PKD","TAL_Cont","TAL_PKD",
                     "DCT_Cont","DCT_PKD","CNT_PC_Cont","CNT_PC_PKD",
                     "IC_Cont","IC_PKD","ENDO_Cont","ENDO_PKD",
                     "FIB_Cont","FIB_PKD","LEUK_Cont","LEUK_PKD")

features <- c("SLC34A1","SLC36A2","SLC5A2","SLC5A1","SLC6A18","VCAM1","CRYAB","CFH","NPHS1","NPHS2",
              "SLC12A1","UMOD","PAPPA2","TCF24","SLC12A3","SALL3","AQP2","SCNN1G","SCNN1B","FAM20A","SLC26A4","SLC26A7",
              "KIT","FLT1","EMCN","ACTA2","MYH11","PDGFRA","GLI2","GLI1","PTPRC","SLC11A1","CSF2RA")

DefaultAssay(sub_atac) <- "RNA"
levels(sub_atac) <- rev(levels(sub_atac))
fig2c <- DotPlot(sub_atac, features = features, cols = c("lightyellow","#992FA0")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
levels(sub_atac) <- rev(levels(sub_atac)) #1400x750



sub_atac@assays[["peaks"]]@fragments[[1]]@path <- "/path/to/fragments.tsv.gz"
DefaultAssay(sub_atac) <- "peaks"
fig2d_1 <- CoveragePlot(
  object = sub_atac,
  region = "chr2-169361500-169363500",
  annotation = T,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0,
) #320x700 #LRP2 prom

fig2d_2 <- CoveragePlot(
  object = sub_atac,
  region = "chr22-32042000-32044000",
  annotation = T,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0
) #320x700 #SLC5A1 prom

fig2d_3 <- CoveragePlot(
  object = sub_atac,
  region = "chr11-32434800-32436800",
  annotation = T,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0
) #320x700 #WT1 prom

fig2d_4 <- CoveragePlot(
  object = sub_atac,
  region = "chr16-20355500-20357500",
  annotation = T,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0
) #320x700 #UMOD prom

fig2d_5 <- CoveragePlot(
  object = sub_atac,
  region = "chr16-56864000-56866000",
  annotation = T,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0
) #320x700 #SLC12A3 prom

fig2d_6 <- CoveragePlot(
  object = sub_atac,
  region = "chr12-49949500-49951500",
  annotation = T,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0
) #320x700 #AQP2 prom

fig2d_7 <- CoveragePlot(
  object = sub_atac,
  region = "chr4-54657000-54659000",
  annotation = T,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0
) #320x700 #KIT prom

fig2d_8 <- CoveragePlot(
  object = sub_atac,
  region = "chr13-28494500-28496500",
  annotation = T,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0
) #320x700 #FLT1 prom

fig2d_9 <- CoveragePlot(
  object = sub_atac,
  region = "chr11-117198000-117200000",
  annotation = T,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0
) #320x700 #TAGLN

fig2d_10 <- CoveragePlot(
  object = sub_atac,
  region = "chr1-198638000-198640000",
  annotation = T,
  peaks = F,
  extend.upstream = 0,
  extend.downstream = 0
) #320x700 #PTPRC

###########################################################
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(org.Hs.eg.db)

Idents(sub_atac) <- "disease"
cont <- subset(sub_atac,idents = "Cont")
pkd <- subset(sub_atac,idents = "PKD")

Idents(cont) <- "celltype"
Idents(pkd) <- "celltype"
idents_cont <- levels(cont)
idents_pkd <- levels(pkd)

Idents(sub_atac) <- "celltype"

darfile <- "/path/to/Supplementary_Data6.xlsx"
idents <- getSheetNames(darfile)

#ADPKD 

list.dar <- lapply(idents, function(x) {
  df <- read.xlsx(darfile, sheet = x, rowNames = T) %>%
    rownames_to_column(var = "coord") %>%
    dplyr::filter(avg_log2FC > 0) %>%
    dplyr::mutate(celltype = x) # annotate each region with its corresponding celltype
})

# convert the DAR to GRanges objects to annotate
list.dar.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$coord, sep = c(":","-"))
  return(gr)
})
names(list.dar.gr) <- idents

# annotate the list of GRanges DAR for each cell type
list.peakAnno <- lapply(list.dar.gr, annotatePeak, TxDb = txdb,
                        tssRegion = c(-3000, 3000), verbose = FALSE)
fig2e_1 <- plotAnnoBar(list.peakAnno) #660x430

#Control
list.dar <- lapply(idents, function(x) {
  df <- read.xlsx(darfile, sheet = x, rowNames = T) %>%
    rownames_to_column(var = "coord") %>%
    dplyr::filter(avg_log2FC < 0) %>%
    dplyr::mutate(celltype = x) 
})


# convert the DAR to GRanges objects to annotate
all_dar.gr <- StringToGRanges(all_dar$coord, sep = c(":","-"))
list.dar.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$coord, sep = c(":","-"))
  return(gr)
})
names(list.dar.gr) <- idents

# annotate the list of GRanges DAR for each cell type
list.peakAnno <- lapply(list.dar.gr, annotatePeak, TxDb = txdb,
                        tssRegion = c(-3000, 3000), verbose = FALSE)

fig2e_2 <- plotAnnoBar(list.peakAnno) #660x430

