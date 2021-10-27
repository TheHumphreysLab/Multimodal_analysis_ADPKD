library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(harmony)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(tibble)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(TFBSTools)
library(patchwork)
set.seed(1234)

# snATACseq data obtained with cellrangerATAC v1.2.0
counts <- Read10X_h5(here("outs_atac/filtered_peak_bc_matrix.h5"))
metadata <- read.csv(here("outs_atac/singlecell.csv"), header = TRUE, row.names = 1)
aggcsv <- read.csv(here("outs_atac/aggregation_csv.csv"),header = TRUE, row.names = 1)

# create Seurat object
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = "outs_atac/fragments.tsv.gz",
  min.cells = 5,
  min.features = 200
)

atacAggr <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# Add the patient information and disease status to the metadata of the Seurat object
gemgroup <- sapply(strsplit(rownames(atacAggr@meta.data), split="-"), "[[", 2) 
current.gemgroups <- seq(length(rownames(aggcsv))) # no. gemgroups is no. samples
orig.ident <- rownames(aggcsv)
sampleID <- plyr::mapvalues(gemgroup, from = current.gemgroups, to = orig.ident)
atacAggr <- AddMetaData(object=atacAggr, metadata=data.frame(orig.ident=sampleID, 
  row.names=rownames(atacAggr@meta.data)))

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(atacAggr) <- annotations

# compute nucleosome signal score per cell
atacAggr <- NucleosomeSignal(object = atacAggr)

# compute TSS enrichment score per cell
atacAggr <- TSSEnrichment(object = atacAggr, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
atacAggr$pct_reads_in_peaks <- atacAggr$peak_region_fragments / atacAggr$passed_filters * 100
atacAggr$blacklist_ratio <- atacAggr$blacklist_region_fragments / atacAggr$peak_region_fragments

# filter the aggregated snATACseq object using empirically-determined QC parameters
atacAggr <- subset(
  x = atacAggr,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 12000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.005 &
    nucleosome_signal < 3 &
    TSS.enrichment > 2
)

atacAggr <- RunTFIDF(atacAggr)
atacAggr <- FindTopFeatures(atacAggr, min.cutoff = 'q0')
atacAggr <- RunSVD(
  object = atacAggr
)

#Batch correction among samples
atacAggr <- RunHarmony(
  object = atacAggr,
  group.by.vars = 'orig.ident',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)

#clustering and visualization
atacAggr <- RunUMAP(atacAggr, dims = 2:20, reduction = 'harmony')
atacAggr <- FindNeighbors(object = atacAggr, reduction = "harmony", dims = 2:20)
atacAggr <- FindClusters(object = atacAggr, verbose = FALSE, algorithm = 3,resolution = 0.3)

#Gene activities
gene.activities <- GeneActivity(atacAggr)

atacAggr[['RNA']] <- CreateAssayObject(counts = gene.activities)
atacAggr <- NormalizeData(
  object = atacAggr,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atacAggr$nCount_RNA)
)

#chromvar score
DefaultAssay(atacAggr) <- 'peaks'
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = T)
)


atacAggr <- AddMotifs(
  object = atacAggr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

#Chromvar score
atacAggr <- RunChromVAR(
  object = atacAggr,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
