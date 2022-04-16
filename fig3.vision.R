library(Seurat)
library(VISION)
library(BuenColors)
set.seed(1234)

rnaAggr <- readRDS("PKDContAggr.rds")
signatures <- "h.all.v6.2.symbols.gmt"
vision.obj <- Vision(rnaAggr, signatures = signatures,pool =F)
options(mc.cores=10)
vision.obj <- analyze(vision.obj)

sigScores <- getSignatureScores(vision.obj)
sigScores <- as.data.frame(sigScores)

rnaAggr[['GSEA']] <- CreateAssayObject(counts = t(sigScores))
DefaultAssay(rnaAggr) <- "GSEA"
Idents(rnaAggr) <- "celltype_all"
new.cluster.ids <- c("PT/FR-PTC","PT/FR-PTC","PEC","ATL","TAL","DCT","CNT_PC","ICA","ICB","PODO","ENDO","FIB","LEUK","URO1","URO2")
names(new.cluster.ids) <- levels(rnaAggr)
rnaAggr <- RenameIdents(rnaAggr, new.cluster.ids)
rnaAggr@meta.data[["celltype_all2"]] <- rnaAggr@active.ident
rnaAggr <- subset(rnaAggr,idents = c("URO1","URO2"),invert = T)
rnaAggr@meta.data$celltype_disease <- paste0(rnaAggr@meta.data$celltype_all2,"_",rnaAggr@meta.data$disease)
Idents(rnaAggr) <- "celltype_disease"

levels(rnaAggr) <- c("PT/FR-PTC_control","PEC_control","ATL_control","TAL_control",
                     "DCT_control", "CNT_PC_control","ICA_control","ICB_control","PODO_control","ENDO_control",
                     "FIB_control","LEUK_control","PT/FR-PTC_PKD","PEC_PKD","ATL_PKD","TAL_PKD",
                     "DCT_PKD", "CNT_PC_PKD","ICA_PKD","ICB_PKD","PODO_PKD","ENDO_PKD",
                     "FIB_PKD","LEUK_PKD")

test <- AverageExpression(rnaAggr,assays = "GSEA")
fig3a <- pheatmap::pheatmap(test[["GSEA"]],scale = "row",cluster_cols = F)

DefaultAssay(rnaAggr) <- "GSEA"
fig3b_1 <- FeaturePlot(rnaAggr,"HALLMARK-IL6-JAK-STAT3-SIGNALING",cols =c("#E4F0F0","#F33D0C"),max.cutoff = "q99",min.cutoff = "q1",split.by = "disease") #920x450
fig3b_2 <- FeaturePlot(rnaAggr,"HALLMARK-TNFA-SIGNALING-VIA-NFKB",cols =c("#E4F0F0","#F33D0C"),max.cutoff = "q99",min.cutoff = "q1",split.by = "disease")
fig3b_3 <- FeaturePlot(rnaAggr,"HALLMARK-TGF-BETA-SIGNALING",cols =c("#E4F0F0","#F33D0C"),max.cutoff = "q99",min.cutoff = "q1",split.by = "disease")
