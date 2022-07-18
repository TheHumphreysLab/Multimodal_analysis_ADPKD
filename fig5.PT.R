library(Seurat)
library(ggplot2)
library(harmony)
library(Rcpp)
library(dplyr)
library(BuenColors)
set.seed(1234)

rnaAggr <- readRDS("PKDContAggr.rds") #Whole dataset

#The target cell type subset using the annotations on an integrated dataset (Fig. 1b) was then further subset with the annotations on each dataset to extract cell type with high confidence. 
#(Supplementary Fig. 1d, Supplementary Fig. 2d) 

Idents(rnaAggr) <- "celltype"
PTC <- subset(rnaAggr,idents = c("PKD-PT","Cont-PT1","Cont-PT2"))
Idents(PTC) <- "celltype_all"
PTC <- subset(PTC,idents = c("PT","FR-PTC"))

#Control PT subclustering
Idents(PTC) <- "disease"
PTC1 <- subset(PTC,idents = "control")
PTC1 <- NormalizeData(PTC1)
PTC1 <- FindVariableFeatures(PTC1, selection.method = "vst", nfeatures = 3000)
PTC1 <- ScaleData(PTC1, verbose = FALSE)
PTC1 <- RunPCA(PTC1, verbose = FALSE)
PTC1 <- RunHarmony(PTC1,group.by.vars = "orig.ident")
PTC1 <- RunUMAP(PTC1, reduction = "harmony", dims = 1:10)
PTC1 <- FindNeighbors(PTC1, reduction = "harmony", dims = 1:10)
PTC1 <- FindClusters(PTC1, resolution = 0.1)

new.cluster.ids <- c("cont_N-PTC","cont_N-PTC","cont_FR-PTC","cont_N-PTC")
names(new.cluster.ids) <- levels(PTC1)
PTC1 <- RenameIdents(PTC1, new.cluster.ids)
levels(PTC1) <- c("cont_N-PTC","cont_FR-PTC")
PTC1@meta.data[["subtype"]] <- PTC1@active.ident
fig5a1 <- DimPlot(PTC1,cols = c("#00BFC4","#F8766D")) #600x450
fig5a2 <- FeaturePlot(PTC1,"VCAM1") #480x450

#ADPKD PT
PTC2 <- subset(PTC,idents = "PKD")
PTC2@meta.data[["subtype"]] <- "ADPKD"

#Clustering on control FR-PTC and ADPKD PT cells (without control N-PTC)
PTC1subtype <- PTC1@meta.data[,c(1,18)]
PTC2subtype <- PTC2@meta.data[,c(1,17)]
PTClist <- rbind(PTC1subtype,PTC2subtype )
PTC <- AddMetaData(PTC,PTClist)
Idents(PTC) <- "subtype"
PTC3 <- subset(PTC,idents = c("ADPKD","cont_FR-PTC"))
PTC3 <- NormalizeData(PTC3)
PTC3 <- FindVariableFeatures(PTC3, selection.method = "vst", nfeatures = 3000)
PTC3 <- ScaleData(PTC3, verbose = FALSE)
PTC3 <- RunPCA(PTC3, verbose = FALSE)
PTC3 <- RunHarmony(PTC3,group.by.vars = "orig.ident")
PTC3 <- RunUMAP(PTC3, reduction = "harmony", dims = 1:10)
PTC3 <- FindNeighbors(PTC3, reduction = "harmony", dims = 1:10)
PTC3 <- FindClusters(PTC3, resolution = 0.2)
new.cluster.ids <- c("PT3","PT1","PT4","PT2")
names(new.cluster.ids) <- levels(PTC3)
PTC3 <- RenameIdents(PTC3, new.cluster.ids)
levels(PTC3) <- c("PT1","PT2","PT3","PT4")
PTC3@meta.data[["subtype2"]] <- PTC3@active.ident

fig5b1 <- DimPlot(PTC3) #550x450
fig5b2 <- DimPlot(PTC3,group.by = "subtype",cols = c("#F8766D","#0090c4"),pt.size = 0.3) #620x450

#Gene expressions in ADPKD subtypes
Idents(PTC3) <- "disease"
PTC4 <- subset(PTC3,idents = "PKD")
Idents(PTC4) <- "subtype2"

features <- c("LRP2","CUBN","VCAM1","DCC","CCL2","CD24","GPX3","FOS","MT2A",
              "CREB5","TPM1","PRICKLE1","PROM1","TGFB2","HAVCR1","CD44")

levels(PTC4) <- rev(levels(PTC4))
fig5c <- DotPlot(PTC4, features = features, cols = c("lightyellow","royalblue")) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #700x450
levels(PTC4) <- rev(levels(PTC4))

fig5e1 <- VlnPlot(PTC4,"CUBN")+NoLegend()
fig5e2 <- VlnPlot(PTC4,"LRP2")+NoLegend()
fig5e3 <- VlnPlot(PTC4,"VCAM1")+NoLegend()

#Transfer the PT subtype identity (N-PTC/FR-PTC in control and PT1/2/3/4 in ADPKD) to whole PT dataset
PTC1subtype <- PTC1@meta.data[,c(1,18)] #N-PTC/FR-PTC in control
PTC4subtype <- PTC4@meta.data[,c(1,19)] #PT1/2/3/4 in ADPKD
colnames(PTC4subtype)[2] <- "subtype" 
PTClist2 <- rbind(PTC1subtype,PTC4subtype )
PTC <- AddMetaData(PTC,PTClist2)
Idents(PTC) <- "subtype"
DefaultAssay(PTC) <- "RNA"
PTC_marker <- FindAllMarkers(PTC,only.pos = T)
write.csv(PTC_marker, "Supplementary_Data13.csv")

############## Correlation of gene expressions with IRI mice dataset ################

library(biomaRt)
library(Matrix)

#IRI.all = IRI mice data

#lift over (mouse to human)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mouse_gene <- rownames(IRI.all)
converted_genelist = getLDS(attributes = "mgi_symbol", filters = "mgi_symbol", values = mouse_gene , 
                            mart = mouse, attributesL = "hgnc_symbol", martL = human, uniqueRows=T)
converted_genelist <- converted_genelist %>% distinct(HGNC.symbol,.keep_all=TRUE)

#Mouse data
mousePT <- subset(IRI.all,idents = c("PTS1","PTS2","PTS3","NewPT1","NewPT2"))
Idents(mousePT) <- "name"
new.cluster.ids <- c("PCT","PCT","PST","Injured","Failed-repair") #PCT=S1/S2, PST=S3
names(new.cluster.ids) <- levels(mousePT)
mousePT <- RenameIdents(mousePT, new.cluster.ids)
levels(mousePT) <- c("PCT","PST","Injured","Failed-repair")
mousePT@meta.data[["subtype"]] <- mousePT@active.ident
mousePT <- FindVariableFeatures(mousePT, selection.method = "vst", nfeatures = 2000) #Highly variable genes that define PT subtypes in IRI mice dataset

#IRI mice PT average gene expression
Avg_mousePT <- AverageExpression(mousePT)
Avg_mousePT <- Avg_mousePT[["RNA"]]
Avg_mousePT <- as.data.frame(Avg_mousePT)
Avg_mousePT$MGI.symbol <- rownames(Avg_mousePT)
Avg_mousePT <-  merge(Avg_mousePT,converted_genelist,by = "MGI.symbol") #27133
Avg_mousePT <- Avg_mousePT[Avg_mousePT$MGI.symbol %in% mousePT@assays[["RNA"]]@var.features, ] #select only variable genes

#Human PT average gene expression
Avg_humanPT <- AverageExpression(PTC)
Avg_humanPT <- Avg_humanPT[["RNA"]]
Avg_humanPT <- as.data.frame(Avg_humanPT)
Avg_humanPT$HGNC.symbol <- rownames(Avg_humanPT)

merged <- merge(Avg_mousePT,Avg_humanPT,by = "HGNC.symbol") #1648 genes left
rownames(merged) <- merged$HGNC.symbol
merged <- merged[,3:12]
cor_merge <- cor(merged) #Pearson correlation
fig5d <- pheatmap::pheatmap(cor_merge[1:4,5:10],cluster_rows = F,cluster_cols = F,border_color = "black") #700x380

######################################## Gene set enrichment analysis(VISION) ################################################
library(VISION)
signatures <- "h.all.v6.2.symbols.gmt" #hallmark gene set (for GSEA)
 
#Run Vision on the whole data (control + ADPKD)
vision.obj <- Vision(PTC, signatures = signatures)
vision.obj <- analyze(vision.obj)

sigScores <- getSignatureScores(vision.obj)
sigScores <- as.data.frame(sigScores)

PTC[['GSEA']] <- CreateAssayObject(counts = t(sigScores))
DefaultAssay(PTC) <- "GSEA"
gsea <- AverageExpression(PTC,assays = "GSEA")
gsea <- as.data.frame((gsea[["GSEA"]]))
# target pathway for analysis
target <- c("HALLMARK-INFLAMMATORY-RESPONSE","HALLMARK-TNFA-SIGNALING-VIA-NFKB","HALLMARK-IL6-JAK-STAT3-SIGNALING",
           "HALLMARK-TGF-BETA-SIGNALING","HALLMARK-OXIDATIVE-PHOSPHORYLATION")
gsea <- gsea[rownames(gsea) %in% target, ] 

#heatmap
fig5g <- pheatmap::pheatmap(gsea,scale = "row", cluster_cols = F,cluster_rows = T) #700x350

######################################## Pathway analysis(PROGENy) ################################################
library(progeny)
library(tidyr)
library(tibble)
PTC <- progeny(PTC, scale=FALSE, organism="Human", top=500, return_assay = TRUE)
PTC <- Seurat::ScaleData(PTC, assay = "progeny")

## Transofm Progeny scores into a data frame to better handling the results
progeny_scores_df <- as.data.frame(t(GetAssayData(PTC, slot = "scale.data", assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

## Progeny scores with the PT subtype
CellsClusters <- PTClist2
CellsClusters$Cell <- rownames(PTClist2)
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

##  Summarize the Progeny scores by subtype
summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, subtype) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## Prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%  
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0,
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength,
                      max(summarized_progeny_scores_df),
                      length.out=floor(paletteLength/2)))
## Plot
fig5h = pheatmap::pheatmap(t(summarized_progeny_scores_df[,]),fontsize=14,
                        fontsize_row = 10,
                        color=myColor, breaks = progenyBreaks,angle_col = 45,
                        treeheight_col = 0,  border_color = NA,cluster_cols = F) #630x500
