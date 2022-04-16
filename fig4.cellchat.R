library(CellChat)
library(patchwork)
library(Seurat)

rnaAggr <- readRDS("PKDContAggr.rds")

new.cluster.ids <- c("PT/FR-PTC","PT/FR-PTC","PEC","TAL2","TAL1","DCT","CNT_PC","ICA","ICB","PODO","ENDO","FIB","LEUK","URO1","URO2")
names(new.cluster.ids) <- levels(rnaAggr)
rnaAggr <- RenameIdents(rnaAggr, new.cluster.ids)
rnaAggr@meta.data[["celltype_all"]] <- rnaAggr@active.ident

rnaAggr <- subset(rnaAggr,idents = c("URO1","URO2"),invert=T)

Idents(rnaAggr) <- "disease"
pkd <- subset(rnaAggr,idents = "PKD")
Idents(pkd) <- "celltype_all"
levels(pkd) <- c("PT/FR-PTC","PEC","TAL1","TAL2","DCT","CNT_PC","ICA","ICB","PODO","ENDO","FIB","LEUK")

data.input <- GetAssayData(pkd, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(pkd)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

CellChatDB.use <- subsetDB(CellChatDB, search = c("IL6","TNF","TGFb","HH"),key = "pathway_name")# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat,population.size=T,type = "truncatedMean",trim = 0.001)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

heat_pkd_TNF <- netVisual_heatmap(cellchat_pkd, signaling = "TNF", color.heatmap = "Reds") #500x400
circ_pkd_TNF <- netVisual_aggregate(cellchat_pkd, signaling = "TNF", layout = "circle") #640x570
trib_pkd_TNF <- netAnalysis_contribution(cellchat_pkd, signaling = "TNF") #640x570

heat_pkd_IL6 <- netVisual_heatmap(cellchat_pkd, signaling = "IL6", color.heatmap = "Reds") 
circ_pkd_IL6 <- netVisual_aggregate(cellchat_pkd, signaling = "IL6", layout = "circle")
trib_pkd_IL6 <- netAnalysis_contribution(cellchat_pkd, signaling = "IL6")

heat_pkd_TGFb <- netVisual_heatmap(cellchat_pkd, signaling = "TGFb", color.heatmap = "Reds") 
circ_pkd_TGFb <- netVisual_aggregate(cellchat_pkd, signaling = "TGFb", layout = "circle")
trib_pkd_TGFb <- netAnalysis_contribution(cellchat_pkd, signaling = "TGFb")