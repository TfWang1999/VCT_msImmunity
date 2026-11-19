# Original CellChat inference for P21 and P12; distinct from drawing cached objects.

sc_compute_cellchat <- function(seurat_obj) {
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(CellChat))
  seurat_obj <- subset(seurat_obj, subset = !Type %in% c("DC", "CD4_Tn","CD8_Tn","Fibro","Endo"))
  seurat_obj <- subset(seurat_obj, subset = !SubCellType %in% c("c90_Epi_MKI67"))
  seurat_obj$Type[seurat_obj$Type %in% c("CD4_Th")] <- "CD4T"
  seurat_obj$Type[seurat_obj$Type %in% c("CD8_em", "CD8_ex")] <- "CD8T"
  table(seurat_obj$Type)
  table(seurat_obj@meta.data[["SubCellType"]])
  seurat_obj$Type[seurat_obj$SubCellType %in% c("c91_Epi_Tumor")] <- "tumor"
  High <- subset(seurat_obj, subset = Patient %in% c("P21"))
  Low <- subset(seurat_obj, subset = Patient %in% c("P12"))
  table(Low$Type)
  table(High$Type)
  cellchat_High <- createCellChat(object = High, meta = High@meta.data, group.by = "Type")
  CellChatDB <- CellChatDB.human
  CellChatDB.use = CellChatDB
  cellchat_High@DB <- CellChatDB.use
  cellchat_High <- subsetData(cellchat_High)
  cellchat_High <- identifyOverExpressedGenes(cellchat_High)
  cellchat_High <- subsetData(cellchat_High)
  cellchat_High <- identifyOverExpressedInteractions(cellchat_High)
  cellchat_High <- computeCommunProb(cellchat_High, raw.use = TRUE, population.size = TRUE)
  cellchat_High <- filterCommunication(cellchat_High, min.cells = 10)
  cellchat_High <- computeCommunProbPathway(cellchat_High)
  cellchat_High <- aggregateNet(cellchat_High)
  cellchat_Low <- createCellChat(object = Low, meta = Low@meta.data, group.by = "Type")
  CellChatDB <- CellChatDB.human
  CellChatDB.use = CellChatDB
  cellchat_Low@DB <- CellChatDB.use
  cellchat_Low <- subsetData(cellchat_Low)
  cellchat_Low <- identifyOverExpressedGenes(cellchat_Low)
  cellchat_Low <- subsetData(cellchat_Low)
  cellchat_Low <- identifyOverExpressedInteractions(cellchat_Low)
  cellchat_Low <- computeCommunProb(cellchat_Low, raw.use = TRUE, population.size = TRUE)
  cellchat_Low <- filterCommunication(cellchat_Low, min.cells = 10)
  cellchat_Low <- computeCommunProbPathway(cellchat_Low)
  cellchat_Low <- aggregateNet(cellchat_Low)
  list(high = cellchat_High, low = cellchat_Low)
}

