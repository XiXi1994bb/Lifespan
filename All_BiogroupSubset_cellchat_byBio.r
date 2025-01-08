
#Args <- commandArgs()
#print(Args)

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(ggalluvial)
library(NMF)
#library(SeuratData)
options(stringsAsFactors = FALSE)

OutDir <- '/home/nifang/mmyy/project/LifeSpan/seurat/sample/S4varGenes3000_PCA25_Res0.6_K15/cellchat/BioGroupSubset/'

All <- readRDS('/home/nifang/mmyy/project/LifeSpan/seurat/sample/S4varGenes3000_PCA25_Res0.6_K15/rds_and_meta/IntegImmune_remove_rename_PrepSCTFindMarkers.rds')
DefaultAssay(object = All) <- "SCT"
All <- subset(All, subset=(predicted_doublets!='True'))
All <- subset(All, subset=(sample!=2))

DF <- read.table(paste0(OutDir,'Merge4celltype_subset.MetaData.removeDou.removeS2_ForCells.txt'),header=T,sep='\t')
TM <- subset(All, cells = DF$CellID)
#Idents(TM) <- DF[, "SubCluster"]
TM@meta.data$SubCluster <- DF$SubCluster
TM$BioSubset <- paste(TM$BioGroup, TM$SubCluster, sep = "_")
Idents(TM) <- 'BioSubset'
print('read RDS done!')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
g <- 'B3'
All_ci <- subset(TM, subset=BioGroup == g)
cellchat_ci <- createCellChat(object = All_ci, meta=All_ci@meta.data, group.by = "BioSubset")
CellChatDB <- CellChatDB.human
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
#使用ECM,Secreted,Cell-Cell所有库
CellChatDB.use <- CellChatDB
cellchat_ci@DB <- CellChatDB.use
cellchat_ci <- subsetData(cellchat_ci)
future::plan("multicore", workers = 4)
cellchat_ci <- identifyOverExpressedGenes(cellchat_ci)
cellchat_ci <- identifyOverExpressedInteractions(cellchat_ci)
# Calculate communication probability and network
cellchat_ci <- computeCommunProb(cellchat_ci)
cellchat_ci <- computeCommunProbPathway(cellchat_ci)
cellchat_ci <- aggregateNet(cellchat_ci)
cellchat_ci <- netAnalysis_computeCentrality(cellchat_ci, slot.name = "netP")
# Make interaction circle plot for BioGroup
#pdf(paste0(OutDir, "Cellchat_number_of_interactions_SubCluster",paste0(g),".pdf"))
#netVisual_circle(cellchat_ci@net$count, weight.scale = T, label.edge= F, title.name = "Cell-Cell Interactions")
#dev.off()
df.net <- subsetCommunication(cellchat_ci)
write.csv(df.net, paste0(OutDir, "Cellchat_net_BiogroupSubset",paste0(g),".txt"),sep='\t')
df.netp <- subsetCommunication(cellchat_ci, slot.name = 'netP')
write.csv(df.netp, paste0(OutDir, "Cellchat_netp_BiogroupSubset",paste0(g),".txt"),sep='\t')
saveRDS(cellchat_ci, file = paste0(OutDir, "Cellchat_BiogroupSubset",paste0(g),".rds"))
#
