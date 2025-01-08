Args <- commandArgs()
print(Args)

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(clustree)
library(stringr)
library(glmGamPoi)
#library(metap)

Dir <- '/home/nifang/mmyy/project/LifeSpan/seurat/sample/S4varGenes3000_PCA25_Res0.6_K15'
outDir <- '/home/nifang/mmyy/project/LifeSpan/seurat/sample/S4varGenes3000_PCA25_Res0.6_K15/GenAge/'

########################################
#1.read object (16 samples:s1 - s16)
########################################
immune.combined.sct <- readRDS(paste0(Dir,'/IntegImmune_remove_rename.rds'))
DefaultAssay(immune.combined.sct) <- "SCT"
#ImmunueCells_obj <- subset(immune.combined.sct, subset = (celltype2!='Granulocyte')&(celltype3!='HSPC')&(celltype3!='Platelet'))
#Idents(ImmunueCells_obj) <- "BioGroup"
#ImmunueCells_obj_sort <- ImmunueCells_obj
#my_levels <- c('B1','B2','B3','B4','B5','B6','B7','B8')
#Idents(ImmunueCells_obj_sort) <- factor(Idents(ImmunueCells_obj_sort), levels= my_levels)

#----------------------------------------------------------------------------
# GenAge_signature
#----------------------------------------------------------------------------
DF <- read.table(paste0(outDir,'GenAge_signature.txt',header=T,sep='\t')

#for (g in colnames(DF)){
#    ImmunueCells_obj_sort <- AddModuleScore(ImmunueCells_obj_sort,features = as.list(DF[g]),name = g)
#    pdf(file=paste0(outDir,paste0(g,'1'),".BioGroup.pdf"),width=2.3, height=2)
#    g1<-VlnPlot(ImmunueCells_obj_sort, features = paste0(g,'1'),pt.size = 0,cols=c("#E71F17", "#F18A1A","#EDE93B", "#81BF25","#19873B", "#1FA1D6","#831182","#143E97"))+NoLegend()+theme(text = element_text(size = 8))
#    print(g1)
#    dev.off()
#}
#MetaData=as.data.frame(ImmunueCells_obj_sort@meta.data)
#write.table(MetaData,file=(paste0(outDir,'IntegImmune.remove.rename.MetabolismSignature.MetaData.txt')),quote =F,sep='\t')

immune.combined.sct <- AddModuleScore(immune.combined.sct,features = as.list(DF['symbol']),name = 'GeneAge')
Idents(immune.combined.sct) <- "BioGroup"
pdf(file=paste0(outDir,"GenAge_signature.BioGroup.pdf"),width=2.3, height=2)
VlnPlot(immune.combined.sct, features = 'GeneAge1',pt.size = 0,cols=c("#E71F17", "#F18A1A","#EDE93B", "#81BF25","#19873B", "#1FA1D6","#831182","#143E97"))+NoLegend()+theme(text = element_text(size = 8))
dev.off()

#----------------------------------------------------------------------------
# cell cycle analysis
#----------------------------------------------------------------------------
cc.genes  <- readLines(con = paste0(outDir,'regev_lab_cell_cycle_genes.txt"))
s.genes   <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
# Assign Cell-Cycle Scores
immune.combined.sct <- CellCycleScoring(object = immune.combined.sct, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

MetaData=as.data.frame(immune.combined.sct@meta.data)
write.table(MetaData,file=(paste0(outDir,'IntegImmune.remove.rename.GenAge.MetaData.txt')),quote =F,sep='\t')
