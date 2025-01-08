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

varGenes <- as.integer(Args[6])
PCAdim <- as.integer(Args[7])
#Resolution:
ResN <- Args[8]
Knum <- Args[9]
Snum <- Args[10]

PATH <- '/home/nifang/mmyy/project/LifeSpan/seurat/sample'

outDir <- paste0(PATH,'/S',Args[10],'varGenes',Args[6],'_PCA',Args[7],'_Res',Args[8],'_K',Args[9])
if (!dir.exists(outDir)){
  dir.create(outDir)}
FilePrefix=paste0(outDir,'/Integ')

########################################
#1.read object (16 samples:s1 - s16)
########################################
sample_list = c(basename(list.dirs("/home/nifang/mmyy/project/LifeSpan/cellranger/out_persample/",recursive = F)))
Object_list = list()

for (sample in sample_list){
  filedir = str_c("/home/nifang/mmyy/project/LifeSpan/cellranger/out_persample/",sample)
  scrna_data <- Read10X(filedir)
  Seurat_object <- CreateSeuratObject(
    counts = scrna_data$`Gene Expression`,
    assay = 'RNA',min.cells = 3)

  Seurat_object[['CMO']] = CreateAssayObject(counts = scrna_data$`Multiplexing Capture`)
  Seurat_object[['percent.mt']] <- PercentageFeatureSet(object=Seurat_object, pattern='^MT-')
  Seurat_object[["sample"]] = sample

  Seurat_object <- subset(Seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

  Object_list[[sample]] = Seurat_object
}

#############################################################
#3.normalize-SCTransform for each dataset independently
#############################################################
#IntegObj.list <- lapply(X = Object_list, FUN = function(x) {
#    x <- NormalizeData(x)
#    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = varGenes)
#})
#features <- SelectIntegrationFeatures(object.list = IntegObj.list)
#IntegObj.list <- lapply(X = IntiegObj.list, FUN = function(x) {
#    x <- ScaleData(x, features = features, verbose = FALSE)
#    x <- RunPCA(x, features = features, verbose = FALSE)
#})
IntegObj.list <- lapply(X = Object_list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = IntegObj.list, nfeatures = varGenes)
print(length(features))
IMMUNOGLOBULIN_Gene=c('IGHA1',	'IGHA2',	'IGHD',	'IGHD1-1',	'IGHE',	'IGHG1',	'IGHG2',	'IGHG3',	'IGHG4',	'IGHJ1',	'IGHM',	'IGHV1-18',	'IGHV1-24',	'IGHV1-3',	'IGHV1-45',	'IGHV1-58',	'IGHV1-69',	'IGHV1-69-2',	'IGHV1-69D',	'IGHV1OR15-1',	'IGHV2-26',	'IGHV2-5',	'IGHV2-70',	'IGHV2-70D',	'IGHV3-11',	'IGHV3-13',	'IGHV3-15',	'IGHV3-16',	'IGHV3-20',	'IGHV3-21',	'IGHV3-23',	'IGHV3-30',	'IGHV3-33',	'IGHV3-35',	'IGHV3-38',	'IGHV3-43',	'IGHV3-48',	'IGHV3-49',	'IGHV3-53',	'IGHV3-64',	'IGHV3-64D',	'IGHV3-66',	'IGHV3-7',	'IGHV3-72',	'IGHV3-73',	'IGHV3-74',	'IGHV4-28',	'IGHV4-31',	'IGHV4-34',	'IGHV4-39',	'IGHV4-4',	'IGHV4-59',	'IGHV4-61',	'IGHV5-10-1',	'IGHV5-51',	'IGHV6-1',	'IGHV7-4-1',	'IGHV7-81',	'IGHV8-51-1',	'IGKC',	'IGKJ1',	'IGKV1-12',	'IGKV1-13',	'IGKV1-16',	'IGKV1-17',	'IGKV1-27',	'IGKV1-37',	'IGKV1-39',	'IGKV1-5',	'IGKV1-6',	'IGKV1-8',	'IGKV1-9',	'IGKV1D-12',	'IGKV1D-13',	'IGKV1D-17',	'IGKV1D-33',	'IGKV1D-37',	'IGKV1D-39',	'IGKV1D-42',	'IGKV1D-43',	'IGKV1D-8',	'IGKV2-24',	'IGKV2-28',	'IGKV2-29',	'IGKV2-30',	'IGKV2-40',	'IGKV2D-24',	'IGKV2D-26',	'IGKV2D-28',	'IGKV2D-29',	'IGKV2D-30',	'IGKV3-15',	'IGKV3-20',	'IGKV3-7',	'IGKV3D-11',	'IGKV3D-15',	'IGKV3D-20',	'IGKV3D-7',	'IGKV4-1',	'IGKV5-2',	'IGKV6-21',	'IGKV6D-21',	'IGKV6D-41',	'IGLC1',	'IGLC2',	'IGLC3',	'IGLC6',	'IGLC7',	'IGLJ1',	'IGLL1',	'IGLL5',	'IGLV1-36',	'IGLV1-40',	'IGLV1-44',	'IGLV1-47',	'IGLV1-50',	'IGLV1-51',	'IGLV10-54',	'IGLV11-55',	'IGLV2-11',	'IGLV2-14',	'IGLV2-18',	'IGLV2-23',	'IGLV2-33',	'IGLV2-8',	'IGLV3-1',	'IGLV3-10',	'IGLV3-12',	'IGLV3-16',	'IGLV3-19',	'IGLV3-21',	'IGLV3-22',	'IGLV3-25',	'IGLV3-27',	'IGLV3-32',	'IGLV3-9',	'IGLV4-3',	'IGLV4-60',	'IGLV4-69',	'IGLV5-37',	'IGLV5-45',	'IGLV5-48',	'IGLV5-52',	'IGLV6-57',	'IGLV7-43',	'IGLV7-46',	'IGLV8-61',	'IGLV9-49',	'JCHAIN',	'LIME1',	'PIGR',	'SYK',	'TRBC1',	'TRBC2',	'TRDC')
print('Remove IMMUNOGLOBULIN_Gene!')
for (g in features){
    if (g %in% IMMUNOGLOBULIN_Gene){
        print(g)
        features=features[features!=g]
    }
}
print(length(features))

IntegObj.list <- PrepSCTIntegration(object.list = IntegObj.list, anchor.features = features)
IntegObj.list <- lapply(X = IntegObj.list, FUN = RunPCA, features = features)

##############################
#4.Perform integration: RPCA
##############################
immune.anchors <- FindIntegrationAnchors(object.list = IntegObj.list, anchor.features = features,
  normalization.method = "SCT", dims = 1:PCAdim,reference = as.numeric(Snum), reduction = "rpca", k.anchor = as.numeric(Knum))
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:PCAdim)
pdf(file=paste0(outDir,'/QC.pdf'))
VlnPlot(immune.combined.sct,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
dev.off()
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
pdf(file=paste0(outDir,'/ElbowPlot.PCA.pdf'))
ElbowPlot(object = immune.combined.sct, ndims = 50)
dev.off()
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:PCAdim)
immune.combined.sct  <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:PCAdim)
immune.combined.sct  <- FindClusters(immune.combined.sct, resolution = as.numeric(ResN))

label_names <- read.csv('/home/nifang/mmyy/project/LifeSpan/seurat/name.txt', header=T,sep='\t')
labers = label_names[match(as.numeric(as.character(immune.combined.sct$sample)), label_names[, 1]), 2]
immune.combined.sct$BioGroup <- labers
B1 <- subset(immune.combined.sct, subset=BioGroup=='B1')
B2 <- subset(immune.combined.sct, subset=BioGroup=='B2')
B3 <- subset(immune.combined.sct, subset=BioGroup=='B3')
B4 <- subset(immune.combined.sct, subset=BioGroup=='B4')
B5 <- subset(immune.combined.sct, subset=BioGroup=='B5')
B6 <- subset(immune.combined.sct, subset=BioGroup=='B6')
B7 <- subset(immune.combined.sct, subset=BioGroup=='B7')
B8 <- subset(immune.combined.sct, subset=BioGroup=='B8')

##############################
#5.Visualization
##############################
pdf(file=paste0(outDir,'/Dimplot_','sample1.pdf'))
DimPlot(immune.combined.sct, reduction = "umap", group.by = "sample")
dev.off()
pdf(file=paste0(outDir,'/Dimplot_','cluster.pdf'))
DimPlot(immune.combined.sct, reduction = "umap", label = TRUE)
dev.off()


DefaultAssay(immune.combined.sct) <- "SCT"
pdf(file=paste0(outDir,'/FeaturePlot_gene1.pdf'))
FeaturePlot(immune.combined.sct, features = c("CD3D", "CD8A",  "IL7R", "CCR7", "S100A4", "KIT"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()

pdf(file=paste0(outDir,'/FeaturePlot_gene2.pdf'))
FeaturePlot(immune.combined.sct, features = c("MS4A1", "CD14", "LYZ", "FCGR3A", "MS4A7", "GNLY"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()

pdf(file=paste0(outDir,'/FeaturePlot_gene3.pdf'))
FeaturePlot(immune.combined.sct, features = c("IL2RB", "CST3", "FCER1A", "PPBP", "KLRC2", "XCL1"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()

pdf(file=paste0(outDir,'/FeaturePlot_gene4.pdf'))
FeaturePlot(immune.combined.sct, features = c("ITM2C", "HBB", "SPINK2", "MKI67", "TRDC", "KLRB1"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()

pdf(file=paste0(outDir,'/FeaturePlot_gene5.pdf'))
FeaturePlot(immune.combined.sct, features = c("PITX1", "COL1A2", "ANK2", "LTBP1", "THBS1", "STAR"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()
pdf(file=paste0(outDir,'/FeaturePlot_gene6.pdf'))
FeaturePlot(immune.combined.sct, features = c("CCL21", "TFPI", "CD74", "SRGN", "LDB2", "CLCN5"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()
pdf(file=paste0(outDir,'/FeaturePlot_gene7.pdf'))
FeaturePlot(immune.combined.sct, features = c("NPAS3", "CHGA", "PTPRC", "RAG1", "IGHA2", "FOXP3"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()
pdf(file=paste0(outDir,'/FeaturePlot_gene8.pdf'))
FeaturePlot(immune.combined.sct, features = c("PPP1R14A", "CLEC9A", "ITM2C", "RAKLRF1G1", "XCL2", "FOXP3"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()
pdf(file=paste0(outDir,'/FeaturePlot_NK1.pdf'))
FeaturePlot(immune.combined.sct, features = c("FCER1G", "MYOM2", "SPON2", "CLIC3", "IL2RB", "GZMB"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()
pdf(file=paste0(outDir,'/FeaturePlot_NK2.pdf'))
FeaturePlot(immune.combined.sct, features = c("KLRC2", "GZMH", "TYROBP", "FCGR3A", "GZMB", "FGFBP2"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()

pdf(file=paste0(outDir,'/FeaturePlot_Myeloid.pdf'))
FeaturePlot(immune.combined.sct, features = c("CD14", "FCRG3A", "CDKN1C", "CD68", "FCN1", "FCER1A","IL1B"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()


pdf(file=paste0(outDir,'/FeaturePlot_Mcro.pdf'))
FeaturePlot(immune.combined.sct, features = c("CD68","MSR1","MRC1","CD14","HLA-DPB1","UGCG"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()

pdf(file=paste0(outDir,'/FeaturePlot_Bnaivemem1.pdf'))
FeaturePlot(immune.combined.sct, features = c("COCH", "AIM2", "SSPN", "IGHM", "IGHD", "CXCR4"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()

pdf(file=paste0(outDir,'/FeaturePlot_Bnaivemem2.pdf'))
FeaturePlot(immune.combined.sct, features = c("IGKC", "IGLC2", "IGLC3", "CCL5", "CD27", "KLRF1"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()

pdf(file=paste0(outDir,'/FeaturePlot_TCD8naive.pdf'))
FeaturePlot(immune.combined.sct, features = c("LEF1", "CD3D", "CCR7", "CD8A", "CD8B", "CCL5"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()

pdf(file=paste0(outDir,'/FeaturePlot_TCD8mem.pdf'))
FeaturePlot(immune.combined.sct, features = c("CCL5", "CD8A", "GZMA", "GZMK", "CD8B", "LTB"), max.cutoff = 3,
    cols = c("grey", "red"))
dev.off()


i<-"COCH"
pdf(file=paste0(outDir,"/",i,".pdf"),width=5, height=2.5)
VlnPlot(immune.combined.sct, features = c(i), pt.size = 0)+NoLegend()+theme(text = element_text(size = 8))
dev.off()
i<-"IGHM"
pdf(file=paste0(outDir,"/",i,".pdf"),width=5, height=2.5)
VlnPlot(immune.combined.sct, features = c(i), pt.size = 0)+NoLegend()+theme(text = element_text(size = 8))
dev.off()

i<-"CCR7"
pdf(file=paste0(outDir,"/",i,".pdf"),width=5, height=2.5)
VlnPlot(immune.combined.sct, features = c(i), pt.size = 0)+NoLegend()+theme(text = element_text(size = 8))
dev.off()
i<-"LEF1"
pdf(file=paste0(outDir,"/",i,".pdf"),width=5, height=2.5)
VlnPlot(immune.combined.sct, features = c(i), pt.size = 0)+NoLegend()+theme(text = element_text(size = 8))
dev.off()
i<-"CCL5"
pdf(file=paste0(outDir,"/",i,".pdf"),width=5, height=2.5)
VlnPlot(immune.combined.sct, features = c(i), pt.size = 0)+NoLegend()+theme(text = element_text(size = 8))
dev.off()

i<-"GZMA"
pdf(file=paste0(outDir,"/",i,".pdf"),width=5, height=2.5)
VlnPlot(immune.combined.sct, features = c(i), pt.size = 0)+NoLegend()+theme(text = element_text(size = 8))
dev.off()
i<-"GZMK"
pdf(file=paste0(outDir,"/",i,".pdf"),width=5, height=2.5)
VlnPlot(immune.combined.sct, features = c(i), pt.size = 0)+NoLegend()+theme(text = element_text(size = 8))
dev.off()
i<-"KLRC2"
pdf(file=paste0(outDir,"/",i,".pdf"),width=5, height=2.5)
VlnPlot(immune.combined.sct, features = c(i), pt.size = 0)+NoLegend()+theme(text = element_text(size = 8))
dev.off()
i<-"CD8A"
pdf(file=paste0(outDir,"/",i,".pdf"),width=5, height=2.5)
VlnPlot(immune.combined.sct, features = c(i), pt.size = 0)+NoLegend()+theme(text = element_text(size = 8))
dev.off()
i<-"CD8B"
pdf(file=paste0(outDir,"/",i,".pdf"),width=5, height=2.5)
VlnPlot(immune.combined.sct, features = c(i), pt.size = 0)+NoLegend()+theme(text = element_text(size = 8))
dev.off()
i<-"CD3D"
pdf(file=paste0(outDir,"/",i,".pdf"),width=5, height=2.5)
VlnPlot(immune.combined.sct, features = c(i), pt.size = 0)+NoLegend()+theme(text = element_text(size = 8))
dev.off()
i<-"CD3E"
pdf(file=paste0(outDir,"/",i,".pdf"),width=5, height=2.5)
VlnPlot(immune.combined.sct, features = c(i), pt.size = 0)+NoLegend()+theme(text = element_text(size = 8))
dev.off()


#MetaData=as.data.frame(immune.combined.sct@meta.data)
#write.table(MetaData,file=paste0(outDir,'/ImmnueInteg.MetaData.txt'),sep='\t',quote=F)
#TSNE=as.data.frame(Embeddings(object = immune.combined.sct, reduction = "umap"))
#write.table(TSNE,file=paste0(outDir,'/ImmnueInteg.umap.txt'),sep='\t',quote=F)

#saveRDS(immune.combined.sct,paste0(outDir,'/IntegImmune.rds'))

#DefaultAssay(immune.combined.sct) <- "SCT"
#immune.combined.sct <- PrepSCTFindMarkers(immune.combined.sct)

#immune.markers <- FindAllMarkers(immune.combined.sct, assay = "SCT", test.use = "wilcox",only.pos=TRUE,min.pct=0.1,return.thresh=0.01,logfc.threshold=0.2)
#MarkerGene=as.data.frame(immune.markers %>% group_by(cluster) %>% top_n(50,avg_log2FC))
#head(MarkerGene)
#write.table(MarkerGene,file=paste0(outDir,'/IntegImmune_Top50MarkerGene_wilcox_SCT.txt'),quote =F,sep='\t')


#DefaultAssay(immune.combined.sct) <- "RNA"
#immune.markers <- FindAllMarkers(immune.combined.sct, test.use = "wilcox",only.pos=TRUE,min.pct=0.1,return.thresh=0.01,logfc.threshold=0.2)
#MarkerGene=as.data.frame(immune.markers %>% group_by(cluster) %>% top_n(50,avg_log2FC))
#head(MarkerGene)
#write.table(MarkerGene,file=paste0(outDir,'/IntegImmune_Top50MarkerGene_wilcox_RNA.txt'),quote =F,sep='\t')

#immnue.combine.sct
#Save RDS
#saveRDS(immune.combined.sct,paste0(outDir,'/IntegImmune.rds'))
