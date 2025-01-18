
library(Signac)
library(Seurat)
library(rtracklayer)
library(GenomicRanges)
library(future)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v86)
Args <- commandArgs()
print(Args)


plan("multicore", workers = 2)
options(future.globals.maxSize = 50000 * 1024^2)

agegroup <- Args[6]
outDir <- paste0('/scATAC/MergeObject/',Args[6])

#------------------------创建对象：create seurat object------------------------#
counts <- Read10X_h5(filename = paste0('cellranger_atac/',Args[6],'/outs/filtered_peak_bc_matrix.h5'))
metadata <- read.csv(file = paste0('cellranger_atac/',Args[6],'/outs/singlecell.csv'),header = TRUE,row.names = 1)


chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = paste0('cellranger_atac/',Args[6],'/outs/fragments.tsv.gz'),
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)
pbmc
#pbmc[['peaks']]

#------------------------添加hg38基因组注释------------------------#
# extract gene annotations from EnsDb
hg38 <- import("gencode.v32.annotation.gtf.gz")
genome(hg38) <- "hg38"
seqlevelsStyle(hg38) <- "UCSC"
hg38$gene_biotype <- hg38$gene_type
Annotation(pbmc) <- hg38

#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg19
#seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
#seqlevelsStyle(annotations) <- 'UCSC'
#genome(annotations) <- "hg38"
# add the gene information to the object
#Annotation(pbmc) <- annotations

#------------------------质控和筛选-QC------------------------#
# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)
# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)
# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
pbmc@meta.data$TSS.enrichment
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 3, 'High', 'Low')

#pdf(file=paste0(outDir,'/QC_TSSPlot.pdf'))
#TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
#dev.off()

#pdf(file=paste0(outDir,'/QC_FragmentHistogram.pdf'))
#FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
#dev.off()

#pdf(file=paste0(outDir,'/QC_VlnPlot.pdf'))
#VlnPlot(object = pbmc,
#  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
#  pt.size = 0.1,
#  ncol = 5
#)
#dev.off()

pbmc <- subset(x = pbmc,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 30000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)
pbmc

#------------------------标准化和降维------------------------#
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 20)
#pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

pdf(file=paste0(outDir,'/QC_DepthCor_LSI.pdf'))
DepthCor(pbmc)
dev.off()

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

pdf(file=paste0(outDir,'/DimPlot_cluster.pdf'))
DimPlot(object = pbmc, label = TRUE) + NoLegend()
dev.off()



#------------------------基因表达及活性矩阵------------------------#
gene.activities <- GeneActivity(pbmc)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

DefaultAssay(pbmc) <- 'RNA'


pdf(file=paste0(outDir,'/FeaturePlot_GeneExp_1.pdf'),width=10, height=6)
FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  cols = c("#DFDFDE", "#A10035")
)
dev.off()


pdf(file=paste0(outDir,'/FeaturePlot_GeneExp_2.pdf'),width=10, height=6)
FeaturePlot(object = pbmc,features = c('CD34', 'KLRF1', 'CD79A', 'S100A9', 'GZMK', 'IL7R'),pt.size = 0.1,max.cutoff = 'q95',ncol = 3,cols = c("#DFDFDE", "#A10035"))
dev.off()

pdf(file=paste0(outDir,'/FeaturePlot_GeneExp_3.pdf'),width=10, height=6)
FeaturePlot(object = pbmc,features = c('LGALS1', 'LTB', 'KLRB1', 'GZMB', 'CCR7', 'CD8A'),pt.size = 0.1,max.cutoff = 'q95',ncol = 3,cols = c("#DFDFDE", "#A10035"))
dev.off()

#------------------------整合 scRNA，按scRNA分群投影scATAC------------------------#
# Load the pre-processed scRNA-seq data for PBMCs
print('------------------test-----------------------')
pbmc_rna <- readRDS("/home/nifang/mmyy/project/LifeSpan/seurat/sample/S4varGenes3000_PCA25_Res0.6_K15/rds_and_meta/IntegImmune_remove_rename_PrepSCTFindMarkers.rds")
DefaultAssay(pbmc_rna) <- 'RNA'
pbmc_rna <- FindVariableFeatures(object = pbmc_rna, nfeatures = 2000)
pbmc_rna
DefaultAssay(pbmc) <- 'RNA'
#pbmc <- FindVariableFeatures(object = pbmc, nfeatures = 2000)
pbmc

#pbmc_rna <- UpdateSeuratObject(pbmc_rna)
transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca')

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype3,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype3',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

pdf(file=paste0(outDir,'/DimPlot_merge_scRNA_scATAC.pdf'))
plot1 + plot2
dev.off()
