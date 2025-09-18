library(Seurat)
library(tidyverse)
library(here)

cluster <- read_csv("Result/highReadRemoved_noBC.cluster.csv")
noBC_umap <- read_csv("Result/highReadRemoved_noBC.umap.csv")
noBC_umap$ident <- cluster$RNA_snn_res.0.2

# load non batch corrected seurat
highReadRemoved <- readRDS("Result/highReadRemoved.rds")
Idents(highReadRemoved) <- "RNA_snn_res.0.2"

cts <- subset(highReadRemoved, subset = RNA_snn_res.0.2!=6)@assays$RNA@counts
cur_seurat <- CreateSeuratObject(counts = cts)
cur_seurat$log10GenesPerUMI <- log10(cur_seurat$nFeature_RNA) / log10(cur_seurat$nCount_RNA)
cur_seurat$percent.mt <- PercentageFeatureSet(cur_seurat,"^mt:")
cur_seurat$percent.rb <- PercentageFeatureSet(cur_seurat,"Rp[LS]")
allSample_normalized <- NormalizeData(cur_seurat)
allSample_normalized <- FindVariableFeatures(allSample_normalized, 
                                             selection.method = "vst", nfeatures = 2000)
allSample_normalized <- ScaleData(allSample_normalized)
allSample_normalized <- RunPCA(allSample_normalized)

min.pc <- 18 # same as before 
allSample_normalized <- RunUMAP(allSample_normalized, dims = 1:min.pc,
                                umap.method = "umap-learn")
allSample_normalized <- RunTSNE(allSample_normalized, 
                                dims = 1:min.pc, reduction = 'pca')
allSample_normalized <- FindNeighbors(allSample_normalized, dims = 1:min.pc)  
allSample_normalized <- FindClusters(allSample_normalized, 
                                     resolution = c(0.1,0.15),
                                     algorithm = 4)
allSample_normalized <- FindClusters(allSample_normalized, 
                                     resolution = seq(.2,.5,0.1),
                                     algorithm = 4)
saveRDS(allSample_normalized, "Result/noCluster6.rds")

rm(cts, cur_seurat)
gc()

allSample_normalized <- readRDS("Result/noCluster6.rds")

Idents(allSample_normalized) <- "RNA_snn_res.0.2"
p1 <- UMAPPlot(allSample_normalized, label=T)
p0 <- UMAPPlot(highReadRemoved, label=T)
p0|p1

# Save clustering results with cluster 6 removed
a <- grep("RNA_snn",names(allSample_normalized@meta.data))
tmp <- (allSample_normalized@meta.data[,a])
a <- grep('^yki',row.names(tmp))
tmp$barcode <- ""
tmp[a,]$barcode <- gsub('-1','-2',gsub('yki_','',row.names(tmp)[a]))
tmp[-a,]$barcode <-gsub('wt_','',row.names(tmp)[-a])
write.csv(tmp, 'Result/noCluster6.cluster.csv', quote = F, row.names = F)
noCluster6_cluster <- read_csv("Result/noCluster6.cluster.csv")

# save umap with cluster6 removed
tmp <- p1$data[,1:2]
a <- grep('^yki',row.names(tmp))
tmp$barcode <- ""
tmp[a,]$barcode <- gsub('-1','-2',gsub('yki_','',row.names(tmp)[a]))
tmp[-a,]$barcode <-gsub('wt_','',row.names(tmp)[-a])
write.csv(tmp, 'Result/noCluster6.umap.csv',quote = F,row.names = F)
