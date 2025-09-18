library(Seurat)
library(dplyr)

# Read in data
Seurat_all_default <- readRDS('tmp_result/Seurat_all_default.rds')

# Find cells that has high nCount & nFeature 
meta.data <- Seurat_all_default@meta.data
potential_doublet <- row.names(meta.data[meta.data$nFeature_RNA>median((meta.data$nFeature_RNA))+2.5*sd(meta.data$nFeature_RNA) | 
                                           meta.data$nCount_RNA > median((meta.data$nCount_RNA))+2.5*sd(meta.data$nCount_RNA),])

counts <- Seurat_all_default[,!colnames(Seurat_all_default)%in%potential_doublet]@assays$RNA@counts

allSample_normalized <- NormalizeData(cur_seurat)
allSample_normalized <- FindVariableFeatures(allSample_normalized, 
                                             selection.method = "vst", nfeatures = 2000)
allSample_normalized <- ScaleData(allSample_normalized)
allSample_normalized <- RunPCA(allSample_normalized)

# determine minimum PC needed
stdv <- allSample_normalized[["pca"]]@stdev
sum.stdv <- sum(allSample_normalized[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                     percent.stdv[2:length(percent.stdv)]) > 0.1), 
            decreasing = T)[1] + 1
min.pc <- min(co1, co2)
min.pc # 18

# finish pre-processing
allSample_normalized <- RunUMAP(allSample_normalized, dims = 1:min.pc,
                                umap.method = "umap-learn")
allSample_normalized <- FindNeighbors(allSample_normalized, dims = 1:min.pc)  
saveRDS(allSample_normalized,'tmp_result/highReadRemoved.rds')

for (rs in seq(.1,.3,.1)){
  allSample_normalized <- FindClusters(allSample_normalized, 
                                       resolution = rs,
                                       algorithm = 4)
  saveRDS(allSample_normalized,'tmp_result/highReadRemoved.rds')
  gc()
}

# cluster at different resolutions
tmp <- (allSample_normalized@meta.data[,c(7,9,10)])
a <- grep('^yki',row.names(tmp))
tmp$barcode <- ""
tmp[a,]$barcode <- gsub('-1','-2',gsub('yki_','',row.names(tmp)[a]))
tmp[-a,]$barcode <-gsub('wt_','',row.names(tmp)[-a])
write.csv(tmp, 'loupe_files/highReadRemoved_noBC.cluster.csv',quote = F,row.names = F)

# umap coordinates
p <- UMAPPlot(allSample_normalized,label=T,split.by='orig.ident')
tmp <- p$data[,1:2]
a <- grep('^yki',row.names(tmp))
tmp$barcode <- ""
tmp[a,]$barcode <- gsub('-1','-2',gsub('yki_','',row.names(tmp)[a]))
tmp[-a,]$barcode <-gsub('wt_','',row.names(tmp)[-a])
write.csv(tmp, 'loupe_files/highReadRemoved_noBC.umap.csv',quote = F,row.names = F)