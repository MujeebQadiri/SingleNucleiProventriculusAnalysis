library(Seurat)
library(tidyverse)

# load 
# clusters annotated using marker gene analysis at resolution 0.15
cluster <- read_csv("Result/cluster.csv")
seuratObj <- readRDS("Result/noCluster6.rds")

# merge
metadata <- (seuratObj@meta.data)
metadata$orig.ident <- row.names(metadata)
tmp <- merge(cluster, metadata[,c(1:6,8)],
             by.y="RNA_snn_res.0.15", by.x="Cluster")
tmp2 <- (tmp[match(row.names(seuratObj@meta.data), tmp$orig.ident),])
row.names(tmp2) <- tmp2$orig.ident
tmp2$condition <- "wt"
tmp2$condition[grep('yki', tmp2$orig.ident)]<- "yki"

seuratObj@meta.data <- tmp2[, setdiff(names(tmp2),"orig.ident")]

saveRDS(seuratObj, "tmp_result/seurat_annotated.rds")
