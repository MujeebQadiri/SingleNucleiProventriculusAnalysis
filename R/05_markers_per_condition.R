library(Seurat)
library(tidyverse)
library(foreach)
library(doParallel)

seuratObj <- readRDS("Result/seurat_annotated.rds")
Idents(seuratObj) <- "celltype"

registerDoParallel(6)

m1 <- foreach(cluster=unique(Idents(seuratObj)),
              .combine = rbind)%dopar%{
                tmp <- FindMarkers(subset(seuratObj, condition=='wt'),
                                   ident.1 = cluster)
                tmp$gene <- row.names(tmp)
                tmp$cluster <- cluster
                
                return(tmp)
              }
m1$condition <- "wt"


m2 <- foreach(cluster=unique(Idents(seuratObj)), .combine = rbind)%dopar%{
  tmp <- FindMarkers(subset(seuratObj, condition!='wt'),
                     ident.1 = cluster)
  tmp$gene <- row.names(tmp)
  tmp$cluster <- cluster
  
  return(tmp)
}
m2$condition <- "yki"
m2%>%arrange(substring(cluster,1,2), p_val_adj, desc(pct.1))
stopImplicitCluster()

output <- rbind(m1%>%arrange(substring(cluster,1,2), p_val_adj, desc(pct.1)),
                m2%>%arrange(substring(cluster,1,2), p_val_adj, desc(pct.1)))

output %>% arrange(condition, substring(cluster,1,2), p_val_adj, desc(pct.1-pct.2))%>%
  write_tsv("Result/FC_08073.marker.txt")
