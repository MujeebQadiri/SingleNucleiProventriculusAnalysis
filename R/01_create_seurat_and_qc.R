# load package 
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(harmony)
  library(SCopeLoomR)
  library(patchwork)
  library(here)
  library(SoupX)
  library(Signac)
  library(gridExtra)
})

samples <- c("BEC1","BEC2")
sam_grp <- c('wt','yki')
in_data_dir <- './Data'

# in_data_dir <- 'Data'
seurat_list <- lapply(1:2, function(i){
  smpl<-samples[i]
  cur_data <- Read10X(file.path(in_data_dir,smpl,"outs/raw_feature_bc_matrix"))
  cur_seurat <- CreateSeuratObject(
    counts = cur_data,
    min.cells=3
  )
  cur_seurat$SampleID <- smpl
  
  # Novelty score
  cur_seurat$log10GenesPerUMI <- log10(cur_seurat$nFeature_RNA) / log10(cur_seurat$nCount_RNA)
  cur_seurat[["percent.mt"]]<-PercentageFeatureSet(cur_seurat, pattern="^mt:")
  cur_seurat[["percent.rb"]]<-PercentageFeatureSet(cur_seurat,pattern="^Rp[LS]")
  
  cur_seurat$condition <- sam_grp[i]
  print("---- Plot per dataset before QC and batch correction ----")
  # library size
  
  
  
  p1<- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(x=condition, fill=condition)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells")
  
  # plot UMI counts (transcripts) per cell
  p2 <- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(color=condition, fill= condition, x=nCount_RNA)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  
  # plot Genes detected per cell
  p3 <- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(color=condition, fill= condition, x=nFeature_RNA)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  
  # Complexity
  p4<- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(x=log10GenesPerUMI,color=condition, fill= condition)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
  
  
  p5<- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA,color=log10GenesPerUMI)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~condition)
  plot_list <- list(p1, p2, p3, p4, p5) 
  svg(here('tmp_result',
           paste0(smpl,".raw.prefilterQC.svg")))
  do.call("grid.arrange", c(plot_list, ncol = 1))  
  dev.off()
  return(cur_seurat)
})


allSample_raw <-  merge(x=seurat_list[[1]], 
                        y = seurat_list[[2]], 
                        add.cell.ids = sam_grp,  
                        project = "FC_08073")

allSample_normalized <- NormalizeData(allSample_raw)


allSample_normalized <- FindVariableFeatures(allSample_normalized, 
                                             selection.method = "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(allSample_normalized), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(allSample_normalized)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- row.names(allSample_normalized)


allSample_normalized <- ScaleData(allSample_normalized, 
                                  features = all.genes,
                                  vars.to.regress = c("nCount_RNA",
                                                      "percent.mt",
                                                      "percent.rb"))

allSample_normalized <- RunPCA(allSample_normalized, 
                               features = VariableFeatures(object = allSample_normalized))
allSample_normalized <- JackStraw(allSample_normalized, 
                                  num.replicate = 100,dims=50)
allSample_normalized <- ScoreJackStraw(allSample_normalized, 
                                       dims = 1:50)
JackStrawPlot(allSample_normalized, dims = 1:50)
ElbowPlot(allSample_normalized,ndims = 50)


allSample_normalized <- FindNeighbors(allSample_normalized, 
                                      dims = 1:40) %>%
  FindClusters(resolution = seq(0.1,1,0.1))


allSample_normalized <- RunUMAP(allSample_normalized, dims = 1:40)
allSample_normalized <- RunTSNE(allSample_normalized, dims = 1:40)

p2 <- UMAPPlot(allSample_normalized,label=T)
p1 <- TSNEPlot(allSample_normalized,label=T)

saveRDS(allSample_normalized,
        file = 'tmp_result/Seurat_all_default.rds')

