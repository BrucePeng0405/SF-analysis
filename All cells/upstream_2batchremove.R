#!/share/home/ouyanghongwei/bin/Rscript
#BSUB -J pz
#BSUB -n 5
#BSUB -R "span[ptile=24]"
#BSUB -o output_%Js
#BSUB -e errput_%J
#BSUB -q normal
rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)


load("mergedraw.Rdata")
#pbmc <- subset(all, downsample = 1000)
head(colnames(all))
all$group <- sapply(X = strsplit(colnames(all), split = "_"), FUN = "[", 1)
all$timepoint <- sapply(X = strsplit(colnames(all), split = "_"), FUN = "[", 2)
table(all$orig.ident)
all <- subset(all, downsample = 8300)
table(all$orig.ident)
gc()

###########harmony 速度快、内存少################
all <-  all%>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData()
all <- RunPCA(all, npcs = 50, verbose = FALSE)

all <- RunTSNE(all, dims = 1:20)
p3 <- DimPlot(all, reduction = "tsne", group.by = "orig.ident", pt.size=0.5)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)

#if(!require(harmony))devtools::install_github("immunogenomics/harmony")



#####run 到PCA再进行harmony，相当于降维########
all=all %>% RunHarmony("orig.ident", plot_convergence = TRUE)


all <- all %>% 
  RunTSNE(reduction = "harmony", dims = 1:20)
save(all,file="batchremoved_raw.Rdata")

p4 <- DimPlot(all, reduction = "tsne", group.by = "orig.ident", pt.size=0.5)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
p3|p4
ggsave("batch_remove.pdf",dpi = 300)

all <- FindNeighbors(all, dims = 1:20)
all <- FindClusters(all, resolution = 0.3)
DimPlot(all, reduction = "tsne",label = T)
ggsave("tsne.clus.pdf")
DimPlot(all, reduction = "tsne",group.by = "orig.ident")
ggsave("tsne.group.pdf")
DimPlot(all, reduction = "tsne",label = T, split.by = "orig.ident")
ggsave("tsne.split.group.pdf",width = 40,height = 6.28)
