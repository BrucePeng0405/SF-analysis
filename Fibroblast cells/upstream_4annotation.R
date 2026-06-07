#!/share/home/ouyanghongwei/bin/Rscript
#BSUB -J pz
#BSUB -n 5
#BSUB -R "span[ptile=24]"
#BSUB -o output_%Js
#BSUB -e errput_%J
#BSUB -q normal
rm(list = ls())
#.libPaths("~/lib/R")
print(.libPaths())

#setwd("/share/home/ouyanghongwei/workspace/PZ/")
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr,tibble)
load("fibro_recluster_4markers_0.2.Rdata")
#all <- subset(all, downsample = 2000)
#gc()
#all1 <- subset(all, downsample = 2000)
#count <- as.data.frame(all1[["RNA"]]@counts)
#write.csv(count, "Countmatrix.csv")
#rm(count,all1)
#gc()

#myofibroblast cells-0,1,2
FeaturePlot(all,reduction = "umap", 
            features = c("ACTA2"))
ggsave("annotation/myofibro.cells.pdf",height = 4.91,width = 6.65)
VlnPlot(all,c("ACTA2"))
ggsave("annotation/myofibro.cells_v.pdf",height = 4.91,width = 6.65)

#Migration cells-0
FeaturePlot(all,reduction = "umap", 
            features = c("ACTA2","ACTG2"))
ggsave("annotation/Migration.Myofibro.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("ACTA2","ACTG2"))
ggsave("annotation/Migration.Myofibro_v.pdf",height = 4.91,width = 6.65*2)

#Adhesion cells-1
FeaturePlot(all,reduction = "umap", 
            features = c("COL12A1","TGFBI"))
ggsave("annotation/Adhesion.Myofibro.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("COL12A1","TGFBI"))
ggsave("annotation/Adhesion.Myofibro_v.pdf",height = 4.91,width = 6.65*2)

#Cycle cells-2
FeaturePlot(all,reduction = "umap", 
            features = c("CCNB2","CDK1"))
ggsave("annotation/Cycle.cells.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("CCNB2","CDK1"))
ggsave("annotation/Cycle.cells_v.pdf",height = 4.91,width = 6.65*2)

#adipose cells-3
FeaturePlot(all,reduction = "umap", 
            features = c("APOD","CD34"))
ggsave("annotation/adipose.cells.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("APOD","CD34"))
ggsave("annotation/adipose.cells_v.pdf",height = 4.91,width = 6.65*2)

#Inflammatory cells-4
FeaturePlot(all,reduction = "umap", 
            features = c("CD14","CHI3L1")) 
ggsave("annotation/Inflammatory.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("CD14","CHI3L1"))
ggsave("annotation/Inflammatory_v.pdf",height = 4.91,width = 6.65*2)

#Adhesion cells-5
FeaturePlot(all,reduction = "umap", 
            features = c("COL8A1","THBS4"))
ggsave("annotation/Adhesion.Fibro.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("COL8A1","THBS4"))
ggsave("annotation/Adhesion.Fibro_v.pdf",height = 4.91,width = 6.65*2)

#Antigen-presenting cells-6
FeaturePlot(all,reduction = "umap", 
            features = c("CD74","ENSOCUG00000013412"))
ggsave("annotation/Antigen-presenting.cells.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("CD74","ENSOCUG00000013412"))
ggsave("annotation/Antigen-presenting.cells_v.pdf",height = 4.91,width = 6.65*2)

# ALL marker
FeaturePlot(all,reduction = "umap", 
            features = c("ACTA2","ACTG2","COL12A1","TGFBI","CCNB2","CDK1",
                         "APOD","CD34","CD14","CHI3L1","COL8A1","THBS4",
                         "CD74","ENSOCUG00000013412"))
ggsave("annotation/ALL.cells.pdf",height = 4.91*4,width = 6.65*4)
VlnPlot(all,c("ACTA2","ACTG2","COL12A1","TGFBI","CCNB2","CDK1",
              "APOD","CD34","CD14","CHI3L1","COL8A1","THBS4",
              "CD74","ENSOCUG00000013412"))
ggsave("annotation/ALL.cells_v.pdf",height = 4.91*4,width = 6.65*4)



cluster_name <- c("Migration Myofibroblasts","Adhesion Myofibroblasts",
                  "Cycle Myofibroblasts","Adipose Fibroblasts",
                  "Inflammatory Fibroblasts","Adhesion Fibroblasts",
                  "Antigen-presenting Fibroblasts")


names(cluster_name) <- c(0:6)
all <- RenameIdents(all, cluster_name)
DimPlot(all, reduction = "umap", label = TRUE, pt.size = 0.5)+NoLegend()
ggsave("annotation/umap.clus.rename.pdf",height = 4.91,width = 6.65)
DimPlot(all, reduction = "umap",label = F, split.by = "orig.ident", ncol = 4)
ggsave("annotation/umap.split.group.pdf",width = 6.65*4,height = 4.91*2)

#Cell type注释
Cell_type <- data.frame(ClusterID=c(0:6), Cell_type=cluster_name, stringsAsFactors = F)
all@meta.data$Cell_type <- "NA"
for (i in 1:nrow(Cell_type)) {
  all@meta.data[which(all@meta.data$seurat_clusters == Cell_type$ClusterID[i]),"Cell_type"] <- Cell_type$Cell_type[i]
}
all@meta.data$Cell_type
DimPlot(all,group.by = "Cell_type",label = F,reduction = "umap")
ggsave("annotation/umap.clus.celltype.pdf",height = 4.91,width = 6.65)
save(all, file = "fibro_recluster_4markers_annotated_0.2.Rdata")
