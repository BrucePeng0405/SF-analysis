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
load("immune_recluster_4markers_0.3.Rdata")
#all <- subset(all, downsample = 2000)
#gc()
#all1 <- subset(all, downsample = 2000)
#count <- as.data.frame(all1[["RNA"]]@counts)
#write.csv(count, "Countmatrix.csv")
#rm(count,all1)
#gc()

#immune cells
FeaturePlot(all,reduction = "umap", 
            features = c("PTPRC")) 
ggsave("annotation/Immune.cells.pdf")
VlnPlot(all,c("PTPRC"))
ggsave("annotation/Immune.cells_v.pdf")
#Myeloid cells
FeaturePlot(all,reduction = "umap", 
            features = c("CSF1R","AIF1"))
ggsave("annotation/Myeloid.cells.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("CSF1R","AIF1"))
ggsave("annotation/Myeloid.cells_v.pdf",height = 4.91,width = 6.65*2)

#Chemotaxis Macrophages-0
FeaturePlot(all,reduction = "umap", 
            features = c("CCL2")) 
ggsave("annotation/Chemotaxis.Macrophages.pdf",height = 4.91,width = 6.65)
VlnPlot(all,c("CCL2"))
ggsave("annotation/Chemotaxis.Macrophages_v.pdf",height = 4.91,width = 6.65)


#Carbohydrate metabolic Macrophages-1
FeaturePlot(all,reduction = "umap", 
            features = c("CHIT1","CHI3L1")) 
ggsave("annotation/Carbohydrate.Macrophages.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("CHIT1","CHI3L1"))
ggsave("annotation/Carbohydrate.Macrophages_v.pdf",height = 4.91,width = 6.65*2)

#Neutrophil-2
FeaturePlot(all,reduction = "umap", 
            features = c("S100A8","CXCL8")) 
ggsave("annotation/Neutrophils.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("S100A8","CXCL8"))
ggsave("annotation/Neutrophils_v.pdf",height = 4.91,width = 6.65*2)

#CCL24 Macrophages-3
FeaturePlot(all,reduction = "umap", 
            features = c("ENSOCUG00000025061")) 
ggsave("annotation/CCL24.Macrophages.pdf",height = 4.91,width = 6.65)
VlnPlot(all,c("ENSOCUG00000025061"))
ggsave("annotation/CCL24.Macrophages_v.pdf",height = 4.91,width = 6.65)

#Monocytes-4
FeaturePlot(all,reduction = "umap", 
            features = c("CD14","TPX2","SMC4")) 
ggsave("annotation/Monocytes.pdf",height = 4.91*2,width = 6.65*2)
VlnPlot(all,c("CD14","TPX2","SMC4"))
ggsave("annotation/Monocytes_v.pdf",height = 4.91,width = 6.65*3)

#T cells-5
FeaturePlot(all,reduction = "umap", 
            features = c("GATA3","CD69")) 
ggsave("annotation/T.cells.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("GATA3","CD69"))
ggsave("annotation/T.cells_v.pdf",height = 4.91,width = 6.65*2)

#DCs-6
FeaturePlot(all,reduction = "umap", 
            features = c("HLA-DOB","HLA-DRA")) 
ggsave("annotation/DCs.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("HLA-DOB","HLA-DRA"))
ggsave("annotation/DCs_v.pdf",height = 4.91,width = 6.65*2)

#B cells-7
FeaturePlot(all,reduction = "umap", 
            features = c("CD79B")) 
ggsave("annotation/B.cells.pdf",height = 4.91,width = 6.65)
VlnPlot(all,c("CD79B"))
ggsave("annotation/B.cells_v.pdf",height = 4.91,width = 6.65)



# ALL marker
FeaturePlot(all,reduction = "umap", 
            features = c("PTPRC","CSF1R","AIF1","CCL2","CHIT1","CHI3L1","S100A8","CXCL8",
                         "ENSOCUG00000025061","CD14","TPX2","GATA3","CD69","HLA-DOB","HLA-DRA","CD79B"))
ggsave("annotation/ALL.cells.pdf",height = 4.91*4,width = 6.65*4)
VlnPlot(all,c("PTPRC","CSF1R","AIF1","CCL2","CHIT1","CHI3L1","S100A8","CXCL8",
              "ENSOCUG00000025061","CD14","TPX2","GATA3","CD69","HLA-DOB","HLA-DRA","CD79B"))
ggsave("annotation/ALL.cells_v.pdf",height = 4.91*4,width = 6.65*4)



cluster_name <- c("Chemotaxis Macrophages","Carbohydrate metabolic Macrophages",
                  "Neutrophils","CCL24 Macrophages","Monocytes",
                  "T cells", "DCs","B cells")


names(cluster_name) <- c(0:7)
all <- RenameIdents(all, cluster_name)
DimPlot(all, reduction = "umap", label = TRUE, pt.size = 0.5)+NoLegend()
ggsave("annotation/umap.clus.rename.pdf",height = 4.91,width = 6.65)
DimPlot(all, reduction = "umap",label = F, split.by = "orig.ident", ncol = 4)
ggsave("annotation/umap.split.group.pdf",width = 6.65*4,height = 4.91*2)

#Cell type注释
Cell_type <- data.frame(ClusterID=c(0:7), Cell_type=cluster_name, stringsAsFactors = F)
all@meta.data$Cell_type <- "NA"
for (i in 1:nrow(Cell_type)) {
  all@meta.data[which(all@meta.data$seurat_clusters == Cell_type$ClusterID[i]),"Cell_type"] <- Cell_type$Cell_type[i]
}
all@meta.data$Cell_type
DimPlot(all,group.by = "Cell_type",label = F,reduction = "umap")
ggsave("annotation/umap.clus.celltype.pdf",height = 4.91,width = 6.65)
save(all, file = "immune_recluster_4markers_annotated_0.3.Rdata")
