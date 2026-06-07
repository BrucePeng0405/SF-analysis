rm(list = ls())

install.packages('D:/velocyto.R', repos = NULL, type = 'source')

library(velocyto.R)
library(Seurat)
library(dplyr)
library(ggplot2)

ldat <- read.loom.matrices("loom_file/cellsorted_A_14D_name_sorted_5J72C.loom")
load("F:/scrna_seq_raw_2'/fibro_recluster_4markers_annotated_0.3.Rdata")

celltypes <- as.data.frame(all@active.ident)
celltypes$barcodes <- row.names(celltypes)
celltypes <- celltypes[grep("A_14D",row.names(celltypes)),]
celltypes$barcodes <- gsub("A_14D_", "", celltypes$barcodes)

colnames(ldat) <- gsub("cellsorted_A_14D_name_sorted_5J72C:", "", colnames(ldat))


emat <- ldat$spliced
emat <- emat[, colSums(emat) >= 1e3]
nmat <- ldat$unspliced

all <- CreateSeuratObject(emat, assay = "spliced", project= "10x")
table(all$orig.ident)

all <-  all%>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData()
all <- RunPCA(all, npcs = 50, verbose = FALSE)

ElbowPlot(all,ndims = 50 )
ggsave("velocity/ElbowPlot.pdf")

all <- FindNeighbors(all, dims = 1:20)
all <- FindClusters(all, resolution = 0.3)
all <- RunUMAP(all, dims = 1:20)
pdf("velocity/draft_umap.pdf", width = 10, height = 10)
UMAPPlot(object = all)
dev.off()

# ALL marker
FeaturePlot(all,reduction = "umap", 
            features = c("LUM","DCN","COL1A2"))
ggsave("velocity/fibro.cells.pdf",height = 4.91*2,width = 6.65*2)
VlnPlot(all,c("LUM","DCN","COL1A2"))
ggsave("velocity/fibro.cells_v.pdf",height = 4.91*2,width = 6.65*2)

cluster_name <- c("M2 Macro","VSMCs","Monocytes",
                  "Fibro-myofibro-1","Neutrophils","Fibroblast-1",
                  "Fibro-myofibro-2","Fibroblast-2","Fibro-myofibro-3",
                  "M1 Macro","Endothelial","Spp1+ monocytes",
                  "T cells","DCs","B cells",
                  "Satellite cells","Lymphatic Endothelial cells","Germinal center B cell")
names(cluster_name) <- c(0:17)
all <- RenameIdents(all, cluster_name)
#Cell type注释
Cell_type <- data.frame(ClusterID=c(0:17), Cell_type=cluster_name, stringsAsFactors = F)
all@meta.data$Cell_type <- "NA"
for (i in 1:nrow(Cell_type)) {
  all@meta.data[which(all@meta.data$seurat_clusters == Cell_type$ClusterID[i]),"Cell_type"] <- Cell_type$Cell_type[i]
}
all@meta.data$Cell_type

all <- subset(all, idents = c("Fibro-myofibro-1","Fibroblast-1",
                              "Fibro-myofibro-2","Fibroblast-2","Fibro-myofibro-3"))
gc()

all <-  all%>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData()
all <- RunPCA(all, npcs = 50, verbose = FALSE)
ElbowPlot(all,ndims = 50 )
ggsave("velocity/fibro/ElbowPlot.pdf")

all <- FindNeighbors(all, dims = 1:20)
all <- FindClusters(all, resolution = 0.3)
all <- RunUMAP(all, dims = 1:20)
pdf("velocity/fibro/draft_umap.pdf", width = 10, height = 10)
UMAPPlot(object = all)
dev.off()

FeaturePlot(all,reduction = "umap", 
            features = c("ACTA2","FN1","ACTG2","APOD","ND4",
                         "CDK1","TNMD","CHI3L1",
                         "S100A4","AIF1"))
ggsave("velocity/fibro/ALL.cells.pdf",height = 4.91*3,width = 6.65*4)
VlnPlot(all,c("ACTA2","FN1","ACTG2","APOD","ND4",
              "CDK1","TNMD","CHI3L1",
              "S100A4","AIF1"))
ggsave("velocity/fibro/ALL.cells_v.pdf",height = 4.91*3,width = 6.65*4)

cluster_name <- c("Migration Myofibroblast","Cytoskeleton Myofibroblast",
                  "Adipose Fibroblast","Mitochondrial Myofibroblast","Cycle Myofibroblast1",
                  "Angiogenesis Fibroblast","Inflammatory Angiogenesis Fibroblast","Cycle Myofibroblast2","NFkb Fibroblast","Inflammatory Fibroblast")


names(cluster_name) <- c(0:9)
all <- RenameIdents(all, cluster_name)
DimPlot(all, reduction = "umap", label = TRUE, pt.size = 0.5)+NoLegend()
ggsave("velocity/fibro/umap.clus.rename.pdf",height = 4.91,width = 6.65)
Cell_type <- data.frame(ClusterID=c(0:9), Cell_type=cluster_name, stringsAsFactors = F)
all@meta.data$Cell_type <- "NA"
for (i in 1:nrow(Cell_type)) {
  all@meta.data[which(all@meta.data$seurat_clusters == Cell_type$ClusterID[i]),"Cell_type"] <- Cell_type$Cell_type[i]
}
all@meta.data$Cell_type

cluster.label = Idents(all)
cell.colors <- rainbow(10)[cluster.label]
names(cell.colors) <- colnames(all)

#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(all@reduction$pca@cell.embeddings)))
#计算RNA速率
rvel <- gene.relative.velocity.estimates(
  emat,nmat,
  deltaT = 1, kCells = 10,
  cell.dist = cell.dist,
  fit.quantile = 0.02, n.cores = 1)

#获取每个细胞在图上的坐标位置
emb <- all@reduction$umap@cell.embeddings
#全局速率可视化
pdf(file = "velocity/fibro/A_14D_velocity.pdf", width = 10, height = 8)
show.velocity.on.embedding.cor(
  emb,rvel,
  cell.colors = cell.colors,
  n = 100,
  scale = 'sqrt',
  cex = 1.2,
  arrow.scale = 0.8,
  show.grid.flow = T,
  min.grid.cell.mass = 0.5,
  grid.n = 20,
  arrow.lwd =1,
  do.par =T,
  cell.border.alpha = 0.5
)
dev.off()




