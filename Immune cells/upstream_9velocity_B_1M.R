rm(list = ls())
#library(xlsx)


library(Seurat)
library("velocyto.R")
library(loomR)
library(dplyr)
library(ggplot2)
library("stringr")

load("immune_recluster_4markers_annotated_0.3.Rdata")

all <- subset(all, subset = orig.ident %in% c("B_1M"))
#all <- RenameCells(all, new.names=gsub("A_7D_", "", colnames(all)))
#将剪切信息筛选为要做velocity的细胞
all <- subset(all, idents = c("Carbohydrate metabolic Macrophages",
                              "Chemotaxis Macrophages","Monocytes","Neutrophils"))
#all <- subset(all, downsample = 1000)
cell.chose <- as.data.frame(all$Cell_type)

cell.chose$barcodes <- row.names(cell.chose)



ldat.A14D <- read.loom.matrices("loom_file/cellsorted_B_1M_name_sorted_GJTKJ.loom")
emat.A14D <- ldat.A14D$spliced
colnames(emat.A14D) <- gsub("cellsorted_B_1M_name_sorted_GJTKJ:", "B_1M_", colnames(emat.A14D))
colnames(emat.A14D) <- gsub("x", "", colnames(emat.A14D))


emat <- emat.A14D


emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]



nmat.A14D <- ldat.A14D$unspliced
colnames(nmat.A14D) <- gsub("cellsorted_B_1M_name_sorted_GJTKJ:", "B_1M_", colnames(nmat.A14D))
colnames(nmat.A14D) <- gsub("x", "", colnames(nmat.A14D))

nmat <- nmat.A14D


nmat <- nmat[, colnames(nmat) %in% colnames(emat)]

cell.type <- cell.chose[colnames(emat),1]
#cell.type <- unique(all@meta.data[["Cell_type"]])

all.colors <- data.frame(
  cluster_name = c("B cells","Carbohydrate metabolic Macrophages","CCL24 Macrophages",
                   "Chemotaxis Macrophages","Monocytes","DCs",
                   "Neutrophils","T cells"),
  alphabet2 = c("#4991c1", "#afc2d9", "#c6307c", "#435b95", "#79b99d", "#e69f84", "#89558d","#d0afc4"))

#cell.type <- as.factor(cell.type)
cell.colors <- c()
for (i in 1:length(cell.type)) {
  cell.colors[i] <- all.colors$alphabet2[which(all.colors$cluster_name == cell.type[i])]
}

names(cell.colors) <- colnames(emat)



#确保细胞数和emat，nmat一致
all <- subset(all, cells = colnames(emat))
#计算细胞间的距离
cell.dist <- as.dist(1 - armaCor(t(all@reductions$pca@cell.embeddings)))




#计算RNA速率
rvel <- gene.relative.velocity.estimates(
  emat,nmat,
  deltaT = 1, kCells = 10,
  cell.dist = cell.dist,
  fit.quantile = 0.02, n.cores = 1)
save(rvel,file = "velocity_all/B_1M_rvel.Rdata")

#load("velocity_all/A_1M_rvel.Rdata")
#获取每个细胞在图上的坐标位置
emb <- all@reductions$umap@cell.embeddings





#pca降维并绘图
pdf(file = "velocity_all/B_1M_velocity.pdf", width = 10, height = 10)

show.velocity.on.embedding.cor(
  emb,rvel,n=300,
  scale='sqrt',
  cell.colors=ac(cell.colors,alpha=1),
  cex=0.8,
  arrow.scale=3,
  show.grid.flow=TRUE,
  min.grid.cell.mass=0.5,
  grid.n=40,
  arrow.lwd=2
)
dev.off()
