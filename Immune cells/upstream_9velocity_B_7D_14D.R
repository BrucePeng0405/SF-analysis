rm(list = ls())
#library(xlsx)


library(Seurat)
library("velocyto.R")
library(loomR)
library(dplyr)
library(ggplot2)
library("stringr")

load("immune_recluster_4markers_annotated_0.3.Rdata")

all <- subset(all, subset = orig.ident %in% c("B_7D","B_14D"))
#all <- RenameCells(all, new.names=gsub("A_7D_", "", colnames(all)))
#将剪切信息筛选为要做velocity的细胞
all <- subset(all, idents = c("Carbohydrate metabolic Macrophages","CCL24 Macrophages",
                              "Chemotaxis Macrophages","Cycle Macrophages","Neutrophils"))
#all <- subset(all, downsample = 1000)
cell.chose <- as.data.frame(all$Cell_type)

cell.chose$barcodes <- row.names(cell.chose)

ldat.A7D <- read.loom.matrices("loom_file/cellsorted_B_7D_name_sorted_U188K.loom")
emat.A7D <- ldat.A7D$spliced
colnames(emat.A7D) <- gsub("cellsorted_B_7D_name_sorted_U188K:", "B_7D_", colnames(emat.A7D))
colnames(emat.A7D) <- gsub("x", "", colnames(emat.A7D))

ldat.A14D <- read.loom.matrices("loom_file/cellsorted_B_14D_name_sorted_LFSUB.loom")
emat.A14D <- ldat.A14D$spliced
colnames(emat.A14D) <- gsub("cellsorted_B_14D_name_sorted_LFSUB:", "B_14D_", colnames(emat.A14D))
colnames(emat.A14D) <- gsub("x", "", colnames(emat.A14D))

emat <- cbind(emat.A7D,emat.A14D)
#emat <- emat.A14D
#emat <- emat.A7D

emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]



nmat.A7D <- ldat.A7D$unspliced
colnames(nmat.A7D) <- gsub("cellsorted_B_7D_name_sorted_U188K:", "B_7D_", colnames(nmat.A7D))
colnames(nmat.A7D) <- gsub("x", "", colnames(nmat.A7D))

nmat.A14D <- ldat.A14D$unspliced
colnames(nmat.A14D) <- gsub("cellsorted_B_14D_name_sorted_LFSUB:", "B_14D_", colnames(nmat.A14D))
colnames(nmat.A14D) <- gsub("x", "", colnames(nmat.A14D))
nmat <- cbind(nmat.A7D,nmat.A14D)
#nmat <- nmat.A14D
#nmat <- nmat.A7D

nmat <- nmat[, colnames(nmat) %in% colnames(emat)]

cell.type <- cell.chose[colnames(emat),1]

all.colors <- data.frame(
  cluster_name = c("B cells","Carbohydrate metabolic Macrophages","CCL24 Macrophages",
                   "Chemotaxis Macrophages","Cycle Macrophages","DCs",
                   "Neutrophils","T cells"),
  alphabet2 = c("#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356", "#16FF32"
                , "#F7E1A0"))

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
save(rvel,file = "velocity_all/B_7D_14D_rvel.Rdata")

#获取每个细胞在图上的坐标位置
emb <- all@reductions$umap@cell.embeddings
#全局速率可视化
pdf(file = "velocity_all/B_7D_14D_velocity.pdf", width = 10, height = 10)

show.velocity.on.embedding.cor(
  emb,rvel,n=300,
  scale='sqrt',
  cell.colors=ac(cell.colors,alpha=0.8),
  cex=0.8,
  arrow.scale=3,
  show.grid.flow=TRUE,
  min.grid.cell.mass=0.5,
  grid.n=40,
  arrow.lwd=2
)
dev.off()


#pca降维并绘图
#pdf(file = "velocity_all/B_14D_velocity_PCA.pdf", width = 10, height = 10)
#pca.velocity.plot(rvel,nPcs=5,
#                  plot.cols=2,
#                  cell.colors=ac(cell.colors,alpha=0.5),
#                  cex=0.8,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1),
#                  arrow.scale = 3,
#                  arrow.lwd = 2)
#dev.off()


