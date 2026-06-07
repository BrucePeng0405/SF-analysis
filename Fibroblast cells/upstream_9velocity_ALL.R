rm(list = ls())
#library(xlsx)


library(Seurat)
library("velocyto.R")
library(loomR)
library(dplyr)
library(ggplot2)
library("stringr")

load("fibro_recluster_4markers_annotated_0.3.Rdata")

all <- subset(all, subset = orig.ident %in% "B_7D")
all <- RenameCells(all, new.names=gsub("B_7D_", "", colnames(all)))

#将剪切信息筛选为要做velocity的细胞
cell.chose <- as.data.frame(all$Cell_type)
cell.chose$barcodes <- row.names(cell.chose)
#cell.chose <- cell.chose[str_detect(row.names(cell.chose),"A_14D"),]
row.names(cell.chose) <- gsub("B_7D_", "", row.names(cell.chose))

ldat <- read.loom.matrices("loom_file/cellsorted_B_7D_name_sorted_U188K.loom")
emat <- ldat$spliced
colnames(emat) <- sapply(X = strsplit(colnames(emat), split = ":"), FUN = "[", 2)
colnames(emat) <- gsub("x", "", colnames(emat))

emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
nmat <- ldat$unspliced
colnames(nmat) <- sapply(X = strsplit(colnames(nmat), split = ":"), FUN = "[", 2)
colnames(nmat) <- gsub("x", "", colnames(nmat))
nmat <- nmat[, colnames(nmat) %in% colnames(emat)]

cell.type <- cell.chose[colnames(emat),1]

all.colors <- data.frame(
  cluster_name = c("Migration Myofibroblast","Adhesion Myofibroblast",
                   "Cycle Myofibroblast","Adipose Fibroblast",
                   "Inflammatory Fibroblast","Adhesion Fibroblast",
                   "Antigen-presenting Fibroblast"),
  alphabet2 = c("#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356", "#16FF32"))

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

#获取每个细胞在图上的坐标位置
emb <- all@reductions$umap@cell.embeddings
#全局速率可视化
pdf(file = "velocity_all/B_7D_velocity.pdf", width = 10, height = 8)
show.velocity.on.embedding.cor(
  emb,rvel,
  cell.colors = cell.colors,
  n = 300,
  scale = 'sqrt',
  cex = 0.8,
  arrow.scale = 3,
  show.grid.flow = T,
  min.grid.cell.mass = 0.5,
  grid.n = 40,
  arrow.lwd =1,
  do.par =F,
  cell.border.alpha = 0.1
)
dev.off()



load("fibro_recluster_4markers_annotated_0.3.Rdata")

all <- subset(all, subset = orig.ident %in% "B_14D")
all <- RenameCells(all, new.names=gsub("B_14D_", "", colnames(all)))

#将剪切信息筛选为要做velocity的细胞
cell.chose <- as.data.frame(all$Cell_type)
cell.chose$barcodes <- row.names(cell.chose)
#cell.chose <- cell.chose[str_detect(row.names(cell.chose),"A_14D"),]
row.names(cell.chose) <- gsub("B_14D_", "", row.names(cell.chose))

ldat <- read.loom.matrices("loom_file/cellsorted_B_14D_name_sorted_LFSUB.loom")
emat <- ldat$spliced
colnames(emat) <- sapply(X = strsplit(colnames(emat), split = ":"), FUN = "[", 2)
colnames(emat) <- gsub("x", "", colnames(emat))

emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
nmat <- ldat$unspliced
colnames(nmat) <- sapply(X = strsplit(colnames(nmat), split = ":"), FUN = "[", 2)
colnames(nmat) <- gsub("x", "", colnames(nmat))
nmat <- nmat[, colnames(nmat) %in% colnames(emat)]

cell.type <- cell.chose[colnames(emat),1]

all.colors <- data.frame(
  cluster_name = c("Migration Myofibroblast","Adhesion Myofibroblast",
                   "Cycle Myofibroblast","Adipose Fibroblast",
                   "Inflammatory Fibroblast","Adhesion Fibroblast",
                   "Antigen-presenting Fibroblast"),
  alphabet2 = c("#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356", "#16FF32"))

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

#获取每个细胞在图上的坐标位置
emb <- all@reductions$umap@cell.embeddings
#全局速率可视化
pdf(file = "velocity_all/B_14D_velocity.pdf", width = 10, height = 8)
show.velocity.on.embedding.cor(
  emb,rvel,
  cell.colors = cell.colors,
  n = 300,
  scale = 'sqrt',
  cex = 0.8,
  arrow.scale = 3,
  show.grid.flow = T,
  min.grid.cell.mass = 0.5,
  grid.n = 40,
  arrow.lwd =1,
  do.par =F,
  cell.border.alpha = 0.1
)
dev.off()


load("fibro_recluster_4markers_annotated_0.3.Rdata")

all <- subset(all, subset = orig.ident %in% "A_14D")
all <- RenameCells(all, new.names=gsub("A_14D_", "", colnames(all)))

#将剪切信息筛选为要做velocity的细胞
cell.chose <- as.data.frame(all$Cell_type)
cell.chose$barcodes <- row.names(cell.chose)
#cell.chose <- cell.chose[str_detect(row.names(cell.chose),"A_14D"),]
row.names(cell.chose) <- gsub("A_14D_", "", row.names(cell.chose))

ldat <- read.loom.matrices("loom_file/cellsorted_A_14D_name_sorted_5J72C.loom")
emat <- ldat$spliced
colnames(emat) <- sapply(X = strsplit(colnames(emat), split = ":"), FUN = "[", 2)
colnames(emat) <- gsub("x", "", colnames(emat))

emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
nmat <- ldat$unspliced
colnames(nmat) <- sapply(X = strsplit(colnames(nmat), split = ":"), FUN = "[", 2)
colnames(nmat) <- gsub("x", "", colnames(nmat))
nmat <- nmat[, colnames(nmat) %in% colnames(emat)]

cell.type <- cell.chose[colnames(emat),1]



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

#获取每个细胞在图上的坐标位置
emb <- all@reductions$umap@cell.embeddings
#全局速率可视化
pdf(file = "velocity_all/A_14D_velocity.pdf", width = 10, height = 8)
show.velocity.on.embedding.cor(
  emb,rvel,
  cell.colors = cell.colors,
  n = 300,
  scale = 'sqrt',
  cex = 0.8,
  arrow.scale = 3,
  show.grid.flow = T,
  min.grid.cell.mass = 0.5,
  grid.n = 40,
  arrow.lwd =1,
  do.par =F,
  cell.border.alpha = 0.1
)
dev.off()


load("fibro_recluster_4markers_annotated_0.3.Rdata")

all <- subset(all, subset = orig.ident %in% "A_7D")
all <- RenameCells(all, new.names=gsub("A_7D_", "", colnames(all)))

#将剪切信息筛选为要做velocity的细胞
cell.chose <- as.data.frame(all$Cell_type)
cell.chose$barcodes <- row.names(cell.chose)
#cell.chose <- cell.chose[str_detect(row.names(cell.chose),"A_7D"),]
row.names(cell.chose) <- gsub("A_7D_", "", row.names(cell.chose))

ldat <- read.loom.matrices("loom_file/cellsorted_A_7D_name_sorted_QH0U0.loom")
emat <- ldat$spliced
colnames(emat) <- sapply(X = strsplit(colnames(emat), split = ":"), FUN = "[", 2)
colnames(emat) <- gsub("x", "", colnames(emat))

emat <- emat[, colSums(emat) >= 1e3]
emat <- emat[, colnames(emat) %in% row.names(cell.chose)]
nmat <- ldat$unspliced
colnames(nmat) <- sapply(X = strsplit(colnames(nmat), split = ":"), FUN = "[", 2)
colnames(nmat) <- gsub("x", "", colnames(nmat))
nmat <- nmat[, colnames(nmat) %in% colnames(emat)]

cell.type <- cell.chose[colnames(emat),1]

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

#获取每个细胞在图上的坐标位置
emb <- all@reductions$umap@cell.embeddings
#全局速率可视化
pdf(file = "velocity_all/A_7D_velocity.pdf", width = 10, height = 8)
show.velocity.on.embedding.cor(
  emb,rvel,
  cell.colors = cell.colors,
  n = 300,
  scale = 'sqrt',
  cex = 0.8,
  arrow.scale = 3,
  show.grid.flow = T,
  min.grid.cell.mass = 0.5,
  grid.n = 40,
  arrow.lwd =1,
  do.par =F,
  cell.border.alpha = 0.1
)
dev.off()
