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


options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr,tibble)
load("batchremoved_raw.Rdata")
table(all$orig.ident)
gc()

MT <- c("COX1","COX2","COX3","ND1","ND2","ND3","ND4L","ND4",
        "ND5","ND6","CYTB","ATP6","ATP8",
        "CO1","CO2","CO3")
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = MT)
all[["log10GenesPerUMI"]] <- log10(all$nFeature_RNA) / log10(all$nCount_RNA)
VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4)
ggsave("rawdataFeature.pdf")
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
ggsave("rawdataFeature2.pdf")

# 为元数据添加细胞ID
metadata <- all@meta.data
metadata$cells <- rownames(metadata)
# 重命名列
metadata <- metadata %>%
  dplyr::rename(sample = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
metadata %>%
  ggplot(aes(color=sample, x=nUMI, fill= sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
ggsave("nUMI.pdf")

# 通过频数图可视化每个细胞检测出的基因数分布
metadata %>%
  ggplot(aes(color=sample, x=nGene, fill= sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)
ggsave("nGene.pdf")
# 通过箱线图可视化每个细胞检测到的基因的分布
metadata %>%
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
ggsave("nGene1.pdf")
# 可视化每个细胞检测到的线粒体基因表达分布
metadata %>%
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 5)
ggsave("percent.mt5.pdf")
all <- subset(all, subset = nFeature_RNA > 300  & nCount_RNA >500  & percent.mt < 5)
dim(all)
table(all$orig.ident)
gc()



ElbowPlot(all,ndims = 50 )
ggsave("ElbowPlot.pdf")


all <- FindNeighbors(all, dims = 1:20)
all <- FindClusters(all, resolution = 0.3)
all <- RunUMAP(all, dims = 1:20)
pdf("plots1/draft_umap_group.pdf", width = 10, height = 10)
UMAPPlot(object = all, group.by = "group")
dev.off()
pdf("plots1/draft_umap_timepoint.pdf", width = 10, height = 10)
UMAPPlot(object = all, group.by = "timepoint")
dev.off()
pdf("plots1/draft_umap.pdf", width = 10, height = 10)
UMAPPlot(object = all)
dev.off()

all.markers <- FindAllMarkers(all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(markers_top10,file = "cluster_markers10.csv")
save(all,file="muscle_regeneration_4markers_0.3.Rdata")
#myGOmarker(all.markers,"allclus15")

#proportion among group
#table(Idents(all))
#table(all$orig.ident)
#cell.prop<-as.data.frame(prop.table(table(Idents(all), all$orig.ident)))
#colnames(cell.prop)<-c("cluster","group","proportion")
#pdf("allProp_anno.pdf",width=12)
#Plotcellprop(cell.prop)
#dev.off()
