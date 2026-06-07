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

library(Seurat)
library(ggplot2)
library(dplyr,tibble)
load("muscle_regeneration_4markers_annotated_0.3.Rdata")
table(all$orig.ident)
table(all@meta.data$Cell_type)


all <- subset(all, idents = c("B cells","DCs","Monocytes","Myeloid-1",
                              "Myeloid-2","Myeloid-3","Myeloid-4","T cells"))

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
save(all,file="immune_recluster_4markers_0.3.Rdata")
#myGOmarker(all.markers,"allclus15")

library(pheatmap)
library(cols4all)
library(ComplexHeatmap)
library(circlize)

load("immune_recluster_4markers_annotated_0.3.Rdata")
sig_markers <- read.csv("cluster_markers5.csv")
#自定义配色：
mycol2<- colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))

#各亚群平均表达量提取
genes<- unique(sig_markers$gene)

aver_dt<- AverageExpression(all,
                            features= genes,
                            group.by = 'Cell_type',
                            slot= 'data')
aver_dt<- as.data.frame(aver_dt$RNA)
#aver_dt <- aver_dt[,c(4,2,7,3,5,8,6,1)]

#归一化：
aver_dtt<- t(scale(t(aver_dt)))
#添加行列注释：
cell_anno <- data.frame(cell_anno = colnames(aver_dt),
           row.names = colnames(aver_dt))

cols<- c("#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356", 
         "#16FF32", "#F7E1A0")
#cols<- c("#782AB6", "#3283FE", "#16FF32", "#85660D", "#565656", "#F7E1A0", 
#         "#1C8356", "#AA0DFE")

names(cols) <- cell_anno$cell_anno

#列注释：
cell<- data.frame(colnames(aver_dtt))
colnames(cell) <- 'cell'
col_anno<- HeatmapAnnotation(df = cell,
                             show_annotation_name= F,
                             gp= gpar(col = 'white', lwd = 2),
                             col= list(cell = cols))

#行注释：
row_cols<- setNames(rep(cols, each = 5), rownames(aver_dtt))
row_cols
row_anno<- rowAnnotation(foo = anno_text(rownames(aver_dtt),
                                         location= 0,
                                         just= "left",
                                         gp= gpar(fill = row_cols,
                                                  col= "black",
                                                  fontface= 'italic'),
                                         width= max_text_width(rownames(aver_dtt))*1.1))

#热图绘制：
aver_dtt <- aver_dtt[,c(4,2,7,3,5,8,6,1)]
pdf("Marker.top5.pdf", width = 8, height = 8)
Heatmap(aver_dtt,
        name= 'Expression',
        col= mycol2,
        cluster_columns= F,
        cluster_rows= F,
        show_column_names = F,
#        column_names_side= c('top'),
#        column_names_rot= 60,
        row_names_gp= gpar(fontsize = 12, fontface = 'italic'),
        rect_gp= gpar(col = "white", lwd = 1.5),
        top_annotation= col_anno) + row_anno
dev.off()
