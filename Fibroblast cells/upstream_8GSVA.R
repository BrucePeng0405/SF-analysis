rm(list = ls())

#读取mdigdir文件，并整理为列表形式
library(Seurat)
library(GSEABase)
library(clusterProfiler)
#library(GSEA)
library(GSVA)
library(msigdbr)
library(tidyverse)
library(patchwork)

load("fibro_recluster_4markers_annotated_0.2.Rdata")
setwd("GSVA")
dir.create("avg_hallmark")
setwd("avg_hallmark/")

#选择基因集
msigdbr_species()
msigdbr_collections()
m_df <- msigdbr(species = "Homo sapiens",
                category = "C5",
                subcategory = "BP")
#按照gs_name给gene_symbol分组，并组成基因集列表形式
fgsea_set <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#计算cluster平均表达量并转换为矩阵
x <- AverageExpression(all)
exp1 <- as.matrix(x[["RNA"]])

#按cluster进行GSVA分析及可视化
es.max <- gsva(exp1, fgsea_set, mx.diff = F, kcdf = "Poisson")
write.csv(es.max, "GSVA_results.csv")



#可视化
library(pheatmap)
library(cols4all)
library(ComplexHeatmap)
library(circlize)


sig_markers <- read.csv("selected_GSVA_results.csv")
sig_markers$X <- gsub("GOBP_","",sig_markers$X)
sig_markers$X <- gsub("_"," ",sig_markers$X)
#自定义配色：
mycol2<- colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))

#各亚群平均表达量提取
genes<- unique(sig_markers$X)

aver_dt <- sig_markers[,2:8]
row.names(aver_dt) <- sig_markers$X


#aver_dt <- aver_dt[,c(4,2,7,3,5,8,6,1)]

#归一化：
aver_dtt<- t(scale(t(aver_dt)))
#添加行列注释：
cell_anno <- data.frame(cell_anno = colnames(aver_dt),
                        row.names = colnames(aver_dt))

cols<- c("#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356", 
         "#16FF32")
#cols<- c("#782AB6", "#3283FE", "#16FF32", "#85660D", "#565656", "#F7E1A0", 
#         "#1C8356", "#AA0DFE")

names(cols) <- cell_anno$cell_anno

#列注释：
cell<- data.frame(colnames(aver_dtt))
colnames(cell) <- 'cell'
col_anno<- HeatmapAnnotation(df = cell,
                             show_annotation_name= F,
                             show_legend = F,
                             gp= gpar(col = 'white', lwd = 2),
                             col= list(cell = cols))

#行注释：
row_cols<- setNames(rep(cols, each = 10), rownames(aver_dtt))
row_cols
row_anno<- rowAnnotation(foo = anno_text(rownames(aver_dtt),
                                         location= 0,
                                         just= "left",
                                         gp= gpar(fill = row_cols,
                                                  col= "black",
                                                  fontface= 'italic'),
                                         width= max_text_width(rownames(aver_dtt))*1.1))



#热图绘制：
#aver_dtt <- aver_dtt[,c(4,2,7,3,5,8,6,1)]
pdf("GSVA.pdf", width = 10, height = 6)
Heatmap(aver_dtt,
        name= 'Pathway activity',
        col= mycol2,
        cluster_columns= F,
        cluster_rows= F,
        show_column_names = F,
        row_names_side = c("left"),
        #        column_names_side= c('top'),
        #        column_names_rot= 60,
        row_names_gp= gpar(fontsize = 4.4, fontface = 'italic'),
        rect_gp= gpar(col = "white", lwd = 1),
        top_annotation= col_anno,
        width = ncol(aver_dtt)*unit(8, "mm"),
        height = ncol(aver_dtt)*unit(18, "mm"),
        heatmap_legend_param = list(title_position = "leftcenter-rot")
        #row_names_max_width = unit(1, "mm")
        ) 
#+ row_anno
dev.off()



