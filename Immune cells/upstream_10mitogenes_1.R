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
library(tidyverse)
library(cowplot)


load("immune_recluster_4markers_annotated_0.3.Rdata")
all <- subset(all, idents = c("Carbohydrate metabolic Macrophages",
                              "Chemotaxis Macrophages",
                              "Monocytes"))
gc()


### Cell type
CELLS <- c("Carbohydrate metabolic Macrophages",
           "Chemotaxis Macrophages",
           "Monocytes")

all@meta.data$Cell_type <- factor(all@meta.data$Cell_type,  levels=CELLS, ordered=T)

all@meta.data$group <- str_replace(all@meta.data$group, "A","SF")
all@meta.data$group <- str_replace(all@meta.data$group, "B","ECM")
all@meta.data$group <- factor(all@meta.data$group,  levels=c("SF","ECM"), ordered=T)

cellGroup <- paste0(all@meta.data$Cell_type, "_", all@meta.data$group)
m <- data.frame("Celltype_Group" = cellGroup)
rownames(m) <- rownames(all@meta.data)
all <- AddMetaData(object = all, metadata = m)

# Assign as main identity
Idents(all) <- "Celltype_Group"





modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
    theme_void() +
    ylab(feature) +
    theme(legend.position = 'none',
          plot.margin = plot.margin,
          title = element_blank(),
          axis.title.y = element_text(hjust = 0.5, angle = 0)
    )
  return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1, size = 15), axis.ticks.x = element_line())
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

"#435b95", "#afc2d9", "#89558d", "#c6307c", "#79b99d", "#d0afc4", 
"#e69f84","#4991c1"
#配色方案(需要注意顺序，从0-7)
#all.colors <- c("#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356", 
#                "#16FF32", "#F7E1A0")
all.colors<- c("#4991c1", "#afc2d9", "#c6307c", "#435b95", "#79b99d", "#e69f84",
               "#89558d", "#d0afc4")
Idents(all) <- factor(all$Cell_type,levels = c("B cells","Carbohydrate metabolic Macrophages",
                                               "CCL24 Macrophages","Chemotaxis Macrophages",
                                               "Monocytes","DCs","Neutrophils","T cells"))
StackedVlnPlot(all, 
               c("ATP6","ND3","CYTB","COX3","ND4","COX2","ND5"), 
               pt.size=0, cols=all.colors)
ggsave("prettyumap/mitognes_Dotplot.pdf",height = 14,width = 10)

markers <- c("ATP6","ND3","CYTB","COX3","ND4","COX2","ND5")
sc_marker <- data.frame(gene = markers)
exp <- GetAssayData(all, slot = "counts")
exp <- log10(exp+1)
head(all$Cell_type)
new_celltype <- sort(all$Cell_type)
group <- sort(all$group)
cs_data <- as.matrix(exp[sc_marker$gene,names(new_celltype)])
ac <- data.frame(cluster = new_celltype,group = group)
rownames(ac) = colnames(cs_data)

library(pheatmap)
pheatmap(cs_data,show_colnames = F,show_rownames = T,
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = ac,
         border_color = "black")


