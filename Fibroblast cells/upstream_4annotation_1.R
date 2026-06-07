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
load("fibro_recluster_4markers_annotated_0.2.Rdata")
#all <- subset(all, downsample = 2000)
#gc()
#all1 <- subset(all, downsample = 2000)
#count <- as.data.frame(all1[["RNA"]]@counts)
#write.csv(count, "Countmatrix.csv")
#rm(count,all1)
#gc()


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

#配色方案(需要注意顺序，从0-7)
all.colors <- c("#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356", 
                "#16FF32")
Idents(all) <- factor(all$Cell_type,levels = c("Migration Myofibroblasts","Adhesion Myofibroblasts",
                                               "Cycle Myofibroblasts","Adipose Fibroblasts",
                                               "Inflammatory Fibroblasts","Adhesion Fibroblasts",
                                               "Antigen-presenting Fibroblasts"))



StackedVlnPlot(all, 
               c("ACTA2","ACTG2","COL12A1","TGFBI","CCNB2","CDK1","APOD","CD34",
                 "CD14","CHI3L1","COL8A1","THBS4","CD74","ENSOCUG00000013412"), 
               pt.size=0, cols=all.colors)
ggsave("annotation/vlolin_all_1.pdf",height = 14,width = 10)

