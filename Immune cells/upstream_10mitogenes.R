# Make pretty tSNE's for all 6 combined 10x replicates.
# Figures 1ab and SupFig 1c
# Matthew Buckley
rm(list = ls())
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(dplyr,tibble)
#library(plyr)

# Load Data
#setwd("~/10X/Seurat/")
load("immune_recluster_4markers_annotated_0.3.Rdata")
#d <- read.csv(file = "F:/scrna_seq_raw_1/10X/Seurat/all_All6_Filtered_meta.data&umap2022-11-18.csv", header=T)
#d <- d[,-1]
#d[1:3,1:11]

# Randomize rows to mix up plotting order
#d <- d[sample(nrow(d)),]

all <- subset(all, idents = c("Carbohydrate metabolic Macrophages",
                              "Chemotaxis Macrophages",
                              "Monocytes"))
gc()


### Cell type
CELLS <- c("Carbohydrate metabolic Macrophages",
           "Chemotaxis Macrophages",
           "Monocytes")

all@meta.data$Cell_type <- factor(all@meta.data$Cell_type,  levels=CELLS, ordered=T)
#DimPlot(object = all, reduction = "umap", cols = "alphabet2", 
#        group.by = "Cell_type",label = T)
#ggsave(paste0("prettyumap/customeUMAP_celltype_axis_nolegend", Sys.Date(), ".pdf"), height = 7, width = 10)
###  DiscretePalette(18, palette = "alphabet2", shuffle = FALSE)
#[1] "#AA0DFE" "#3283FE" "#85660D" "#782AB6" "#565656" "#1C8356" "#16FF32" "#F7E1A0" "#E2E2E2"
#[10] "#1CBE4F" "#C4451C" "#DEA0FD" "#FE00FA" "#325A9B" "#FEAF16" "#F8A19F" "#90AD1C" "#F6222E"
alphabet2<- c("#afc2d9", "#435b95", "#79b99d")

#"#435b95", "#afc2d9", "#89558d", "#c6307c", "#79b99d", "#d0afc4", 
#"#e69f84","#4991c1"


all@meta.data$group <- str_replace(all@meta.data$group, "A","SF")
all@meta.data$group <- str_replace(all@meta.data$group, "B","ECM")
all@meta.data$group <- factor(all@meta.data$group,  levels=c("SF","ECM"), ordered=T)

VlnPlot(all,c("ATP6"),
        group.by = "group",split.by = "Cell_type",
        cols = alphabet2)
ggsave(paste0("prettyumap/mitognes_vlnplot_ATP6.", Sys.Date(), ".pdf"), height = 6, width = 8)

VlnPlot(all,c("ND3"),
        group.by = "group",split.by = "Cell_type",
        cols = alphabet2)
ggsave(paste0("prettyumap/mitognes_vlnplot_ND3.", Sys.Date(), ".pdf"), height = 6, width = 8)

VlnPlot(all,c("CYTB"),
        group.by = "group",split.by = "Cell_type",
        cols = alphabet2)
ggsave(paste0("prettyumap/mitognes_vlnplot_CYTB.", Sys.Date(), ".pdf"), height = 6, width = 8)

VlnPlot(all,c("ND3"),
        group.by = "group",split.by = "Cell_type",
        cols = alphabet2)
ggsave(paste0("prettyumap/mitognes_vlnplot_ND3.", Sys.Date(), ".pdf"), height = 6, width = 8)

VlnPlot(all,c("COX2"),
        group.by = "group",split.by = "Cell_type",
        cols = alphabet2)
ggsave(paste0("prettyumap/mitognes_vlnplot_COX2.", Sys.Date(), ".pdf"), height = 6, width = 8)


DotPlot(all,features = c("ATP6","ND3","CYTB","COX3","ND4","COX2","ND5"),
        assay = "RNA",cols = c("#AFC2D9","#C6307C"),split.by = "group")



ggsave(paste0("prettyumap/mitognes_Dotplot.", Sys.Date(), ".pdf"), height = 6, width = 9)
