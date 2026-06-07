# Make pretty tSNE's for all 6 combined 10x replicates.
# Figures 1ab and SupFig 1c
# Matthew Buckley
rm(list = ls())
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(dplyr,tibble)

# Load Data
#setwd("~/10X/Seurat/")
load("muscle_regeneration_4markers_annotated_0.3.Rdata")
#d <- read.csv(file = "F:/scrna_seq_raw_1/10X/Seurat/all_All6_Filtered_meta.data&umap2022-11-18.csv", header=T)
#d <- d[,-1]
#d[1:3,1:11]

# Randomize rows to mix up plotting order
#d <- d[sample(nrow(d)),]

# Plot parameters
a <-  0.7
s <- .8
sampleColors <- c("#f41f1f", "darkorange", "firebrick", "deepskyblue", "#ad42f4","slateblue")
groupColors <- c("firebrick", "deepskyblue")
ben_colors <- c("#03c03c","#0054b4","#966fd6","#aec6cf","#ffdf00",
	"#ffb347","#e5aa70","#db7093","#e8000d", "#555f6d", "brown","deepskyblue")
timeColors <- c("#f41f1f", "darkorange", "firebrick")


### Cell type
CELLS <- c("Myeloid-1","Myeloid-2","VSMCs","Fibro-myofibro-1",
           "Myeloid-3","Fibro-myofibro-2","Fibro-myofibro-3","Fibro-myofibro-4",
           "Fibro-myofibro-5","Myeloid-4","Endothelial cells","T cells",
           "Monocytes","DCs","B cells","Satellite cells",
           "Lymphatic Endothelial cells","Germinal center B cell")
#all@meta.data$Cell_type <- factor(all@meta.data$Cell_type,  levels=CELLS, ordered=T)
DimPlot(object = all, reduction = "umap", cols = "alphabet2", 
        group.by = "Cell_type",label = F)
ggsave(paste0("prettyumap/customeUMAP_celltype_axis_nolabel", Sys.Date(), ".pdf"), height = 7, width = 10)

##三种细胞类型
Three_CELLS <- c("Immune cells","Immune cells","Peri-vascular cells","Fibroblast lineage cells",
           "Immune cells","Fibroblast lineage cells","Fibroblast lineage cells","Fibroblast lineage cells",
           "Fibroblast lineage cells","Immune cells","Peri-vascular cells","Immune cells",
           "Immune cells","Immune cells","Immune cells","Fibroblast lineage cells",
           "Peri-vascular cells","Immune cells")
names(Three_CELLS) <- CELLS
all <- RenameIdents(all, Three_CELLS)
DimPlot(object = all, reduction = "umap", cols = c("#d5151e","#00008e","#f19149"),label = T, label.size = 7)
ggsave(paste0("prettyumap/customeUMAP_3celltype_axis_nolegend", Sys.Date(), ".pdf"), height = 7, width = 10)


###  DiscretePalette(18, palette = "alphabet2", shuffle = FALSE)
#[1] "#AA0DFE" "#3283FE" "#85660D" "#782AB6" "#565656" "#1C8356" "#16FF32" "#F7E1A0" "#E2E2E2"
#[10] "#1CBE4F" "#C4451C" "#DEA0FD" "#FE00FA" "#325A9B" "#FEAF16" "#F8A19F" "#90AD1C" "#F6222E"
alphabet2 <- c("#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356", "#16FF32",
               "#F7E1A0","#E2E2E2","#1CBE4F", "#C4451C", "#DEA0FD", "#FE00FA", "#325A9B",
               "#FEAF16", "#F8A19F", "#90AD1C","#F6222E")

### Sample
all$orig.ident = factor(all$orig.ident, levels=c('A_7D','A_14D','A_1M','A_2M','B_7D','B_14D','B_1M','B_2M'))

DimPlot(object = all, reduction = "umap", cols = "alphabet2", 
        group.by = "orig.ident",label = F)
ggsave(paste0("prettyumap/customeUMAP_LEG_sample_size.8_", Sys.Date(), ".pdf"), height = 7, width = 8)
DimPlot(object = all, reduction = "umap", cols = "alphabet2",
        group.by = "Cell_type",label = F, split.by = "orig.ident", ncol = 4)
ggsave(paste0("prettyumap/customeUMAP_LEG_splite_sample_", Sys.Date(), ".pdf"), height = 7*2, width = 8*4)


### time
DimPlot(object = all, reduction = "umap", 
             group.by = "timepoint",label = F)
ggsave(paste0("prettyumap/customeUMAP_LEG_time", Sys.Date(), ".pdf"), height = 7, width = 8)

### Group
DimPlot(object = all, reduction = "umap", 
        group.by = "group",label = F)
ggsave(paste0("prettyumap/customeUMAP_group_axis_legend", Sys.Date(), ".pdf"), height = 7, width = 8)

DimPlot(object = all, reduction = "umap", 
        group.by = "group",label = F, split.by = "group")
ggsave(paste0("prettyumap/customeUMAP_group_axis_legend", Sys.Date(), ".pdf"), height = 7, width = 16)


#================================画堆积柱状图================================
cellnum <- table(all$Cell_type, all$orig.ident)
cell.prop<-as.data.frame(prop.table(cellnum,2))
colnames(cell.prop)<-c("Celltype","Group","Proportion")

cell.prop$Group = factor(cell.prop$Group, levels=c('A_7D','B_7D','A_14D','B_14D','A_1M','B_1M','A_2M','B_2M'))
p.bar <- ggplot(cell.prop,aes(Group,Proportion,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=alphabet2)+#自定义fill的颜色
  ggtitle("Cell proportion")+
  theme_test()+
  theme(axis.ticks.length=unit(0.4,'cm'))+
  guides(fill=guide_legend(title=NULL))
p.bar
#ggsave(paste0("prettyumap/cellsummary_group_axis_legend", Sys.Date(), ".pdf"), height = 7, width = 10)
ggsave("prettyumap/cellsummary.pdf",height = 4.91,width = 6.65)
dev.off()

