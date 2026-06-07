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
load("muscle_regeneration_4markers_0.3.Rdata")
#all <- subset(all, downsample = 2000)
#gc()
#all1 <- subset(all, downsample = 2000)
#count <- as.data.frame(all1[["RNA"]]@counts)
#write.csv(count, "Countmatrix.csv")
#rm(count,all1)
#gc()
FeaturePlot(all,reduction = "umap", 
            features = c("PTPRC","ACTA2","COL1A1","AIF1",
                         "CD68","SPARCL1","VCAM1","CD34","MYF5","CD74",
                         "CD14","PECAM1","MYH11"))


#immune cells
FeaturePlot(all,reduction = "umap", 
            features = c("PTPRC")) 
ggsave("annotation/Immune.cells.pdf")
VlnPlot(all,c("PTPRC"))
ggsave("annotation/Immune.cells_v.pdf")
#Myeloid cells
FeaturePlot(all,reduction = "umap", 
            features = c("CSF1R","AIF1"))
ggsave("annotation/Myeloid.cells.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("CSF1R","AIF1"))
ggsave("annotation/Myeloid.cells_v.pdf",height = 4.91,width = 6.65*2)
#Macrophages
FeaturePlot(all,reduction = "umap", 
            features = c("CD68","IGSF6","IL1B","MRC1","SPP1")) 
ggsave("annotation/Macrophages.pdf",height = 4.91*3,width = 6.65*2)
VlnPlot(all,c("CD68","IGSF6","IL1B","MRC1","SPP1"))
ggsave("annotation/Macrophages_v.pdf",height = 4.91*2,width = 6.65*3)
#MONOcytes
FeaturePlot(all,reduction = "umap", 
            features = c("CD14")) 
ggsave("annotation/Monocytes.pdf")
VlnPlot(all,c("CD14"))
ggsave("annotation/Monocytes_v.pdf")
#DC
FeaturePlot(all,reduction = "umap", 
            features = c("CD83","CD74")) 
ggsave("annotation/DC.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("CD83","CD74","HLA-DRA"))
ggsave("annotation/DC_v.pdf",height = 4.91,width = 6.65*2)
#T cells 
FeaturePlot(all,reduction = "umap", 
            features = c("GATA3","CD2")) 
ggsave("annotation/T.cells.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("GATA3","CD2"))
ggsave("annotation/T.cells_v.pdf",height = 4.91,width = 6.65*2)
#B cells
FeaturePlot(all,reduction = "umap", 
            features = c("CD79B","MS4A1"))
ggsave("annotation/B.cells.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("CD79B","MS4A1"))
ggsave("annotation/B.cells_v.pdf",height = 4.91,width = 6.65*2)
#endothelial cells
FeaturePlot(all,reduction = "umap", 
            features = c("PECAM1","ITGA6","VWF","CCL21","PROX1")) 
ggsave("annotation/Endothelial.cells.pdf",height = 4.91*3,width = 6.65*2)
VlnPlot(all,c("PECAM1","ITGA6","VWF","CCL21","PROX1"))
ggsave("annotation/Endothelial.cells_v.pdf",height = 4.91*2,width = 6.65*3)
#Vascular Smooth muscle cells
FeaturePlot(all,reduction = "umap", 
            features = c("ACTA2","MYH11")) 
ggsave("annotation/Vascular.smooth.muscle.cells.pdf",height = 4.91,width = 6.65*2)
VlnPlot(all,c("ACTA2","MYH11"))
ggsave("annotation/Vascular.smooth.muscle.cells_v.pdf",height = 4.91,width = 6.65*2)
# FAP cells
FeaturePlot(all,reduction = "umap", 
            features = c("LUM","DCN","COL1A2","PDGFRL")) 
ggsave("annotation/Fibroblast.myofibroblast.pdf",height = 4.91*2,width = 6.65*2)
VlnPlot(all,c("LUM","DCN","COL1A2","PDGFRL"))
ggsave("annotation/Fibroblast.myofibroblast_v.pdf",height = 4.91*2,width = 6.65*3)
# Muscle satellite cells
FeaturePlot(all,reduction = "umap", 
            features = c("SPARCL1","CADM1","SPATS2L","MYF5"))
ggsave("annotation/Muscle.satellite.cells.pdf",height = 4.91*2,width = 6.65*2)
VlnPlot(all,c("SPARCL1","CADM1","SPATS2L","MYF5"))
ggsave("annotation/Muscle.satellite.cells_v.pdf",height = 4.91*2,width = 6.65*3)

# ALL marker
FeaturePlot(all,reduction = "umap", 
            features = c("PTPRC","AIF1","CD14","CD68",
                         "CD74","CD79B","GATA3",
                         "COL1A2","ACTA2","PECAM1","SPATS2L"))
ggsave("annotation/ALL.cells.pdf",height = 4.91*3,width = 6.65*4)
VlnPlot(all,c("PTPRC","AIF1","CD14","CD68",
              "CD74","CD79B","GATA3",
              "COL1A2","ACTA2","PECAM1","SPATS2L"))
ggsave("annotation/ALL.cells_v.pdf",height = 4.91*3,width = 6.65*4)



cluster_name <- c("Myeloid-1","Myeloid-2","VSMCs","Fibro-myofibro-1",
                  "Myeloid-3","Fibro-myofibro-2","Fibro-myofibro-3","Fibro-myofibro-4",
                  "Fibro-myofibro-5","Myeloid-4","Endothelial cells","T cells",
                  "Monocytes","DCs","B cells","Satellite cells",
                  "Lymphatic Endothelial cells","Germinal center B cell")

names(cluster_name) <- c(0:17)
all <- RenameIdents(all, cluster_name)
DimPlot(all, reduction = "umap", label = TRUE, pt.size = 0.5)+NoLegend()
ggsave("annotation/umap.clus.rename.pdf",height = 4.91,width = 6.65)
DimPlot(all, reduction = "umap",label = F, split.by = "orig.ident", ncol = 4)
ggsave("annotation/umap.split.group.pdf",width = 6.65*4,height = 4.91*2)

#Cell type注释
Cell_type <- data.frame(ClusterID=c(0:17), Cell_type=cluster_name, stringsAsFactors = F)
all@meta.data$Cell_type <- "NA"
for (i in 1:nrow(Cell_type)) {
  all@meta.data[which(all@meta.data$seurat_clusters == Cell_type$ClusterID[i]),"Cell_type"] <- Cell_type$Cell_type[i]
}
all@meta.data$Cell_type
DimPlot(all,group.by = "Cell_type",label = F,reduction = "umap")
ggsave("annotation/umap.clus.celltype.pdf",height = 4.91,width = 6.65)
save(all, file = "muscle_regeneration_4markers_annotated_0.3.Rdata")

# FAP <- subset(all,idents=c("FAP-1","FAP-2","FAP-3","FAP-4","FAP-5"))
# FAP <- RunPCA(FAP, features = VariableFeatures(object = FAP))
# VizDimLoadings(FAP, dims = 1:2, reduction = "pca")
# DimPlot(FAP, reduction = "pca")
# ElbowPlot(FAP,ndims = 50 )
# FAP.2 <- FindNeighbors(FAP, dims = 1:8,k.param = 20)
# FAP.2 <- FindClusters(FAP.2, resolution = 0.4)
# FAP.2 <- RunTSNE(FAP.2, dims = 1:8)
# DimPlot(FAP.2, reduction = "tsne",label = T)
# #M2
# VlnPlot(FAP.2,c("TCF4","MYL9","ACTA2","VCAM1","SERPING1","CDK1","SPARCL1"))
# #M1
# VlnPlot(macrophage.2,c("Il1b","Ccr7","Tnf","Cd80","RT1-Da"))
# ##Initiation/Recruitment of endothelial cells
# VlnPlot(macrophage.2,c("Vegfa","Fgf2","Il8","Ccl5"))
# ##Maturation/Recruitment of endothelial cells
# VlnPlot(macrophage.2,c("Pdgfb","Hbegf")) #M2a #M1
# VlnPlot(macrophage.2,c("Mmp9","Timp3")) #M1,M2c, #M2a

# VlnPlot(macrophage.2,c("Cd163","Socs3","Nos2","Cxcl10","Il10","Ccl2","Ccl17","Ccl22","Ccl24","Ccl9"),pt.size = 0)
# FeaturePlot(macrophage,"Cd163")
# ggsave("Marcophage_umap.clus.pdf")
# DimPlot(macrophage, reduction = "umap")
# ggsave("Marcophage_umap.group.pdf")
# macrophage <- RunTSNE(macrophage, dims = 1:6)
# DimPlot(macrophage, reduction = "tsne",label = T)
# ggsave("Marcophage_tsne.clus.pdf")
# DimPlot(macrophage, reduction = "tsne",group.by = "orig.ident")
# ggsave("Marcophage_tsne.group.pdf")
# save(macrophage,file="macrophage.Rdata")
# FAP.2.markers <- FindAllMarkers(FAP.2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# FAP.2_markers_top10 <- FAP.2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# write.csv(FAP.2_markers_top10,file = "FAP.2_markers_top10.csv")
# save(FAP.2_markers_top10,file="FAP.2.allmarkers.Rdata")
# #cluster_name <- c("M1","M1","M2","M1","M1","M2","M2","Unclassified")
# VlnPlot(macrophage.2,c("Cd163"),pt.size = 0)
# FeaturePlot(macrophage.2,c("Ccl18","Ccl22","Mrc1","Pdgfb","Timp3"),pt.size = 0)
# pdf("Macrophage.heatmap.orig.pdf")
# DoHeatmap(macrophage,features = macrophage.markers$gene)
# dev.off()
# pdf("macrophageProp.orig.pdf")
# #Plotcellprop(cell.prop,"macro")
# PlotPiefacet(cell.prop,"macro")
# dev.off()
# cluster_name <- c("M2a","M1","M2c","M1","M2a","M2c","M2c","Unclassified")
# names(cluster_name) <- c(0:7)
# macrophage.2 <- RenameIdents(macrophage.2, cluster_name)
# # DimPlot(macrophage, reduction = "umap", label = TRUE, pt.size = 0.5)+NoLegend()
# # ggsave("marcophage_M1M2.pdf")
# # m1m2.markers <- FindAllMarkers(macrophage, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# # m1m2_markers_top10 <- m1m2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# # write.csv(m1m2_markers_top10,file="m1m2.markers_top10.csv")
# macro6v5 <- FindMarkers(macrophage.2,ident.1 = "6",idents.2="5",only.pos = F,min.pct = 0.25, logfc.threshold = 0.25)
# #myGOmarker(all.markers,"all_clus15")
# t <- macrophage.markers %>% filter(avg_logFC>0.5)
# table(t$cluster)
# Exportmarkers4GO(t,"macrophage.markers4GO.csv")
# DoHeatmap(macrophage ,features = t$gene,angle = 45,label=F,draw.lines =T)
# 
# myGOmarker(macrophage.markers,"macrophage")
# #proportion among group
# table(Idents(macrophage))
# table(macrophage$group)
# cell.prop<-as.data.frame(table(Idents(macrophage.2), macrophage.2$group))
# colnames(cell.prop)<-c("cluster","group","proportion")
# 
# pdf("macrophageProp.pdf")
# Plotcellprop(cell.prop,"macro")
# PlotPiefacet(cell.prop,"macro")
# dev.off()
# save(macrophage,file="Macrophage.Rdata")
# DEmarkers <- FindMarkers(macrophage,only.pos = F,assay = "RNA",slot = "data",ident.1=rownames(subset(macrophage@meta.data,group=="SilkFibroin")),ident.2 = rownames(subset(macrophage@meta.data,group=="Purilon")))
# myGO(rownames(DEmarkers[DEmarkers$avg_logFC>0,]))
# ggsave("DE.Macrophge.SFvPurilon.pdf",width=10)
# 
# Fibroblast <- subset(all,idents=c("0","13"))
# Fibroblast$all.ident <- Fibroblast@active.ident
# Fibroblast <- RunPCA(Fibroblast, features = VariableFeatures(object = Fibroblast))
# #VizDimLoadings(Fibroblast, dims = 1:2, reduction = "pca")
# DimPlot(Fibroblast, reduction = "pca")
# ElbowPlot(Fibroblast,ndims = 30 )
# Fibroblast <- FindNeighbors(Fibroblast, dims = 1:12,k.param = 20)
# Fibroblast <- FindClusters(Fibroblast, resolution = 0.1)
# Fibroblast <- RunUMAP(Fibroblast, dims = 1:12)
# FeaturePlot(Fibroblast,c("Blimp1", "Itga8", "Cd26", "Lrig1", "EphB2", "Trps1","Crabp1", "AlkPhos"))
# DimPlot(Fibroblast, reduction = "umap",label = T)
# DimPlot(Fibroblast,group.by = "all.ident")
# ggsave("Fibroblast_umap.clus.pdf")
# DimPlot(Fibroblast, reduction = "umap",group.by = "orig.ident")
# ggsave("Fibroblast_umap.group.pdf")
# Fibroblast <- RunTSNE(Fibroblast, dims = 1:12)
# DimPlot(Fibroblast, reduction = "tsne",label = T)
# ggsave("Fibroblast_tsne.clus.pdf")
# DimPlot(Fibroblast, reduction = "tsne",group.by = "orig.ident")
# ggsave("Fibroblast_tsne.group.pdf")
# save(Fibroblast,file="Fibroblast.Rdata")
# table(Idents(Fibroblast))
# cell.prop<-as.data.frame(prop.table(table(Idents(Fibroblast), Fibroblast$group)))
# colnames(cell.prop)<-c("cluster","group","proportion")
# pdf("FibroblastProp_anno.pdf")
# Plotcellprop(cell.prop,"Fibroblast")
# PlotPiefacet(cell.prop,"Fibroblast")
# dev.off()
# Fibroblast.markers <- FindAllMarkers(Fibroblast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Fibroblast_markers_top10 <- Fibroblast.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# write.csv(Fibroblast_markers_top10,file = "Fibroblast_clus4_markers10.csv")
# write.csv(Fibroblast.markers,file="Fibroblast_clus4.allmarkers.csv")
# 
# t <- Fibroblast.markers %>% filter(avg_logFC>0.4)
# table(t$cluster)
# myGOmarker(Fibroblast.markers,"Fibroblast")
# 
# Exportmarkers4GO(t,"fibroblast_markers4GO.csv")
# #对成纤维细胞的功能分析
# #ECM分泌，炎症因子分泌，生长因子分泌
# Fibroblast <- subset(all,idents=c("0","13"))
# DEmarkers <- FindMarkers(Fibroblast,only.pos = F,assay = "RNA",slot = "data",ident.1=rownames(subset(Fibroblast@meta.data,group=="SilkFibroin")),ident.2 = rownames(subset(Fibroblast@meta.data,group=="Purilon")))
# DEmarkers2  <- DEmarkers[DEmarkers$avg_logFC>0,]
# myGO(rownames(DEmarkers2))
# ggsave("DE.fibroblast.SFvPurilonGO_neg.pdf",width = 8)
# write.csv(DEmarkers ,file="DE.fibroblast.SFvPurilon.csv")
# DEclus <-data.frame(class="ECM&Reconstruction-related",gene=c("Pfn1","Myl6","Tgfbi","Vim","Serpinh1","Col12a1","Fbln2"))
# DEclus <- rbind(DEclus,data.frame(class="GrowthFactor-related",gene=c("Vim","Ssr2","Sec61b","Spp1","Krt17","Fbln2","Selenom","Ccn3")))
# DEclus <- rbind(DEclus,data.frame(class="Inflammation-related",gene=c("Tuba1b","S100a9","Fcer1g","Crip1","S100a8","Cd74","Tyrobp","RT1-Da","Trem2","Cd68","RT1-Ba","RT1-Db1","Il1rn","Il1b","Cd14","Mmp9","Pf4","Cxcl2","Lyz2","Aif1","Ccl3","Ddit4","Ccl4","Acp5")))
# DEclus$gene <- as.character(DEclus$gene)
# DoHeatmap(Fibroblast,features = DEclus$gene,cells = c(rownames(subset(Fibroblast@meta.data,group=="SilkFibroin")), rownames(subset(Fibroblast@meta.data,group=="Purilon"))),group.by = "orig.ident")
# 
# 
# endo <- subset(all,idents=c("Ramp2+VEC","Il6+VEC"))
# endo <- RunPCA(endo, features = VariableFeatures(object = endo))
# #VizDimLoadings(endo, dims = 1:2, reduction = "pca")
# DimPlot(endo, reduction = "pca")
# ElbowPlot(endo,ndims = 50 )
# endo <- FindNeighbors(endo, dims = 1:8,k.param = 10)
# endo <- FindClusters(endo, resolution = 0.1)
# endo <- RunUMAP(endo, dims = 1:8)
# DimPlot(endo, reduction = "umap",label = T)
# DimPlot(endo, reduction = "umap",group.by = "orig.ident")
# save(endo,file="endo.Rdata")
# table(Idents(endo))
# table(endo$group)
# cell.prop<-as.data.frame(table(Idents(endo), endo$group))
# colnames(cell.prop)<-c("cluster","group","proportion")
# pdf("endoProp_anno.pdf")
# Plotcellprop(cell.prop,"endo")
# PlotPiefacet(cell.prop,"endo")
# dev.off()
# 
# VlnPlot(endo,c("Cd34","Cxcr4", "Ptprc", "Cdh5","Pecam1","Prom1","Kit","Kdr","Csf1r"),pt.size=0)
# FeaturePlot(endo,c("Cd34","Cxcr4", "Cdh5","Kit","Kdr"),pt.size=0)
# endo.markers <- FindAllMarkers(endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# endo_markers_top10 <- endo.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# write.csv(endo_markers_top10,file = "endo_markers10.csv")
# save(endo.markers,file="endo.allmarkers.Rdata")
# 
# t <- endo.markers %>% filter(avg_logFC>0.4)
# table(t$cluster)
# myGO(t$gene[t$cluster==4])
# myGOmarker(t,"endo")
# DEmarkers <- FindMarkers(endo,only.pos = F,ident.1=rownames(subset(endo@meta.data,group=="SilkFibroin")),ident.2 = rownames(subset(endo@meta.data,group=="Purilon")))
# diff_GO <- myGO(rownames(DEmarkers))
# angiogenesis <- unlist(strsplit(diff_GO[["data"]][grep("angio",diff_GO[["data"]]$Description),"geneID"],"/"))
# t <- endo.markers[which(endo.markers$gene %in% angiogenesis),]
# write.csv(DEmarkers ,file="DE.endo.SFvPurilon.csv")
# DEmarkers <- FindMarkers(endo,only.pos = F,ident.1=rownames(subset(endo@meta.data,group=="SilkFibroin")),ident.2 = rownames(subset(endo@meta.data,group=="Control")))
# angiogenesis <- diff_GO@geneSets[["GO:0001525"]]
# angiogenesis.vControl <- DEmarkers[rownames(DEmarkers) %in% angiogenesis,]
# t <- endo.markers[which(endo.markers$gene %in% angiogenesis),]
# write.csv(DEmarkers ,file="DE.endo.SFvPurilon.csv")
# 
# 
# 
# epi <- subset(all,idents=c("Stfa3+Epithelial","Krt79+Epithelial","Krt15+Epithelial","Top2a+Epithelial"))
# epi <- RunPCA(epi, features = VariableFeatures(object = epi))
# #VizDimLoadings(epi, dims = 1:2, reduction = "pca")
# DimPlot(epi, reduction = "pca")
# ElbowPlot(epi,ndims = 20 )
# epi <- FindNeighbors(epi, dims = 1:10)
# epi <- FindClusters(epi, resolution = 0.)
# epi <- RunUMAP(epi, dims = 1:10)
# DimPlot(epi, reduction = "umap",label = T)
# ggsave("epi_umap.clus.pdf")
# DimPlot(epi, reduction = "umap",group.by = "orig.ident")
# ggsave("epi_umap.group.pdf")
# epi <- RunTSNE(epi, dims = 1:10)
# DimPlot(epi, reduction = "tsne",label = T)
# ggsave("epi_tsne.clus.pdf")
# DimPlot(epi, reduction = "tsne",group.by = "orig.ident")
# ggsave("epi_tsne.group.pdf")
# t <- as.data.frame(Idents(all))
# epi$all.ident <- t[colnames(epi),]
# save(epi,file="epi.Rdata")
# 
# table(Idents(epi))
# table(epi$group)
# cell.prop<-as.data.frame(prop.table(table(Idents(epi), epi$group)))
# colnames(cell.prop)<-c("cluster","group","proportion")
# pdf("epiProp_anno.clus7.pdf")
# Plotcellprop(cell.prop,"epi")
# PlotPie(cell.prop,"epi")
# dev.off()
# 
# pdf("epi.heatmap.pdf")
# DoHeatmap(epi)
# dev.off()
# 
# epi.markers.clu7<- FindAllMarkers(epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# epi_markers_top10 <- epi.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# write.csv(epi_markers_top10,file = "epi_clus5_markers10.csv")
# write.csv(epi.markers,file = "epi_clus5.allmarkers.csv")
# save(epi.markers,file="epi_clus5.allmarkers.Rdata")
# VlnPlot(epi,c("Krt71","Krt25","Krt73"),pt.size = 0) #是inner root sheath
# FeaturePlot(epi,c("Krt71","Krt25","Krt73"),pt.size = 0) #是inner root sheath
# FeaturePlot(epi,c("Cd200","Cd34","Itga6","Krt15","Lhx2","Sox9"),pt.size = 0)
# ggsave("bulge stem cells_Vln_marker.pdf",width = 9)
# #认为第2群是Basal cells 基底细胞，表达分化和转录因子
# #ref https://discovery.lifemapsc.com/in-vivo-development/hair/bulge/bulge-stem-cells
# #myGOmarker(epi.markers,"epi")
# 
# t <- epi.markers.clu7 %>% filter(avg_logFC>0.4)
# table(t$cluster)
# myGO(t[t$cluster==4,"gene"])
# myGOmarker(t,"epi.clus7")
# epi_markers4GO <- Exportmarkers4GO(t,"epi_markers4GO.csv")
# DoHeatmap(epi,features = t$gene,angle = 45,label=F,draw.lines =T)
# 
# 
# DC <- subset(all,idents="10")
# DC <- RunPCA(DC, features = VariableFeatures(object = DC))
# #VizDimLoadings(DC, dims = 1:2, reduction = "pca")
# DimPlot(DC, reduction = "pca")
# ElbowPlot(DC,ndims = 20 )
# DC <- FindNeighbors(DC, dims = 1:5,k.param = 5)
# DC <- FinDClusters(DC, resolution = 0.1)
# DC <- RunUMAP(DC, dims = 1:5)
# DimPlot(DC, reduction = "umap",label = T)
# ggsave("DC_umap.clus.pdf")
# DimPlot(DC, reduction = "umap",group.by = "orig.ident")
# ggsave("DC_umap.group.pdf")
# DC <- RunTSNE(DC, dims = 1:6)
# DimPlot(DC, reduction = "tsne",label = T)
# ggsave("DC_tsne.clus.pdf")
# DimPlot(DC, reduction = "tsne",group.by = "orig.ident")
# ggsave("DC_tsne.group.pdf")
# save(DC,file="DC.Rdata")
# FeaturePlot(DC,"Cd207") #Langerhans cells (immature dendritic cells)
# table(Idents(DC))
# table(DC$group)
# cell.prop<-as.data.frame(prop.table(table(Idents(DC), DC$group)))
# colnames(cell.prop)<-c("cluster","group","proportion")
# pdf("DCProp_anno.pdf")
# Plotcellprop(cell.prop,"DC")
# PlotPie(cell.prop,"DC")
# dev.off()
# DC.markers <- FindAllMarkers(DC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# DC_markers_top10 <- DC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# write.csv(DC_markers_top10,file = "DC_clus5_markers10.csv")
# write.csv(DC.markers,file = "DC_clus5.allmarkers.csv")
# save(DC.markers,file="DC_clus5.allmarkers.Rdata")
# t <- DC.markers %>% filter(avg_logFC>0.4)
# table(t$cluster)
# myGOmarker(t,"DC")
# DC_markers4GO <- Exportmarkers4GO(t,"DC_markers4GO.csv")
# DoHeatmap(DC,features = t$gene,angle = 45,label=F,draw.lines =T)
# 
# #neutrophil
# setwd("./ploting/")
# dir.create("Neutrophil")
# setwd("./Neutrophil/")
# Neutrophil <- subset(all,idents="Neutrophil")
# Neutrophil <- RunPCA(Neutrophil, features = VariableFeatures(object = Neutrophil))
# #VizDimLoadings(Neutrophil, dims = 1:2, reduction = "pca")
# DimPlot(Neutrophil, reduction = "pca")
# ElbowPlot(Neutrophil,ndims = 50 )
# Neutrophil <- FindNeighbors(Neutrophil, dims = 1:5,k.param = 20)
# Neutrophil <- FindClusters(Neutrophil, resolution = 0.2)
# Neutrophil <- RunUMAP(Neutrophil, dims = 1:5)
# Neutrophil <- RunTSNE(Neutrophil, dims = 1:5)
# DimPlot(Neutrophil, reduction = "umap", label = TRUE,repel = T, pt.size = 0.5)+NoLegend()
# ggsave("Neutrophil.umap.clus.pdf",width=5)
# DimPlot(Neutrophil, reduction = "umap",group.by = "orig.ident", pt.size = 0.5,cols=palette_all)
# ggsave("Neutrophil.umap.group.pdf",width=6.5)
# DimPlot(Neutrophil, reduction = "tsne",label = T,repel = T, pt.size = 0.5)+NoLegend()
# ggsave("Neutrophil.tsne.clus.pdf",width=5)
# DimPlot(Neutrophil, reduction = "tsne",group.by = "orig.ident", pt.size = 0.5,cols=palette_all)
# ggsave("Neutrophil.tsne.group.pdf",width=6.5)
# save(Neutrophil,file="../../data/Neutrophil.Rdata")
# FeaturePlot(Neutrophil,c("Vegfr1","Vegfr2","Vegfa", "Cxcr4"))
# table(Idents(Neutrophil))
# table(Neutrophil$group)
# cell.prop<-as.data.frame((table(Idents(Neutrophil), Neutrophil$group)))
# colnames(cell.prop)<-c("cluster","group","proportion")
# pdf("NeutrophilProp_anno.pdf")
# PlotPiefacet(cell.prop,"Neutrophil")
# dev.off()
# Neutrophil.markers <- FindAllMarkers(Neutrophil, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Neutrophil.markers.top10 <- Neutrophil.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# write.csv(Neutrophil.markers.top10,file = "Neutrophil.markers10.csv")
# write.csv(Neutrophil.markers,file = "Neutrophil.allmarkers.csv")
# t <- Neutrophil.markers %>% filter(avg_logFC>0.4)
# table(t$cluster)
# myGOmarker(t,"Neutrophil")
# DoHeatmap(Neutrophil,features = t$gene,angle = 45,label=F,draw.lines =T)
# 
# {
#   rm(macrophage)
#   rm(endo)
#   rm(epi)
#   rm(all)
#   rm(Neutrophil)
#   rm(Fibroblast)
#   rm(diff_GO)
# }
# 
# 
# 
# 
# 
