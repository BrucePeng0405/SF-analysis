rm(list = ls())
#library(xlsx)
library(readxl)
library(clusterProfiler)
#library(dbplyr)
library(DOSE)


Adhesion_Myo_SF <- read.csv("GO_KEGG/Adhesion Myofibroblast-SF-UP-GO.csv")
Adhesion_Myo_SF <- Adhesion_Myo_SF[which(Adhesion_Myo_SF$pvalue<0.05),]
Adhesion_Myo_SF$Group <- c(rep("Adhesion Myofibroblasts",153))

Adhesion_Fibro_SF <- read.csv("GO_KEGG/Adhesion Fibroblast-SF-UP-GO.csv")
Adhesion_Fibro_SF <- Adhesion_Fibro_SF[which(Adhesion_Fibro_SF$pvalue<0.05),]
Adhesion_Fibro_SF$Group <- c(rep("Adhesion Fibroblasts",44))

Migra_Myo_SF <- read.csv("GO_KEGG/Migration Myofibroblast-SF-UP-GO.csv")
Migra_Myo_SF$Group <- c(rep("Migration Myofibroblasts",43))

Infla_Fibro_SF <- read.csv("GO_KEGG/Inflammatory Fibroblast-SF-UP-GO.csv")
Infla_Fibro_SF$Group <- c(rep("Inflammatory Fibroblasts",20))

Cycle_Myo_SF <- read.csv("GO_KEGG/Cycle Myofibroblast-SF-UP-GO.csv")
Cycle_Myo_SF$Group <- c(rep("Cycle Myofibroblasts",150))


kk <- rbind(Adhesion_Myo_SF,Adhesion_Fibro_SF,Migra_Myo_SF,Infla_Fibro_SF,Cycle_Myo_SF)


kk$GeneRatio <- parse_ratio(kk$GeneRatio)
a <- as.data.frame(table(kk$Description))
a <- a[order(a$Freq,decreasing = T),]
write.csv(a, "GO_KEGG/heatmap/SF-GO-frequncy.csv")

selected_pathway <- as.data.frame(a[1:20,])

colnames(selected_pathway) <- c("Description","Freq")
tmp <- merge(kk,selected_pathway,by.x = "Description", by.y = "Description")
tmp$Description <- factor(tmp$Description,levels = rev(selected_pathway$Description))
write.csv(tmp,"GO_KEGG/heatmap/SF-all-GO.csv")

library(stringr)
library(ggplot2)

#tiff(filename = "pp_overlap_GO-2.tiff",width = 430*2, height = 185)

p = ggplot(tmp,aes(Group,Description))
#p=p + geom_point(aes(size=GeneRatio))
p = p+ geom_point(aes(size=GeneRatio,color=pvalue))
p = p+scale_color_gradient(low="#C6307C",high = "#AFC2D9")
p = p+labs(color=expression(pvalue),size="GeneRatio",x = "", y="")
p = p+scale_y_discrete(labels=function(x) str_wrap(x, width=55))
p+theme(axis.text.y =element_text(size=12,face="bold",colour = "black"),
        axis.text.x =element_text(size=12,face="bold",colour = "black", angle = 15,hjust = 1),
        panel.background = element_rect(color = "black", size = 2,fill = "NA"),
        panel.grid = element_blank())
#dev.off()
ggsave("GO_KEGG/heatmap/SF-UP-1.pdf",width = 12, height = 6.5)

##################SIS#######################################
rm(list = ls())
#library(xlsx)
library(readxl)
library(clusterProfiler)
#library(dbplyr)
library(DOSE)


Adhesion_Myo_SF <- read.csv("GO_KEGG/Adhesion Myofibroblast-SIS-UP-GO.csv")
Adhesion_Myo_SF <- Adhesion_Myo_SF[which(Adhesion_Myo_SF$pvalue<0.05),]
Adhesion_Myo_SF$Group <- c(rep("Adhesion Myofibroblasts",409))

Adhesion_Fibro_SF <- read.csv("GO_KEGG/Adhesion Fibroblast-SIS-UP-GO.csv")
Adhesion_Fibro_SF <- Adhesion_Fibro_SF[which(Adhesion_Fibro_SF$pvalue<0.05),]
Adhesion_Fibro_SF$Group <- c(rep("Adhesion Fibroblasts",161))

Migra_Myo_SF <- read.csv("GO_KEGG/Migration Myofibroblast-SIS-UP-GO.csv")
Migra_Myo_SF$Group <- c(rep("Migration Myofibroblasts",92))

Infla_Fibro_SF <- read.csv("GO_KEGG/Inflammatory Fibroblast-SIS-UP-GO.csv")
Infla_Fibro_SF$Group <- c(rep("Inflammatory Fibroblasts",256))

Cycle_Myo_SF <- read.csv("GO_KEGG/Cycle Myofibroblast-SIS-UP-GO.csv")
Cycle_Myo_SF$Group <- c(rep("Cycle Myofibroblasts",225))


kk <- rbind(Adhesion_Myo_SF,Adhesion_Fibro_SF,Migra_Myo_SF,Infla_Fibro_SF,Cycle_Myo_SF)


kk$GeneRatio <- parse_ratio(kk$GeneRatio)
a <- as.data.frame(table(kk$Description))
a <- a[order(a$Freq,decreasing = T),]
write.csv(a, "GO_KEGG/heatmap/SIS-GO-frequncy.csv")

selected_pathway <- as.data.frame(a[1:20,])

colnames(selected_pathway) <- c("Description","Freq")
tmp <- merge(kk,selected_pathway,by.x = "Description", by.y = "Description")
tmp$Description <- factor(tmp$Description,levels = rev(selected_pathway$Description))
write.csv(tmp,"GO_KEGG/heatmap/SIS-all-GO.csv")

library(stringr)
library(ggplot2)

#tiff(filename = "pp_overlap_GO-2.tiff",width = 430*2, height = 185)

p = ggplot(tmp,aes(Group,Description))
#p=p + geom_point(aes(size=GeneRatio))
p = p+ geom_point(aes(size=GeneRatio,color=pvalue))
p = p+scale_color_gradient(low="#C6307C",high = "#AFC2D9")
p = p+labs(color=expression(pvalue),size="GeneRatio",x = "", y="")
p = p+scale_y_discrete(labels=function(x) str_wrap(x, width=55))
p+theme(axis.text.y =element_text(size=12,face="bold",colour = "black"),
        axis.text.x =element_text(size=12,face="bold",colour = "black", angle = 15,hjust = 1),
        panel.background = element_rect(color = "black", size = 2,fill = "NA"),
        panel.grid = element_blank())
#dev.off()
ggsave("GO_KEGG/heatmap/SIS-UP-1.pdf",width = 12, height = 6.5)


