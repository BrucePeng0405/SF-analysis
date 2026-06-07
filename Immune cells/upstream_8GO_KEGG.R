rm(list = ls())
#library(xlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

#FDR<0.05 & |FC|>2
#load("all.markers.mast_All6.rda")
all.markers <- read.csv("mast/mast_df.csv")
all.markers.sig <- all.markers[which(all.markers$p_val_adj<0.05),]
all.markers.sig <- all.markers.sig[which(all.markers.sig$avg_log2FC < 0),]
ord <- data.frame(celltype = c("B cells","Carbohydrate metabolic Macrophages","CCL24 Macrophages","Chemotaxis Macrophages",
                               "Cycle Macrophages","DCs","Neutrophils","T cells"))
for (x in 0:7) {
  print(x)
  cell <- ord[x+1,1]
  print(cell)

genelist <- all.markers.sig[which(all.markers.sig[,8]==cell),2]
genelist <- as.list(genelist)
ids <-  bitr(genelist, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ego <- enrichGO(gene = ids$ENTREZID,
                keyType = "ENTREZID",
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = T)

write.csv(ego, file = paste0("GO_KEGG/", cell, "-SIS-UP-GO.csv"))
barplot(ego, showCategory=10, font.size = 14)
ggsave(filename = paste0("GO_KEGG/", cell, "-SIS-UP-GO.pdf"), width = 8, height = 7)


#  kk <- enrichKEGG(gene = ids$ENTREZID,
#                 organism = 'hsa',
 #                pvalueCutoff = 0.05)

#  write.csv(kk,file = paste0("GO_KEGG/", cell, "-KEGG.csv"))
#  dotplot(kk,showCategory=10, font.size = 14)
#  ggsave(filename = paste0("GO_KEGG/", cell, "-KEGG.pdf"), width = 8, height = 7)
}








Fibro2_go <- read.table("./GO_KEGG/david/SIS-Fibro2-GO.txt",sep = "\t",header = T)
Fibro2_go <- Fibro2_go[1:10,]

pdf("./GO_KEGG/david/SIS-Fibro2-GO.pdf",width = 8, height = 7)
ggplot(data=Fibro2_go,aes(x=Term,y=Count, fill = PValue))+
  geom_bar(stat = "identity", width = 0.8)+
  coord_flip()+
  xlab("GO Terms")+
  ylab("Gene Numbers")+
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()
dev.off()

Fibro2_KEGG <- read.table("./GO_KEGG/david/SIS-Fibro2-KEGG.txt",sep = "\t",header = T)
Fibro2_KEGG <- Fibro2_KEGG[1:10,]
pdf("./GO_KEGG/david/SIS-Fibro2-KEGG.pdf",width = 8, height = 7)
ggplot(data=Fibro2_KEGG,aes(x=Term,y=Count, fill = PValue))+
  geom_bar(stat = "identity", width = 0.8)+
  coord_flip()+
  xlab("GO Terms")+
  ylab("Gene Numbers")+
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()
dev.off()

