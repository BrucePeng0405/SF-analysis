rm(list = ls())
library(Seurat)
library(tibble)

load("immune_recluster_4markers_annotated_0.3.Rdata")
#=======================================================================================
# Age-associated Differential Expression for each cell type
#=======================================================================================

## Append age to celltype label and add to metadata
cellGroup <- paste0(all@meta.data$Cell_type, "_", all@meta.data$group)
m <- data.frame("Celltype_Group" = cellGroup)
rownames(m) <- rownames(all@meta.data)
all <- AddMetaData(object = all, metadata = m)

# Assign as main identity
Idents(all) <- "Celltype_Group"

CELLTYPES <- unique(all@meta.data$Cell_type)

mast_list <- vector(mode="list", length = 18)
names(mast_list) = CELLTYPES

# Compare young/old for each celltype, save each matrix in a list. 
for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  # Find  cluster marker genes
  all.mast.de <- FindMarkers(object = all,
                             ident.1 = paste0(CELLTYPE, "_A"),
                             ident.2 = paste0(CELLTYPE, "_B"),
                             only.pos = FALSE, 
                             min.pct = 0,
                             logfc.threshold = 0,
                             test.use = "MAST")
  all.mast.de <- as.data.frame(all.mast.de)
  all.mast.de <- rownames_to_column(all.mast.de, var = "gene")
  all.mast.de$celltype <- rep(CELLTYPE, dim(all.mast.de)[1])
  mast_list[[CELLTYPE]] <- all.mast.de
}

# Combine all matrices into one dataframe
mast_df <- data.frame()
for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  mast_df  <- rbind(mast_df, mast_list[[CELLTYPE]])
}
dim(mast_df); head(mast_df) # 200976 by 7

# Save differential expression results
write.csv(mast_df, "mast/mast_df.csv")
#save(all, file = "annotated_cluster.Rdata")

#=======================================================================================
sessionInfo()


#==================================================================================
#mast results pretty plots
#==================================================================================

library(tidyverse)

library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(scales)
library(NCmisc) #v1.1.5

# See sessionInfo() at end of script for version information.

# Modify path below. All subsequent paths are relative.
#setwd("~/Desktop/Dropbox/DBN2019/10x/mast")

#==================================================================================================
df <- mast_df
df <- read.csv("mast/mast_df.csv")

# Add z-score based on two sided null hypothesis.
df$z <- p.to.Z(df$p_val) * sign(df$avg_log2FC)
df$z.adj <- p.to.Z(df$p_val_adj) * sign(df$avg_log2FC)

# Randomize rows to reduce overplotting issues.
df <- df[sample(nrow(df)), ]

# Adjust factors
CELLS <- c("B cells","Carbohydrate metabolic Macrophages","CCL24 Macrophages","Chemotaxis Macrophages",
           "Monocytes","DCs","Neutrophils","T cells")
df$celltype_factor <- factor(df$celltype,  levels=CELLS, ordered=T)

#==================================================================================================

# Adjust colors  alphabet2
new_colors <- c("#4991c1", "#afc2d9", "#c6307c", "#435b95", "#79b99d", "#e69f84", 
                "#89558d","#d0afc4")
names(new_colors) <- levels(df$celltype_factor)

# Make custom color column to facilitate grey coloring by threshold.
col <- new_colors[df$celltype_factor]
col[df$p_val_adj > 0.05] <- "#D3D3D3" # grey
df$col <- as.factor(col)

q <- ggplot(df, aes(x = celltype_factor, y = z, color = col)) +
  geom_jitter(width = 0.40, alpha = .55, size = 1) +theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) +
  ggtitle("Changes in gene expression with group") +
  theme(axis.title.y = element_text(size = 20, face = "plain")) +
  theme(axis.text.y = element_text(size = 20,colour = "black")) +
  theme(axis.text.x = element_text(size = 11,colour = "black")) +
  theme(plot.title = element_text(size=20, face = "plain")) +
  theme(panel.grid = element_blank())+
  labs(y = "Z-score") +
  theme(legend.position="none") +
  scale_color_manual(values = levels(df$col)) +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")
q
ggsave(paste0("mast/group.dependent.degs", Sys.Date(), ".pdf"), height = 5, width = 7)
dev.off()

#==================================================================================================
# Celltype age enrichment for each replicate
#==================================================================================================

#load(file = "../seurat/all_All6_Filtered_meta.data&umap2022-11-21.rda")
m <- as.matrix(all[["RNA"]]@counts)

save(m, file = "mast/10x_RawExpressionMatrix.rdata")
#write.csv(m, file="10x_RawExpressionMatrix.csv")
#m <- read.csv(file = "../seurat/10x_RawExpressionMatrix.csv",header = T)

embeds <- all@reductions$umap@cell.embeddings
uMAP_1 <- embeds[,1] * -1 # Flip x-axis for aesthetics; does not affect any clustering/inferences.
uMAP_2 <- embeds[,2]
all_alldata <- cbind(all@meta.data, uMAP_1, uMAP_2)
all_alldata <- all_alldata[, !grepl("res.0*", colnames(all_alldata))]
rownames(all_alldata) <- NULL

#d <- cbind(all_alldata, t(m))
d <- cbind(all_alldata, colnames(m))
d <- d[, !grepl("res.0*", colnames(d))]
rownames(d) <- NULL

# Group by cell type and sample, then count cells.
# Then divide sample cell type counts by total sample cell counts. 
keep = c("Cell_type", "orig.ident")
counts <- d %>% select(one_of(keep)) %>% group_by(orig.ident, Cell_type) %>% dplyr::mutate(count = n()) %>% distinct()
tmp <- data.frame(Cell_type="CCL24 Macrophages",orig.ident=c("A_14D","A_7D"),count=0)
counts <- rbind(tmp, counts)
totals <- aggregate(counts$count, by=list(orig.ident=counts$orig.ident), FUN=sum)
colnames(totals) <- c("orig.ident", "sample.sum")
prop_df <- merge(counts, totals)
#prop_df <- rbind(prop_df,c("B_14D","NFkb Fibroblast","1","4470"),c("A_7D","NFkb Fibroblast","1","3969"))
prop_df$prop <- prop_df$count / prop_df$sample.sum
write.csv(prop_df, file = "mast/CelltypeCounts.csv")

enrich <- c(); rep <- c(); cell <- c(); y_prop <- c(); o_prop <- c()
for (celltype in unique(prop_df$Cell_type)) {
  for (r in c("7D","14D","1M","2M")) {
    old_prop <- prop_df[prop_df$orig.ident==paste0("A_", r) & prop_df$Cell_type==celltype, "prop"]
    young_prop <- prop_df[prop_df$orig.ident==paste0("B_", r) & prop_df$Cell_type==celltype, "prop"]
    y_prop <- c(y_prop, young_prop)
    o_prop <- c(o_prop, old_prop)
    e <- old_prop/young_prop
    enrich <- c(enrich, e)
    rep <- c(rep, r)
    cell <- c(cell, celltype)
  }
}

enrich_df <- data.frame("Rep" = rep, "Celltype" = cell, "Y_prop" = y_prop, "O_prop" = o_prop, "Log2Enrich" = log2(enrich))

# Supplementary Table 3, sheet 2:
write.csv(enrich_df, file = "mast/ProportionEnrichments.csv" )

#==================================================================================================
enrich_df <- read.csv(file = "mast/ProportionEnrichments.csv")
library(plyr)
# Calculate standard error
sumr <- plyr::ddply(enrich_df, c("Celltype"), summarise,
                    N    = length(Log2Enrich),
                    mean = mean(Log2Enrich),
                    sd   = sd(Log2Enrich),
                    se   = sd / sqrt(N)
)
enrich_df <- merge(enrich_df, sumr, by = "Celltype")
head(enrich_df)

# Adjust colors
enrich_df$celltype_factor <- factor(enrich_df$Celltype,  levels=CELLS, ordered=T)
names(new_colors) <- levels(enrich_df$celltype_factor)
enrich_df[9:10,5] <- c(0,0)


p <- ggplot(enrich_df, aes(x = celltype_factor, y = Log2Enrich)) +
  geom_bar(aes(fill = celltype_factor), stat = "summary", fun.y = "mean",  alpha = .5) +
  geom_errorbar(aes(x=celltype_factor, ymin=mean-se, ymax=mean+se), width=.15) +
  #geom_jitter(width = 0.40, alpha = .55, size = 1) +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) +
  ggtitle("Changes in gene expression with group") +
  theme(axis.title.y = element_text(size = 20, face = "plain")) +
  theme(axis.text.y = element_text(size = 20,colour = "black")) +
  theme(axis.text.x = element_text(size = 11,colour = "black")) +
  theme(plot.title = element_text(size=20, face = "plain")) +
  theme(panel.grid = element_blank())+
  labs(y = "Log2 Enrichment") +
  theme(legend.position="none") +  
  scale_fill_manual(values = new_colors) +
  scale_color_manual(values = new_colors) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.3, binwidth = 1)+
  ggtitle("Changes in cell number with group") +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")
p
ggsave(paste0("mast/group.enrich", Sys.Date(), ".pdf"), height = 5, width = 7)
dev.off()

#==================================================================================================
# Combine into one plot

pq <- plot_grid(p, NULL, q, ncol = 1, rel_heights = c(1.5, 0, 1.5))
ggsave("mast/ProportionsPlusDE.png", pq, height = 10, width = 7, dpi = 300)

#==================================================================================================
sessionInfo()
