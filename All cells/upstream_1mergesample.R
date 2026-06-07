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

setwd("F:/scrna_seq_raw_2/")
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr,tibble)
#source("./scfunctions.R")
A_7D <- CreateSeuratObject(Read10X('./raw/A_7D/'),
                           "A_7D",min.cells=3,min.features=200)
A_14D <- CreateSeuratObject(Read10X('./raw/A_14D/'),
                             "A_14D",min.cells=3,min.features=200)
A_1M <- CreateSeuratObject(Read10X('./raw/A_1M/'),"A_1M",min.cells=3,min.features=200)
A_2M <- CreateSeuratObject(Read10X('./raw/A_2M/'),"A_2M",min.cells=3,min.features=200)
B_7D <- CreateSeuratObject(Read10X('./raw/B_7D/'),"B_7D",min.cells=3,min.features=200)
B_14D <- CreateSeuratObject(Read10X('./raw/B_14D/'),"B_14D",min.cells=3,min.features=200)
B_1M <- CreateSeuratObject(Read10X('./raw/B_1M/'),"B_1M",min.cells=3,min.features=200)
B_2M <- CreateSeuratObject(Read10X('./raw/B_2M/'),"B_2M",min.cells=3,min.features=200)

all <- merge(A_7D,
                 y = c(A_14D,A_1M,A_2M,B_7D,B_14D,B_1M,B_2M),
                 add.cell.ids = c("A_7D","A_14D","A_1M","A_2M","B_7D","B_14D","B_1M","B_2M"),
                 project = "MuscleRegeneration")
save(all,file="mergedraw.Rdata")
                 