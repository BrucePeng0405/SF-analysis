#!/share/home/ouyanghongwei/bin/Rscript
#BSUB -J pz
#BSUB -n 5
#BSUB -R "span[ptile=24]"
#BSUB -o output_%Js
#BSUB -e errput_%J
#BSUB -q normal
rm(list = ls())


library(Seurat)
library(ggplot2)
library(dplyr,tibble)
load("muscle_regeneration_4markers_annotated_0.3.Rdata")
#all <- subset(all, downsample = 2000)
gc()

unique(all$orig.ident)
unique(all$group)
table(all$Cell_type)

#画甜甜圈图
library(plotrix)
library(ggsci)
library(RColorBrewer)
cols <-c(brewer.pal(10,"Paired"),brewer.pal(8,"Accent"))
doughnut <- function (x, labels = names(x), edges = 200, outer.radius = 0.8,
                      inner.radius=0.6, clockwise = FALSE,
                      init.angle = if (clockwise) 90 else 0, density = NULL,
                      angle = 45, col = NULL, border = FALSE, lty = NULL,
                      main = NULL, ...)
{
  if (!is.numeric(x) || any(is.na(x) | x < 0))
    stop("'x' values must be positive.")
  if (is.null(labels))
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L])
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col))
    col <- if (is.null(density))
      palette()
  else par("fg")
  col <- rep(col, length.out = nx)
  border <- rep(border, length.out = nx)
  lty <- rep(lty, length.out = nx)
  angle <- rep(angle, length.out = nx)
  density <- rep(density, length.out = nx)
  twopi <- if (clockwise)
    -2 * pi
  else 2 * pi
  t2xy <- function(t, radius) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p),
         y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n),
              outer.radius)
    polygon(c(P$x, 0), c(P$y, 0), density = density[i],
            angle = angle[i], border = border[i],
            col = col[i], lty = lty[i])
    Pout <- t2xy(mean(x[i + 0:1]), outer.radius)
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      lines(c(1, 1.05) * Pout$x, c(1, 1.05) * Pout$y)
      text(1.1 * Pout$x, 1.1 * Pout$y, labels[i],
           xpd = TRUE, adj = ifelse(Pout$x < 0, 1, 0),
           ...)
    }      
    Pin <- t2xy(seq.int(0, 1, length.out = n*nx),
                inner.radius)
    polygon(Pin$x, Pin$y, density = density[i],
            angle = angle[i], border = border[i],
            col = "white", lty = lty[i])
  }
  
  title(main = main, ...)
  invisible(NULL)
}
df <- table(all$Cell_type) %>% as.data.frame()
labs <- paste0(df$Var1,"(",round(df$Freq/sum(df$Freq)*100,2),"%)")
p.circle <- doughnut(
  df$Freq,
  labels = labs,
  init.angle = 90,#初始角度
  col = cols,#颜色
  border = "white",#边框颜色
  inner.radius = 0.4,#内环大小
  cex = 0.8#字体大小
)

#画堆积柱状图
cellnum <- table(all$Cell_type, all$orig.ident)
cell.prop<-as.data.frame(prop.table(cellnum))
colnames(cell.prop)<-c("Celltype","Group","Proportion")

cell.prop$Group = factor(cell.prop$Group, levels=c('A_7D','B_7D','A_14D','B_14D','A_1M','B_1M','A_2M','B_2M'))
p.bar <- ggplot(cell.prop,aes(Group,Proportion,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=cols)+#自定义fill的颜色
  ggtitle("Cell proportion")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))
p.bar
