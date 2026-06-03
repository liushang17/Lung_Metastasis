library(ggplot2)
library(Seurat)
library(ggsci)
library(ggrepel)
library(dplyr)
library(plyr)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(sfsmisc)
col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")
#########
infile <- args[1] # The Spatial seurat rds for tumor cells

#########
Combine <- readRDS(infile)

ann <- Combine@meta.data
ann$immune <- rowSums(ann[,11+c(3:7,9,13:18,20:22)])

require(splines)

colnames(ann)
pos <- which(colnames(ann) %in% "immune")
annt <- ann[,c(11,pos)]
annt$Celltype <- feas[i]
colnames(annt) <- c("Pseudotime","Value","Celltype")

model <- loess(annt$Value~annt$Pseudotime,data=annt,span = 0.3)
y <- predict(model)
x <- annt$Pseudotime
#plot(x,y,main="Original fit") + abline(h = 0, v = c(5,24.3,28.4,27.5,29.22474))

anns <- data.frame(x= x,y = y)
anns <- anns[order(anns$x),]
anns$f <- D1ss(anns$x,anns$y)
#plot(anns$x,anns$f,main="first fit") + abline(h = 0, v = c(5,24.3,27.3,28.4,30.7))

tes <- D2ss(anns$x,anns$y)
anns$s <- tes$y
#plot(anns$x,anns$s,main="Der fit") + abline(h = 0, v = c(5,24.3,27.3,28.4,30.7))

