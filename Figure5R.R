library(Seurat)
library(ggplot2)



library(ggplot2)
library(Seurat)
library(ggsci)
library(ggrepel)
library(dplyr)
library(plyr)
library(ggsci)
library(ggplot2)
library(ggpubr)
col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")

infile <- args[1] # The Spatial seurat rds for tumor cells
infile1 <- args[2] # The LR score for TME
obj <- readRDS(infile1)
Combine <- readRDS(infile)

ann <- Combine@meta.data
ann$immune <- rowSums(ann[,11+c(3:7,9,13:18,20:22)])

ann$point <- gsub("Casp8_","",ann$time)
ann$point <- gsub("_.*","",ann$point)

ann$stage <- "Stage1"
pos <- which(ann$pseudotime > 5)
ann$stage[pos] <- "Stage2"
pos <- which(ann$pseudotime > 24.3)
ann$stage[pos] <- "Stage3"
pos <- which(ann$pseudotime > 27.4)
ann$stage[pos] <- "Stage5"
pos <- which(ann$pseudotime > 31)
ann$stage[pos] <- "Stage4"

ann1 <- ann[colnames(obj),]

Idents(obj) <- factor(ann1$stage,levels = c("Stage1","Stage2","Stage3","Stage5","Stage4"))
obj$TME_Cluster <-factor(ann1$stage,levels = c("Stage1","Stage2","Stage3","Stage5","Stage4"))

ligrec <- rownames(obj@assays$RNA@counts)
chemo <- ligrec[grepl("Ccl|Cxcl|Cx3|Xcl|Cx3cl|Csf",ligrec)]

DotPlot(obj,features=chemo,cols = c("#fffbd5","#b20a2c"), group.by="TME_Cluster")+RotatedAxis()

chemo_sub <- c("Ccr1-Ccl7","Ccr1-Ccl8","Cxcr6-Cxcl16","Cxcl2-Cxcr2","Cx3cl1-Cx3cr1","Cx3cr1-Cx3cl1")
chemo_sub <- c("Cmklr1-Rarres2","Podxl-Sell","Nrp1-Vegfa","Sema3e-Plxnd1","Sell-Cd34","Cxcl10-Dpp4","Igf1-Igf1r")
pos <- grep("Cd274",rownames(obj))
chemo_sub <- c(rownames(obj)[pos])
pos <- grep("Cd47",rownames(obj))
chemo_sub1 <- c(rownames(obj)[pos])
chemo_sub <- c(chemo_sub,chemo_sub1)
DotPlot(obj,features=rev(chemo_sub),cols = c("#fffbd5","#b20a2c"), group.by="TME_Cluster")+RotatedAxis()+coord_flip()
DotPlot(obj,features="Igf1r-Igf1",cols = c("#fffbd5","#b20a2c"), group.by="TME_Cluster")+RotatedAxis()+coord_flip()

VlnPlot(obj,features = c("Igf1-Igf1r"),pt.size = 0)+geom_boxplot()

obj1 <- subset(obj,idents = paste0("Stage",c(1,2,3,4,5)))

DotPlot(obj1,features=c("Igf1-Igf1r"),cols = c("#fffbd5","#b20a2c"), group.by="TME_Cluster",dot.scale= 10)+RotatedAxis()
######
clus.mak <- FindMarkers(obj,ident.1 = "Stage5",group.by = "TME_Cluster",logfc.threshold = 0,min.pct = 0.05)

clus.mak$pro <- abs(clus.mak$pct.1-clus.mak$pct.2)
pos <- which(clus.mak$pro > 0.1)
clus.mak <- clus.mak[pos,]

genename <- c("Cd44","Nrp1","Ror1","Igf1r","Fgfr2")
pos1 <- grep(genename[1],rownames(clus.mak))
pos2 <- grep(genename[2],rownames(clus.mak))
pos3 <- grep(genename[3],rownames(clus.mak))
pos4 <- grep(genename[4],rownames(clus.mak))
pos5 <- grep(genename[5],rownames(clus.mak))
pos <- c(pos1,pos2,pos3,pos4,pos5)
lrc <- rownames(clus.mak)[pos]

clus.mak$type <- "Other"
pos <- which(clus.mak$p_val_adj < 0.01 & clus.mak$avg_log2FC > 0.1)
clus.mak$type[pos] <- "Up"
pos <- which(clus.mak$p_val_adj < 0.01 & clus.mak$avg_log2FC < -0.1)
clus.mak$type[pos] <- "Dn"
clus.mak$P <- -log10(clus.mak$p_val_adj)
clus.mak$name <- NA
pos <- which(rownames(clus.mak) %in% lrc)
clus.mak$name[pos] <- rownames(clus.mak)[pos]
ggplot(clus.mak,aes(x=avg_log2FC,y=P))+geom_point(size = 2,aes(color = type))+
  geom_text_repel(data = clus.mak,aes(avg_log2FC, P, label = name))+
  theme_classic()+scale_colour_manual(values = c("grey","red"))

###
clus.mak <- FindMarkers(obj,ident.1 = "Stage5",group.by = "TME_Cluster",logfc.threshold = 0,min.pct = 0.05)
pos <- which(clus.mak$avg_log2FC > 0 & clus.mak$p_val_adj < 0.05)
clus.mak1 <- clus.mak[pos,]
lr <- rownames(clus.mak1)
lr1 <- gsub(".*-","",lr)
lr2 <- gsub("-.*","",lr)

mit <- as.matrix(obj@assays$RNA@counts)
clus.mak1 <- data.frame(lr = rownames(mit))
lr <- rownames(mit)
lr1 <- gsub(".*-","",lr)
lr2 <- gsub("-.*","",lr)
clus.mak1$l <- lr1
clus.mak1$r <- lr2

pos <- which(rownames(de1) %in% c(lr1,lr2) )
de1[pos,]
pos <- which(clus.mak1$l %in% rownames(de1) | clus.mak1$r %in% rownames(de1))
clus.mak2 <- clus.mak1[pos,]

pos <- which(poi$avg_log2FC > 0 & poi$p_val_adj < 0.01)
poi <- poi[pos,]
pos <- which(clus.mak2$l %in% rownames(poi) | clus.mak2$r %in% rownames(poi))
clus.mak3 <- clus.mak2[pos,]

pos1 <- which(obj$TME_Cluster %in% "Stage5")

for(i in 1:nrow(clus.mak3)){
  pos <- which(rownames(mit) %in% clus.mak3$lr[i])
  sed1 <- as.numeric(as.character(mit[pos,pos1]))
  sed2 <- as.numeric(as.character(mit[pos,-pos1]))
  tmp1 <- length(which(sed1 > 0))
  tmp2 <- length(sed1) - tmp1
  tmp3 <- length(which(sed2 > 0))
  tmp4 <- length(sed2) - tmp3
  if(tmp1 > 10 & tmp2 > 10){
    data <- matrix(c(tmp1,tmp2,tmp3,tmp4),nrow = 2)
    tmp4 <- fisher.test(data)
    clus.mak3$pvd[i]<- tmp4$p.value
    clus.mak3$odd[i] <- tmp4$estimate
  }else{
    data <- matrix(c(tmp1,tmp2,tmp3,tmp4),nrow = 2)
    tmp4 <- fisher.test(data)
    clus.mak3$pvd[i]<- 1
    clus.mak3$odd[i]<- tmp4$estimate
  }
}

clus.mak4 <- clus.mak3[order(-clus.mak3$odd),]
pos <- which(clus.mak4$pvd > 0.5)
clus.mak4 <- clus.mak4[-pos,]
terms <- as.character(clus.mak4$lr)
ggplot(clus.mak4,aes(x=factor(lr,levels = terms),y=odd))+geom_bar(stat = "identity")+
  theme_classic()+theme(axis.text.x = element_text(angle = 60))

