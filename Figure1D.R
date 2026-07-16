library(Seurat)
library(ggsci)
library(ggrepel)
library(dplyr)
library(harmony)
library(umap)
col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")

infile <- args[1]
outdir <- args[2]

Combine <- readRDS(infile)


ann <- Combine@meta.data
ann$time <- gsub("Casp8_","",ann$time)
ann$time <- gsub("_.*","",ann$time)
ann$slide <- gsub(":.*","",ann$allbin)
ann$time <- factor(ann$time,levels = paste0("G",c(1,2,4,5,6,7,9,11,12)))

ann <- Combine@meta.data
ann1 <- ann[,c(12:33)]
#tms <- FindNeighbors(as.matrix(ann1[,c(3:7,9,12:18,20:22)]))
tms <- FindNeighbors(as.matrix(ann1))
tms1 <- FindClusters(tms$snn,resolution = 0.1)
colnames(tms1) <- "cluster"
ann$new_cluster <- tms1$cluster

table(ann$new_cluster,ann$seurat_clusters)

tsne_out <- umap((ann1))
umapinfo <- data.frame(tsne_out$layout)
umapinfo$cluster <- ann$new_cluster

ggplot(umapinfo,aes(x=as.numeric(X1),y=as.numeric(X2),color = cluster))+geom_point()+
  scale_color_manual(values = col) + labs(x="UMAP_1",y="UMAP_2") + theme_classic()

ann$stage <- "Stage1"
pos <- which(ann$pseu > 5)
ann$stage[pos] <- "Stage2"
pos <- which(ann$pseu > 24.3)
ann$stage[pos] <- "Stage3"
pos <- which(ann$pseu > 28.4)
ann$stage[pos] <- "Stage4"
table(ann$stage,ann$new_cluster)
ann3 <- ann
ann3$batch <- gsub("Casp8_","",ann3$time)
ann3$batch <- gsub("_.*","",ann3$batch)
clus <- as.character(c("Stage1","Stage2","Stage3","Stage4"))
clus <- c("G1","G2","G4","G5","G6","G7","G9","G11","G12")
mergeclus <- unique(as.character((ann3$new_cluster)))
mit <- matrix(nrow = length(clus),ncol = length(mergeclus))
for(i in 1:length(clus)){
  for(j in 1:length(mergeclus)){
    pos1 <- which(ann3$batch %in% clus[i] & ann3$new_cluster %in% mergeclus[j])
    pos2 <- which(ann3$new_cluster %in% mergeclus[j])
    pos3 <- which(ann3$batch %in% clus[i])
    mit[i,j] <- length(pos1) / length(pos2) / length(pos3) * 10000
  }
}
rownames(mit) <- clus
colnames(mit) <- mergeclus
pheatmap::pheatmap(mit,cluster_rows = F)
mit <- mit[,mergeclus[c(1,2,3,5,4)]]
pheatmap::pheatmap(mit,cluster_rows = F,cluster_cols = F)
mit[mit > 4] <- 4
pheatmap::pheatmap((mit),cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("lightblue","white","red"))(100))

saveRDS(umapinfo,file = paste0(outdir,"/TME.cluster.rds"))
