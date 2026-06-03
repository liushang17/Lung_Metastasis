args = commandArgs(T)

library("Seurat")
lung.merge2 <- readRDS(args[1])

names <- basename(args[1])
names1 <- gsub(".around.seurat.rds","",names)
names2 <- gsub("Casp8\\w+.","",names1)

lung.merge2$clus <- "Other"
pos <- grep(names2,colnames(lung.merge2))
lung.merge2$clus[pos] <- names2

clusters <- unique(as.character(lung.merge2$seurat_clusters))

sites <- NULL
for(i in 1:length(clusters)){
	pos1 <- which(lung.merge2$seurat_clusters %in% clusters[i])
	pos2 <- which(lung.merge2$seurat_clusters %in% clusters[i] & lung.merge2$clus %in% names2)
	tmp <- length(pos2) / length(pos1)
	sites <- c(sites,tmp)
}

pos <- which(sites > 0.5)
clus_chose <- clusters[pos]

if(length(clus_chose) > 0){
	de_markers <- FindAllMarkers(lung.merge2,only.pos = TRUE,assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
	pos <- which(de_markers$p_val_adj < 0.05)
	de_markers1 <- de_markers[pos,]
	pos <- which(de_markers1$cluster %in% clus_chose)
	if(length(pos) > 0){de_markers2 <- de_markers1[-pos,]}else{de_markers2 <- de_markers1}
	de_markers2$slide <- names2
	filename <- paste0(args[2],"/",names1,".blacklist.xls")
	write.table(de_markers2,file = filename,sep = "\t",quote=F,row.names=F,col.names=F)
}
