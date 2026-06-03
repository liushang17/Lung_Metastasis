library(dplyr)
library(Seurat)
library(data.table)
library(Matrix)
library(ggplot2)
library(ggsci)
library(plotly)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(cowplot)
library(ggthemes)
library(SeuratWrappers)


args <- commandArgs(T)
seuratSpaObj <- readRDS(args[1])
bs <- as.numeric(args[2])
flitG <- as.numeric(args[3])
outdir <- args[4]

dir.create(outdir)
setwd(outdir)

pro <- tail(unlist(strsplit(outdir,"/")),1)
pro <- gsub(".txt|.gz|_filterd|.bin1.Lasso.gene","",pro)

##############################  3. Spatial Analyse  ##############################
Q1 <- quantile(seuratSpaObj$nFeature_Spatial)[2]
Q3 <- quantile(seuratSpaObj$nFeature_Spatial)[4]
upper <- as.numeric(Q3+1.5*(Q3-Q1))
lower <- as.numeric(Q1-1.5*(Q3-Q1))

seuratSpaObj <- subset(seuratSpaObj,subset = nFeature_Spatial>flitG &  percent.mt < 30)
pdf(paste0("figures/",pro,"_QC.Feature.pdf"))
VlnPlot(seuratSpaObj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol=3,pt.size = 0) +theme(axis.text.x=element_text(angle=20,size=9))+labs(x=paste0("nGene: ",dim(seuratSpaObj)[1],"; ","nBIN: ",dim(seuratSpaObj)[2]))
plot(density(seuratSpaObj$nFeature_Spatial))+abline(v=c(500,650,800),col="grey",lwd=2,lty=6)
SpatialFeaturePlot(seuratSpaObj, features="nFeature_Spatial") + theme(legend.position = "right")
SpatialFeaturePlot(seuratSpaObj, features="nCount_Spatial") + theme(legend.position = "right")
dev.off()
seuratSpaObj

###############  SCTransform
set.seed(6)

seuratSpaObj <- SCTransform(seuratSpaObj,assay = "Spatial", verbose = FALSE)
seuratSpaObj <- RunPCA(seuratSpaObj,assay = "SCT",verbose = F)
seuratSpaObj <- FindNeighbors(seuratSpaObj, reduction = "pca", dims = 1:15)
seuratSpaObj <- RunUMAP(seuratSpaObj, reduction = "pca", dims = 1:15)
seuratSpaObj <- FindClusters(seuratSpaObj, verbose = FALSE,resolution = 0.5)

col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")

if(length(levels(seuratSpaObj))>length(col)){
  col <- colorRampPalette(col)(length(levels(seuratSpaObj)))}else{
    col <- col[1:length(levels(seuratSpaObj))]
  }

plot2 <- ElbowPlot(seuratSpaObj, ndims=50, reduction="pca")
pdf(paste0("figures/",pro,"_SCTran_UMAP.pdf"))
DimPlot(seuratSpaObj, reduction="umap", label=TRUE,cols = col,size=3,label.size = 6)
SpatialDimPlot(seuratSpaObj,cols=col,stroke = 0)
FeaturePlot(seuratSpaObj, features="nFeature_Spatial") + theme(legend.position = "right")
print(plot2)
dev.off()

png(paste0("figures/",pro,"_SCTran_Cluster.png"), width=1000,height=1000,res=100)
SpatialDimPlot(seuratSpaObj, cells.highlight=CellsByIdentities(object=seuratSpaObj),pt.size.factor = 0.1,stroke = 0,cols.highlight = c("#DE2D26", "grey90"), facet.highlight=TRUE)
dev.off()

outinfor <- data.frame(seuratSpaObj@meta.data,fig.X = seuratSpaObj@images$slice1@coordinates$col, fig.Y= (max(seuratSpaObj@images$slice1@coordinates$row)-seuratSpaObj@images$slice1@coordinates$row),bin.x=seuratSpaObj@images$slice1@coordinates$col,bin.y=seuratSpaObj@images$slice1@coordinates$row)
outinfor <- outinfor[,-c(1,7)]
fwrite(outinfor,paste0(pro,"_cellcluster_Axis.txt"),col.names=T,row.names=T,sep="\t",quote=F)


size <- max(diff(range(outinfor$fig.X)),diff(range(outinfor$fig.Y)))
outinfor$seurat_clusters <- as.factor(outinfor$seurat_clusters)
pdf(paste0("figures/",pro,"_cluter.pdf"),0.025*size,0.025*size)
p <- ggplot(outinfor,aes(x=fig.X,y=fig.Y,color=seurat_clusters))+
  geom_point(size=0.0006*size,shape=15)+
  scale_color_manual(values = col)+
  coord_fixed(ratio = 1)+
  theme_classic()
print(p)
dev.off()

p <- ggplot(outinfor,aes(x=fig.X,y=fig.Y,color=seurat_clusters))+
  geom_point(size=0.5,shape=15)+
  scale_color_manual(values = col)+
  coord_fixed(ratio = 1)+
  theme_classic()
p <- ggplotly(p)
htmlwidgets::saveWidget(as_widget(p), paste0(outdir,"/figures/",pro,"_cluter.html"))

system( paste0("rm -rf ",outdir,"/figures/",pro,"_cluter_files"))

####################### FeaturePlot ###########################
pdf(paste0("figures/",pro,"_FeaturePlot.pdf"), 10, 10)
for (i in 1:ceiling(length(genes)/12)) {
  if(i<ceiling(length(genes)/12)){
    plot(VlnPlot(seuratSpaObj,features=genes[(12*(i-1)+1):(12*i)], assay="Spatial",slot="data", pt.size=0,ncol = 3,log=T,cols=col))
  }else{
    plot(VlnPlot(seuratSpaObj,features=genes[(12*(i-1)+1):length(genes)], assay="Spatial",slot="data", pt.size=0,ncol = 3,log=T,cols=col))
  }
}
dev.off()

outinfor <- cbind(outinfor,t(as.matrix(seuratSpaObj@assays$Spatial@counts[genes,])))
cols <- c("#FFFFFF00","#0C3383","#005EA3","#0A88BA","#00C199","#F2D338","#F6B132","#F28F38","#E76124","#D91E1E","#9E0142")

xylim <- c(range(outinfor$fig.X),range(outinfor$fig.Y))

plot_list <- list()
for (i in 1:length(genes)) {
  data_temp <- outinfor[,c("fig.X","fig.Y",genes[i])]
  data_temp <- data_temp[data_temp[,3]>0,]
  plot_data <- unlist(apply(data_temp, 1, function(x) rep(c(x[1:2]),x[3])))
  plot_data <- as.data.frame(matrix(plot_data,ncol = 2,byrow = T))
  colnames(plot_data) <- c("X","Y")
  
  plot_list[[i]] <- ggplot(plot_data,aes(x=X,y=Y)) + 
    stat_density_2d(mapping = aes(fill = ..density.. ), geom = "raster", contour = FALSE,h = 2,n = 500) +
    scale_fill_gradientn(colours = cols) + xlim(xylim[1], xylim[2]) + ylim(xylim[3], xylim[4])+
    coord_fixed(ratio = 1)+
    labs(title = genes[i],x="",y="")+
    theme_few() + theme(legend.position='none',plot.title=element_text(face="italic",size=20, hjust=0.5),axis.text = element_blank(),axis.ticks = element_blank())
}

for(i in 1:ceiling(length(genes)/9)){
  png(paste0("figures/",pro,"_SpatialFeaturePlot_",i,".png"), width=1000,height=1000,res=100)
  if(i<ceiling(length(genes)/9)){
    p1 <- cowplot::plot_grid(plotlist=plot_list[(9*(i-1)+1):(9*i)])
  }else{
    p1 <- cowplot::plot_grid(plotlist=plot_list[(9*(i-1)+1):length(genes)])}
  print(p1)
  dev.off()
}

DefaultAssay(seuratSpaObj) <- "Spatial"
p <- FeaturePlot(seuratSpaObj,features = genes,label=0.02,label.size=0.02,order=T,combine = FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes()
}

for(i in 1:ceiling(length(genes)/9)){
  png(paste0("figures/",pro,"_SCTran_Markers_",i,".png"), width=1000,height=1000,res=100)
  if(i<ceiling(length(genes)/9)){
    p1 <- cowplot::plot_grid(plotlist=p[(9*(i-1)+1):(9*i)])
  }else{
    p1 <- cowplot::plot_grid(plotlist=p[(9*(i-1)+1):length(genes)])}
  print(p1)
  dev.off()
}

de_markers <- RunPrestoAll(seuratSpaObj,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,max.cells.per.ident=1000)
de_markers <- de_markers[de_markers$p_val_adj<0.05,]
top10 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
seuratSpaObj <- SCTransform(seuratSpaObj,assay = "Spatial", residual.features=rownames(seuratSpaObj))

pdf(paste0("figures/",pro,"_SCTran_top10Heat.pdf"),10,10)
DoHeatmap(seuratSpaObj, features=top10$gene) + NoLegend()
dev.off()
write.table(de_markers,file=paste0(pro,"Spatial_markers.txt"),sep = "\t",quote = F)

saveRDS(seuratSpaObj@assays$Spatial@counts,file=paste0(pro,"_counts_mat.rds"))

if(dim(seuratSpaObj@assays$Spatial@scale.data)[1]>2){
    seuratSpaObj@assays$Spatial@scale.data <- seuratSpaObj@assays$Spatial@scale.data[1:2,1:2]
}
if(dim(seuratSpaObj@assays$SCT@scale.data)[1]>2){
    seuratSpaObj@assays$SCT@scale.data <- seuratSpaObj@assays$SCT@scale.data[1:2,1:2]
}

saveRDS(seuratSpaObj,file=paste0(pro,"_Spatial_Seurat.rds"))

