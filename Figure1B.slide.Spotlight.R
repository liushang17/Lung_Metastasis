library(SPOTlight)
library(data.table)
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggsci)

args <- commandArgs(T)
infile <- args[1]
outdir <- args[2]

dir.create(paste0(outdir,"/figures"),recursive = T)
setwd(outdir)

pro <- basename(args[1])
pro <- gsub("_Spatial_Seurat.rds","",pro)

lung_st <- readRDS(infile)


celltype_markers_all <- read.table(args[3],header = T,stringsAsFactors = F,sep = "\t")
celltype_markers_all <- celltype_markers_all[order(celltype_markers_all$p_val),]
celltype_markers_all <- celltype_markers_all[order(celltype_markers_all$cluster),]

Idents(lung_sc) <- lung_sc@meta.data$Cell_Type
set.seed(123)

counts_mat <- as.matrix(lung_st@assays$Spatial@counts)
gc()

spotlight_ls <- spotlight_deconvolution(se_sc = lung_sc,
                                        counts_spatial = counts_mat, clust_vr = "free_annotation",
                                        cluster_markers = celltype_markers_all, cl_n = 100,
                                        hvg = 0,
                                        ntop = 100,
                                        transf = "uv",
                                        method = "nsNMF",
                                        min_cont = 0.1)

saveRDS(object = spotlight_ls, file = paste0(pro,"_spotlight_ls.rds"))

nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

pdf(paste0("figures/",pro,"_SPOTlight_predict_result.pdf"),17,15)
h <- NMF::coef(nmf_mod[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(h = h,train_cell_clust = nmf_mod[[2]])

topic_profile_plts[[2]] + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 90),
  axis.text = ggplot2::element_text(size = 12))

cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]

lung_st@meta.data <- cbind(lung_st@meta.data, decon_mtrx)
predict_CellType = apply(decon_mtrx[,which(colnames(decon_mtrx) != "res_ss")], 1,
                         function(x){
                           index = which.max(x)
                           a = colnames(decon_mtrx)[index]
                           return(a)}
)

lung_st@meta.data$predict_CellType = predict_CellType

col <- c(pal_npg("nrc")(10),"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")


DimPlot(lung_st, reduction="umap",cols = col,group.by = "predict_CellType",pt.size = 1.2)
SpatialDimPlot(lung_st,cols=col,stroke = 0,group.by = "predict_CellType",pt.size = 1.2)
SpatialDimPlot(lung_st,cols=col,stroke = 0,group.by = "seurat_clusters")
dev.off()

outinfor <- data.frame(lung_st@meta.data,fig.X = lung_st@images$slice1@coordinates$col, fig.Y= (max(lung_st@images$slice1@coordinates$row)-lung_st@images$slice1@coordinates$row),check.names=F)
outinfor <- outinfor[,-c(1)]

write.table(lung_st@meta.data,file=paste0(pro,"_predict_CellType.txt"),quote=F,col.name=T,row.names=T)


size <- max(diff(range(outinfor$fig.X)),diff(range(outinfor$fig.Y)))
pdf(paste0("figures/",pro,"_SPOTlight_predict_CellType.pdf"),0.025*size,0.025*size)
p <- ggplot(outinfor,aes(x=fig.X,y=fig.Y,color=predict_CellType))+
  geom_point(size=0.0006*size,shape=15)+
  scale_color_manual(values = col)+
  coord_fixed(ratio = 1)+
  theme_classic()
print(p)
dev.off()

