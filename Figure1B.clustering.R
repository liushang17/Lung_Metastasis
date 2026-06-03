library(Seurat)
library(pheatmap)
library(ggsci)
library(ggrepel)
col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")

args = commandArgs(T)
infile <- args[1]
blackfile <- args[2]
outdir <- args[3]

##########
lung.tumor <- readRDS(infile)
batchs <- sort(unique(as.character(lung.tumor$batch)))

celltype <- "Hepa1.6"
genefilter <- NULL
for(i in 1:length(batchs)){
  filepath <- paste0(blackfile,"/",batchs[i],".de.xls")
  if(file.exists(filepath)){
    genes <- read.table(filepath,sep = "\t")
    pos <- which(genes$V7 %in% celltype)
    genes1 <- genes[pos,]
    genefilter <- c(genefilter,as.character(genes1$V6))
  }
}
genefilter <- unique(genefilter)
genefilter <- c(genefilter)

pos <- which(rownames(lung.tumor@assays$RNA@counts) %in% c(genefilter))
genename1 <- rownames(lung.tumor@assays$RNA@counts)[-pos]
lung.tumor <- subset(lung.tumor,features = genename1)


ob.list <- list() 
for(i in 1:length(batchs)){
  pos <- which(lung.tumor$batch %in% batchs[i])
  cellname <- colnames(lung.tumor)[pos]
  lung.tumor1 <- subset(lung.tumor,cells = cellname)
  Combine.tmp <- CreateSeuratObject(lung.tumor1@assays$RNA@counts)
  Combine.tmp <- NormalizeData(object = Combine.tmp,normalization.method = "LogNormalize", scale.factor = 10000)
  Combine.tmp <- FindVariableFeatures(Combine.tmp, selection.method = "vst", nfeatures = 2000 )
  Combine.tmp@meta.data$time <- lung.tumor1$batch
  tmp1 <- list(Combine.tmp)
  ob.list <- c(ob.list,tmp1)
}

dim.u <- 25
Combine.anchors = FindIntegrationAnchors(object.list = ob.list,dims = 1:dim.u)
Combine = IntegrateData(anchorset = Combine.anchors,dims = 1:dim.u)
double_gene_info<-data.frame(ID=rownames(Combine@assays$RNA@data),gene_short_name=rownames(Combine@assays$RNA@data))
rownames(double_gene_info)<-double_gene_info$ID
DefaultAssay(Combine) <- "integrated"
Combine
Combine = ScaleData(Combine,verbose = FALSE)
Combine = RunPCA(Combine, npcs = 30, verbose = FALSE)
Combine = FindNeighbors(Combine, reduction = "pca",dims = 1:dim.u)
Combine = FindClusters(Combine, resolution = 0.6)
Combine = RunUMAP(Combine, reduction = "pca", dims = 1:25)

saveRDS(Combine,file = paste0(outdir,"/Combine.all.rds"))
