library("ArchR")
library(ggplot2)


addArchRGenome("mm10")

######### The arrow file
setwd("arrow.mouse/") # the arrow directory of ATAC-seq
ArrowFiles <- list.files("/arrow.mouse/")
ArrowFiles <- ArrowFiles[grep("arrow",ArrowFiles)]
ArrowFiles_1=paste("arrow.mouse/",ArrowFiles,sep="/")

addArchRThreads(threads = 1)

setwd("ATAC.res/") # The output directory
proj <- ArchRProject(
  ArrowFiles = ArrowFiles_1, 
  outputDirectory = "Results",
  copyArrows = FALSE
)

plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(1), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  seed=6,
  resolution = 1,
  corCutOff = 0.3
)

proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  corCutOff = 0.3
)

proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.3, 
  metric = "cosine"
)


plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")

proj <- addImputeWeights(proj)

markerGenes  <- c(
  "Afp",  #Early Progenitor
  "Alb", #Erythroid
  "Rhox5", "Apoc4", "Pvt1", #B-Cell Trajectory
  "Myc", "Stra6"
)
markerGenes <- "Myc"
plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(proj)
)
saveArchRProject(ArchRProj = proj, outputDirectory = "Results", load = FALSE)
######save
######
proj <- loadArchRProject("ATAC.res/Results/")


#####
ann <- data.frame(proj@cellColData)

#########
mat <- read.table("merge.xls") # The cut and ATAC cell pair metadata
mat$sample <- gsub("_BC.*","",mat$V2)
mat$cellall <- paste0(mat$sample,"#",mat$V2)
pos <- which(mat$cellall %in% rownames(ann))
mat1 <- mat[pos,]
mat1$sample <- gsub("_CUT_BC.*","",mat1$V3)
mat1$cellall1 <- paste0(mat1$sample,"#",mat1$V3)


cut <- readRDS("CUT.meta.data.rds") # The CUT meta data info
pos <- which(rownames(cut) %in% mat1$cellall1)
cut1 <- cut[pos,]

genescore <- readRDS("CUT.final.rds") # The CUT score 

############## Tumor subset
library(ggsci)
library(Matrix)
library(ggplot2)
library(Matrix)
library(ggpubr)
library(ggsignif)
col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")

ann <- data.frame(proj@cellColData)
pos <- which(ann$Clusters %in% "C3")
ann1 <- ann[pos,]

pos <- which(mat1$cellall %in% rownames(ann1))
mat2 <- mat1[pos,]

pos <- which(rownames(cut) %in% mat2$cellall1)
cut1 <- cut[pos,]
cut1$stage <- gsub("_.*","",cut1$Sample)
library(ggplot2)
ggplot(cut1,aes(x=stage,y=log2(nFrags)))+geom_boxplot()+labs(y="CUT.frags")+theme_classic()+
  theme_classic()+ geom_signif(comparisons =  list(c("AB2", "AB1"),c("AB2", "AB3")))

library(clusterProfiler)
tf <- read.gmt("s3_Hepa1.6.regulons.gmt") # The TF-targets files
pos <- grep("Rela",tf$term)
tf1 <- tf[pos,]

geenname <- tf1$gene

genescore <- readRDS("CUT.final.rds")
genescore <- as.matrix(genescore)
colnames(genescore) <- gsub("\\#",".",colnames(genescore))
pos <- which(rownames(genescore) %in% genename)
genescore1 <- genescore[pos,]
mittmp <- data.frame(cell = colnames(genescore),tf = colMeans(genescore1))

cut1$cellname <- gsub("\\#",".",rownames(cut1))
mittmp1 <- mittmp[as.character(cut1$cellname),]
cut1$tf.CUT <- mittmp1$tf
ggplot(cut1,aes(x=stage,y=as.numeric(as.character(cut1$tf.CUT))))+geom_boxplot()+labs(y="CUT.frags") +
  theme_classic()+ geom_signif(comparisons =  list(c("AB2", "AB1"),c("AB2", "AB3")))

genescore <- readRDS("ATAC.final.rds") # The ATAC score file
genescore <- as.matrix(genescore)
colnames(genescore) <- gsub("\\#",".",colnames(genescore))
pos <- which(rownames(genescore) %in% genename)
genescore1 <- genescore[pos,]
mittmp <- data.frame(cell = colnames(genescore),tf = colMeans(genescore1))

ann1$cellname <- gsub("\\#",".",rownames(ann1))
mittmp1 <- mittmp[as.character(ann1$cellname),]
ann1$tf.CUT <- mittmp1$tf
ann1$stage <- gsub("_.*","",ann1$Sample)
ggplot(ann1,aes(x=stage,y=as.numeric(as.character(ann1$tf.CUT))))+geom_boxplot()+labs(y="CUT.frags") +
  theme_classic()+ geom_signif(comparisons =  list(c("AB2", "AB1"),c("AB2", "AB3")))








