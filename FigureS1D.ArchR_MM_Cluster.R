### Get the parameters
parser = argparse::ArgumentParser(description="Script to Cluster scATAC data by ArchR")
parser$add_argument('-R','--ref', help='reference genome')
parser$add_argument('-I','--inputpath', help='input arrow directory')
parser$add_argument('-P','--id', help='tissue|select ID')
parser$add_argument('-O','--out', help='out directory')
parser$add_argument('-D','--filterRatio', help='filterRatio')
parser$add_argument('-X','--threads', help='threads')
args = parser$parse_args()

###
library("ArchR")
library(ggplot2)

if(args$ref=="mm10"){
addArchRGenome("mm10")}else{
addArchRGenome("mm9")
}

addArchRThreads(threads = as.numeric(args$threads))
######Creating Arrow Files

#################################################################################################################
#####------------------------------------------------ArchR--------------------------------------------------#####
#################################################################################################################

#########-----------------------------Creating Arrow Files---------------------------------#########
setwd(args$inputpath)
ArrowFiles <- list.files(args$inputpath)
ArrowFiles <- ArrowFiles[grep("arrow",ArrowFiles)]
ArrowFiles_1=paste(args$inputpath,ArrowFiles,sep="/")


setwd(args$out)

########------------------------------Creating an ArchRProject----------------------------########
proj <- ArchRProject(
  ArrowFiles = ArrowFiles_1, 
  outputDirectory = "Results",
  copyArrows = FALSE
)
#saveArchRProject(ArchRProj = proj, outputDirectory = "Save-proj", load = FALSE)


p2 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

p4 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
   
plotPDF(p2,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj, addDOC = FALSE, width = 4, height = 4)

#p1 <- plotFragmentSizes(ArchRProj = proj)
#p2 <- plotTSSEnrichment(ArchRProj = proj)
#plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


########-----------------------------Filtering Doublets from an ArchRProject-------------########
#proj <- filterDoublets(proj,filterRatio = as.numeric(args$filterRatio))
########-----------------------------Iterative Latent Semantic Indexing (LSI)------------########
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
########-----------------------------Batch Effect Correction wtih Harmony----------------########
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)
#######-----------------------------Clustering using Seurat’s FindClusters() function-----#######
proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    seed=6,
    resolution = 1,
    corCutOff = 0.3
)
######------------------------------Single-cell Embeddings--------------------------------#######
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

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")

plotPDF(p1,p2,p3,p4, name = paste(args$id,"Plot-UMAP2Harmony-Sample-Clusters.pdf",sep="_"), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
######------------------------------Marker Genes Imputation with MAGIC--------------------------------#######
#proj <- addImputeWeights(proj)
#saveArchRProject(ArchRProj = proj, outputDirectory = "Results", load = FALSE)
######------------------------------Identifying Marker Genes--------------------------------#######
markersGS <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0")

ml=as.data.frame(markerList)
write.table(ml,paste(args$id,"marker_list.xls",sep="_"),sep = "\t",quote = FALSE,row.names = FALSE)
######------------------------------Marker Genes Imputation with MAGIC--------------------------------#######
proj <- addImputeWeights(proj)
saveArchRProject(ArchRProj = proj, outputDirectory = "Results", load = FALSE)
