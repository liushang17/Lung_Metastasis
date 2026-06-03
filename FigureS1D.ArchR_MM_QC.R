### Get the parameters
parser = argparse::ArgumentParser(description="Script to QC scATAC data by ArchR")
parser$add_argument('-R','--ref', help='reference genome')
parser$add_argument('-I','--inputFiles', help='input fragment file')
parser$add_argument('-T','--filtertss', help='tss threhold')
parser$add_argument('-F','--filterFrags', help='fragment number threhold')
parser$add_argument('-O','--out', help='out directory')
args = parser$parse_args()

###
library("ArchR")
library(ggplot2)

set.seed(1)
addArchRThreads(threads = 1)

if(args$ref=="mm10"){
addArchRGenome("mm10")}else{
addArchRGenome("mm9")
}

#################################################################################################################
#####------------------------------------------------ArchR--------------------------------------------------#####
#################################################################################################################

inputFiles <- args$inputFiles
names(inputFiles) <- gsub("\\.fragments\\.tsv\\.gz","",tail(unlist(strsplit(inputFiles,"/")),1))

if(!file.exists(args$out)){dir.create(args$out)}
setwd(args$out)

#########-----------------------------Creating Arrow Files---------------------------------#########
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = as.numeric(args$filtertss), #Dont set this too high because you can always increase later
  minFrags = as.numeric(args$filterFrags), 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

#########-----------------------------Inferring scATAC-seq Doublets with ArchR-------------#########
addArchRThreads(threads = 1)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

