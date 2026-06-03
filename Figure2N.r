library(Seurat)
library(ggplot2)
########
args = commandArgs(T)
infile <- args[1] ### The seurat file of metastasis data
outdir <- args[2]

lung <- readRDS(infile)


#########
quit.m <- read.table("Stage3.DE.xls",sep = "\t",header = T,row.names = 1)
pos <-which(quit.m$avg_log2FC > 0.5 & quit.m$p_val_adj < 1e-10)
quit.m1 <- quit.m[pos,]
ge <- toupper(c(rownames(quit.m1),"CKB","INMT"))

genename1 <- ge
pos <- which(quit.m$p_val_adj < 0.05 & quit.m$avg_log2FC < -0.4)
quit.m <- quit.m[pos,]
genename2 <- toupper(rownames(quit.m))
geneset <- list(up = genename1,dn = genename2)
DefaultAssay(lung) <- "RNA"
lung <- AddModuleScore(lung,features = geneset,name = names(geneset))

DefaultAssay(lung) <- "RNA"
lung <- AddModuleScore(lung,features = geneset,name = names(geneset))

########
#########
library(splines)
library(sfsmisc)
data <- lung$up1
data1 <- density(data)

threshold <- otsu_threshold(data)

plot(density(data))+abline(v = threshold)

data1 <- density(data)
# Find the "elbow" point (point of maximum curvature)
elbow_point <- data1$x[which(diff(sign(diff(data1$y))) == 2) + 1]

# Plot the density curve
plot(data1, main = "Density Plot with Elbow Detection")
abline(v = elbow_point, col = "red", lty = 2)
elbow_point
mean((elbow_point + threshold) /2)
mean((elbow_point ) /2)

data <- lung$up1
data2 <- sort(data)
sit <- length(data) * 0.995
data2[c(floor(sit),ceiling(sit))]
plot(density(data))+abline(v = c(data2[c(floor(sit),ceiling(sit))]))

data <- lung$dn2
data2 <- sort(data)
sit <- length(data) * 0.005
data2[c(floor(sit),ceiling(sit))]
plot(density(data))+abline(v = c(data2[c(floor(sit),ceiling(sit))]))

lung$dif <- lung$up1 - lung$dn2
data <- lung$dif
data2 <- sort(data)
sit <- length(data) * 0.995
data2[c(floor(sit),ceiling(sit))]
plot(density(data))+abline(v = c(data2[c(floor(sit),ceiling(sit))]))


########
pos1 <- which(lung$up1 > 0.5218479)
lung$clus <- "Bak"
lung$clus[pos1] <- "PHGDH"

mit1 <- lung@meta.data
write.table(mit1,file = paste0(outdir,"lung.meta.new.pro.xls"),sep = "\t",quote = F,row.names = F)
