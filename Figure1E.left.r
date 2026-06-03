library(Seurat)
library(pheatmap)
library(monocle3)
library(ggplot2)
#################################################
cluster_rows = TRUE
hclust_method = "ward.D2" 
num_clusters = 21

hmcols = NULL

add_annotation_row = NULL
add_annotation_col = NULL
show_rownames = FALSE
use_gene_short_name = TRUE

norm_method = c("log", "vstExprs")
scale_max=3
scale_min=-3

trend_formula = '~sm.ns(Pseudotime, df=3)'

return_heatmap=FALSE
hmcols <- NULL


#################################################

table.ramp <- function(n, mid = 0.5, sill = 0.5, base = 1, height = 1)
{
  x <- seq(0, 1, length.out = n)
  y <- rep(0, length(x))
  sill.min <- max(c(1, round((n - 1) * (mid - sill / 2)) + 1))
  sill.max <- min(c(n, round((n - 1) * (mid + sill / 2)) + 1))
  y[sill.min:sill.max] <- 1
  base.min <- round((n - 1) * (mid - base / 2)) + 1
  base.max <- round((n - 1) * (mid + base / 2)) + 1
  xi <- base.min:sill.min
  yi <- seq(0, 1, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  xi <- sill.max:base.max
  yi <- seq(1, 0, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  height * y
}

#' @importFrom grDevices rgb
rgb.tables <- function(n,
                       red = c(0.75, 0.25, 1),
                       green = c(0.5, 0.25, 1),
                       blue = c(0.25, 0.25, 1))
{
  rr <- do.call("table.ramp", as.list(c(n, red)))
  gr <- do.call("table.ramp", as.list(c(n, green)))
  br <- do.call("table.ramp", as.list(c(n, blue)))
  rgb(rr, gr, br)
}

matlab.like <- function(n) rgb.tables(n)

matlab.like2 <- function(n)
  rgb.tables(n,
             red = c(0.8, 0.2, 1),
             green = c(0.5, 0.4, 0.8),
             blue = c(0.2, 0.2, 1))

blue2green2red <- matlab.like2


################################################
args = commandArgs(T)
infile <- args[1] # monocle3 rds
outdir <- args[2]

Combine <- readRDS(infile)
degene <- read.table(paste0(outdir,"/Combine.monocle.de.rds"),sep = "\t",header = T)

pos <- which(rownames(Combine@assays$RNA@data) %in% degene$gene_id)
mat <- data.frame(t(Combine@assays$RNA@data[pos,]))
mat$pro <- Combine$pseu
mat <- mat[order(as.numeric(mat$pro)),]

mit <- matrix(nrow = (nrow(mat)),ncol = (ncol(mat)-1))
for(i in 1:(ncol(mat)-1)){
  y <- as.numeric(as.character(as.matrix(mat[,i])))
  x <- as.numeric(mat$pro)
  lo <- loess(y~x,span = 0.5)
  tmp <- predict(lo)
  mit[,i] <- tmp
}
colnames(mit) <- colnames(mat)[1:(ncol(mat)-1)]
rownames(mit) <- rownames(mat)

mit1 <- t(mit)

############
m <- mit1
m=m[!apply(m,1,sum)==0,]
m=m[!apply(m,1,sd)==0,]
m=Matrix::t(scale(Matrix::t(m),center=TRUE))
m=m[is.na(row.names(m)) == FALSE,]
m[is.nan(m)] = 0
m[m>scale_max] = scale_max
m[m<scale_min] = scale_min

heatmap_matrix <- m

### features show
degene1 <- degene[order(degene$Status,degene$GeneName),]
pos <- which(degene1$GeneName %in% rownames(m))
degene2 <- degene1[pos,]
degene3 <- data.frame(cluster = degene2$Status)
rownames(degene3) <- degene2$GeneName
heatmap_matrix1 <- heatmap_matrix[rownames(degene3),]
pdf(paste0(outdir,"/test.pdf"))
pheatmap::pheatmap(heatmap_matrix1,annotation_row = degene3,show_rownames = F,cluster_rows = F,cluster_cols = F,show_colnames = F,gaps_col = c(2562,12635,13484)) # 2421 10934 12491
dev.off()
#### features re-clustering

row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1

if(is.null(hmcols)) {
  bks <- seq(-3.1,3.1, by = 0.1)
  hmcols <- blue2green2red(length(bks) - 1)
}else {
  bks <- seq(-3.1,3.1, length.out = length(hmcols))
} 

ph <- pheatmap(heatmap_matrix, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=cluster_rows, 
               show_rownames=F, 
               show_colnames=F, 
               clustering_distance_rows=row_dist,
               clustering_method = hclust_method,
               cutree_rows=num_clusters,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols)

if(cluster_rows) {
  annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
} else {
  annotation_row <- NULL
}

if(!is.null(add_annotation_row)) {
  old_colnames_length <- ncol(annotation_row)
  annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])  
  colnames(annotation_row)[(old_colnames_length+1):ncol(annotation_row)] <- colnames(add_annotation_row)
}

if(!is.null(add_annotation_col)) {
  if(nrow(add_annotation_col) != 100) {
    stop('add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!')
  }
  annotation_col <- add_annotation_col
} else {
  annotation_col <- NA
}


feature_label <- row.names(heatmap_matrix)
if(!is.null(annotation_row)) row_ann_labels <- row.names(annotation_row)


row.names(heatmap_matrix) <- feature_label
if(!is.null(annotation_row)) row.names(annotation_row) <- row_ann_labels

colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))

ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                   useRaster = T,
                   cluster_cols = FALSE, 
                   cluster_rows = cluster_rows, 
                   show_rownames=T, 
                   show_colnames=F, 
                   #scale="row",
                   clustering_distance_rows=row_dist, #row_dist
                   clustering_method = hclust_method, #ward.D2
                   cutree_rows=num_clusters,
                   # cutree_cols = 2,
                   annotation_row=annotation_row,
                   annotation_col=annotation_col,
                   treeheight_row = 20, 
                   breaks=bks,
                   fontsize = 6,
                   color=hmcols, 
                   border_color = NA,
                   silent=TRUE,
                   filename=NA
)
print(ph_res)

ann1 <- annotation_row
ann1$cluster <- ann1$Cluster
clus <- c(15,11,7,12,4,20,18,17,3,16,19,10,13,5,14,9,8,6,21,1,2)
mittmp <- NULL
for(i in 1:length(clus)){
    pos <- which(ann1$cluster %in% clus[i])
    annt <- ann1[pos,]
    mittmp <- rbind(mittmp,annt)
}
heatmap_matrix1 <- heatmap_matrix[rownames(mittmp),]
annotation_row <- data.frame(Cluster = factor(mittmp$Cluster,levels = clus))
rownames(annotation_row) <- rownames(mittmp)
ph_res <- pheatmap(heatmap_matrix1[, ],
                   useRaster = T,cluster_cols = FALSE, 
                   cluster_rows = F,show_rownames=T, 
                   show_colnames=F, 
                   #scale="row",
                   clustering_distance_rows=row_dist, #row_dist
                   clustering_method = hclust_method, #ward.D2
                   cutree_rows=num_clusters,
                   # cutree_cols = 2,
                   annotation_row=annotation_row,
                   annotation_col=annotation_col,
                   treeheight_row = 20, 
                   breaks=bks,
                   fontsize = 6,
                   color=hmcols, 
                   border_color = NA,
                   silent=TRUE,
                   filename=NA
)
pdf(paste0(outdir,"/Spatial.TF.along.pes.pdf"))
print(ph_res)
dev.off()
write.table(annotation_row,file = paste0(outdir,"/Spatial.TF.along.cluster.xls"),sep = "\t",quote=F)

#####
gene_clus <- read.table(paste0(outdir,"/DE.along.final.xls"),sep = "\t",header = T)
pos <- which(gene_clus$GeneName %in% pbmc.markers2$gene_id)
gene_clus1 <- gene_clus[pos,]
gene_clus1$Status <- "Phase2"
pos <- which(gene_clus1$Cluster %in% c(7,11,15))
gene_clus1$Status[pos] <- "Phase1"
pos <- which(gene_clus1$Cluster %in% c(5,10))
gene_clus1$Status[pos] <- "Phase3"
pos <- which(gene_clus1$Cluster %in% c(1,2,6,8,9,13,14,21))
gene_clus1$Status[pos] <- "Phase4"
table(gene_clus1$Status)
write.table(gene_clus1,file = paste0(outdir,"/Spatial.TF.along.cluster.filter.xls"),sep = "\t",quote=F,row.names = F)




