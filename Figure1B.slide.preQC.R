library(dplyr)
library(Seurat)
library(data.table)
library(Matrix)
library(rjson)
library(ggplot2)


args <- commandArgs(T)
infile <- args[1]
bs <- as.numeric(args[2])
outdir <- args[3]

pro <- tail(unlist(strsplit(infile,"/")),1)
pro <- gsub(".tissue.gem.gz|.txt|_filterd|.bin1.Lasso.gene|.bin1.Gene_Expression_table.tsv|.bin1.lasso.gem.gz|_bin50_information.txt.gz","",pro)

outdir <- paste0(outdir,"/",pro)
dir.create(outdir,recursive = T)
setwd(outdir)
system(paste0("mkdir"," ","figures"))

############################## 1. bin data  ##############################
dat <- fread(file = infile)
if(length(grep("MIDCounts|MIDCount",colnames(dat))>0)){
  colnames(dat) <- gsub("MIDCounts|MIDCount","UMICount",colnames(dat))}
out <- as.data.frame(dat)

dat$x <- trunc((dat$x - min(dat$x)) / bs + 1)
dat$y <- trunc((max(dat$y) - dat$y) / bs + 1)

out <- cbind(dat$y,dat$x,out)
colnames(out)[1:2] <- c(paste0("bin",bs,".y"),paste0("bin",bs,".x"))
fwrite(out,paste0(pro,"_bin",bs,"_information.txt"),col.names=T,row.names=F,sep="\t",quote=F)

dat <- dat[, sum(UMICount), by = .(geneID, x, y)]
dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x
bin.coor <- dat[, sum(V1), by = .(x, y)]

out <- as.data.frame(cbind(paste0("BIN.",unique(dat$bin_ID)),bin.coor$y,bin.coor$x))
colnames(out) <- c(paste0("BIN.",bs),paste0("bin",bs,".y"),paste0("bin",bs,".x"))
rownames(out) <- out[,1]
fwrite(out,paste0(pro,"_bin",bs,"_position.txt"),col.names=T,row.names=F,sep="\t",quote=F)

## geneID
geneID <- seq(length(unique(dat$geneID)))
hash.G <- data.frame(row.names = unique(dat$geneID), values = geneID)
gen <- hash.G[dat$geneID, 'values']
## bin_ID
bin_ID <- unique(dat$bin_ID)
hash.B <- data.frame(row.names = sprintf('%d', bin_ID), values = bin_ID)
bin <- hash.B[sprintf('%d', dat$bin_ID), 'values']
##
cnt <- dat$V1

##
rm(dat)
gc()

##
tissue_lowres_image <- matrix(1, max(bin.coor$y), max(bin.coor$x))

tissue_positions_list <- data.frame(row.names = paste('BIN', rownames(hash.B), sep = '.'),
                                    tissue = 1,
                                    row = bin.coor$y, col = bin.coor$x,
                                    imagerow = bin.coor$y, imagecol = bin.coor$x)

scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1,
                                 tissue_hires_scalef = 1,
                                 tissue_lowres_scalef = 1))


##
mat <- sparseMatrix(i = gen, j = bin, x = cnt)
rownames(mat) <- rownames(hash.G)
colnames(mat) <- paste('BIN', sprintf('%d', seq(max(hash.B[, 'values']))), sep = '.')

##
generate_spatialObj <- function (image, scale.factors, tissue.positions, filter.matrix = TRUE)
{
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
  }
  
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  
  spot.radius <- unnormalized.radius / max(dim(image))
  
  return(new(Class = 'VisiumV1',
             image = image,
             scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                          fiducial = scale.factors$fiducial_diameter_fullres,
                                          hires = scale.factors$tissue_hires_scalef,
                                          lowres = scale.factors$tissue_lowres_scalef),
             coordinates = tissue.positions,
             spot.radius = spot.radius))
}

spatialO <- generate_spatialObj(image = tissue_lowres_image,
                                scale.factors = fromJSON(scalefactors_json),
                                tissue.positions = tissue_positions_list)



############################## 2. creat Spatial Object  ##############################
seuratSpaObj <- CreateSeuratObject(mat, project = 'Spatial', assay = 'Spatial',min.cells=5, min.features=5)

##
spatialO <- spatialO[Cells(seuratSpaObj)]
DefaultAssay(spatialO) <- 'Spatial'
seuratSpaObj[['slice1']] <- spatialO


rm(mat)
rm(bin.coor)
rm(hash.G)
rm(hash.B)
rm(bin)
rm(gen)
rm(cnt)

##############################  3. Spatial Analyse  ##############################
seuratSpaObj[["percent.mt"]] <- PercentageFeatureSet(seuratSpaObj, pattern = "^MT-")

Q1 <- quantile(seuratSpaObj$nFeature_Spatial)[2]
Q3 <- quantile(seuratSpaObj$nFeature_Spatial)[4]
upper <- as.numeric(Q3+1.5*(Q3-Q1))
lower <- as.numeric(Q1-1.5*(Q3-Q1))

pdf(paste0("figures/",pro,"_bin",bs,"_preQC.pdf"),width=8,height=8)
p1 <- VlnPlot(seuratSpaObj, features=c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol=3, pt.size=0)+theme(axis.text.x=element_text(angle=0,size=4),axis.title.x=element_text(angle=30,size=10))+labs(x=paste0("nGene:",dim(seuratSpaObj)[1],"; ","nBIN:",dim(seuratSpaObj)[2]))
print(p1)
plot(density(seuratSpaObj$nFeature_Spatial))+abline(v=c(150,200,300,500),col="grey",lwd=2,lty=6)+abline(v=lower,col="blue",lwd=2,lty=6)+abline(v=upper, col="red",lwd=2,lty=6)
p2 <- ggplot(seuratSpaObj@meta.data,aes(x=nFeature_Spatial)) +geom_density(colour="black") + theme_classic() + theme(plot.title=element_text(hjust=0.5,size=18, face="bold.italic"), legend.position="none",axis.title=element_text(size=15, face="bold.italic"),axis.text.x=element_text(size=12),axis.ticks.x=element_blank()) + geom_vline(aes(xintercept=100,colour="#999999",linetype="twodash")) + geom_vline(aes(xintercept=200,colour="#999999",linetype="twodash"))+geom_vline(aes(xintercept=300,colour="#999999",linetype="twodash"))+geom_vline(aes(xintercept=500,colour="#999999",linetype="twodash")) + geom_vline(aes(xintercept=lower, colour="#377EB8",linetype="twodash")) + geom_vline(aes(xintercept=upper, colour="#E41A1C", linetype="twodash"))+xlim(min(seuratSpaObj@meta.data$nFeature_Spatial),max(seuratSpaObj@meta.data$nFeature_Spatial)) + ggtitle(paste0(pro,".nBIN_",bs,":",dim(seuratSpaObj@meta.data)[1]))
print(p2)
dev.off()

outinfor <- data.frame(seuratSpaObj@meta.data,fig.X = seuratSpaObj@images$slice1@coordinates$col, fig.Y= (max(seuratSpaObj@images$slice1@coordinates$row)-seuratSpaObj@images$slice1@coordinates$row),bin.x=seuratSpaObj@images$slice1@coordinates$col,bin.y=seuratSpaObj@images$slice1@coordinates$row)
outinfor <- outinfor[,-c(1)]
fwrite(outinfor,paste0(pro,"_bin",bs,"_cellcluster_Axis.txt"),col.names=T,row.names=T,sep="\t",quote=F)

saveRDS(seuratSpaObj,file=paste0(pro,"_bin",bs,"_preQC.rds"))
saveRDS(seuratSpaObj@assays$Spatial@counts,file=paste0(pro,"_bin",bs,"_counts_mat.rds"))

