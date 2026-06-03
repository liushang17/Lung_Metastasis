library(pheatmap)
library(Seurat)
library(ggplot2)
#######
args = commandArgs(T)
infile <- args[1] # The scRNA-seq seurat rds for Tumor cells
infile1 <- args[2] # The Scenic TF AUCell txt

######

lung <- readRDS(infile)
tfs <- read.table("s3_scRNA.AUCell.txt",sep= "\t",header = T,row.names = 1)
lung <- subset(lung,cells = rownames(tfs))
tfs <- tfs[colnames(lung),]
lung@meta.data <- cbind(lung@meta.data,tfs)
VlnPlot(lung,features = tfnames)

############
tfnames <- c("Irf1...","Irf7...","Irf3...","Stat1...","Nfkb2...","Rela...","Sp1...","Spi1...","Jun...")
genename.all <- c("Phgdh","Scd1","Rsrp1","Hopx","Vamp8","Filip1l","Cavin1","Txnip","Gsn",
             "Sparcl1","Ckb","Fah","Lbh","Ddx5","Hspg2","Chil3","Nkain2","Hexb","Nrp1","Inmt")
genename.all <- c("Phgdh","Hopx","Inmt","Sparcl1","Hexb","Gm48099","Filip1l","Nkain2","Cavin1","Gm19951","Scd1","Txnip",
                  "Ckb","AY036118","Chil3","Nrp1","Fah","Gsn","Mt-Nd2","H2-K1")
############
mitall <- NULL
for(m in 1:length(tfnames)){
  feas  <- tfnames[m]
  pos <- which(colnames(lung@meta.data) %in% feas)
  tmp <- as.numeric(as.character(lung@meta.data[,pos]))
  mat <- lung@assays$RNA@data
  
  pos <- which(rownames(mat) %in% genename.all)
  mat <- mat[pos,]
  siz <- NULL
  pva <- NULL
  nus <- NULL
  for(i in 1:nrow(mat)){
    tmp1 <- as.numeric(as.character(as.matrix(mat[i,])))
    pos <- which(tmp1 > 0)
    tmp1 <- tmp1
    tmp3 <- tmp
    if(length(pos) > 3){
      tmp2 <- cor.test(tmp1,tmp3,method = "pearson")
      siz <- c(siz,as.numeric(tmp2$estimate))
      pva <- c(pva,as.numeric(tmp2$p.value))
      nus <- c(nus,length(pos))
    }else{
      siz <- c(siz,0)
      pva <- c(pva,1)
      nus <- c(nus,length(pos))
    }
  }
  mit <- data.frame(genename = rownames(mat),cor = siz,pvalue = pva,num = nus)
  pos <- which(mit$cor %in% NA)
  mit1 <- mit
  pos <- which( mit$num > 10)
  mit1 <- mit[pos,]
  mit1$tfname <- tfnames[m]
  mitall <- rbind(mitall,mit1)
}

library(ggplot2)
sui <- matrix(nrow = length(genename.all),ncol = 2)
for(i in 1:length(genename.all)){
  pos <- which(mitall$genename %in% genename.all[i])
  mitall1 <- mitall[pos,]
  sui[i,1] <- genename.all[i]
  sui[i,2] <- median(mitall1$cor)
}
sui1 <- data.frame(sui)
sui1 <- sui1[order(-as.numeric(as.character(sui1$X2))),]
geneall <- as.character(sui1$X1)
geneall <- rev(geneall)
ggplot(mitall,aes(x=factor(genename,geneall),y = cor))+geom_boxplot()+geom_point()+theme_classic()+
  coord_flip()+labs(x="",y="")
mitall.tf <- mitall

