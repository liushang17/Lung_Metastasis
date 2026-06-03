library(Seurat)
lung <- readRDS("/Users/shangliu/01.terms/01.lung_metastasis/11.sub/data13/C04493G4.tumor.rds")

lung <- NormalizeData(lung)

#########
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
c