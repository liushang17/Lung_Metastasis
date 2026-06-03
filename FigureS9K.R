library(Seurat)
library(ggplot2)
lung <- readRDS("/Users/shangliu/01.terms/01.lung_metastasis/05.Other/new/All.tumor.0117.rds")

###########
ann <- lung@meta.data
ann$clus <- "Bak"
tmpann <- read.table("/Users/shangliu/01.terms/01.lung_metastasis/11.sub/data1/liver.meta.new.xls",sep = " ",header = T)
pos <- which(tmpann$clus %in% "PHGDH")
tmpann1 <- tmpann[pos,]
pos <- which(rownames(ann) %in% rownames(tmpann1))
ann$clus[pos] <- "PHGDH"
tmpann <- read.table("/Users/shangliu/01.terms/01.lung_metastasis/11.sub/data1/lung.meta.xls",sep = "\t",header = T)
pos <- which(tmpann$clus %in% "PHGDH")
tmpann1 <- tmpann[pos,]
pos <- which(rownames(ann) %in% rownames(tmpann1))
ann$clus[pos] <- "PHGDH"
tmpann <- read.table("/Users/shangliu/01.terms/01.lung_metastasis/11.sub/data1/ova.meta.xls",sep = "\t",header = T)
pos <- which(tmpann$clus %in% "PHGDH")
tmpann1 <- tmpann[pos,]
pos <- which(rownames(ann) %in% rownames(tmpann1))
ann$clus[pos] <- "PHGDH"

lung@meta.data <- ann
########### Gene
Idents(lung) <- factor(lung$tumor,levels = c("liver","lung","ovarian"))

genename <- c("PHGDH","HOPX","DDX5","VAMP8","TXNIP","GSN","ITGB1","CKB","EFNA1","CAV1","RSRP1","SCD","CCL2","CXCL10","CXCL2","IRF1")
VlnPlot(lung,features = "ZEB1",split.by = "clus",pt.size = 0,ncol = 2)
VlnPlot(lung,features = "STMN1",split.by = "clus",pt.size = 0,log =F)

genename <- c("PHGDH","FAH","CKB","FILIP1L")
genename <- genetmp[-c(9,28,25,31,14,26,1)]
DotPlot(lung, features = genename, split.by = "clus") + RotatedAxis()
########## hallmark
# HALLMARK_INTERFERON_GAMMA_RESPONSE HALLMARK_G2M_CHECKPOINT 
#HALLMARK_TNFA_SIGNALING_VIA_NFKB HALLMARK_INTERFERON_ALPHA_RESPONSE HALLMARK_MYC_TARGETS_V1 HALLMARK_INFLAMMATORY_RESPONSE
ann <- lung@meta.data
genename <- "HALLMARK_G2M_CHECKPOINT"
#pos <- which(lung$HALLMARK_G2M_CHECKPOINT > 0.25)
#lung$HALLMARK_G2M_CHECKPOINT[pos] <- 0.25
#pos <- which(lung$HALLMARK_INTERFERON_GAMMA_RESPONSE > 0.3)
#lung$HALLMARK_INTERFERON_GAMMA_RESPONSE[pos] <- 0.3

#pdf("/Users/shangliu/02.report/01.lung_metastasis/Figure/Figure3/main.new.2/All.feature.pdf",width = 8,height = 20)
VlnPlot(lung,features = "HALLMARK_G2M_CHECKPOINT",split.by = "clus",pt.size = 0)+geom_boxplot()
#dev.off()
pos <- which(lung$tumor %in% "ovarian" & lung$clus %in% "PHGDH")
tmp1 <- lung$HALLMARK_INFLAMMATORY_RESPONSE[pos]
pos <- which(lung$tumor %in% "ovarian" & lung$clus %in% "Bak")
tmp2 <- lung$HALLMARK_INFLAMMATORY_RESPONSE[pos]
t.test(tmp1,tmp2,alternative = "less")
wilcox.test(tmp1,tmp2,alternative = "less")

#####
chomk <- read.table("/Users/shangliu/01.terms/01.lung_metastasis/00.basic/Chemokin_gene.list",sep = "\t")
#genename <- c("CXCL2","CXCL3","CCL2","CCL4","CCL5","CXCL10","CX3CL1","CCL8","CCL7")
genename <- toupper(c("Cxcl2","Cxcl5","Cxcl3","Ccl20","Ccl5","Ccl4","Cxcl10","Cxcl9","Ccl7","Ccl8","Ccl24","Ccl26","Ccl11","Ccl3","Cxcl11","Cxcl16","Ccl22","Ccl17"))

genename <- toupper(chomk$V1)
pos <- which(genename %in% "CCL28")
genename <- genename[-pos]
tmp <- list(chemokine = genename)
DefaultAssay(lung) <- "RNA"
lung <- AddModuleScore(lung,features = tmp,name = "chemokine")

pos <- which(rownames(lung@assays$RNA@data) %in% genename)
lung$chemokine2 <- rowMeans(t(as.matrix(lung@assays$RNA@data[pos,])))

pos <- which(lung$chemokine1 > 0.3)
lung$chemokine1[pos] <- 0.3

VlnPlot(lung,features = "chemokine2",group.by = "tumor",split.by = "clus",pt.size = 0,ncol = 1,log = T) + geom_boxplot()
pos <- which(lung$tumor %in% "ovarian" & lung$clus %in% "PHGDH")
tmp1 <- lung$chemokine2[pos]
pos <- which(lung$tumor %in% "ovarian" & lung$clus %in% "Bak")
tmp2 <- lung$chemokine2[pos]
t.test(tmp1,tmp2,alternative = "less")
wilcox.test(tmp1,tmp2,alternative = "less")


#####
library(clusterProfiler)
genes <- read.gmt("/Users/shangliu/01.terms/01.lung_metastasis/03.TME/14.dormancy/KEGG_GLYCINE_SERINE_AND_THREONINE_METABOLISM.v2023.1.Hs.gmt")
genes <- c("PHGDH","PSAT1","PSPH","SHMT1","SHMT2","MTHFD1","MTHFR","MTR","MT2A")
genest <- list(ser = genes)
lung <- AddModuleScore(lung,features = genest,name = names(genest))
VlnPlot(lung,features = "ser1",split.by = "clus",pt.size = 0) +geom_boxplot()
VlnPlot(lung,features = "STAT1",group.by = "tumor",split.by = "clus",pt.size = 0) 

######
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggprism)

ann <- lung@meta.data
mitall <- ann[,c("tumor","clus","HALLMARK_G2M_CHECKPOINT")]
colnames(mitall) <- c("type","condition","exp")
yst <- min(mitall$exp)
mitall$exp <- mitall$exp - yst
ggplot(mitall,aes(x=type,y = exp,color = condition))+
  geom_bar(stat="summary",fun=mean,position = position_dodge2(preserve = 'single', padding = 0.2),fill = "white")+ #绘制柱状图
  stat_summary(geom = "errorbar",fun.data  = 'mean_se', width = 0.3,position=position_dodge(0.95))+#误差棒
  labs(x="",y=NULL)+#标题
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45)+
  geom_signif(data=mitall,aes(xmin=0.75, xmax=1.25, annotations="**", y_position=0.08),
              textsize = 5, vjust = 0.05, tip_length = c(0, 0.0),
              manual=TRUE, color = "black")+
  geom_signif(data=mitall,aes(xmin=1.75, xmax=2.25, annotations="**", y_position=0.08),
              textsize = 5, vjust = 0.05, tip_length = c(0, 0.0),
              manual=TRUE, color = "black")+
  geom_signif(data=mitall,aes(xmin=2.75, xmax=3.25, annotations="**", y_position=0.08),
              textsize = 5, vjust = 0.05, tip_length = c(0, 0.0),
              manual=TRUE, color = "black")

pos <- which(lung$clus %in% "PHGDH")
ann1 <- lung@meta.data[pos,]
sui <- data.frame(table(ann1$tumor,ann1$Phase))
sui$Var2 <- factor(sui$Var2,levels = c("G1","S","G2M"))
ggplot(sui,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat = "identity",position = "fill")+
  theme_classic()

######
#########
de <- read.table("/Users/shangliu/01.terms/01.lung_metastasis/11.sub/data/Stage3.DE.xls",sep = "\t",header = T,row.names = 1)
go <- read.table("/Users/shangliu/01.terms/01.lung_metastasis/11.sub/data1/Stage3.go.txt",sep = "\t")

geneall <- NULL
for(i in c(1,5,7)){
  genename <- as.character(strsplit(go$V4[5],split = ",")[[1]])
  geneall <- c(geneall,genename)
}
geneall <- toupper(unique(geneall))

genest <- list(sti = toupper(geneall))
lung <- AddModuleScore(lung,features = genest,name = names(genest))
VlnPlot(lung,features = "sti1",group.by = "tumor",split.by = "clus",pt.size = 0,log = F) +geom_boxplot()
# 1,8,5,7
genename1 <- c("SUZ12","RBBP7","RBBP4","EED","EZH2","EZH1")
#VlnPlot(lung,features = genename1,group.by = "tumor",split.by = "clus",pt.size = 0)+geom_boxplot()
genename2 <- c("SETDB1","SUV39H1","SUV39H2","EHMT1","EHMT2","ESET")
#VlnPlot(lung,features = genename2,group.by = "tumor",split.by = "clus",pt.size = 0)+geom_boxplot()

genest <- list(k27me3 = genename1,k9me3 = genename2)
lung <- AddModuleScore(lung,features = genest,name = names(genest))
colnames(lung@meta.data)[(ncol(lung@meta.data)-1): ncol(lung@meta.data)] <- c("H3K27me3","H3K9me3")
VlnPlot(lung,features = "H3K9me3",split.by = "clus",pt.size = 0)+geom_boxplot()

pos <- which(lung$tumor %in% "liver" & lung$clus %in% "PHGDH")
tmp1 <- lung$H3K27me3[pos]
pos <- which(lung$tumor %in% "liver" & lung$clus %in% "Bak")
tmp2 <- lung$H3K27me3[pos]
t.test(tmp1,tmp2,alternative = "greater")
wilcox.test(tmp1,tmp2,alternative = "greater")


pos <- which(rownames(lung@assays$RNA@data) %in% genename1)
lung$H3K27me3.exp <- colMeans(as.matrix(lung@assays$RNA@data[pos,]))
pos <- which(rownames(lung@assays$RNA@data) %in% genename2)
lung$H3K9me3.exp <- colMeans(as.matrix(lung@assays$RNA@data[pos,]))


######
pos1 <- which(rownames(lung@assays$RNA@data) %in% "INMT")
pos2 <- which(lung$clus %in% "PHGDH" & lung$tumor %in% "lung")
tmp1 <- as.numeric(as.character(as.matrix(lung@assays$RNA@data[pos1,pos2])))
pos2 <- which(lung$clus %in% "Bak" & lung$tumor %in% "lung")
tmp2 <- as.numeric(as.character(as.matrix(lung@assays$RNA@data[pos1,pos2])))
t.test(tmp1,tmp2)
wilcox.test(tmp1,tmp2)

pos2 <- which(lung$clus %in% "PHGDH" & lung$tumor %in% "ovarian")
tmp1 <- as.numeric(as.character(as.matrix(lung@assays$RNA@data[pos1,pos2])))
pos2 <- which(lung$clus %in% "Bak" & lung$tumor %in% "ovarian")
tmp2 <- as.numeric(as.character(as.matrix(lung@assays$RNA@data[pos1,pos2])))
t.test(tmp1,tmp2)
wilcox.test(tmp1,tmp2)

pos2 <- which(lung$clus %in% "PHGDH" & lung$tumor %in% "liver")
tmp1 <- as.numeric(as.character(as.matrix(lung@assays$RNA@data[pos1,pos2])))
pos2 <- which(lung$clus %in% "Bak" & lung$tumor %in% "liver")
tmp2 <- as.numeric(as.character(as.matrix(lung@assays$RNA@data[pos1,pos2])))
t.test(tmp1,tmp2)
wilcox.test(tmp1,tmp2)

######
pos2 <- which(lung$clus %in% "PHGDH" & lung$tumor %in% "lung")
tmp1 <- as.numeric(as.character(as.matrix(lung$sti1[pos2])))
pos2 <- which(lung$clus %in% "Bak" & lung$tumor %in% "lung")
tmp2 <- as.numeric(as.character(as.matrix(lung$sti1[pos2])))
t.test(tmp1,tmp2,alternative = "less")
wilcox.test(tmp1,tmp2,alternative = "less")

pos2 <- which(lung$clus %in% "PHGDH" & lung$tumor %in% "ovarian")
tmp1 <- as.numeric(as.character(as.matrix(lung$sti1[pos2])))
pos2 <- which(lung$clus %in% "Bak" & lung$tumor %in% "ovarian")
tmp2 <- as.numeric(as.character(as.matrix(lung$sti1[pos2])))
t.test(tmp1,tmp2,alternative = "less")
wilcox.test(tmp1,tmp2,alternative = "less")


##########
ann <- lung@meta.data
pas <- unique(as.character(ann$patient_id))
mit <- matrix(nrow = length(pas),ncol = 2)
for(i in 1:length(pas)){
  pos1 <- which(ann$patient_id %in% pas[i] & ann$clus %in% "PHGDH")
  pos2 <- which(ann$patient_id %in% pas[i])
  mit[i,1] <- pas[i]
  mit[i,2] <- length(pos1) / length(pos2)
}

#######
library(clusterProfiler)
library(GSVA)
library(GSEABase)

geneset <- getGmt("/Users/shangliu/01.terms/01.lung_metastasis/00.basic/h.all.v7.2.symbols.gmt")
geneset <- geneIds(geneset)
names(geneset)<- gsub("HALLMARK_","",names(geneset))

DefaultAssay(lung) <- "RNA"
lung <- AddModuleScore(lung,features = geneset,name = names(geneset))
colnames(lung@meta.data)[(ncol(lung@meta.data)-49):ncol(lung@meta.data)] <- names(geneset)
genenames <- geneset$INFLAMMATORY_RESPONSE

#VlnPlot(lung,features = c("SRI","CDKN1A","ITGA5"),split.by = "clus",pt.size = 0,ncol = 4)+geom_boxplot()
#geneall <- strsplit("IL1R1 ITGB8 NAMPT TNFRSF10 IFNGR2 BST2 CD55 MYC RELA CD82 TIMP1 SRI NFKBIA ATP2B1 KLF6 TAPBP HIF1A NMI IFNAR1 DCBLD2 ADRM1 ATP2A2 BTG2 EIF2AK2 KIF1B SELENOS P2RX4 ABI1",split =" ")[[1]]
pos <- which(rownames(lung@assays$RNA@data) %in% genenames)
tmp <- colMeans(as.matrix(lung@assays$RNA@data[pos,]))
lung$exp <- tmp

VlnPlot(lung,features = "INFLAMMATORY_RESPONSE",split.by = "clus",group.by = "tumor",pt.size = 0) +geom_boxplot()

pos <- which(lung$tumor %in% "ovarian")
cellnames <- colnames(lung)[pos]
pbmc <-subset(lung,cells = cellnames)
clus.mak <- FindMarkers(pbmc,group.by = "clus",ident.1 = "PHGDH")
pos <- which(rownames(clus.mak) %in% geneall)
clus.mak[pos,]

pos <- which(lung$tumor %in% "ovarian" & lung$clus %in% "PHGDH")
tmp1 <- lung$exp[pos]
pos <- which(lung$tumor %in% "ovarian" & lung$clus %in% "Bak")
tmp2 <- lung$exp[pos]
t.test(tmp1,tmp2)
wilcox.test(tmp1,tmp2)
