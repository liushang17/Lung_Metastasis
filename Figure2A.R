library(Seurat)
library(ggsci)
library(ggrepel)
library(dplyr)
library(harmony)
library(ggsignif)
library(ggprism)
col <- c(pal_npg("nrc")(10)[1:6],"#FAFD7CFF","#FF6F00FF",pal_lancet("lanonc")(9)[c(1,3,7)],"#660099FF","#B5CF6BFF","#B24745FF","#CCFF00FF",
         "#FFCD00FF","#800000FF","#20854EFF","#616530FF","#FF410DFF","#EE4C97FF","#FF1463FF","#00FF00FF","#990080FF","#00FFFFFF",
         "#666666FF","#CC33FFFF","#00D68FFF","#4775FFFF","#C5B0D5FF","#FDAE6BFF","#79CC3DFF","#996600FF","#FFCCCCFF","#0000CCFF",
         "#7A65A5FF","#1A5354FF","#24325FFF")
args = commandArgs(T)

infile <- args[1]

Combine <- readRDS(infile)
Combine$time <- gsub("Casp8_","",Combine$time)
Combine$time <- gsub("_.*","",Combine$time)
Combine$slide <- gsub(":.*","",Combine$allbin)
Combine$time <- factor(Combine$time,levels = paste0("G",c(1,2,4,5,6,7,9,11,12)))

ann <- Combine@meta.data
ann$immune <- rowSums(ann[,11+c(3:7,9,13:18,20:22)])

ann$stage <- "Phase1"
pos <- which(ann$pseudotime > 5)
ann$stage[pos] <- "Phase2"
pos <- which(ann$pseudotime > 24.3)
ann$stage[pos] <- "Phase3"
pos <- which(ann$pseudotime > 28.4)
ann$stage[pos] <- "Phase4"

Combine$stage <- ann$stage
Combine$immune <- ann$immune

DimPlot(Combine,group.by = "stage")
DimPlot(Combine,group.by = "seurat_clusters")

mit <- data.frame(table(ann$slide,ann$time))
pos <- which(mit$Freq == 0)
mit <- mit[-pos,]
mit$num <- 0
mit$pro <- 0

for(i in 1:nrow(mit)){
  pos <- which(ann$stage %in% "Phase3" & ann$slide %in% mit$Var1[i])
  mit$num[i] <- length(pos)
  mit$pro[i] <- length(pos) / mit$Freq[i]
}

mit <- read.table("/Users/shangliu/01.terms/01.lung_metastasis/11.sub/data/Phase3_cellpro.csv",sep = ",",header = T)

mit$Var2 <- factor(mit$Var2,levels = paste0("G",c(1,2,4,5,6,7,9,11,12)))
ggplot(mit,aes(x=Var2,y = pro,color = "white"))+
  geom_bar(stat="summary",fun=mean,position = position_dodge2(preserve = 'single', padding = 0.2),fill = "red")+ #绘制柱状图
  labs(x="",y=NULL)+#标题
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45)

posn.d <- position_dodge(width=0.7)
# Function for median and IQR
median_IQR <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}

mead_sd <- function(x) {
  data.frame(y = mean(x),
             ymin = mean(x)-sd(x), # 1st quartile
             ymax = mean(x)+sd(x))  # 3rd quartile
}
mead_sem <- function(x) {
  data.frame(y = mean(x),
             ymin = mean(x)-sd(x)/sqrt(length(x)), # 1st quartile
             ymax = mean(x)+sd(x)/sqrt(length(x)))  # 3rd quartile
}

ggplot(mit,aes(x=Var2,y = Total_pct,color = "white"))+
  #geom_jitter(width = 0.1,shape=2, color="#C44E52") +
  geom_smooth(method="loess",se=FALSE)+
  stat_summary(geom = "linerange",fun.data = mead_sem, position = posn.d, color="#C44E52") +
  stat_summary(fun = mean, geom = 'line', color="#C44E52")+
  stat_summary(fun = mean, size = 2, geom = "point",  color="#C44E52")+
#  geom_line(size=1) +
  labs(y="Num of Phase3",x="")+
  theme_classic()+theme( panel.grid= element_blank(),plot.title = element_text(size = 15, hjust=0.5))

ggplot(mit,aes(x=Var2,y = pro,color = "white"))+
  geom_bar(stat="summary",fun.data=mead_sem,position = posn.d,fill = "red")+ #绘制柱状图
  stat_summary(geom = "linerange",fun.data = mead_sd, position = posn.d, color="#C44E52") +
  stat_summary(fun = mean, geom = 'line', color="#C44E52")+
  stat_summary(fun = mean, size = 2, geom = "point",  color="#C44E52")+
  #geom_line(size=1) +
  labs(y="Num of Phase3",x="")+
  theme_classic()+theme( panel.grid= element_blank(),plot.title = element_text(size = 15, hjust=0.5))

ggplot(mit,aes(x=Var2,y = pro,color = "white"))+
  geom_bar(stat="summary",fun=mean,position = position_dodge2(preserve = 'single', padding = 0.2),fill = "white")+ #绘制柱状图
  stat_summary(geom = "errorbar",fun.data = 'mead_sem', width = 0.3,position=position_dodge(0.95))+#误差棒
  labs(x="",y=NULL)+#标题
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45)

ggplot(mit, aes(x=Var2, y=Total_pct,group = 1)) +
  #geom_jitter(width = 0.1,shape=0, color="#3C5488") +
  # geom_smooth(method="loess",se=FALSE)+
  stat_summary(geom = "linerange",fun.data = mead_sem, position = posn.d,color="#3C5488") +
  stat_summary(fun = median, geom = 'line', color="#3C5488")+
  stat_summary(fun = median, size = 2, geom = "point",  color="#3C5488")+
  #geom_line(size=1) +
  labs(y="Pct of tumor cells",x="")+
  theme_classic()+theme( panel.grid= element_blank(),plot.title = element_text(size = 15, hjust=0.5))+ylim(1,5.5)

############
