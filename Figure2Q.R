library(clusterProfiler)
library(org.Hs.eg.db)
library(survival)
library(survminer)
#######
args = commandArgs(T)
infile <- args[1] # The RNA count matrix
infile1 <- args[2] # The clinical matrix


########

RNAMats <- readRDS(infile)
clinicalMat <- readRDS(infile1)

RNAMat <- RNAMats[,rownames(clinicalMat)]

OS.p.best <- c()
OS.p.median <- c()
RFS.p.best <- c()
RFS.p.median <- c()

###################### a. Overall Survival
quit.m <- read.table("Stage3.DE.xls",sep = "\t",header = T,row.names = 1)
pos <-which(quit.m$avg_log2FC > 0.5 & quit.m$p_val_adj < 1e-10)
quit.m1 <- quit.m[pos,]
ge <- toupper(c(rownames(quit.m1),"CKB","INMT"))


pos <- which(ge %in% rownames(RNAMat))
ge <- ge[pos]
ov <- cbind(colMeans(RNAMat[ge,]),as.numeric(clinicalMat[,8]),as.numeric(clinicalMat[,6]),colMeans(RNAMat[c("ALB","APOC4","TTR"),]))
ov <- as.data.frame(ov)
names(ov) <- c("Gene","OS","month","ALB")
ov$RNA <- ov$Gene / ov$ALB
res.cut  <-  surv_cutpoint(ov, time = "month", event = "OS",variables = c("RNA"))
summary(res.cut)
res.cat <- surv_categorize(res.cut)
head(res.cat)

#### best
fit  <-  survfit(Surv(month, OS) ~RNA, data = res.cat)
OS.p.b <- surv_pvalue(fit, method = "survdiff")$pval
OS.p.best <- c(OS.p.best, OS.p.b)

print(ggsurvplot(fit,
                 pval = TRUE,
                 conf.int = F,
                 risk.table = TRUE, # Add risk table
                 #risk.table.col = "strata", # Change risk table color by groups
                 xlab="Time (months)", ylab="Overall survival probability",
                 title=paste0("HCC RNA of ",ge," in ",Type)
                 #linetype = "strata", # Change line type by groups
                 #surv.median.line = "hv", # Specify median survival
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 #palette = c("red","#2E9FDF")
))

######################
dat <- ov
dat <- dat[order(dat$RNA),]

dat3 <- dat

med2 <- median(dat3$RNA)

dat3$type <- "high"
pos <- which(dat3$RNA <= med2)
dat3$type[pos] <- "low"

surv_rnaseq.cat <- dat3[,c("month","OS","type")]
colnames(surv_rnaseq.cat)[3] <- "gene_score_norm"

fit <- survfit(Surv(month, OS) ~ gene_score_norm, data = surv_rnaseq.cat)

ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for
  # point estimaes of survival curves.
  #xlim = c(0,24),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 12,    # break X axis in time intervals by 500.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE, # show bars instead of names in text annotations
  xlab="Month (log-rank)",
  ylab="Overall survival (%)"
)

