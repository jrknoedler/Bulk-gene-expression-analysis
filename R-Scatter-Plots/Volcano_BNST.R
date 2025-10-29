library(EnhancedVolcano)
library(tidyverse)
MvP <- read.table("Tac1cluster.txt", header=TRUE, sep="\t", row.names=1) 
head(MvP)
data2 <- read.table("BNST_MvU_Volcano.txt", header=TRUE)
data3 <- read.table("BNST_PvU_Volcano.txt", header=TRUE)

keyvals1 <- ifelse(MvP$log2fc < 0 & MvP$padj < 0.05, 'red', ifelse(MvP$log2fc > 0 & MvP$padj < 0.05, 'blue', 'gray')) 
keyvals1[is.na(keyvals1)] <- 'gray'
names(keyvals1)[keyvals1=='red'] <- 'Fr' 
names(keyvals1)[keyvals1=='blue'] <- 'M' 
names(keyvals1)[keyvals1=='gray'] <- 'NC' 
head(keyvals1)
selectLab <- c("Tac1","Etl4","Socs2","Rfx3")
pdf(file="BNSTcluster16_MvP_volcano.pdf")
EnhancedVolcano(MvP, lab =rownames(MvP), x='log2fc', y='padj', ylim=c(0,8), selectLab=selectLab, gridlines.major=FALSE, labSize=6, labvjust=-3, labhjust=-4, labCol="black", gridlines.minor=FALSE, colCustom=keyvals1, pCutoff=NA, FCcutoff=NA, drawConnectors=TRUE, widthConnectors=0.75, colConnectors="black")
dev.off()


keyvals2 <- ifelse(data2$log2FoldChange < -0.584962501 & data2$padj < 0.05, 'green', ifelse(data2$log2FoldChange > 0.584962501 & data2$padj < 0.05, 'blue', 'gray')) 
keyvals2[is.na(keyvals2)] <- 'gray'
names(keyvals2)[keyvals2=='green'] <- 'Fu' 
names(keyvals2)[keyvals2=='blue'] <- 'M' 
names(keyvals2)[keyvals2=='gray'] <- 'NC' 
head(keyvals2)


keyvals3 <- ifelse(data3$log2FoldChange < -0.584962501 & data3$padj < 0.05, 'green', ifelse(data3$log2FoldChange > 0.584962501 & data3$padj < 0.05, 'red', 'gray')) 
keyvals3[is.na(keyvals3)] <- 'gray'
names(keyvals3)[keyvals3=='green'] <- 'Fu' 
names(keyvals3)[keyvals3=='red'] <- 'Fr' 
names(keyvals3)[keyvals3=='gray'] <- 'NC' 
head(keyvals3)

pdf(file="BNST_MvP_volcano.pdf")
EnhancedVolcano(MvP, lab = NA, x='log2FoldChange', y='padj', gridlines.major=FALSE, gridlines.minor=FALSE, colCustom=keyvals1, selectLab=c(), pCutoff=NA, FCcutoff=0.584962501)
dev.off()

pdf(file="BNST_MvU_volcano.pdf")
EnhancedVolcano(data2, lab = NA, x='log2FoldChange', y='padj', gridlines.major=FALSE, gridlines.minor=FALSE, colCustom=keyvals2, selectLab=c(), pCutoff=NA, FCcutoff=0.584962501)
dev.off()

pdf(file="BNST_PvU_volcano.pdf")
EnhancedVolcano(data3, lab = NA, x='log2FoldChange', y='padj', gridlines.major=FALSE, gridlines.minor=FALSE, colCustom=keyvals3, selectLab=c(), pCutoff=NA, FCcutoff=0.584962501)
dev.off()