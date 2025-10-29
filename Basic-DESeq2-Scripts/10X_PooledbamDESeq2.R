#!/usr/bin/env Rscript
# variables
cond1 <- "Male"
cond2 <- "UnprimedFemale"

output <- "DESeq2Scripts/Esr1Adult_Pooled10XBNSTMalevPrimed_"
siglfc <- 1		# |Log2| fold change to be used as significant. Set at 1 by default, meaning a 2 fold change (log2(2)=1)
sigpadj <- 0.05		# Adjusted p-value used as upper threshold for significance


# This script uses DESeq2 in R to run differential-expression analysis on a reads count matrix of genes and automatically generate some graphs.


# Load libraries
library(tximport)
library(rhdf5)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(calibrate)
library(pheatmap)

counts <- read.table("BNST_MalevPrimed_features.txt", row.names=1, header=TRUE)
counts <- as.matrix(counts, rownames=TRUE)
counts
cols <- data.frame(condition=factor(c(rep("Primed"),rep("Male"))))
cols
dds=DESeqDataSetFromMatrix(countData=counts, colData=cols, design =~condition)
dds <- DESeq(dds)

res <- results(dds)

write.csv(res, file=paste0(output, "results.csv"))

