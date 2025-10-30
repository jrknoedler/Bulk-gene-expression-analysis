#!/usr/bin/env Rscript
# variables
cond1 <- "Primed"
cond2 <- "Unprimed"

output <- "DESeq2Scripts/Esr1Adult_VMH_GLM_4way_dropY_droplow_dropmt_updatedcastrate"
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

dir <- "/scratch/users/knoedler/Kallisto_ProteinRef/"
samples <- read.table("Kallistosheets/AVPV.tsv", header=TRUE)
samples
files <- file.path(dir, samples$sample, "abundance.h5")
files
names(files) <- paste0(samples$name, 1:4)
tx2gene <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/vM21_Tx2Gene_Updated.txt")
head(tx2gene)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

ychromgenes <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt")
ychromgenes <- unlist(ychromgenes)
head(ychromgenes)

sampleTable <- data.frame(condition = factor(c(rep(cond1, 2), rep(cond2, 2))))
rownames(sampleTable) <- colnames(txi.kallisto$counts)

write.csv(txi.kallisto, file=paste0(output,"_import.csv"))
# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromTximport(txi.kallisto, sampleTable, design =~condition)
dds <- dds[!(rownames(dds) %in% ychromgenes),]
keep <- rowSums(counts(dds) >= 5) >=3
dds <- dds[keep,]

dds <- DESeq(dds)
# run DESeq2
resultsNames(dds)
# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)


# Save differential expression results
resMalevPrimed <- results(dds, contrast=c("condition",cond1,cond2))
table(resMalevPrimed$padj<sigpadj)
## Order by adjusted p-value
resMalevPrimed <- resMalevPrimed[order(resMalevPrimed$padj), ]
## Merge with normalized count data
resMalevPrimeddata <- merge(as.data.frame(resMalevPrimed), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resMalevPrimeddata)[1] <- "Gene"
## Write results
write.csv(resMalevPrimeddata, file=paste0(output, "_PrimedvUnprimed-ExpDiff.csv"))


