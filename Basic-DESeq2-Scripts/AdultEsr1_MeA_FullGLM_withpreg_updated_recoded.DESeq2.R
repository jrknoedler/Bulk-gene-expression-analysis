#!/usr/bin/env Rscript
# variables
cond1 <- "Male"
cond2 <- "Virgin"
cond3 <- "Pregnant"
output <- "DESeq2Scripts/Esr1Adult_MeA_GLM_Preg4way_droplowdropYdropmt_GLM_updated_recoded"
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

dir <- "/scratch/users/knoedler/Joe_Sherlock_Backup_17/Joe_Sherlock_Backup/Kallisto_Output/Kallisto_ProteinRef/"
samples <- read.table("Joe_Sherlock_Backup_17/Joe_Sherlock_Backup/Data_TextFiles/Kallistosheets/MeA_Preg_FULL.tsv", header=TRUE)
samples
files <- file.path(dir, samples$sample, "abundance.tsv")
files
names(files) <- paste0(samples$name, 1:12)
tx2gene <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/vM21_Tx2Gene_Updated.txt")
head(tx2gene)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

write.csv(txi.kallisto, file=paste0(output,"_import.csv"))
ychromgenes <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt")
ychromgenes <- unlist(ychromgenes)
head(ychromgenes)

sampleTable <- data.frame(condition = factor(c(rep(cond1, 3), rep(cond2, 6), rep(cond3, 3))))
rownames(sampleTable) <- colnames(txi.kallisto$counts)


# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromTximport(txi.kallisto, sampleTable, design =~condition)
dds <- dds[!(rownames(dds) %in% ychromgenes),]
keep <- rowSums(counts(dds) >= 5) >=9
dds <- dds[keep,]

dds <- DESeq(dds)
# run DESeq2
resultsNames(dds)
# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)


# Save differential expression results
resAll <- results(dds)
table(resAll$padj<sigpadj)
## Order by adjusted p-value
resAll <- resAll[order(resAll$padj), ]
## Merge with normalized count data
resAlldata <- merge(as.data.frame(resAll), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resAlldata)[1] <- "Gene"
## Write results
write.csv(resAlldata, file=paste0(output, "_All-ExpDiff.csv"))
resPregvVirgin <- results(dds, contrast=c("condition",cond2,cond3))
table(resPregvVirgin$padj<sigpadj)
## Order by adjusted p-value
resPregvVirgin <- resPregvVirgin[order(resPregvVirgin$padj), ]
## Merge with normalized count data
resPregvVirgindata <- merge(as.data.frame(resPregvVirgin), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resPregvVirgin)[1] <- "Gene"
## Write results
write.csv(resPregvVirgindata, file=paste0(output, "_VirginvPreg-ExpDiff.csv"))
