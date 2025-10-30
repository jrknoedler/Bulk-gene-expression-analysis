#!/usr/bin/env Rscript
# variables
cond1 <- "Male"
cond2 <- "Primed"
cond3 <- "Unprimed"
cond4 <- "CastrateMale"
cond5 <- "PregnantFemale"
output <- "DESeq2Scripts/Esr1Adult_POA_GLM_4way_dropY_droplow_dropmt_updatedcastratewithpreg_gonadsexsplit"
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
samples <- read.table("Kallistosheets/AdultEsr1Kallisto_POAFULLcaspregupdated.tsv", header=TRUE)
samples
files <- file.path(dir, samples$sample, "abundance.h5")
files
names(files) <- paste0(samples$name, 1:15)
tx2gene <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/vM21_Tx2Gene_Updated.txt")
head(tx2gene)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

ychromgenes <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt")
ychromgenes <- unlist(ychromgenes)
head(ychromgenes)

sampleTable <- data.frame(condition = factor(c(rep(cond1, 3), rep(cond2, 3), rep(cond3, 3), rep(cond4, 3), rep(cond5, 3))), batch=factor(c(rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("4",5), rep("5",1))), sex=factor(c(rep("male",3), rep("female", 6), rep("male",3), rep("female",3))), gonad=factor(c(rep("hormonal", 6), rep("castrate",6), rep("hormonal",3))), gravid=factor(c(rep("no",12), rep("yes",3))))
rownames(sampleTable) <- colnames(txi.kallisto$counts)

write.csv(txi.kallisto, file=paste0(output,"_import.csv"))
# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromTximport(txi.kallisto, sampleTable, design =~condition)
dds <- dds[!(rownames(dds) %in% ychromgenes),]
keep <- rowSums(counts(dds) >= 5) >=9
dds <- dds[keep,]

dds <- DESeq(dds)
# run DESeq2
resultsNames(dds)
names <- resultsNames(dds)
write.csv(names, file=paste0(output,"ResultsNames.csv"))
# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)


resCastratevPreg <- results(dds, contrast=c("condition",cond4,cond5))
table(resCastratevPreg$padj<sigpadj)
## Order by adjusted p-value
resCastratevPreg <- resCastratevPreg[order(resCastratevPreg$padj), ]
## Merge with normalized count data
resCastratevPregdata <- merge(as.data.frame(resCastratevPreg), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resCastratevPreg)[1] <- "Gene"
## Write resuls
write.csv(resCastratevPregdata, file=paste0(output, "CasvPreg-ExpDiff.csv"))