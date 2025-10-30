#!/usr/bin/env Rscript
# variables
cond1 <- "BNSTMale"
cond2 <- "BNSTPrimed"
cond3 <- "BNSTUnprimed"
cond4 <- "BNSTCastrateMale"
cond5 <- "BNSTPregnantFemale"
cond6 <- "MeAMale"
cond7 <- "MeAPrimed"
cond8 <- "MeAUnprimed"
cond9 <- "MeACastrateMale"
cond10 <- "MeAPregnantFemale"
cond11 <- "POAMale"
cond12 <- "POAPrimed"
cond13 <- "POAUnprimed"
cond14 <- "POACastrateMale"
cond15 <- "POAPregnantFemale"
cond16 <- "VMHMale"
cond17 <- "VMHPrimed"
cond18 <- "VMHUnprimed"
cond19 <- "VMHCastrateMale"
cond20 <- "VMHPregnantFemale"
output <- "DESeq2Scripts/Esr1Adult_AllSamples_GLMtests_nestedTest"
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
samples <- read.table("Kallistosheets/AllAdultEsr1preg_newsamples.tsv", header=TRUE)
samples
files <- file.path(dir, samples$sample, "abundance.h5")
files
names(files) <- paste0(samples$name, 1:60)
tx2gene <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/vM21_Tx2Gene_Updated.txt")
head(tx2gene)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

ychromgenes <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt")
ychromgenes <- unlist(ychromgenes)
head(ychromgenes)

sampleTable <- data.frame(condition = factor(c(rep(cond1, 3), rep(cond2, 3), rep(cond3, 3), rep(cond4, 3), rep(cond5, 3), rep(cond6, 3), rep(cond7, 3), rep(cond8, 3), rep(cond9, 3), rep(cond10, 3), rep(cond11, 3), rep(cond12, 3), rep(cond13, 3), rep(cond14, 3), rep(cond15, 3), rep(cond16, 3), rep(cond17,3), rep(cond18, 3), rep(cond19, 3), rep(cond20,3))), region=factor(c(rep("BNST", 15), rep("MeA", 15), rep("POA", 15), rep("VMH",15))), batch=factor(c(rep("1"),rep("2"),rep("3"),rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"), rep("4",6), rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"), rep("4",5),rep("5",1), rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"), rep("4",5),rep("5",1), rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"), rep("4",6))), gravid=factor(c(rep("no",12), rep("yes",3),rep("no",12), rep("yes",3),rep("no",12), rep("yes",3),rep("no",12), rep("yes",3))), state=factor(c(rep("Male",3), rep("Primed",3), rep("Unprimed",3), rep("Castrate",3), rep("Pregnant",3), rep("Male",3), rep("Primed",3), rep("Unprimed",3), rep("Castrate",3), rep("Pregnant",3), rep("Male",3), rep("Primed",3), rep("Unprimed",3), rep("Castrate",3), rep("Pregnant",3), rep("Male",3), rep("Primed",3), rep("Unprimed",3), rep("Castrate",3), rep("Pregnant",3))))
rownames(sampleTable) <- colnames(txi.kallisto$counts)

write.csv(txi.kallisto, file=paste0(output,"_import.csv"))
# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromTximport(txi.kallisto, sampleTable, design =~batch+region+region:state)
dds <- dds[!(rownames(dds) %in% ychromgenes),]
keep <- rowSums(counts(dds) >= 5) >=9
dds <- dds[keep,]

dds <- DESeq(dds)
# run DESeq2
names <- resultsNames(dds)
write.csv(names, file=paste0(output,"Resultnames.csv"))

# Regularized log transformation for clustering/heatmaps, etc
#rld <- varianTransformation(dds)
#Save DE results
resPOAPrimedvPreg <- results(dds, contrast=list("regionPOA.statePrimed","regionPOA.statePregnant"))
#resPOAPrimedvPreg <- results(dds, name="statePregnant.regionPOA")
table(resPOAPrimedvPreg$padj<sigpadj)
## Order by adjusted p-value
resPOAPrimedvPreg <- resPOAPrimedvPreg[order(resPOAPrimedvPreg$padj), ]
## Merge with normalized count data
resPOAPrimedvPregdata <- merge(as.data.frame(resPOAPrimedvPreg), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resPOAPrimedvPregdata)[1] <- "Gene"
## Write results
write.csv(resPOAPrimedvPregdata, file=paste0(output, "_POA_PrimedvPreg-ExpDiff.csv"))

