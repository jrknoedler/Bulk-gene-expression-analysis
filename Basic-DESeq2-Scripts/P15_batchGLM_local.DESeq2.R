#!/usr/bin/env Rscript
# variables
cond1 <- "FemaleBNST"
cond2 <- "MaleBNST"
cond3 <- "FemaleMeA"
cond4 <- "MaleMeA"
cond5 <- "FemalePOA"
cond6 <- "MalePOA"
cond7 <- "FemaleVMH"
cond8 <- "MaleVMH"
output <- "DESeq2Scripts/P15_GlobalBatchCorrection_local_MaleVMHupdated"
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

dir <- "/scratch/users/knoedler/Trio_Coding_Redo2/"
samples <- read.table("Kallistosheets_cp/AllP15Esr1IPKallisto.txt", header=TRUE)
samples
files <- file.path(dir, samples$sample, "abundance.h5")
files
names(files) <- paste0(samples$name, 1:16)
tx2gene <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/vM21_Tx2Gene_Final.txt")
head(tx2gene)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

sampleTable <- data.frame(condition=factor(c(rep(cond1, 2), rep(cond2, 2), rep(cond3, 2), rep(cond4, 2), rep(cond5, 2), rep(cond6, 2), rep(cond7, 2), rep(cond8, 2))), batch = c(rep("1"), rep("2",3), rep("1"), rep("2",3), rep("1"), rep("2",3), rep("1"), rep("2",3)))
rownames(sampleTable) <- colnames(txi.kallisto$counts)


# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromTximport(txi.kallisto, sampleTable, design=~batch+condition)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, fitType="local")																													# run DESeq2
resultsNames(dds)

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)

# Save differential expression results
resBNST <- results(dds, contrast=c("condition", cond1, cond2))
## Order by adjusted p-value
resBNST <- resBNST[order(resBNST$padj), ]
## Merge with normalized count data
resBNSTdata <- merge(as.data.frame(resBNST), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resBNSTdata)[1] <- "Gene"
## Write results
write.csv(resBNSTdata, file=paste0(output, "BNST-ExpDiff.csv"))
resMeA <- results(dds, contrast=c("condition",cond3,cond4))
table(resMeA$padj<sigpadj)
## Order by adjusted p-value
resMeA <- resMeA[order(resMeA$padj), ]
## Merge with normalized count data
resMeAdata <- merge(as.data.frame(resMeA), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resMeAdata)[1] <- "Gene"
## Write results
write.csv(resMeAdata, file=paste0(output, "MeA-ExpDiff.csv"))
resPOA <- results(dds, contrast=c("condition",cond5,cond6))
table(resPOA$padj<sigpadj)
## Order by adjusted p-value
resPOA <- resPOA[order(resPOA$padj), ]
## Merge with normalized count data
resPOAdata <- merge(as.data.frame(resPOA), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resPOAdata)[1] <- "Gene"
## Write results
write.csv(resPOAdata, file=paste0(output, "POA-ExpDiff.csv"))
resVMH <- results(dds, contrast=c("condition",cond7,cond8))
table(resVMH$padj<sigpadj)
## Order by adjusted p-value
resVMH <- resVMH[order(resVMH$padj), ]
## Merge with normalized count data
resVMHdata <- merge(as.data.frame(resVMH), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resVMHdata)[1] <- "Gene"
## Write results
write.csv(resVMHdata, file=paste0(output, "VMH-ExpDiff.csv"))
##write transformed counts if you want to do 3D PCA or what have you
write.csv(assay(rld), file=paste0(output, "_transformed_counts.csv"))
# Principal Components Analysis
pdf(paste0(output, "all-PCA.pdf"))
plotPCA(rld, intgroup="condition",ntop=2000)
garbage <- dev.off() # Save to file


# Sample Distance Matrix
dmcol <- colorRampPalette(brewer.pal(12, 'RdBu'))(100)
sampleDists <- as.matrix(dist(t(assay(rld))))

pdf(paste0(output, "all-DistanceMatrix.pdf"))
heatmap.2(as.matrix(sampleDists), key=FALSE, symm=TRUE, trace="none", col=rev(dmcol), margin=c(15, 15), main="all samples Matrix")
garbage <- dev.off() # Save to file


# Heatmap of significant hits ( padj<sigpadj and |log2FoldChange|>=siglfc )
hmcol <- colorRampPalette(brewer.pal(9, 'RdYlBu'))(100)
counts <- counts(dds,normalized=TRUE)																								# Get normalized read counts (sorted by padj)
counts <- counts[apply(counts, 1, function(row) all(row !=0 )),]																	# remove genes with zero reads
## You will probably need to change conditions to perform a meaningful comparison here
respairwise = results(dds)																		# Select comparison to perform (for log2 changes)
sig <- rownames(respairwise[!is.na(respairwise$padj) & respairwise$padj<sigpadj & abs(respairwise$log2FoldChange)>=siglfc,])[1:100]	# Select the first 30 significant hits (sorted by padj)
sig <- sig[!is.na(sig)]					
sigcounts <- counts(dds,normalized=TRUE)[sig,]																						# Get normalized read counts (sorted by padj) for significan genes
try({																																# May fail if only 0 or 1 gene is significant
	sigcounts <- sigcounts[apply(sigcounts, 1, function(row) all(row !=0 )),]														# remove genes with zero reads
})

# By defaults hits are not clustered and thus stay sorted by their padj value
pdf(paste0(output, "all-HeatmapSig.pdf"), onefile=FALSE)
try({																																# May fail if only 0 or 1 gene is significant
	pheatmap(log2(sigcounts), col=rev(hmcol), scale='row', cluster_rows=FALSE, cluster_cols=FALSE, main="all Top Hits")
})
garbage <- dev.off() # Save to file

# But we can cluster them using the following command
pdf(paste0(output, "all-HeatmapSigClust.pdf"), onefile=FALSE)
try({																																# May fail if only 0 or 1 gene is significant
	pheatmap(log2(sigcounts), col=rev(hmcol), scale='row', cluster_rows=TRUE, cluster_cols=TRUE, main="all Clustered Top Hits")
})
garbage <- dev.off() # Save to file

# And even plot everything, this time using the old school Red/Green color palette
pdf(paste0(output, "all-HeatmapAllClust.pdf"), onefile=FALSE)
pheatmap(log2(counts), col=redgreen(75), scale='row', cluster_rows=TRUE, cluster_cols=TRUE, main="all Clustered Reads Count")
garbage <- dev.off() # Save to file

