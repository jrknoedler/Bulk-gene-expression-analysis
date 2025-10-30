#!/usr/bin/env Rscript
# variables
cond1 <- "PrimedBNST"
cond2 <- "UnprimedBNST"
cond3 <- "MaleBNST"
cond4 <- "MaleMeA"
cond5 <- "MalePOA"
cond6 <- "MaleVMH"
cond7 <- "PrimedMeA"
cond8 <- "UnprimedMeA"
cond9 <- "PrimedPOA"
cond10 <- "UnprimedPOA"
cond11 <- "PrimedVMH"
cond12 <- "UnprimedVMH"
output <- "DESeq2Scripts/Esr1Adult_FULL_GLM_batch_adjust"
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
samples <- read.table("Kallistosheets_cp/AdultEsr1Kallisto_GLM.tsv", header=TRUE)
samples
files <- file.path(dir, samples$sample, "abundance.h5")
files
names(files) <- paste0(samples$name, 1:36)
tx2gene <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/vM21_Tx2Gene_Final.txt")
head(tx2gene)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

sampleTable <- data.frame(condition = factor(c(rep(cond1, 3), rep(cond2, 3), rep(cond3, 3), rep(cond4, 3), rep(cond5, 3), rep(cond6, 3), rep(cond7, 3), rep(cond8, 3), rep(cond9, 3), rep(cond10, 3), rep(cond11, 3), rep(cond12, 3))), batch= factor(c(rep("1"), rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("3"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"),rep("1"),rep("1"),rep("2"))))
rownames(sampleTable) <- colnames(txi.kallisto$counts)


# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromTximport(txi.kallisto, sampleTable, design =~batch+condition)
dds <- DESeq(dds)																													# run DESeq2
resultsNames(dds)
# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)

# Save differential expression results
resBNST1 <- results(dds, contrast=c("condition",cond1,cond3))
table(resBNST1$padj<sigpadj)
## Order by adjusted p-value
resBNST1 <- resBNST1[order(resBNST1$padj), ]
## Merge with normalized count data
resBNST1data <- merge(as.data.frame(resBNST1), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resBNST1data)[1] <- "Gene"
## Write results
write.csv(resBNST1data, file=paste0(output, "BNST_MvP-ExpDiff.csv"))
resBNST2 <- results(dds, contrast=c("condition",cond2,cond3))
table(resBNST2$padj<sigpadj)
## Order by adjusted p-value
resBNST2 <- resBNST2[order(resBNST2$padj), ]
## Merge with normalized count data
resBNST2data <- merge(as.data.frame(resBNST2), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resBNST2data)[1] <- "Gene"
## Write results
write.csv(resBNST2data, file=paste0(output, "BNST_MvU-ExpDiff.csv"))
resBNST3 <- results(dds, contrast=c("condition",cond1,cond2))
table(resBNST3$padj<sigpadj)
## Order by adjusted p-value
resBNST3 <- resBNST3[order(resBNST3$padj), ]
## Merge with normalized count data
resBNST3data <- merge(as.data.frame(resBNST3), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resBNST1data)[1] <- "Gene"
## Write results
write.csv(resBNST3data, file=paste0(output, "BNST_PvU-ExpDiff.csv"))
resMeA1 <- results(dds, contrast=c("condition",cond4,cond7))
table(resMeA1$padj<sigpadj)
## Order by adjusted p-value
resMeA1 <- resMeA1[order(resMeA1$padj), ]
## Merge with normalized count data
resMeA1data <- merge(as.data.frame(resMeA1), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resMeA1data)[1] <- "Gene"
## Write results
write.csv(resMeA1data, file=paste0(output, "MeA_MvP-ExpDiff.csv"))
resMeA2 <- results(dds, contrast=c("condition",cond4,cond8))
table(resMeA2$padj<sigpadj)
## Order by adjusted p-value
resMeA2 <- resMeA2[order(resMeA2$padj), ]
## Merge with normalized count data
resMeA2data <- merge(as.data.frame(resMeA2), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resMeA2data)[1] <- "Gene"
## Write results
write.csv(resMeA2data, file=paste0(output, "MeA_MvU-ExpDiff.csv"))
resMeA3 <- results(dds, contrast=c("condition",cond7,cond8))
table(resMeA3$padj<sigpadj)
## Order by adjusted p-value
resMeA3 <- resMeA3[order(resMeA1$padj), ]
## Merge with normalized count data
resMeA3data <- merge(as.data.frame(resMeA3), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resMeA3data)[1] <- "Gene"
## Write results
write.csv(resMeA3data, file=paste0(output, "MeA_PvU-ExpDiff.csv"))
resPOA1 <- results(dds, contrast=c("condition",cond5,cond9))
table(resPOA1$padj<sigpadj)
## Order by adjusted p-value
resPOA1 <- resPOA1[order(resPOA1$padj), ]
## Merge with normalized count data
resPOA1data <- merge(as.data.frame(resPOA1), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resPOA1data)[1] <- "Gene"
## Write results
write.csv(resPOA1data, file=paste0(output, "POA_MvP-ExpDiff.csv"))
resPOA2 <- results(dds, contrast=c("condition",cond5,cond10))
table(resPOA2$padj<sigpadj)
## Order by adjusted p-value
resPOA2 <- resPOA2[order(resPOA2$padj), ]
## Merge with normalized count data
resPOA2data <- merge(as.data.frame(resPOA2), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resPOA2data)[1] <- "Gene"
## Write results
write.csv(resPOA2data, file=paste0(output, "POA_MvU-ExpDiff.csv"))
resPOA3 <- results(dds, contrast=c("condition",cond9,cond10))
table(resPOA3$padj<sigpadj)
## Order by adjusted p-value
resPOA3 <- resPOA3[order(resPOA3$padj), ]
## Merge with normalized count data
resPOA3data <- merge(as.data.frame(resPOA3), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resPOA3data)[1] <- "Gene"
## Write results
write.csv(resPOA3data, file=paste0(output, "POA_PvU-ExpDiff.csv"))
resVMH1 <- results(dds, contrast=c("condition",cond6,cond11))
table(resVMH1$padj<sigpadj)
## Order by adjusted p-value
resVMH1 <- resVMH1[order(resVMH1$padj), ]
## Merge with normalized count data
resVMH1data <- merge(as.data.frame(resVMH1), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resVMH1data)[1] <- "Gene"
## Write results
write.csv(resVMH1data, file=paste0(output, "VMH_MvP-ExpDiff.csv"))
resVMH2 <- results(dds, contrast=c("condition",cond6,cond12))
table(resVMH2$padj<sigpadj)
## Order by adjusted p-value
resVMH2 <- resVMH2[order(resVMH2$padj), ]
## Merge with normalized count data
resVMH2data <- merge(as.data.frame(resVMH2), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resVMH2data)[1] <- "Gene"
## Write results
write.csv(resVMH2data, file=paste0(output, "VMH_MvU-ExpDiff.csv"))
resVMH3 <- results(dds, contrast=c("condition",cond11,cond12))
table(resVMH3$padj<sigpadj)
## Order by adjusted p-value
resVMH3 <- resVMH3[order(resVMH3$padj), ]
## Merge with normalized count data
resVMH3data <- merge(as.data.frame(resVMH3), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resVMH3data)[1] <- "Gene"
## Write results
write.csv(resVMH3data, file=paste0(output, "VMH_PvU-ExpDiff.csv"))
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


# Adapted from:
# http://www.bioconductor.org/help/workflows/rnaseqGene/
# https://gist.github.com/stephenturner/f60c1934405c127f09a6
# https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/

