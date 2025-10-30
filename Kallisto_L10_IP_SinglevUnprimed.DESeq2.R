#!/usr/bin/env Rscript
# variables
cond1 <- "PrimedIP"
cond2 <- "SingleHousedMale"
output <- "All_RNASeq_Analyses/DESeq2/Kallisto_UnprimedvSingle_"
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

dir <- "/scratch/users/knoedler/L10_Kallisto_vm17_Deopt_Tile_adjacent/"
samples <- read.table("NewKallistoSheets/L10VMH_IP_UnprimedvSingle.txt", header=TRUE)
samples
files <- file.path(dir, samples$sample, "abundance.h5")
files
names(files) <- paste0(samples$name, 1:6)
tx2gene <- read.table("Bioinformatics3/Bioinformatics/vM17.tx2gene.FINAL.txt")
head(tx2gene)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

sampleTable <- data.frame(condition = factor(c(rep(cond1, 3), rep(cond2, 3))))
rownames(sampleTable) <- colnames(txi.kallisto$counts)


# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromTximport(txi.kallisto, sampleTable, design =~condition)
dds <- DESeq(dds)																													# run DESeq2

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)


# Save differential expression results
res <- results(dds)
table(res$padj<sigpadj)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
## Write results
write.csv(resdata, file=paste0(output, "all-ExpDiff.csv"))

# Plot dispersions
pdf(paste0(output, cond1, "_vs_", cond2, "-DispersionPlot.pdf"))
plotDispEsts(dds, main=paste0(cond1, " vs ", cond2, " Dispersion Plot"))
garbage <- dev.off() # Save to file

# Principal Components Analysis
pdf(paste0(output, "all-PCA.pdf"))
plotPCA(rld, intgroup="condition",ntop=250)
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
respairwise = results(dds, contrast=c("condition",cond1,cond2))																		# Select comparison to perform (for log2 changes)
sig <- rownames(respairwise[!is.na(respairwise$padj) & respairwise$padj<sigpadj & abs(respairwise$log2FoldChange)>=siglfc,])[1:30]	# Select the first 30 significant hits (sorted by padj)
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

