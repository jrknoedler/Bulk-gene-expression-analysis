#!/usr/bin/env Rscript
# variables
cond1 <- "Male"
cond2 <- "Etrus"
cond3 <- "Diestrus"
cond4 <- "Pregnant"
output <- "DESeq2Scripts/Esr1Adult_VMH_GLM_Preg4way_droplowdropYdropmt_GLMdegpca"
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
Pregdets <- read.table("VMH_Pregdegs2.txt", header=FALSE)
#Pregdets2 <- read.table("VMH_Pregdegs2.txt", header=FALSE)
#Pregdets2 <- unlist(Pregdets2)
#Pregdets2.df <- data.frame(Pregdets2.df)
Pregdets <- unlist(Pregdets)
Pregdets.df <- data.frame(Pregdets)
#Preddets.wtf <- cbind(Pregdets.df, Pregdets2.df)
dir <- "/scratch/users/knoedler/Kallisto_ProteinRef/"
samples <- read.table("Kallistosheets/AdultEsr1Kallisto_VMH_FULLpreg.tsv", header=TRUE)
samples
files <- file.path(dir, samples$sample, "abundance.h5")
files
names(files) <- paste0(samples$name, 1:12)
tx2gene <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/vM21_Tx2Gene_Updated.txt")
head(tx2gene)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

write.csv(txi.kallisto, file=paste0(output,"_import.csv"))
ychromgenes <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt")
ychromgenes <- unlist(ychromgenes)
head(ychromgenes)

sampleTable <- data.frame(condition = factor(c(rep(cond1, 3), rep(cond2, 3), rep(cond3, 3), rep(cond4, 3))))
rownames(sampleTable) <- colnames(txi.kallisto$counts)
sampleTable <- sampleTable[Pregdets,]

# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromTximport(txi.kallisto, sampleTable, design =~condition)
dds <- dds[!(rownames(dds) %in% ychromgenes),]
keep <- rowSums(counts(dds) >= 5) >=9


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
resMalevPrimed <- results(dds, contrast=c("condition",cond1,cond2))
table(resMalevPrimed$padj<sigpadj)
## Order by adjusted p-value
resMalevPrimed <- resMalevPrimed[order(resMalevPrimed$padj), ]
## Merge with normalized count data
resMalevPrimeddata <- merge(as.data.frame(resMalevPrimed), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resMalevPrimeddata)[1] <- "Gene"
## Write results
write.csv(resMalevPrimeddata, file=paste0(output, "_MvP-ExpDiff.csv"))
resPrimedvUnprimed <- results(dds, contrast=c("condition",cond2,cond3))
table(resPrimedvUnprimed$padj<sigpadj)
## Order by adjusted p-value
resPrimedvUnprimed <- resPrimedvUnprimed[order(resPrimedvUnprimed$padj), ]
## Merge with normalized count data
resPrimedvUnprimeddata <- merge(as.data.frame(resPrimedvUnprimed), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resPrimedvUnprimed)[1] <- "Gene"
## Write results
write.csv(resPrimedvUnprimeddata, file=paste0(output, "PvU-ExpDiff.csv"))
resMalevUnprimed <- results(dds, contrast=c("condition",cond1,cond3))
table(resMalevUnprimed$padj<sigpadj)
## Order by adjusted p-value
resMalevUnprimed <- resMalevUnprimed[order(resMalevUnprimed$padj), ]
## Merge with normalized count data
resMalevUnprimeddata <- merge(as.data.frame(resMalevUnprimed), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resMalevUnprimed)[1] <- "Gene"
## Write results
write.csv(resMalevUnprimeddata, file=paste0(output, "MvU-ExpDiff.csv"))
resMalevPreg <- results(dds, contrast=c("condition",cond1,cond4))
table(resMalevPreg$padj<sigpadj)
## Order by adjusted p-value
resMalevPreg <- resMalevPreg[order(resMalevPreg$padj), ]
## Merge with normalized count data
resMalevPregdata <- merge(as.data.frame(resMalevPreg), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resMalevPreg)[1] <- "Gene"
## Write resuls
write.csv(resMalevPregdata, file=paste0(output, "MvPreg-ExpDiff.csv"))
resPrimedvPreg <- results(dds, contrast=c("condition",cond2,cond4))
table(resPrimedvPreg$padj<sigpadj)
## Order by adjusted p-value
resPrimedvPreg <- resPrimedvPreg[order(resPrimedvPreg$padj), ]
## Merge with normalized count data
resPrimedvPregdata <- merge(as.data.frame(resPrimedvPreg), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resPrimedvPreg)[1] <- "Gene"
## Write resuls
write.csv(resPrimedvPregdata, file=paste0(output, "PvPreg-ExpDiff.csv"))
resUnprimedvPreg <- results(dds, contrast=c("condition",cond3,cond4))
table(resUnprimedvPreg$padj<sigpadj)
## Order by adjusted p-value
resUnprimedvPreg <- resUnprimedvPreg[order(resUnprimedvPreg$padj), ]
## Merge with normalized count data
resUnprimedvPregdata <- merge(as.data.frame(resUnprimedvPreg), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resPrimedvPreg)[1] <- "Gene"
## Write resuls
write.csv(resUnprimedvPregdata, file=paste0(output, "UvPreg-ExpDiff.csv"))
##write transformed counts if you want to do 3D PCA or what have you
write.csv(assay(rld), file=paste0(output, "_transformed_counts.csv"))
# Principal Components Analysis
pdf(paste0(output, "all-PCA.pdf"))



plotPCA(rld, intgroup="condition",ntop=500)
garbage <- dev.off() # Save to file


# Sample Distance Matrix
dmcol <- colorRampPalette(brewer.pal(12, 'RdBu'))(100)
sampleDists <- as.matrix(dist(t(assay(rld))))

pdf(paste0(output, "all-DistanceMatrix.pdf"))
heatmap.2(as.matrix(sampleDists), key=FALSE, symm=TRUE, trace="none", col=rev(dmcol), margin=c(15, 15), main="all samples Matrix")
garbage <- dev.off() # Save to file
pdf(paste0(output,"Disp.pdf"))
plotDispEsts(dds)
dev.off()

# Heatmap of significant hits ( padj<sigpadj and |log2FoldChange|>=siglfc )
hmcol <- colorRampPalette(brewer.pal(9, 'RdYlBu'))(100)
counts <- counts(dds,normalized=TRUE)																								# Get normalized read counts (sorted by padj)
counts <- counts[apply(counts, 1, function(row) all(row !=0 )),]																	# remove genes with zero reads
## You will probably need to change conditions to perform a meaningful comparison here
respairwise = results(dds, contrast=c("condition",cond1,cond2))																		# Select comparison to perform (for log2 changes)
sig <- rownames(respairwise[!is.na(respairwise$padj) & respairwise$padj<sigpadj,])	# Select the first 30 significant hits (sorted by padj)
sig <- sig[!is.na(sig)]					
sigcounts <- counts(dds,normalized=TRUE)[sig,]	
sigcounts <- sigcounts[,c(1:6)]																				# Get normalized read counts (sorted by padj) for significan genes
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

