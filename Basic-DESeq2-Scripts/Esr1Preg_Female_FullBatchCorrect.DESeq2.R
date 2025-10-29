#!/usr/bin/env Rscript
# variables
cond1 <- "PrimedPOA"
cond2 <- "UnprimedPOA"
cond3 <- "PrimedBNST"
cond4 <- "UnprimedBNST"
cond5 <- "PrimedMeA"
cond6 <- "UnprimedMeA"
cond7 <- "PrimedVMH"
cond8 <- "UnprimedVMH"
cond9 <- "BNSTpreg"
cond10 <- "MeApreg"
cond11 <- "POApreg"
cond12 <- "VMHpreg"

output <- "DESeq2Scripts/Esr1Adult_PregGLM_Femaleonly_"
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
samples <- read.table("Kallistosheets/AllAdultEsr1preg_FemaleOnly.tsv", header=TRUE)
samples
files <- file.path(dir, samples$sample, "abundance.h5")
files
names(files) <- paste0(samples$name, 1:34)
tx2gene <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/vM21_Tx2Gene_Updated.txt")
head(tx2gene)
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

sampleTable <- data.frame(condition = factor(c(rep(cond1, 3), rep(cond2, 3), rep(cond3, 3), rep(cond4, 3),rep(cond5, 3), rep(cond6, 3), rep(cond7, 3), rep(cond8, 3), rep(cond9, 3), rep(cond10, 2), rep(cond11, 2), rep(cond12, 3))), region=factor(c(rep("POA",6), rep("BNST",6), rep("MeA", 6), rep("VMH",6), rep("BNST",3),rep("MeA",2),rep("POA",2), rep("VMH",3))), batch = factor(c(rep("1"),("1"), ("2"),rep("1"),rep("1"), ("2"),rep("1"),rep("1"), ("2"),rep("1"),rep("1"), ("2"),rep("1"),rep("1"), ("2"),rep("1"),rep("1"), ("2"),rep("1"),rep("1"), ("2"),rep("1"),rep("1"), ("2"), rep("3",10))), state=factor(c(rep("Primed",3), rep("Unprimed",3), rep("Primed",3), rep("Unprimed",3), rep("Primed",3), rep("Unprimed",3), rep("Primed",3), rep("Unprimed",3),rep("Pregnant",10)))) 
rownames(sampleTable) <- colnames(txi.kallisto$counts)

ychromgenes <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt")
ychromgenes <- unlist(ychromgenes)
head(ychromgenes)
write.csv(txi.kallisto, file=paste0(output,"_import.csv"))
# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromTximport(txi.kallisto, sampleTable, design =~0+state:region)
keep <- rowSums(counts(dds) >= 5) >=7
dds <- dds[keep,]
dds <- DESeq(dds)
# run DESeq2
resultsNames(dds)
# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)


# Save differential expression results
resEstvPregVMH <- results(dds, contrast=list("statePrimed.regionVMH","statePregnant.regionVMH"))
table(resEstvPregVMH$padj<sigpadj)
## Order by adjusted p-value
resEstvPregVMH <- resEstvPregVMH[order(resEstvPregVMH$padj), ]
## Merge with normalized count data
resEstvPregVMHdata <- merge(as.data.frame(resEstvPregVMH), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resEstvPregVMHdata)[1] <- "Gene"
## Write results
write.csv(resEstvPregVMHdata, file=paste0(output, "_VMH_EstrusvPreg-ExpDiff.csv"))

##write transformed counts if you want to do 3D PCA or what have you
write.csv(assay(rld), file=paste0(output, "_transformed_counts.csv"))
# Principal Components Analysis
pdf(paste0(output, "all-PCA.pdf"))
plotPCA(rld, intgroup="condition",ntop=1500)
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

