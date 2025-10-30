#!/usr/bin/env Rscript
# variables
cond1 <- "Male"
cond2 <- "Primed"
cond3 <- "Unprimed"
cond4 <- "PregnantFemale"
cond5 <- "P15Female"
cond6 <- "P15Male"
output <- "DESeq2Scripts/Esr1Adult_POA_GLM_4way_dropY_droplow_dropmt_withP151stpass"
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
samples <- read.table("Kallistosheets/AdultEsr1Kallisto_POAFULLpreg_updated_withP15.tsv", header=TRUE)
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

sampleTable <- data.frame(condition = factor(c(rep(cond1, 3), rep(cond2, 3), rep(cond3, 3), rep(cond4, 2), rep(cond5,2), rep(cond6,2))))
rownames(sampleTable) <- colnames(txi.kallisto$counts)

write.csv(txi.kallisto, file=paste0(output,"_import.csv"))
# Run DESeq2 to get the result for differential-expression
dds=DESeqDataSetFromTximport(txi.kallisto, sampleTable, design =~condition)
dds <- dds[!(rownames(dds) %in% ychromgenes),]
keep <- rowSums(counts(dds) >= 5) >=13
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

resMale_P15vadult <- results(dds, contrast=c("condition",cond1,cond6))
table(resMale_P15vadult$padj<sigpadj)
## Order by adjusted p-value
resMale_P15vadult <- resMale_P15vadult[order(resMale_P15vadult$padj), ]
## Merge with normalized count data
resMale_P15vadultdata <- merge(as.data.frame(resMale_P15vadult), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resMale_P15vadult)[1] <- "Gene"
## Write results
write.csv(resMale_P15vadultdata, file=paste0(output, "Male_P15vadult-ExpDiff.csv"))

resFemale_P15vPrimed <- results(dds, contrast=c("condition",cond2,cond5))
table(resFemale_P15vPrimed$padj<sigpadj)
## Order by adjusted p-value
resFemale_P15vPrimed <- resFemale_P15vPrimed[order(resFemale_P15vPrimed$padj), ]
## Merge with normalized count data
resFemale_P15vPrimeddata <- merge(as.data.frame(resFemale_P15vPrimed), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resFemale_P15vPrimed)[1] <- "Gene"
## Write results
write.csv(resFemale_P15vPrimeddata, file=paste0(output, "Female_P15vPrimed-ExpDiff.csv"))

resFemale_P15vUnprimed <- results(dds, contrast=c("condition",cond3,cond5))
table(resFemale_P15vUnprimed$padj<sigpadj)
## Order by adjusted p-value
resFemale_P15vUnprimed <- resFemale_P15vUnprimed[order(resFemale_P15vUnprimed$padj), ]
## Merge with normalized count data
resFemale_P15vUnprimeddata <- merge(as.data.frame(resFemale_P15vUnprimed), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resFemale_P15vUnprimed)[1] <- "Gene"
## Write results
write.csv(resFemale_P15vUnprimeddata, file=paste0(output, "Female_P15vUnprimed-ExpDiff.csv"))

resFemale_P15vPreg <- results(dds, contrast=c("condition",cond4,cond5))
table(resFemale_P15vPreg$padj<sigpadj)
## Order by adjusted p-value
resFemale_P15vPreg <- resFemale_P15vPreg[order(resFemale_P15vPreg$padj), ]
## Merge with normalized count data
resFemale_P15vPregdata <- merge(as.data.frame(resFemale_P15vPreg), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resFemale_P15vPreg)[1] <- "Gene"
## Write results
write.csv(resFemale_P15vPregdata, file=paste0(output, "Female_P15vPreg-ExpDiff.csv"))


## Order by adjusted p-value

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
pdf(paste0(output,"Disp.pdf"))
plotDispEsts(dds)
dev.off()

# Heatmap of significant hits ( padj<sigpadj and |log2FoldChange|>=siglfc )
hmcol <- colorRampPalette(brewer.pal(9, 'RdYlBu'))(100)
counts <- counts(dds,normalized=TRUE)																								# Get normalized read counts (sorted by padj)
counts <- counts[apply(counts, 1, function(row) all(row !=0 )),]																	# remove genes with zero reads
## You will probably need to change conditions to perform a meaningful comparison here
respairwise = results(dds)																		# Select comparison to perform (for log2 changes)
siggenes <- read.table("DESeq2Scripts/VMH_allsig.txt")
siggenes <- unlist(siggenes)
				
sigcounts <- counts(dds,normalized=TRUE)	
sigcounts													# Get normalized read counts (sorted by padj) for significan genes
try({																																# May fail if only 0 or 1 gene is significant
	sigcounts <- sigcounts[apply(sigcounts, 1, function(row) all(row !=0 )),]														# remove genes with zero reads
})
ntd <- normTransform(dds)
# By defaults hits are not clustered and thus stay sorted by their padj value
pdf(paste0(output, "all-HeatmapSig.pdf"), onefile=FALSE)
try({																																# May fail if only 0 or 1 gene is significant
	pheatmap(assay(ntd)[siggenes,], col=rev(hmcol), scale='row', cluster_rows=FALSE, cluster_cols=FALSE, main="all Top Hits")
})
garbage <- dev.off() # Save to file

# But we can cluster them using the following command
pdf(paste0(output, "all-HeatmapSigClust.pdf"), onefile=FALSE)
try({																																# May fail if only 0 or 1 gene is significant
	pheatmap(assay(ntd)[siggenes,], col=hmcol, scale='row', cluster_rows=TRUE, cluster_cols=TRUE, main="all Clustered Top Hits")
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

