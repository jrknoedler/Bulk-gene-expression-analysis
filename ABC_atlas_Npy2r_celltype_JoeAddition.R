rm(list=ls())

library(Seurat)
library(ggplot2)
library(patchwork)

# Load POA scRNAseq data
mySeurat <- readRDS('C:/Users/Joe/Downloads/Merged_WMB-HY_mPOA_subset_adata.rds')
mySeurat
# An object of class Seurat 
# 32285 features across 23383 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 2 layers present: counts, data

# log-normalize expression levels (RNA)
mySeurat <- NormalizeData(mySeurat)

unique(mySeurat$subclass)

# Set idents as subclass
Idents(mySeurat) <- mySeurat$subclass

VlnPlot(mySeurat, features = c('ENSMUSG00000028004','ENSMUSG00000030500','ENSMUSG00000070880'), ncol = 1, pt.size = 0) &
  theme(axis.text.x = element_blank())

# ! The gene names are endembl IDs not gene symbols
# ENSMUSG00000028004: Npy2r
# ENSMUSG00000019768: Esr1
# ENSMUSG00000030500: Slc17a6 (Vglut2)
# ENSMUSG00000070880: Gad1

VlnPlot(mySeurat, features = c('ENSMUSG00000028004'), ncol = 1, pt.size = 0) 

# There are non-POA cells - remove them
remove <- unique(mySeurat$subclass[grepl("^(177|185|190|277|285|287)", mySeurat$subclass)])
subset.mySeurat <- subset(mySeurat, idents = remove, invert = TRUE)
subset.mySeurat
# An object of class Seurat 
# 32285 features across 23347 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 2 layers present: counts, data
unique(subset.mySeurat$subclass[grepl("^(177|185|190|277|285|287)", subset.mySeurat$subclass)]) # Check

Idents(subset.mySeurat) <- subset.mySeurat$subclass
VlnPlot(subset.mySeurat, features = c('ENSMUSG00000028004','ENSMUSG00000030500','ENSMUSG00000070880'), ncol = 1, pt.size = 1) &
  theme(axis.text.x = element_blank())

Idents(subset.mySeurat) <- subset.mySeurat$supertype
plots <- VlnPlot(subset.mySeurat, features = c('ENSMUSG00000028004','ENSMUSG00000019768','ENSMUSG00000030500','ENSMUSG00000070880'), ncol = 1, pt.size = 0) &
  theme(axis.text.x = element_blank())
new_labels <- c('Npy2r', 'Esr1', 'Vglut2', 'Gad1')

# Replace titles of each ggplot object
for (i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + ggtitle(new_labels[i])
}

wrap_plots(plots)

# Make it nicer 
Vlnplot.data.Npy2r <- plots[[1]]$data
Vlnplot.data.Esr1 <- plots[[2]]$data
Vlnplot.data.Slc17a6 <- plots[[3]]$data
Vlnplot.data.Gad1 <- plots[[4]]$data


vln.Npy2r <- ggplot(Vlnplot.data.Npy2r, aes(x = ident, y = ENSMUSG00000028004, fill = ident)) +
  geom_violin(trim = TRUE, scale = "width") +
  theme_classic() +
  theme(
    axis.line.y  = element_line(size = 0.5, color = "darkgray"),
    axis.line.x  = element_line(size = 0.5, color = "darkgray"),
    axis.text.y = element_text(size = 12, colour = 'black'),
    axis.text.x = element_blank(),
    axis.ticks = element_line(size = 0.5, color = "darkgray"),
    
    legend.position = "none"
  ) +
  guides(fill = "none") +
  labs(title = NULL,x = NULL, y = NULL)
vln.Npy2r

vln.Esr1 <- ggplot(Vlnplot.data.Esr1, aes(x = ident, y = ENSMUSG00000019768, fill = ident)) +
  geom_violin(trim = TRUE, scale = "width") +
  theme_classic() +
  theme(
    axis.line.y  = element_line(size = 0.5, color = "darkgray"),
    axis.line.x  = element_line(size = 0.5, color = "darkgray"),
    axis.text.y = element_text(size = 12, colour = 'black'),
    axis.text.x = element_blank(),
    axis.ticks = element_line(size = 0.5, color = "darkgray"),
    
    legend.position = "none"
  ) +
  guides(fill = "none") +
  labs(title = NULL,x = NULL, y = NULL)
vln.Esr1

vln.Slc17a6 <- ggplot(Vlnplot.data.Slc17a6, aes(x = ident, y = ENSMUSG00000030500, fill = ident)) +
  geom_violin(trim = TRUE, scale = "width") +
  theme_classic() +
  theme(
    axis.line.y  = element_line(size = 0.5, color = "darkgray"),
    axis.line.x  = element_line(size = 0.5, color = "darkgray"),
    axis.text.y = element_text(size = 12, colour = 'black'),
    axis.text.x = element_blank(),
    axis.ticks = element_line(size = 0.5, color = "darkgray"),
    
    legend.position = "none"
  ) +
  guides(fill = "none") +
  labs(title = NULL,x = NULL, y = NULL)
vln.Slc17a6

vln.Gad1 <- ggplot(Vlnplot.data.Gad1, aes(x = ident, y = ENSMUSG00000070880, fill = ident)) +
  geom_violin(trim = TRUE, scale = "width") +
  theme_classic() +
  theme(
    axis.line.y  = element_line(size = 0.5, color = "darkgray"),
    axis.line.x  = element_line(size = 0.5, color = "darkgray"),
    axis.text.y = element_text(size = 12, colour = 'black'),
    axis.text.x = element_blank(),
    axis.ticks = element_line(size = 0.5, color = "darkgray"),
    
    legend.position = "none"
  ) +
  guides(fill = "none") +
  labs(title = NULL,x = NULL, y = NULL)
vln.Gad1

vln.Gad1.modified <- vln.Gad1+theme(axis.text.x = element_text(color="black", size=12, angle = 90,
                                                               vjust=0.5,
                                                               hjust=1,))
vln.Gad1.modified

pdf("ABC_Npy2r_Vlnplot_100125.pdf", width = 15, height = 9)
vln.Npy2r+
  vln.Esr1+
  vln.Slc17a6+
  vln.Gad1.modified+
  plot_layout(nrow = 4, heights = c(1,1,1,1))
dev.off()


## Check Npy2r+ (25% cells express Npy2r)
Dotplot <- DotPlot(subset.mySeurat, features = 'ENSMUSG00000028004')
Dotplot.data <- Dotplot$data
Npy2r.pos <- subset(Dotplot.data, pct.exp > 25, select = id)
Npy2r.pos <- Npy2r.pos$id

# # Subset 85 and 106 Npy2r+ POA
# keep <- unique(subset.mySeurat$subclass[grepl("^(106|085)", subset.mySeurat$subclass)]) # Check
# keep
# POA.Npy2r.mySeurat <- subset(subset.mySeurat, idents = keep)
# POA.Npy2r.mySeurat
# 
# Idents(POA.Npy2r.mySeurat) <- POA.Npy2r.mySeurat$supertype
# VlnPlot(POA.Npy2r.mySeurat, features = c('ENSMUSG00000028004'), ncol = 1, pt.size = 0) + theme(legend.position = 'none')
# 
# VlnPlot(POA.Npy2r.mySeurat, features = c('ENSMUSG00000028004','ENSMUSG00000019768'), ncol = 1, pt.size = 0) + theme(legend.position = 'none')
# 
#Find non-Npy2r markers across cell types
DefaultAssay(subset.mySeurat) <- "RNA"
subset.mySeurat <- NormalizeData(subset.mySeurat)

marks <- FindAllMarkers(subset.mySeurat, assay="RNA", logfc.threshold=0.5, min.pct=0.1, min.diff.pct=0.25, only.pos=TRUE)

#Add gene symbols to marker list and save to file
Symbols <- subset.mySeurat[["RNA"]]
GeneSymbols <- Symbols@meta.features
GeneSymbols.df <- data.frame(GeneSymbols)
GeneSymbols.df$gene <- rownames(GeneSymbols.df)
marks.ann <- merge(GeneSymbols.df, marks, by="gene")


#Just for fun recluster it using Joe's parameters
SCT <- SCTransform(subset.mySeurat)
