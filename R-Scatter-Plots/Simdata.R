#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
data <- read.table("Castrate2reps.tsv", header=TRUE)
head(data)
png("Simdata.png")
plot <- ggplot(data, aes(x=Rep1, y=Rep2)) + geom_point(aes(color=gray,alpha=0.1))
plot
dev.off()
