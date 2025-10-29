#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
data <- read.table("Simdata_adjusted.txt", header=TRUE)
head(data)
png("Simdata.png")
plot <- ggplot(data, aes(x=Primed, y=Unprimed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Gray","Green","Orange"))+scale_alpha_manual(values=c("0.1","1","1"))+theme(panel.background = element_rect(fill = "black", color = "black", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+coord_equal(ratio=1)+theme(plot.background= element_rect(fill="black"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="white"))+labs(x="Primed Log2 TPM", y="Unprimed Log2 TPM") + theme(axis.title.x = element_text(color="white", size=20)) + theme(axis.title.y = element_text(color="white", size=20))+theme(axis.text.x=element_text(size=16, color="white"), axis.text.y = element_text(color="white", size=16))+theme(legend.position="none")
plot
dev.off()
