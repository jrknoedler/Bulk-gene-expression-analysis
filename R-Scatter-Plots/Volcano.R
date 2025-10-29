#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
data <- read.table("P15VMH_IPvInput_Volcanodata.txt", header=TRUE)
head(data)
png("P15IPvInputvolcano.png")
plot <- ggplot(data, aes(x=log2FoldChange, y=log10pvalue)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Green","Orange","Gray"))+scale_alpha_manual(values=c("1","1","0.1"))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Log2 Fold Change", y="Log10 FDR-adujusted p value") + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ geom_text_repel(data=subset(data, Gene=="Esr1"), point.padding=1, label="Esr1", size=10, aes(fontface="bold")) + geom_text_repel(data=subset(data, Gene=="Nedd9"), label="Nedd9", point.padding=1, size=10, aes(fontface="bold"))+geom_text_repel(data=subset(data, Gene=="Bcar1"), label="Bcar1", point.padding=1, size=10, aes(fontface="bold"))+geom_text_repel(data=subset(data, Gene=="mt-Co1"), label="mt-Co1", point.padding=1, size=10, aes(fontface="bold")) +geom_text_repel(data=subset(data, Gene=="mt-Co2"), label="mt-Co2", point.padding=1, size=10, aes(fontface="bold")) + geom_text_repel(data=subset(data, Gene=="Cx3cr1"), label="Cx3cr1", size=10, point.padding=1, aes(fontface="bold"))
plot
dev.off()
