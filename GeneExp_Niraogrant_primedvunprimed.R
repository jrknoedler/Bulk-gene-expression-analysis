#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
data <- read.table("L10_Grantplot.txt", header=TRUE)
head(data)
png("VMH_L10_PrimedvUnprimed_publabels.png")
plot <- ggplot(data, aes(x=Input, y=IP)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Gray","Red","Green"))+scale_alpha_manual(values=c("0.1","1","1"))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Log2 Counts Input)", y="Log2 Counts IP") + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+geom_text_repel(data=subset(data, Gene=="Cckar"), point.padding=1, label="Cckar", nudge_y=2, nudge_x=-1, size=10, aes(fontface="bold"))+geom_text_repel(data=subset(data, Gene=="Pgr"), nudge_y=2, nudge_x=-4, point.padding=1, label="PR", size=10, aes(fontface="bold")) +geom_text_repel(data=subset(data, Gene=="Greb1"), point.padding=1, label="Greb1", nudge_y=0, nudge_x=-4, size=10, aes(fontface="bold"))
plot
dev.off()
