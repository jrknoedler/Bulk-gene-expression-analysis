#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(patchwork)
library(dplyr)
data <- read.table("BNST_MvP_TPM_Log10pluseone.txt", header=TRUE)


pdf(file="FigureS6_Tac1.pdf")
ggplot(data, aes(y=Male, x=Female)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("deeppink1","Blue","Gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_blank(), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1) +xlim(0,5)+ylim(0,5) +theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))+geom_text_repel(data=subset(data, Gene=="Tac1"), point.padding=1, color="white", segment.color="black", label="Tac1", nudge_y=1, nudge_x=-1, size=25, aes(fontface="plain")) + geom_point(data=data %>% filter(Gene == "Tac1"), pch=21, size=6, color="black", stroke=2)
dev.off()