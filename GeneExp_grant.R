#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
data <- read.table("BNST_MvP_TPM_Log10pluseone.txt", header=TRUE)
head(data)
pdf("BNST_Esr1Adult_Malevprimed_TPMlog10.pdf")
plot <- ggplot(data, aes(x=Male, y=Female)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Red","Blue","Gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Male TPM (Log10+1)", y="Female TPM (Log10+1)") + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ coord_equal(ratio=1)+geom_text_repel(data=subset(data, Gene=="Rps6ka6"), point.padding=1, label="Rps6ka6", nudge_y=-1, nudge_x=5, size=10, aes(fontface="bold"))+geom_text_repel(data=subset(data, Gene=="Cckar"), point.padding=1, label="Cckar", nudge_y=0, nudge_x=2, size=10, aes(fontface="bold"))+geom_text_repel(data=subset(data, Gene=="Greb1"), nudge_x=2, nudge_y=-1.5, point.padding=1, label="Greb1", size=10, aes(fontface="bold"))+geom_text_repel(data=subset(data, Gene=="Cyp19a1"), point.padding=1, label="Cyp19a1", size=10, nudge_x=1, nudge_y=-3, aes(fontface="bold"))+geom_text_repel(data=subset(data, Gene=="Glra3"), point.padding=1, label="Glra3", size=10, nudge_x=2, nudge_y=-0.75, aes(fontface="bold")) + geom_text_repel(data=subset(data, Gene=="Sytl4"), point.padding=1, label="Sytl4", size=10, nudge_x=2, nudge_y=0.5, aes(fontface="bold")) 
plot
dev.off()
