#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
data <- read.table("L10_IPvInput_Primed.txt", header=TRUE)
head(data)
png("L10_Primed_InputvIP.png")
plot <- ggplot(data, aes(x=Input, y=IP)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Red","Blue","Gray"))+scale_alpha_manual(values=c("1","1","0.1"))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+coord_equal(ratio=1)+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Input Normalized counts", y="IP Normalized counts") + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ geom_text(data=subset(data, Gene=="Esr1"), label="Esr1", size=10, nudge_x = -2, nudge_y=2) + geom_text(data=subset(data, Gene=="Gfap"), label="Gfap") + geom_text(data=subset(data, Gene=="Olig1"), label="Olig1") + geom_text(data=subset(data, Gene=="Pgr"), label="Pgr") + geom_text(data=subset(data, Gene=="Mobp"), label="Mobp") 
plot
dev.off()
