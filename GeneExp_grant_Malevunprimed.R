#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
data <- read.table("VMH_Esr1Adult_MalevUnprimed.txt", header=TRUE)
head(data)
png("VMH_Esr1Adult_MalevUnprimed.png")
plot <- ggplot(data, aes(x=Male, y=Unprimed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Blue","Gray","Green"))+scale_alpha_manual(values=c("1","0.1","1"))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Log2 Counts (Male)", y="Log2 Counts (Unprimed Female)") + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ coord_equal(ratio=1)+geom_text_repel(data=subset(data, Gene=="Gpr88"), point.padding=1, label="Gpr88", nudge_y=2, nudge_x=-3, size=10, aes(fontface="bold"))+geom_text_repel(data=subset(data, Gene=="Moxd1"), nudge_y=-2, nudge_x=1, point.padding=1, label="Moxd1", size=10, aes(fontface="bold"))+geom_text_repel(data=subset(data, Gene=="Greb1"), point.padding=1, label="Greb1", size=10, nudge_x=3, nudge_y=1, aes(fontface="bold"))+geom_text_repel(data=subset(data, Gene=="Cyp19a1"), point.padding=1, label="Cyp19a1", size=10, nudge_x=2, aes(fontface="bold"))+geom_text_repel(data=subset(data, Gene=="Rprm"), point.padding=1, label="Rprm", size=10, nudge_x=-2, aes(fontface="bold"))+geom_text_repel(data=subset(data, Gene=="Tac2"), point.padding=1, label="Tac2", size=10, nudge_y=2, aes(fontface="bold"))
plot
dev.off()
