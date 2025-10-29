#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
data <- read.table("VMH_MvU_1.5draft.txt", header=TRUE)
head(data)
pdf("VMH_Esr1Adult_MalevUnprimed1.5_brandeis.pdf")
plot <- ggplot(data, aes(x=Male, y=Diestrus)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Green","Gray","Orange"))+scale_alpha_manual(values=c(1,0.1,1))+theme(panel.background = element_rect(fill = "black", color = "black", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="black"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="white"))+labs(x="Log2 Counts (Male)", y="Log2 Counts (Unprimed Female)") + theme(axis.title.x = element_text(color="white", size=20)) + theme(axis.title.y = element_text(color="white", size=20))+theme(axis.text.x=element_text(size=16, color="white"), axis.text.y = element_text(color="white", size=16))+theme(legend.position="none")+ coord_equal(ratio=1)+ geom_text_repel(data=subset(data, Row.names=="Cckar"), label="Cckar", color="white", nudge_x=0.5, nudge_y=-0.5)
plot
dev.off()
