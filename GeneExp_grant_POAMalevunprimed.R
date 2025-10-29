#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
data <- read.table("POA_MvU_draft1.5.txt", header=TRUE)
head(data)
pdf("POA_Esr1Adult_MalevUnprimed1.5.pdf")
plot <- ggplot(data, aes(x=Male, y=Unprimed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Blue","Gray","Green"))+scale_alpha_manual(values=c(1,0.1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Log2 Counts (Male)", y="Log2 Counts (Unprimed Female)") + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ coord_equal(ratio=1)
plot
dev.off()
