#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(patchwork)

data <- read.table("P15_BNST_July2022.txt", header=TRUE)

BNST <- ggplot(data, aes(y=Female, x=Male)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("blue","orange","gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_blank(), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1) + ylim(min=0,5)  +theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))


data <- read.table("P15_MeA_July2022.txt", header=TRUE)
MeA <- ggplot(data, aes(y=Female, x=Male)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("blue","orange","gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_blank(), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1) + ylim(min=0,5)  +theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))
data <- read.table("P15_POA_July2022.txt", header=TRUE)
POA <- ggplot(data, aes(y=Female, x=Male)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("blue","orange","gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_blank(), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1) + ylim(min=0,5)  +theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))

data <- read.table("P15_VMH_July2022.txt", header=TRUE)
blank <- ggplot() +theme_void()
VMH <- ggplot(data, aes(y=Female, x=Male)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("blue","orange","gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_blank(), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1) + ylim(min=0,5)  +theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))
Row1 <- (BNST |MeA| POA|VMH)
Row2 <- (BNST |MeA|POA|VMH) 
pdf(file="P15TPMPlots.pdf", width=22, height=5.5)
Row1 & theme(axis.text.y=element_text(color="white", size=11))
dev.off()
