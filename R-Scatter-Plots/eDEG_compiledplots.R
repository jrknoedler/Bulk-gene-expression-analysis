#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(patchwork)

data <- read.table("BNST_PvUTpm.txt", header=TRUE)

BNST <- ggplot(data, aes(y=Primed, x=Unprimed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Gray","deeppink1","Green"))+scale_alpha_manual(values=c(0.1,1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_text(size=30, color="black"), axis.text.y = element_text(size=30, color="black"))+theme(legend.position="none")+ coord_equal(ratio=1)+xlim(0,5)+ylim(0,5)+theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))


data <- read.table("MeAPvUTPMedited.txt", header=TRUE)
MeA <- ggplot(data, aes(y=Primed, x=Unprimed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Gray","deeppink1","Green"))+scale_alpha_manual(values=c(0.1,1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_text(color="black", size=30), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1)+xlim(0,5)+ylim(0,5)+theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))
data <- read.table("POA_PvUTpm.txt", header=TRUE)
POA <- ggplot(data, aes(y=Primed, x=Unprimed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Gray","deeppink1","Green"))+scale_alpha_manual(values=c(0.1,1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_text(color="black", size=30), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1)+xlim(0,5)+ylim(0,5)+theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))

data <- read.table("VMU_PvUTpm.txt", header=TRUE)
blank <- ggplot() +theme_void()
VMH <- ggplot(data, aes(y=Primed, x=Unprimed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Gray","deeppink1","Green"))+scale_alpha_manual(values=c(0.1,1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_text(color="black", size=30), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1)+xlim(0,5)+ylim(0,5)+theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))
pdf(file="Fig2eDEGcompiled.pdf", width=25, height=10)
(BNST |blank| MeA|blank | POA|blank |VMH)/(BNST |blank| MeA|blank | POA|blank |VMH)  + plot_layout(widths=c(1,0.01,1,0.01,1,0.01,1,1,0.01,1,0.01,1,0.01,1))
dev.off()
