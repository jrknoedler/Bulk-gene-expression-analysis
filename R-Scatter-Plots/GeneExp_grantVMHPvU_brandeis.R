#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
data <- read.table("VMU_PvUTpm.txt", header=TRUE)
head(data)
pdf("VMH_Esr1Adult_PrimedvUnprimed_TPMlog10_brandeis.pdf")
plot <- ggplot(data, aes(x=Primed, y=Unprimed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Gray","red","yellow"))+scale_alpha_manual(values=c(0.1,1,1))+theme(panel.background = element_rect(fill = "black", color = "black", linetype="solid")) 
+ theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="black"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="white"))+labs(x="Primed TPM (Log10+1)", y="Unprimed TPM (Log10+1)") + xlim(0,5) + theme(axis.title.x = element_text(color="white", size=20)) + theme(axis.title.y = element_text(color="white", size=20))+theme(axis.text.x=element_text(size=16, color="white"), axis.text.y = element_text(color="white", size=16))+theme(legend.position="none")+ coord_equal(ratio=1) + ylim(0,5)
plot
dev.off()
