#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
data <- read.table("AdultVMH_MvP_Tpm_Log10pluseone.txt", header=TRUE)
head(data)
pdf("VMH_Esr1Adult_Malevprimed_TPMlog10_Brandeis.pdf")
plot <- ggplot(data, aes(x=Male, y=Primed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Red","Green","Gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "black", color = "black", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="black"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="white"))+labs(x="Male TPM (Log10+1)", y="Female TPM (Log10+1)") + xlim(0,5) + theme(axis.title.x = element_text(color="white", size=20)) + theme(axis.title.y = element_text(color="white", size=20))+theme(axis.text.x=element_text(size=16, color="white"), axis.text.y = element_text(color="white", size=16))+theme(legend.position="none")+ coord_equal(ratio=1)+ geom_text_repel(data=subset(data, Gene=="Cckar"), label="Cckar", color="white", nudge_x=-0.5, nudge_y=0.5, size=16)
plot
dev.off()
