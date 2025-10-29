#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
data <- read.table("AdultPOA_MvP_Tpm.csv", header=TRUE)
head(data)
pdf("POA_Esr1Adult_Malevprimed_TPMlog10_nolabels.pdf")
plot <- ggplot(data, aes(x=Male, y=Female)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Red","Blue","Gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Male TPM (Log10+1)", y="Female TPM (Log10+1)") + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ coord_equal(ratio=1)
plot
dev.off()
