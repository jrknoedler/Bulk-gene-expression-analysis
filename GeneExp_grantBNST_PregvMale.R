#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
data <- read.table("BNST_PregvMale_TPM.txt", header=TRUE)
head(data)
pdf("BNST_Esr1Adult_PregvMale_TPMlog10_pink.pdf")
plot <- ggplot(data, aes(x=Male, y=Preg)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("blue","Gray","orange"))+scale_alpha_manual(values=c(1,0.11,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Male (Log10 TPM)", y="Pregnant (Log10 TPM)") + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ coord_equal(ratio=1) + ylim(0,5)
plot
dev.off()
