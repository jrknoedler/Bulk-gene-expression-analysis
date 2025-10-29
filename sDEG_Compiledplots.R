#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(patchwork)
data <- read.table("AdultVMH_MvP_Tpm_Log10pluseone.txt", header=TRUE)

VMHMvP <- ggplot(data, aes(y=Male, x=Primed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("deeppink1","blue","Gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_blank(), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1)+xlim(0,5)+ylim(0,5)+theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))

data <- read.table("BNST_MvU_1.5tpmDraft.txt", header=TRUE)


BNSTMvU <- ggplot(data, aes(y=Male, x=Diestrus)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Blue","Gray","Green"))+scale_alpha_manual(values=c(1,0.1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_text(size=30, color="black"), axis.text.y = element_text(color="black", size=30))+theme(legend.position="none")+ coord_equal(ratio=1) +xlim(0,5)+ylim(0,5) +theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))

data <- read.table("MeA_MvU_draft1.5.txt", header=TRUE)

MeAMvU <- ggplot(data, aes(y=Male, x=Diestrus)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Blue","Gray","Green"))+scale_alpha_manual(values=c(1,0.1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_text(size=30, color="black"), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1) +xlim(0,5)+ylim(0,5) +theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))

data <- read.table("POA_MvU_draft1.5.txt", header=TRUE)

POAMvU <- ggplot(data, aes(y=Male, x=Unprimed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Blue","Gray","Green"))+scale_alpha_manual(values=c(1,0.1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_text(size=30, color="black"), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1) +xlim(0,5)+ylim(0,5) +theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))

data <- read.table("VMH_MvU_1.5draft.txt", header=TRUE)

VMHMvU <- ggplot(data, aes(y=Male, x=Diestrus)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Blue","Gray","Green"))+scale_alpha_manual(values=c(1,0.1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_text(size=30, color="black"), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1)+xlim(0,5)+ylim(0,5) +theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))

data <- read.table("BNST_MvP_TPM_Log10pluseone.txt", header=TRUE)

BNSTMvP <- ggplot(data, aes(y=Male, x=Female)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("deeppink1","Blue","Gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_blank(), axis.text.y = element_text(color="black", size=30))+theme(legend.position="none")+ coord_equal(ratio=1) +xlim(0,5)+ylim(0,5) +theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))

data <- read.table("MeA_MvP_again.txt", header=TRUE)

MeAMvP <- ggplot(data, aes(y=Male, x=Primed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("deeppink1","Blue","Gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_blank(), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1) + ylim(min=0,5)  +theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))

data <- read.table("AdultPOA_MvP_Tpm.csv", header=TRUE)
blank <- ggplot() +theme_void()
POAMvP <- ggplot(data, aes(y=Male, x=Female)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("deeppink1","Blue","Gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + xlim(0,5)+ylim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_blank(), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1) +xlim(0,5)+ylim(0,5) +theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))
pdf(file="sDEG_Row1.pdf", width=25, height=10)
(BNSTMvP |blank|  MeAMvP |blank|  POAMvP|blank| VMHMvP) + plot_layout(widths=c(1,0.2,1,0.2,1,0.2,1))
dev.off()
Row1 <- (BNSTMvP |blank|  MeAMvP |blank|  POAMvP|blank| VMHMvP) + plot_layout(widths=c(1,0.01,1,0.01,1,0.01,1))
Row2 <- (BNSTMvU |blank|MeAMvU | blank|POAMvU |blank| VMHMvU) + plot_layout(widths=c(1,0.01,1,0.01,1,0.01,1))
pdf(file=paste0("sDEG_Fig1compiled.pdf"), width=25, height=10)
(Row1/Row2)
dev.off()

#pdf(file=paste0("sDEG_Fig1compiled.pdf"), width=25, height=10)
#(BNSTMvP | blank| MeAMvP | blank| POAMvP| blank|  VMHMvP)/(BNSTMvU |blank|MeAMvU | blank|POAMvU |blank| VMHMvU) + plot_layout(widths=c(1,0.2,1,0.2,1,0.2,1,1,0.2,1,0.2,1,0.2,1))
#dev.off()