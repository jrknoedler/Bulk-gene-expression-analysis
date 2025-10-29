#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
data <- read.table("MeA_MvP_MA.txt", header=TRUE)
head(data)
neg <- data[which(data$fc < 0),]
pos <- data[which(data$fc > 0),]
pos[,3] <- 2^(pos[,3])
neg[,3] <- abs(neg[,3])
neg[,3] <- 2^(neg[,3])
neg[,3] <- -1*neg[,3]
data2 <- rbind(pos,neg)
#data[,3] <- 2^(data[,3])
#head(data)
pdf("MeA_Esr1Adult_MalevPrimedMAplot.pdf")
plot <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Red","Blue","Gray"))+xlim(0.9,5)+ylim(-10,10)+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank())+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_text(color="black", size=20)) +theme(legend.position="none") + geom_hline(yintercept=-1.5, linetype="dashed", color="black", size=1) + geom_hline(yintercept=1.5, linetype="dashed", color="black",size=1)+ geom_hline(yintercept=2, linetype="dashed", color="dark gray",size=0.5) + geom_hline(yintercept=-2, linetype="dashed", color="dark gray",size=0.5) + geom_hline(yintercept=-4, linetype="dashed", color="dark gray",size=0.5) + geom_hline(yintercept=4, linetype="dashed", color="dark gray",size=0.5) + geom_hline(yintercept=-8, linetype="dashed", color="dark gray",size=0.5) + geom_hline(yintercept=8, linetype="dashed", color="dark gray",size=0.5) 
plot
dev.off()
