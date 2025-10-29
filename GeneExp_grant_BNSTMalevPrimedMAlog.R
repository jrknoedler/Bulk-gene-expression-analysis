#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
data <- read.table("BNST_MvP_MA.txt", header=TRUE)
#expo <- function(x){2^x}
head(data)


head(data)
#data[,3] <- as.data.frame(lapply(data[,3], FUN=function(x){sapply(x, FUN=expo)}))
#neg <- data[which(data$fc < 0),]
#pos <- data[which(data$fc > 0),]
#pos[,3] <- 2^(pos[,3])
#neg[,3] <- abs(neg[,3])
#neg[,3] <- 2^(neg[,3])
#neg[,3] <- -1*neg[,3]
#data2 <- rbind(pos,neg)
#data[,3] <- 2^(data[,3])
#head(data)
data$order <- ifelse(data$Category=="NC", 1,2)
data$level <- ifelse(data$fc > 3 | data$fc < -3,"outlier","regular") 
outliers <- data[data$level=="outlier",]
outliers
head(data)
data[data$fc>3,"fc"] <- 3
min <- -3
data$fc[data$fc < min] <- -3
head(data)
data2 <- data[order(as.numeric(factor(data$order))),]
head(data2)
pdf(file="BNST_MalevPrimedMAplotlog.pdf")
plot <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category, shape=level))+scale_shape_manual(values=c(16,2))+scale_color_manual(values=c("Red","Blue","Gray"))+xlim(0.9,5)+ylim(-3,3)+scale_alpha_manual(values=c(1,1,0.2))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank())+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_text(color="black", size=20)) +theme(legend.position="none") + geom_hline(yintercept=-0.584963, linetype="dashed", color="black",size=0.75) + geom_hline(yintercept=0.584963, linetype="dashed", color="black",size=0.75)
plot
dev.off()
