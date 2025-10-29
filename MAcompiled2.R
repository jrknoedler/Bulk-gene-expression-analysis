#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(patchwork)
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

BNSTMvP <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category, shape=level))+scale_shape_manual(values=c(16,2))+scale_color_manual(values=c("deeppink1","royalblue3","Gray"))+xlim(0.9,5)+coord_cartesian(ylim = c(-3, 3))+scale_alpha_manual(values=c(1,1,0.2))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_text(color="white",size=30),  axis.text.x=element_blank(), axis.text.y=element_text(size=40, color="black"), axis.title.y=element_text(color="white", size=30))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+theme(legend.position="none") + geom_hline(yintercept=-0.584963, linetype="dashed", color="black",size=0.75) + geom_hline(yintercept=0.584963, linetype="dashed", color="black",size=0.75)+theme(axis.ticks.length = unit(0.33, "cm")) + scale_y_continuous(breaks=seq(-3,3,1))

data <- read.table("BNST_MvU_MA.txt", header=TRUE)
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

BNSTMvU <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category, shape=level))+scale_shape_manual(values=c(2,16))+scale_color_manual(values=c("royalblue3","Gray","chartreuse3"))+xlim(0.9,5)+ylim(-3,3)+scale_alpha_manual(values=c(1,0.2,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_blank(), axis.text.y=element_text(color="black", size=40), axis.text.x=element_blank(), axis.title.y=element_blank())+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_blank()) +theme(legend.position="none") + geom_hline(yintercept=-0.584963, linetype="dashed", color="black",size=0.75)  + geom_hline(yintercept=0.584963, linetype="dashed", color="black",size=0.75) +theme(axis.ticks.length = unit(0.33, "cm")) + scale_y_continuous(breaks=seq(-3,3,1)) 

data <- read.table("BNST_PvU_MA.txt", header=TRUE)
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

BNSTPvU <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category, shape=level))+scale_shape_manual(values=c(2,16))+scale_color_manual(values=c("Gray","deeppink1","chartreuse3"))+xlim(0.9,5)+coord_cartesian(ylim = c(-3, 3))+scale_alpha_manual(values=c(0.2,1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_blank(), axis.text.y=element_text(color="black", size=40), axis.text.x=element_text(color="black", size=40), axis.title.y=element_blank())+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_blank()) +theme(legend.position="none") + geom_hline(yintercept=-0.584963, linetype="dashed", color="black",size=0.75) + geom_hline(yintercept=0.584963, linetype="dashed", color="black",size=0.75) + theme(axis.ticks.length = unit(0.33, "cm")) + scale_y_continuous(breaks=seq(-3,3,1))


data <- read.table("MeA_MvP_MA.txt", header=TRUE)

#expo <- function(x){2^x}

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


MeAMvP <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category, shape=level))+scale_shape_manual(values=c(2,16))+scale_color_manual(values=c("deeppink1","royalblue3","Gray"))+xlim(0.9,5)+coord_cartesian(ylim = c(-3, 3))+scale_alpha_manual(values=c(1,1,0.2))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_blank(),  axis.text.x=element_blank(),axis.text.y=element_text(size=40, color="white"), axis.title.y=element_blank())+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_blank()) +theme(legend.position="none") + geom_hline(yintercept=-0.584963, linetype="dashed", color="black",size=0.75) + geom_hline(yintercept=0.584963, linetype="dashed", color="black",size=0.75)+ theme(axis.ticks.length = unit(0.33, "cm")) + scale_y_continuous(breaks=seq(-3,3,1))


data <- read.table("MeA_MvU_MA.txt", header=TRUE)
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

MeAMvU <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category, shape=level))+scale_shape_manual(values=c(2,16))+scale_color_manual(values=c("royalblue3","Gray","chartreuse3"))+xlim(0.9,5)+ylim(-3,3)+scale_alpha_manual(values=c(1,0.2,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_blank(),  axis.text.x=element_blank(), axis.text.y=element_text(size=40, color="white"), axis.title.y=element_blank())+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_blank()) +theme(legend.position="none")  + geom_hline(yintercept=0.584963, linetype="dashed", color="black",size=0.75) + geom_hline(yintercept=-0.584963, linetype="dashed", color="black",size=0.75) + theme(axis.ticks.length = unit(0.33, "cm")) + scale_y_continuous(breaks=seq(-3,3,1))


data <- read.table("MeA_PvU_MA.txt", header=TRUE)
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

MeAPvU <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category, shape=level))+scale_shape_manual(values=c(2,16))+scale_color_manual(values=c("Gray","deeppink1","chartreuse3"))+xlim(0.9,5)+coord_cartesian(ylim = c(-3, 3))+scale_alpha_manual(values=c(0.2,1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_blank(), axis.text.y=element_text(size=40, color="white"), axis.text.x=element_text(color="black", size=40), axis.title.y=element_blank())+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_blank()) +theme(legend.position="none") + geom_hline(yintercept=-0.584963, linetype="dashed", color="black",size=0.75) + geom_hline(yintercept=0.584963, linetype="dashed", color="black",size=0.75) + theme(axis.ticks.length = unit(0.33, "cm")) + scale_y_continuous(breaks=seq(-3,3,1))


data <- read.table("POA_MvP_MA.txt", header=TRUE)
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

POAMvP <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category, shape=level))+scale_shape_manual(values=c(2,16))+scale_color_manual(values=c("deeppink1","royalblue3","Gray"))+xlim(0.9,5)+ylim(-3,3)+scale_alpha_manual(values=c(1,1,0.2))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_blank(),  axis.text.x=element_blank(), axis.text.y=element_text(size=40, color="white"), axis.title.y=element_blank())+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_blank()) +theme(legend.position="none") + geom_hline(yintercept=-0.584963, linetype="dashed", color="black",size=0.75) + geom_hline(yintercept=0.584963, linetype="dashed", color="black",size=0.75) + theme(axis.ticks.length = unit(0.33, "cm")) + scale_y_continuous(breaks=seq(-3,3,1))


data <- read.table("POA_MvU_MA.txt", header=TRUE)
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

POAMvU <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category, shape=level))+scale_shape_manual(values=c(2,16))+scale_color_manual(values=c("royalblue3","Gray","chartreuse3"))+xlim(0.9,5)+ylim(-3,3)+scale_alpha_manual(values=c(1,0.2,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_blank(),  axis.text.x=element_blank(), axis.text.y=element_text(size=40, color="white"), axis.title.y=element_blank())+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_blank()) +theme(legend.position="none")+ geom_hline(yintercept=-0.584963, linetype="dashed", color="black",size=0.75) + geom_hline(yintercept=0.584963, linetype="dashed", color="black",size=0.75) + theme(axis.ticks.length = unit(0.33, "cm")) + scale_y_continuous(breaks=seq(-3,3,1))


data <- read.table("POA_PvU_MA.txt", header=TRUE)
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
POAPvU <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category, shape=level))+scale_shape_manual(values=c(2,16))+scale_color_manual(values=c("Gray","deeppink1","chartreuse3"))+xlim(0.9,5)+ylim(-3,3)+scale_alpha_manual(values=c(0.2,1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_blank(),axis.text.x=element_text(color="black", size=40), axis.text.y=element_text(size=40, color="white"), axis.title.y=element_blank())+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_blank()) +theme(legend.position="none") + geom_hline(yintercept=0.584963, linetype="dashed", color="black",size=0.75) + geom_hline(yintercept=-0.584963, linetype="dashed", color="black",size=0.75) + theme(axis.ticks.length = unit(0.33, "cm")) + scale_y_continuous(breaks=seq(-3,3,1))


data <- read.table("VMH_MvP_MA.txt", header=TRUE)
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

VMHMvP <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category, shape=level))+scale_shape_manual(values=c(2,16))+scale_color_manual(values=c("deeppink1","royalblue3","Gray"))+xlim(0.9,5)+ylim(-3,3)+scale_alpha_manual(values=c(1,1,0.2))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_blank(),  axis.text.x=element_blank(), axis.text.y=element_text(size=40, color="white"), axis.title.y=element_blank())+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_blank()) +theme(legend.position="none") + geom_hline(yintercept=0.584963, linetype="dashed", color="black",size=0.75) + geom_hline(yintercept=-0.584963, linetype="dashed", color="black",size=0.75) + theme(axis.ticks.length = unit(0.33, "cm")) + scale_y_continuous(breaks=seq(-3,3,1))


data <- read.table("VMH_MvU_MA.txt", header=TRUE)
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

VMHMvU <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category, shape=level))+scale_shape_manual(values=c(2,16))+scale_color_manual(values=c("royalblue3","Gray","chartreuse3"))+xlim(0.9,5)+ylim(-3,3)+scale_alpha_manual(values=c(1,0.2,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_blank(),  axis.text.x=element_blank(), axis.text.y=element_text(size=40, color="white"), axis.title.y=element_blank())+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_blank()) +theme(legend.position="none")+ geom_hline(yintercept=-0.584963, linetype="dashed", color="black",size=0.75) + geom_hline(yintercept=0.584963, linetype="dashed", color="black",size=0.75) + theme(axis.ticks.length = unit(0.33, "cm")) + scale_y_continuous(breaks=seq(-3,3,1)) 


data <- read.table("VMH_PvU_MA.txt", header=TRUE)
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

VMHPvU <- ggplot(data2, aes(x=mean, y=fc)) + geom_point(aes(color=Category,alpha=Category, shape=level))+scale_shape_manual(values=c(2,16))+scale_color_manual(values=c("Gray","deeppink1","chartreuse3"))+xlim(0.9,5)+ylim(-3,3)+scale_alpha_manual(values=c(0.2,1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"), axis.title.x=element_blank(), axis.text.x=element_text(color="black", size=40), axis.text.y=element_text(size=40, color="white"), axis.title.y=element_blank())+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+ theme(axis.title.x = element_blank()) +theme(legend.position="none") + geom_hline(yintercept=-0.584963, linetype="dashed", color="black",size=0.75) + geom_hline(yintercept=0.584963, linetype="dashed", color="black",size=0.75) + theme(axis.ticks.length = unit(0.33, "cm")) + scale_y_continuous(breaks=seq(-3,3,1))  


#blank <- ggplot() +theme_void()
#Row1 <- (BNSTMvP |blank| MeAMvP | blank| POAMvP | blank| VMHMvP) + plot_layout(widths=c(1,0.2,1,0.2,1,0.2,1))  
#Row4 <- (blank|blank|blank|blank|blank|blank|blank) + plot_layout(widths=c(1,0.01,1,0.01,1,0.01,1))
#Row2 <- (BNSTMvU | blank| MeAMvU | blank| POAMvU | blank| VMHMvU)+ plot_layout(widths=c(1,0.2,1,0.2,1,0.2,1))  
#Row5 <- (blank|blank|blank|blank|blank|blank|blank) + plot_layout(widths=c(1,0.01,1,0.01,1,0.01,1))
#Row3 <- (BNSTPvU | blank| MeAPvU | blank| POAPvU|blank| VMHPvU) + plot_layout(widths=c(1,0.2,1,0.2,1,0.2,1)) 
Row1 <- (BNSTMvP | MeAMvP |POAMvP |  VMHMvP) 
Row2 <- (BNSTMvU |  MeAMvU | POAMvU | VMHMvU)
Row3 <- (BNSTPvU | MeAPvU | POAPvU| VMHPvU)
full <- (Row1/Row2/Row3) 
pdf(file="MAcompiled.pdf", width=28, height=18)
full & theme(aspect.ratio=5/6) & theme(axis.text.x=element_blank()) & theme(axis.title.y=element_text(color="white", size=30), axis.title.x=element_text(color="white", size=30)) & theme(axis.text.y=element_blank()) 
dev.off()

#(BNSTMvP | MeAMvP | POAMvP | VMHMvP)/(BNSTMvU | MeAMvU | POAMvU | VMHMvU)/(BNSTPvU | MeAPvU | POAPvU|VMHPvU) & theme(aspect.ratio=5/6) + plot_layout(heights=c(1,0.01,1,0.01,1))
#dev.off()