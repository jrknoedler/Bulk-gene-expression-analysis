
#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(patchwork)

data <- read.table("POA_P15vAdult.txt", header=TRUE)
Male <- ggplot(data, aes(y=Adult, x=P15)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Red","Gray","blue"))+scale_alpha_manual(values=c(1,0.1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_blank(), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1) + ylim(min=0,5)  +theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))


data <- read.table("POA_P15vFemale.txt", header=TRUE)
Female <- ggplot(data, aes(y=Adult, x=P15)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("red","gray","blue"))+scale_alpha_manual(values=c(1,0.1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black")) + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="white", size=40))+theme(axis.text.x=element_blank(), axis.text.y = element_blank())+theme(legend.position="none")+ coord_equal(ratio=1) + ylim(min=0,5)  +theme(axis.title.x=element_blank()) +theme(axis.ticks.length = unit(0.33, "cm"))

Row1 <- (Male|Female)

pdf(file="P15vAdult_1stpass.pdf", width=22, height=11)
Row1 & theme(axis.text.y=element_text(color="white", size=11))
dev.off()
