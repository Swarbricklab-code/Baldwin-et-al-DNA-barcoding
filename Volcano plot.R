
library(ggplot2)
library(ggrepel)
timeStamp <- format(Sys.time(), "%Y_%m_%d")
species <- "mouse"
ensVer <- 84


#setwd("~/contrib1/nenbar/projects/simon/scripts")
homedir="/share/ScratchGeneral/nenbar"
#homedir="../../../../"



######## directory structure #######
projectDir=paste0(homedir,"/projects/simon")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures")


data1=read.table("IE1vsBulk.txt",header=T,sep="\t")

#volcano plot
thresholdDE<-data1$FDR<=0.01
data1$threshold <- thresholdDE 
data1_ordered <- data1[order(data1$FDR), ] 
data1_ordered$genelabels <- ""
data1_ordered$genelabels[1:50] <- data1_ordered$Gene[1:50]
#repeats_Lx9CVB4_vs_CVB4_ordered$genelabels[repeats_Lx9CVB4_vs_CVB4_ordered$logFC<0]<-""
pdf(paste0(imageDir,"/volcano_IE1vsBulk_2.pdf"),width=12,height=8)
ggplot(data1_ordered) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  ggtitle("IE1") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_text_repel(size = 6,aes(x = logFC, y = -log10(FDR), label = genelabels)) +
  xlim(-12,12) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

dev.off()


data1=read.table("IE2vsBulk.txt",header=T,sep="\t")

#volcano plot
thresholdDE<-data1$FDR<=0.01
data1$threshold <- thresholdDE 
data1_ordered <- data1[order(data1$FDR), ] 
data1_ordered$genelabels <- ""
data1_ordered$genelabels[1:50] <- data1_ordered$Gene[1:50]
#repeats_Lx9CVB4_vs_CVB4_ordered$genelabels[repeats_Lx9CVB4_vs_CVB4_ordered$logFC<0]<-""
pdf(paste0(imageDir,"/volcano_IE2vsBulk_2.pdf"),width=12,height=8)
ggplot(data1_ordered) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  ggtitle("IE2") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_text_repel(size = 6,aes(x = logFC, y = -log10(FDR), label = genelabels)) +
  xlim(-12,12) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

dev.off()








