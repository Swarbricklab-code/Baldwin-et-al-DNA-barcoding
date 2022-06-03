# Script information

# author L Baldwin, N Bartonicek
# input files: matrix_fraction.csv -> a file from the output of ClonTracer barcode count.py script.
# outputs: pdf ggplots - unique barcode counts per group, shannon diversity per group

# use in paper fig2F, 2G, 3D, 3E, 3F


# ---------------------------------

# Session set up

library(tidyverse)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(Hmisc)
library(corrplot)
library(corrr)
library(pvclust)
library(EntropyExplorer)
library(reshape2)
library(gridExtra)
library(UpSetR)
library(ComplexHeatmap)

figDir=paste0("../figures/")
tabDir=paste0("../tables")

cutoff=0.0001

# ##########
# #data tidying step - only necessary when ClonTracer has been run lane by lane, not sample by sample

# matrix <- read.table("matrix_count.txt")
 
# #need to collapse columns by seq lane read

# colnames <- colnames(matrix)
# colnames <- str_replace(colnames, pattern = "_L001.ds", "")
# colnames <- str_replace(colnames, pattern = "_L002.ds", "")
# colnames <- str_replace(colnames, pattern = "_L003.ds", "")
# colnames <- str_replace(colnames, pattern = "_L004.ds", "")

# colnames(matrix) <- colnames
# #collapse by column name
# test <- pivot_longer(matrix, cols=2:37, names_to = "sample", values_to = "counts")
# test2 <- test %>% group_by(Barcode, sample) %>% summarise_all(sum)
# new_matrix <- pivot_wider(test2, names_from = sample, values_from = counts)

# write.csv(x = new_matrix, file = "./new_matrix_count.csv" )
########


#read in  matrix
dataSample <- read.csv("./matrix_count.csv", header = T, row.names = 1)
dataSample<-as.matrix(dataSample)
dataSampleShort<-dataSample

#get fractions
cellCounts<-colSums(dataSampleShort)

dataSampleShortfrac<-as.data.frame(dataSampleShort %*% diag(1/cellCounts))
colnames(dataSampleShortfrac)<-colnames(dataSampleShort)
dataSampleShortfrac$barcode<-row.names(dataSampleShort)

#melt
dataM<-melt(dataSampleShortfrac)

#merge with IDs
#ids$class<-gsub(".*_","",ids$sampleName)
#merged<-merge(dataM,ids,by.x="variable",by.y="id")
merged<-dataM

#Plot raw counts across samples
raw<-data.frame(id=names(cellCounts),counts=as.numeric(cellCounts))

pdf(paste0(figDir,"LB11_4_raw_counts.pdf"),width=12,height=5)
p<-ggplot(raw,aes(sampleName,counts,fill=sampleName))
p<-p+geom_col()+scale_fill_manual(values=cols)
p<-p+theme(axis.text.x=element_text(angle =- 90, vjust = 0.5))
p<-p+ggtitle(paste0("Raw counts across samples"))
print(p)
dev.off()

dataSampleShortBinary<-dataSampleShort
dataSampleShortBinary[,]<-ifelse(dataSampleShortBinary %in% c("0"),FALSE,TRUE)
uniqueCellCounts<-colSums(dataSampleShortBinary)

raw<-data.frame(id=names(uniqueCellCounts),counts=as.numeric(uniqueCellCounts))
mergedRaw<-merge(raw,ids,by.x="id",by.y="id")
mergedRaw$sampleName<-factor(mergedRaw$sampleName,levels=ids$sampleName[order(ids$class)])

pdf(paste0(figDir,"LB11_4_unique_counts.pdf"),width=12,height=5)
p<-ggplot(mergedRaw,aes(sampleName,counts,fill=sampleName))
p<-p+geom_col()+scale_fill_manual(values=cols)
p<-p+theme(axis.text.x=element_text(angle =- 90, vjust = 0.5))
p<-p+ggtitle(paste0("Unique counts across samples"))
print(p)
dev.off()


## Shannon diveristy
#the paramater shift defines how much is to be added to zero-values for algorithm to work. It influences mostly the 
#margins of the distribution

shift=1e-5
dataAggregateShort1=dataSampleShort[,c(1,2,3,4,5)]
dataAggregateShort2=dataSampleShort[,c(1,2,3,4,5,6)+5]
dataAggregateShort3=dataSampleShort[,c(1,2,3,4,5,6)+11]
entropy1=EntropyExplorer(dataAggregateShort1,dataAggregateShort2,dataAggregateShort3, dmetric="dse",otype="ba",shift=c(shift,shift))
hist(entropy1[,3],xlim=c(-1,1),breaks=40)


#diversity function
shannon<-function(data){
  H<-numeric(nrow(data))
  mult<-numeric(length(data[1,]))
  
  for(i in 1:nrow(data)){
    prop<-data[i,]/sum(data[i,])
    
    for(j in 1:ncol(data)){
      mult[j]<-prop[j]*log(prop[j])
    }
    
    H[i]<--sum(mult, na.rm=T)
  }
  plot.number<-1:nrow(data)
  return(rbind(plot.number,H))
}

t_IT=shannon(t(dataAggregateShort1))[2,]
t_CONT=shannon(t(dataAggregateShort2))[2,]
t_NSG=shannon(t(dataAggregateShort3))[2,]

df<-rbind(t_NSG,t_CONT,t_IT)
df=t(df)
dataM<-melt(df)

ggplot(dataM,aes(Var2,value))+geom_violin()
ggplot(dataM,aes(Var2,value))+geom_point()
ggplot(dataM,aes(Var2,value))+geom_boxplot()

sampleList<-list()
sampleList[["IT"]]<-t_IT
sampleList[["CONT"]]<-t_CONT
sampleList[["NSG"]]<-t_NSG

df<-data.frame(value=unlist(sampleList),sampleName=rep(names(sampleList),times=sapply(sampleList,length)))
df$sampleName<-factor(df$sampleName,levels=c("NSG","CONT","IT"))
ggplot(df,aes(sampleName,value))+geom_violin()


#----