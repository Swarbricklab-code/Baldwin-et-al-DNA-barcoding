
library(GenomicRanges)
library(ShortRead)
library("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome)
library(reshape2)
library(ggplot2)
library(edgeR)
library(rtracklayer)
library(RColorBrewer)
library(org.Mm.eg.db)
library(DESeq)
#library(Go.db)
library(ggrepel)
timeStamp <- format(Sys.time(), "%Y_%m_%d")
species <- "mouse"
ensVer <- 84


homedir="../../.."

inPath=paste0(homedir,"/projects/simon/project_results/RNAseq_1.rsem/")

######## directory structure #######
projectDir=paste0(homedir,"/projects/simon")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures")

chrs=seqlengths(Mmusculus)[!grepl("_",names(seqlengths(Mmusculus)))]
gr<-GRanges(seqnames=names(chrs),IRanges(start=1,end=chrs))

########## 1. Load in the files

geneResults<-list()
FPKMs<-list()
TPMs<-list()
for(file in list.files(inPath,full.names=T,pattern="genes.results")){

  sampleName<-basename(file)
  sampleName<-gsub(".genes.*","",sampleName)

  cat(sampleName)
  cat("\n")
  data<-read.table(file,header=T)  

  geneResults[[sampleName]]<-as.integer(data$expected_count)
  FPKMs[[sampleName]]<-as.integer(data$FPKM)
  TPMs[[sampleName]]<-as.integer(data$TPM)
}

########## 2. Outputs

df<-as.data.frame(geneResults)
colnames(df)<-names(geneResults)
row.names(df)<-as.character(data$gene_id)
#eliminate the outlier group with fewer replicates
df<-df[,!grepl("NT3",colnames(df))]
write.table(df,file="../project_results/all_counts.tsv",sep="\t",quote=F)

FPKMsdf<-as.data.frame(FPKMs)
colnames(FPKMsdf)<-names(FPKMsdf)
row.names(FPKMsdf)<-as.character(data$gene_id)
FPKMsdf<-FPKMsdf[,!grepl("NT3",colnames(FPKMsdf))]
write.table(FPKMsdf,file="../project_results/FPKMS.tsv",sep="\t",quote=F)

TPMsdf<-as.data.frame(TPMs)
colnames(TPMsdf)<-names(TPMsdf)
row.names(TPMsdf)<-as.character(data$gene_id)
TPMsdf<-TPMsdf[,!grepl("NT3",colnames(TPMsdf))]
write.table(TPMsdf,file="../project_results/TPMs.tsv",sep="\t",quote=F)

########## 3. PCA

cols<-rep(brewer.pal(6,"Set1"),each=3)
pdf("../project_results/figures/samples_pca.pdf",width=12,height=8)
pca<-princomp(df)
plot(pca$loading,pch=19, cex=2,col=cols)
text(pca$loading, colnames(df),pos = 1)
dev.off()

pdf("../project_results/figures/samples_pca_components.pdf",width=12,height=8)
plot(pca)
dev.off()

########## 4. Change IDs to symbols

#Annotate with symbols, aggregate
egENSEMBL <- toTable(org.Mm.egENSEMBL)
row.names(df)=gsub("\\..*","",row.names(df))
m <- match(row.names(df), egENSEMBL$ensembl_id)
df$EntrezGene<-egENSEMBL$gene_id[m]

egSYMBOL <- toTable(org.Mm.egSYMBOL)
m <- match(df$EntrezGene, egSYMBOL$gene_id)
df$symbol<-egSYMBOL$symbol[m]

#eliminate duplicated symbols
o <- order(rowSums(df[,1:12]), decreasing=TRUE)
df=df[o,]
d<-duplicated(df$symbol)
df<-df[!d,]

#eliminate lowly expressed
include<-apply(df[,1:20],1,function(x){sum(x>=10)>=3})
df<-df[include,]

df=df[!is.na(df$symbol),]
row.names(df)<-df$symbol

########## 5. Build models and perform DE

#remove IE1 and IE2 outliers
dfS<-df[,!colnames(df) %in% "IE1_3"]
dfS<-dfS[,!colnames(dfS) %in% "IE2_1"]
#rename 4T1 as bulk
colnames(dfS)<-gsub("4T1","Bulk",colnames(dfS))
row.names(dfS)<-dfS$symbol
dfS<-dfS[,!colnames(dfS) %in% c("EntrezGene","symbol")]

#prepare model to control for the batch effect
cellType<-gsub("_.*","",colnames(dfS))[1:18]
batch<-c(1,2,2,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1)

design<-model.matrix(~0+cellType+batch)

my.contrasts <- makeContrasts(
 IE1vsIE2 = cellTypeIE1-cellTypeIE2,
 IE1vsBulk = cellTypeIE1-(cellTypeNT1/2+cellTypeNT2/2-cellTypeBulk),
 IE2vsBulk = cellTypeIE2-(cellTypeNT1/2+cellTypeNT2/2-cellTypeBulk),
 NT1vsBulk = cellTypeNT1-(cellTypeNT1/2+cellTypeNT2/2-cellTypeBulk),
 NT2vsBulk = cellTypeNT2-(cellTypeNT1/2+cellTypeNT2/2-cellTypeBulk),
 levels=design)

########### IE1 vs IE2
expr <- DGEList(counts=dfS)
expr <- calcNormFactors(expr,method="RLE")
expr <- estimateDisp(expr,design)
fit <- glmQLFit(expr,design)

for(contrast in colnames(my.contrasts)){
  qlf <- glmQLFTest(fit, contrast=my.contrasts[,contrast])
  logFCs<-as.data.frame(topTags(qlf,nrow(dfS)))
  write.csv(logFCs,paste0(contrast,"_BatchCorr.tsv"))
}


