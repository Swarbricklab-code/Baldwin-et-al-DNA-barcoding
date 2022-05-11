# Annotation of SvABA variant calls 
# 220511
# Louise Baldwin

#input files: svaba output file called svaba.indel.vcf (or similar)
#output files: sample specific csv files listing indels and their annotations


# R set up -------------------------------------------------

library(GenomicRanges)
library(annotatr)
library(rtracklayer)
library(ShortRead)


projectDir="/share/ScratchGeneral/loubal/projects/WGS"
inPath=paste0(projectDir,"/results_svaba/BALB_cJ_v1")

# Load data and analysis ----------------------------------

samples<-list.files(inPath, pattern="svaba.indel.vcf")
#> samples"IE1.svaba.indel.vcf" "IE2.svaba.indel.vcf" "NT1.svaba.indel.vcf" "NT2.svaba.indel.vcf"

results<-GRangesList()
resultsVCF<-list()
for(sampleName in samples){
	inDir<-paste0(inPath,"/",sampleName)
	inFile<-paste0(inDir)
	data=read.table(inFile)
	gr<-GRanges(seqnames=data$V1,IRanges(start=data$V2,width=nchar(as.character(data$V4))))
	results[[sampleName]]=gr
	resultsVCF[[sampleName]]=data
}

#can find the individual sample results by going results$IE1.svaba.indel.vcf
#do the filtering
grReduce<-reduce(unlist(results))

#compare the subclones
#gr1=IE1, gr2=IE2, gr3=NT1, gr4=NT2
gr1<-results[[1]]
gr2<-results[[2]]
gr3<-results[[3]]
gr4<-results[[4]]

gr1<-results$IE1.svaba.indel.vcf
gr2<-results$IE2.svaba.indel.vcf
gr3<-results$NT1.svaba.indel.vcf
gr4<-results$NT2.svaba.indel.vcf


#datashort is an abbreviated version of the granges object, can just give an overview of the indels
dataShort_1<-resultsVCF[[1]]
dataShort_2<-resultsVCF[[2]]
dataShort_3<-resultsVCF[[3]]
dataShort_4<-resultsVCF[[4]]

#how long is each insert/deletion?
#remember: $V4 = reference genome, $V5 is what was sequenced in the sample
dataShort_1$V11<-nchar(as.character(dataShort_1$V4))
dataShort_1$V12<-nchar(as.character(dataShort_1$V5))

#length of 
dataShort_1$indel<-apply(dataShort_1[,c("V11","V12")],1,function(x){ifelse(x[1]>x[2],"deletion","insertion")})
dataShort_1$width<-apply(dataShort_1[,c("V11","V12")],1,function(x){abs(x[1]-x[2])})
#sample2
dataShort_2$V11<-nchar(as.character(dataShort_2$V4))
dataShort_2$V12<-nchar(as.character(dataShort_2$V5))

#length of 
dataShort_2$indel<-apply(dataShort_2[,c("V11","V12")],1,function(x){ifelse(x[1]>x[2],"deletion","insertion")})
dataShort_2$width<-apply(dataShort_2[,c("V11","V12")],1,function(x){abs(x[1]-x[2])})

#sample3
dataShort_3$V11<-nchar(as.character(dataShort_3$V4))
dataShort_3$V12<-nchar(as.character(dataShort_3$V5))

#length of 
dataShort_3$indel<-apply(dataShort_3[,c("V11","V12")],1,function(x){ifelse(x[1]>x[2],"deletion","insertion")})
dataShort_3$width<-apply(dataShort_3[,c("V11","V12")],1,function(x){abs(x[1]-x[2])})

#sample4
dataShort_4$V11<-nchar(as.character(dataShort_4$V4))
dataShort_4$V12<-nchar(as.character(dataShort_4$V5))

#length of 
dataShort_4$indel<-apply(dataShort_4[,c("V11","V12")],1,function(x){ifelse(x[1]>x[2],"deletion","insertion")})
dataShort_4$width<-apply(dataShort_4[,c("V11","V12")],1,function(x){abs(x[1]-x[2])})


#write these to a csv
write.table(x=dataShort_1, file=paste0(inPath,"/IE1_indels.csv"), sep=",")
write.table(x=dataShort_2, file=paste0(inPath,"/IE2_indels.csv"), sep=",")
write.table(x=dataShort_3, file=paste0(inPath,"/NT1_indels.csv"), sep=",")
write.table(x=dataShort_4, file=paste0(inPath,"/NT2_indels.csv"), sep=",")

