library(cn.mops)
setwd("/home/Users/cnmops/BamFiles_4T1/") #this is the directory of the bam files from the WGS data
BAMFILES<-list.files(pattern='.bam$') #read the list of bam files
#prepare read counts for cnmops of the copy number analysis
Samples_bam<-getReadCountsFromBAM(BAMFILES,refSeqName=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX"),WL=10000,mode='paired')
#save(Samples_bam,file='Samples_bam_Aug1920_chr1_X.Rdata')
#run cnmops analysis on the selected samples
Samples_cnv<-cn.mops(Samples_bam)
#save(Samples_cnv,file='Samples_cnv_Aug1920_chr1_X.Rdata')
#plot the copy number segmentation from the cnmops results
segplot(Samples_cnv,sampleIdx = 1)