#!/bin/bash

#load all modules

module load gi/samtools/1.0
module load nenbar/star/2.4.0d
#number of cores
ncore=8


#genome directory
genome="mm10_2.7.6a"
genomesDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/$genome/"
#annotationDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/sequin"
starLine="/share/ScratchGeneral/nenbar/local/src/STAR-2.7.6a/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir $genomesDir  --genomeFastaFiles $genomesDir/mm10.fa --runThreadN $ncore"

qsub -N mm10 -b y -cwd -j y -R y -pe smp $ncore -V $starLine
