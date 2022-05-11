#!/bin/bash

module load gi/gcc/4.8.2
module load gi/samtools/1.0
module load gi/novosort/precompiled/1.03.08

numcores=30
tag="-P TumourProgression"
#tag=""

homedir="/share/ScratchGeneral/loubal"
projectDir="$homedir/projects/WGS"
resultsDir="$projectDir/results"
scriptsPath="/share/ScratchGeneral/loubal/projects/WGS/scripts"
logDir="$scriptsPath/logs"
mkdir -p $resultsDir
mkdir -p $logDir

genomeFile="/share/ScratchGeneral/loubal/projects/WGS/genomes/mm10.fa" 
projectnames=( "4T1_WGS" )

i=0

samples=( "H7YWTCCX2_5_200424_FD07688460_Mus-musculus__R_200423_SIMJUN_DNA_M001"
	"H7YWTCCX2_6_200424_FD07688475_Mus-musculus__R_200423_SIMJUN_DNA_M001" 
	"H7YWTCCX2_7_200424_FD07688477_Mus-musculus__R_200423_SIMJUN_DNA_M001" 
	"H7YWTCCX2_8_200424_FD07688484_Mus-musculus__R_200423_SIMJUN_DNA_M001" 
	)
for sample in ${samples[@]}; do

	inPath="$projectDir/results/mm10/$sample"
	outPath="$projectDir/results_svaba/mm10/$sample"
        #log and command files for bsub
	logPath=$logDir
    #make the directory structure   
	mkdir -p $outPath
	mkdir -p $logPath

	files=`ls $inPath/*.bam`
	for file in ${files[@]};do
		uniqueID=`basename $file | sed s/.bam//`
		mkdir -p $outPath/$uniqueID
			svaba_line="/share/ScratchGeneral/loubal/local/svaba/bin/svaba run -t $file -p $numcores -L 6 -I -a 4T1_WGS_1 -G $genomeFile"
		qsub -N $uniqueID -b y -wd $outPath/$uniqueID -j y -R y -pe smp $numcores $tag -q short.q -V $svaba_line
	done;
done;
