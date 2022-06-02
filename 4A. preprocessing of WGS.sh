#!/bin/bash


# Script for the preprocessing of raw files from WGS seq

# Input: .fastq
# Output: sorted .bam

# Usage in paper: First step in all WGS analysis

# -------------------------------------------------------

module load gi/bwa/0.7.8
module load gi/gcc/4.8.2
module load gi/samtools/1.2
module load gi/novosort/precompiled/1.03.08

numcores=50
tag="-P TumourProgression"

homedir="/share/ScratchGeneral/loubal"
projectDir="$homedir/projects/WGS"
resultsDir="$projectDir/results"
scriptsPath="/share/ScratchGeneral/loubal/projects/WGS/scripts"
logDir=$scriptsPath"/logs"
mkdir -p $logDir
mkdir -p $scriptsPath
mkdir -p $resultsDir

genomeFile="/share/ScratchGeneral/loubal/projects/WGS/genomes/BALB_cJ_v1.fa"
inExt="fastq.gz"

i=0

	inPath="/share/ScratchGeneral/loubal/projects/WGS/raw_files"
	outPath="$projectDir/results/BALB_cJ_v1.fa"
        #log and command files for bsub
        logPath=$logDir
        commandPath="commands"
        #make the directory structure   
        mkdir -p $outPath
        mkdir -p $logPath

        subs=0

        #get the name of the script for the logs
        scriptName=`basename $0`
        
        echo $inPath
	files=`ls $inPath/*.fastq.gz`
	for file in ${files[@]};do
                        echo The file used is: $file
                        filesTotal[i]=$file;
                        let i++;
        done;

j=0
echo ${#filesTotal[@]}
while [ $j -lt ${#filesTotal[@]} ]; do

        #dir=`echo ${filesTotal[$j]}`
        #files=`ls $inPath/$dir/*.$inExt`
	
        inFile1=${filesTotal[$j]}
        inFile2=${filesTotal[$(($j+1))]}
        uniqueID=`basename $inFile1 | sed s/_R1.fastq.gz//`
        name=$uniqueID
        outDir=$outPath/$name/
	mkdir -p $outDir
        echo $name
	#echo $command_line

	bwaJobName="bwa."$name
	samSortJobName="samSort"$name
	bamJobName="bam."$name
	sortJobName="sort."$name
	
	indexJobName="index."$name
	indexStatsJobName="indexstats."$name
	outSam=$outDir"Aligned.out.sam"
	outSortedSam=$outDir"Aligned.sorted.sam"
	outBam=$outDir"$name.bam"
	outSortedBam=$outDir"Aligned.sortedByCoord.out.bam"

        bwa_line="bwa mem $genomeFile -t $numcores $inFile1 $inFile2 | samtools view -b -u -S - |novosort -m 7G -c $numcores -i -o $outDir$name.sorted.bam -"
        commandline="qsub -N $bwaJobName -hold_jid BALB_cJ_v1.fa -b y -wd $logDir -j y -R y -pe smp $numcores -q short.q -V"
	$commandline $bwa_line

        j=$(($j+2))
done;

# -----------------------------------