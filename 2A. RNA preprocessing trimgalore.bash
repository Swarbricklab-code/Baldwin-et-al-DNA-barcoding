#!/bin/bash

module load gi/trim_galore/0.4.0
module load gi/fastqc/0.11.5

############## directory hierarchy ##############
#raw files directory
homedir="/share/ScratchGeneral/nenbar"
inPath="$homedir/projects/simon/raw_files/RNAseq_1"
#inPath="/share/ClusterScratch/nenbar/simon/raw_files/arrayexpress/"

#extension of the files to be used
inExt="fastq.gz"

#scripts directory
scriptsPath="$homedir/projects/simon/scripts/QC"

#name to append to projectname and create a folder
inType="trimgalore"
projectname="RNAseq_1"

#out directory
outPath="/share/ScratchGeneral/nenbar/projects/simon/project_results/"$projectname.$inType
#outPath="/share/ClusterScratch/nenbar/simon/project_results/"$projectname.$inType

#log and command files for bsub
logDir=$scriptsPath/"logs"
commandPath="commands"

#make the directory structure   
mkdir -p $outPath
mkdir -p $logDir
mkdir -p $commandPath
rm -f $commandFile

############## fetch file names ##############

i=0   
files=( $(ls $inPath/$projectname/*.fastq.gz) )
#files=( $(ls $inPath/*.fastq.gz) )
for file in ${files[@]};do
        echo The file used is: $file
        filesTotal[i]=$file;
        let i++;
done;

############## perform analysis in pairs ##############


j=0
echo -e "The total number of files is:"
echo ${#filesTotal[@]}
echo -e

while [ $j -lt ${#filesTotal[@]} ]; do

        inFile1=${files[$j]}
        inFile2=${files[$(($j+1))]}

        uniqueID=`basename $inFile1 | sed s/_1.fastq.gz//`
        echo $uniqueID
        name=$uniqueID
        outDir=$outPath/$uniqueID/
        mkdir -p $outDir
        echo $name
        #echo $command_line

        command_line="trim_galore $inFile1 $inFile2 --gzip --fastqc --paired --length 16 -o $outDir"
        #echo $command_line
        qsub -b y -wd $logDir -j y -N trimgalore -R y -pe smp 1 -V $command_line
        j=$(($j+2))

done;

