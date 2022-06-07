#!/bin/bash

#module load gi/star/2.3.0e
module load gi/samtools/1.0
module load gi/novosort/precompiled/1.03.08
module load nenbar/star/2.4.0d
module load gi/gcc/4.8.2

numcores=24
tag="-P TumourProgression"


#directory hierarchy
#raw files directory
homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/simon"
resultsDir="$projectDir/project_results/"

genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/mm10"
#genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/mm10_sequin"
#extension of the files to be used
inExt="fastq.gz"

#scripts directory
scriptsPath="/share/ClusterShare/biodata/contrib/nenbar/projects/simon/scripts/rnaseq"
logDir=$scriptsPath"/logs"

#name to append to projectname and create a folder
inType="trimgalore"

projectnames=( "RNAseq_1" )

for projectname in ${projectnames[@]}; do


        #out directory
        
        inPath="$homedir/projects/simon/project_results/$projectname.$inType/"
	outPath="/share/ScratchGeneral/nenbar/projects/simon/project_results/$projectname.star"
        #log and command files for bsub
        logPath="logs"
        commandPath="commands"
        #make the directory structure   
        mkdir -p $outPath
        mkdir -p $logPath
        mkdir -p $commandPath

        rm -f $commandFile

        
        #echo Reading files from $sampleFile
        if [[ ! (-f $sampleFile) ]]; then
                echo The file with the list of samples $sampleFile does not exist
        fi

        subs=0

        #get the name of the script for the logs
        scriptName=`basename $0`
        i=0
        echo $inPath
	files=`ls $inPath`
        for file in ${files[@]};do
                        echo The file used is: $file
                        filesTotal[i]=$file;
                        let i++;
        done 
done;

j=0
echo ${#filesTotal[@]}
while [ $j -lt ${#filesTotal[@]} ]; do

        dir=`echo ${filesTotal[$j]}`
        files=`ls $inPath/$dir/*.$inExt`
	
	inFile1=${files[0]}
	inFile2=${files[1]}
        uniqueID=`basename $dir`
        name=$uniqueID
        outDir=$outPath/$uniqueID/
	mkdir -p $outDir
        echo $name
	#echo $command_line


	starJobName="star."$name
	samSortJobName="samSort"$name
	bamJobName="bam."$name
	sortJobName="sort."$name
	
	indexJobName="index."$name
	indexStatsJobName="indexstats."$name
	outSam=$outDir"Aligned.out.sam"
	outSortedSam=$outDir"Aligned.sorted.sam"
	outBam=$outDir"$name.bam"
	outSortedBam=$outDir"$name.sorted.bam"

	#star_line="/home/nenbar/local/lib/STAR-STAR_2.4.0i/source/STAR --genomeDir $genomeDir --runMode alignReads --readFilesIn $inFile1 $inFile2 --outFileNamePrefix $outDir --runThreadN 4 --outSAMattributes Standard --outSAMstrandField intronMotif --sjdbOverhang 99" 
	
	star_line="STAR \
                --runMode alignReads \
		--genomeDir $genomeDir \
		--readFilesIn $inFile1 $inFile2 \
                --outFileNamePrefix $outDir \
		--runThreadN $numcores \
		--sjdbOverhang 100 \
		--readFilesCommand zcat \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS NM MD\
        --outFilterMultimapNmax 999 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1500000 \
        --alignMatesGapMax 1500000 \
        --alignSJoverhangMin 6 \
        --alignSJDBoverhangMin 1 \
        --quantMode TranscriptomeSAM \
        --outFilterMatchNmin 40 \
        --outSAMtype BAM Unsorted "

	bam_line="samtools view -m 16G -h -S $outSam -b -o $outBam"
	samtools_line="samtools sort -m 16G $outBam $outDir$name.sorted"
        qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp $numcores $tag -V"
	
        $qsubLine -N $starJobName -hold_jid trimgalore  $star_line 
	$qsubLine -N $bamJobName -hold_jid $starJobName $bam_line
	$qsubLine -N $sortJobName -hold_jid $bamJobName $samtools_line 
	$qsubLine -N $indexJobName -hold_jid $sortJobName "samtools index $outSortedBam;"
	$qsubLine -N $indexStatsJobName -hold_jid $indexJobName "samtools idxstats $outSortedBam | head -n1 | cut -f3 >$outDir$name.tmp;"

        j=$(($j+1))


done;

