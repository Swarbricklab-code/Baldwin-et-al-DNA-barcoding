#$ -S /bin/bash

module load centos6.10/gi/boost/1.53.0
module load centos6.10/gi/bowtie/1.0.0 
module load centos6.10/borgue/rsem/1.2.26
module load centos6.10/gi/gcc/4.8.2

numcores=15
tag="-P TumourProgression" 

homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/simon"
resultsDir="$projectDir/project_results/"
projectname="RNAseq_1"

#scripts directory
scriptsPath="$homedir/projects/simon/scripts/QC"
logDir=$scriptsPath"/logs"

mkdir -p $logDir

genome="mm10"
genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/$genome"
annotationFile="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/$genome/gencode.vM23.annotation.gtf"
indexDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/$genome/mm10"
rsemIndex="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/$genome/$genome"


#input/output
inExt="sorted.bam"
inType="star"
inPath="$homedir/projects/simon/project_results/$projectname.$inType/"

outTool="rsem"
outPath="$projectDir/project_results/$projectname.$outTool"
mkdir -p $outPath



#Get the subpath
files=`ls $inPath`
i=0
for file in ${files[@]};do
	echo The file used is: $file
	filesTotal[i]=$file;
	let i++;
done 

files=`ls $inPath/**/*.$inExt`

qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp $numcores $tag -V"

j=0
echo ${#filesTotal[@]}
while [ $j -lt ${#filesTotal[@]} ]; do

    	dir=`echo ${filesTotal[$j]}`
    	file=($(ls $inPath/$dir/*.$inExt))
	
	uniqueID=`basename $dir`
	name=$uniqueID
	echo $name

	sortJobName="sort."$name
	outType="rsem"
	outPath="$homedir/projects/simon/project_results/$projectname.$outType/$name"
	mkdir -p $outPath
	qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp $numcores $tag -V"
	rsem_line="rsem-calculate-expression -p $numcores --bam --no-bam-output --forward-prob 0 --paired-end --bam $file $rsemIndex $outPath"
	echo $rsem_line                

	$qsubLine -N RSEM_count_$name -hold_jid $sortJobName $rsem_line

	j=$(($j+1))

done;

