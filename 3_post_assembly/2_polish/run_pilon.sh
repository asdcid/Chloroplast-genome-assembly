#!/bin/bash


polish()
{
    inputR1=$1
    inputR2=$2
    outputDir=$3
    ref=$4
    threads=$5

    mkdir $outputDir

    echo $ref
    echo $inputFile
    refDir=$outputDir/ref
    bowtie2Dir=$outputDir/bowtie2
    resultDir=$outputDir/result

    mkdir $refDir
    mkdir $bowtie2Dir
    mkdir $resultDir

    coverage=$(basename ${inputR1%%.R1.fastq.gz})

    #index
    ln $ref $refDir
    ref=$refDir/*fasta
    build_index $ref

    #run bowtie2
    outputBowtie2=$bowtie2Dir/$coverage.sort.bam
    run_bowtie2 $inputR1 $inputR2 $bowtie2Dir $ref $threads $coverage $outputBowtie2

    #run nanopolish
    outputPilon=$resultDir/ #$coverage.polished.pilon.fa
    run_pilon $outputBowtie2 $outputPilon $ref $threads
    echo $ref
    echo $outputBowtie2
}

build_index()
{
    indexFile=$1
    bowtie2-build -q $indexFile $indexFile
}

run_bowtie2()
{
    in1=$1
    in2=$2
    outputDir=$3
    ref=$4
    threads=$5
    coverage=$6
    outputFile=$7

    bowtie2 \
        -q \
        -x $ref \
        -1 $in1 \
        -2 $in2 \
        -p $threads \
        -S $outputDir/$coverage.sam

    samtools sort \
        -o $outputFile \
        $outputDir/$coverage.sam

    samtools index \
        $outputFile
}

run_pilon()
{
echo run_pilon
bowtie2_inputFile=$1
outputDir=$2
ref=$3
threads=$4

#should be change
java -jar /home/raymond/devel/polish/pilon/pilon-1.22.jar \
    --genome $ref \
    --frags $bowtie2_inputFile \
    --outdir $outputDir \
    --changes \
    --vcf \
    --threads $threads \
    --fix all 
}


inputR1=$1
inputR2=$2
outputDir_all=$3
ref_ori=$4
threads=$5
checkDiff=$6

cd $outputDir_all


n=1

outputDir=round-$n
echo polish
polish $inputR1 $inputR2 $outputDir $ref_ori $threads
echo "finish 1"


DIFF=$(python $checkDiff $ref_ori round-1/result/pilon.fasta 2>&1)
if [ "$DIFF" == 'same' ]
then
    echo $outputDir_all 'same'
    cp  $ref ../final/$(basename $outputDir_all).fasta
    exit
fi

ref=round-$n/result/pilon.fasta
n=$(( $n + 1 ))
echo $n
outputDir=round-$n
polish $inputR1 $inputR2 $outputDir $ref $threads
echo "finish 2"

LAST=round-$(( $n - 1 ))/result/pilon.fasta 
echo $LAST
CURRENT=round-$n/result/pilon.fasta
echo $CURRENT
DIFF=$(python $checkDiff $LAST $CURRENT 2>&1)
while [ "$DIFF" == 'different' ]
do
    echo "begin round"
    ref=round-$n/result/pilon.fasta
    n=$(( $n + 1 ))
    echo round $n
    outputDir=round-$n
    polish $inputR1 $inputR2 $outputDir $ref $threads

    LAST=round-$(( $n - 1 ))/result/pilon.fasta
    CURRENT=round-$n/result/pilon.fasta
    DIFF=$(diff -q $LAST $CURRENT)
    echo $outputDir_all 'diff'
done

cp round-$n/result/pilon.fasta ../final/$(basename $outputDir_all).fasta
#rm -r round*
