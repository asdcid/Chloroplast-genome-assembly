#!/bin/bash


ulimit -c unlimited

#export PATH='/home/raymond/devel/racon/racon/bin/':$PATH
export PATH='racon/path':$PATH

polish()
{
    inputFile=$1
    outputDir=$2
    ref=$3
    threads=$4

    mkdir $outputDir

    echo $ref
    echo $inputFile
    refDir=$outputDir/ref
    bwaDir=$outputDir/bwa
    resultDir=$outputDir/result

    mkdir $refDir
    mkdir $bwaDir
    mkdir $resultDir

    coverage=$(basename ${inputFile%%.fastq})

    #index
    ln $ref $refDir
    ref=$refDir/*fa
    build_index $ref

    #run bwa
    outputBwa=$bwaDir/$coverage.sam
    run_bwa $inputFile $bwaDir $ref $threads $coverage $outputBwa

    #run nanopolish
    outputRacon=$resultDir/$coverage.polished.racon.fa
    run_racon $inputFile $outputRacon $outputBwa $ref 
    echo $ref
    echo $outputBwa
}

build_index()
{
    indexFile=$1
    bwa index $indexFile
}

run_bwa()
{
    inputFile=$1
    outputDir=$2
    ref=$3
    threads=$4
    coverage=$5
    outputFile=$6

    bwa mem \
        -x ont2d \
        -t $threads \
        $ref \
        $inputFile > $outputFile

}

run_racon()
{
echo run_racon
inputFile=$1
outputFile=$2
bwa_inputFile=$3
ref=$4

racon \
    --sam \
    $inputFile \
    $bwa_inputFile \
    $ref \
    $outputFile


}


inputFile=$1
outputDir=$2
ref=$3
threads=$4


cd $outputDir


n=1

outputDir=round-$n
echo polish
polish $inputFile $outputDir $ref $threads
echo "finish 1"

ref=round-$n/result/*fa
n=$(( $n + 1 ))
echo $n
outputDir=round-$n
polish $inputFile $outputDir $ref $threads
echo "finish 2"

LAST=round-$(( $n - 1 ))/result/*.fa 
echo $LAST
CURRENT=round-$n/result/*.fa
echo $CURRENT
DIFF=$(diff -q $LAST $CURRENT)
while [ "$DIFF" ]
do
    echo "begin round"
    ref=round-$n/result/*fa
    n=$(( $n + 1 ))
    echo round $n
    outputDir=round-$n
    polish $inputFile $outputDir $ref $threads
    if [ $n -gt 10 ]
    then
        cp round-$n/result/* ../final/
        exit
    fi
    LAST=round-$(( $n - 1 ))/result/*.fa 
    CURRENT=round-$n/result/*fa
    DIFF=$(diff -q $LAST $CURRENT)
    echo 'diff' $DIFF
done

cp round-$n/result/* ../final/
#rm -r round*
