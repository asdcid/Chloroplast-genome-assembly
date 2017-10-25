#!/bin/bash


#$ -N trim
#$ -cwd

export PATH='nanoFilt/path':$PATH
#export PATH='/home/raymond/devel/python/install/bin/':$PATH


#inputDir='/home/raymond/work/Eucalyptus_pauciflora/genome/bin/raw_read_check/8_branches/long_read/2_adapter_trim/porechop/result/'
#outputDir='result/'
#quality=9
#minLen=5000

inputDir=$1
outputDir=$2
quality=$3
minLen=$4


for inputFile in $inputDir/*fastq.gz
do
    id=$(basename $inputFile)
    outputFile=$outputDir/$id

    echo $id
    echo $outputFile

    zcat $inputFile | NanoFilt -q $quality -l $minLen | gzip > $outputFile
done
