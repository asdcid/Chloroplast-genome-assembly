#!/bin/bash


ulimit -c unlimited

set -e

#export PATH='/home/raymond/devel/python3/install/bin/':$PATH
#export PATH='/home/raymond/devel/porechop/Porechop/':$PATH

export PATH='porchop/path':$PATH

#inputDir='/home/raymond/work/Eucalyptus_pauciflora/genome/data/nanopore/'
#outputDir='result'
#threads=20

inputDir=$1
outputDir=$2
threads=$3


for inputFile in $inputDir/*fastq.gz
do
    id=$(basename $inputFile)
    outputFile=$outputDir/$id

    echo $id
    echo $outputFile

    porechop-runner.py \
        -i $inputFile \
        -o $outputFile \
        -t $threads 

done
