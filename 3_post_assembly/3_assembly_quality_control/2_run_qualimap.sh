#!/bin/bash

set -e

#export PATH="/home/raymond/devel/qualimap/qualimap_v2.2.1/":$PATH
export PATH='qualimap/path':$PATH

#inputDir='../bowtie2/result/'
#outputDir='result/'
#threads=10

#the dir includes the 100x validation data mapping result, sort.bam
inputDir=$1
outputDir=$2
threads=$3

for inputFile in $inputDir/*.sort.bam
do
    echo "Processing " $inputFile
    qualimap \
            bamqc \
            -bam $inputFile \
            -outdir $outputDir/$(basename ${inputFile/".sort.bam"/""}) \
            -nt $threads \
            -c
done

