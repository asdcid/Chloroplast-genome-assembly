#!/bin/bash

inputFile=$1
outputFile=$2
ref=$3
threads=$4


 ngmlr \
        -t $threads \
        -r $ref \
        -q $inputFile \
        -o $outputFile.sam \
        -x ont \
        --skip-write
 
 
        
 nanopolish variants \
        -r $inputFile \
        -b $outputFile.sort.bam \
        -g $ref \
        -t $threads \
        --min-candidate-frequency 0.1 \
        -o $outputFile.vcf
