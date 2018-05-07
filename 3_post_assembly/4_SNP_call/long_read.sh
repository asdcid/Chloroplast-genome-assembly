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
 
samtools view -b -F 2308 -@ $threads $outputFile > $outputFile.bam
samtools sort -@ $threads $outputFile.bam > $outputFile.sort.bam
samtools index $outputFile.sort.bam
 
        
 nanopolish variants \
        -r $inputFile \
        -b $outputFile.sort.bam \
        -g $ref \
        -t $threads \
        --min-candidate-frequency 0.1 \
        -o $outputFile.vcf
