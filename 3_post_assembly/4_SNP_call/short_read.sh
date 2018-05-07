#!/bin/bash

in1=$1
in2=$2
ref=$3
threads=$4
outputFile=$5

bowtie2 \
        -q \
        -x $ref \
        -1 $in1 \
        -2 $in2 \
        -p $threads \
        -S $outputFile
samtools view -b -F 264 -@ $threads $outputFile > $outputFile.bam
samtools sort -@ $threads $outputFile.bam > $outputFile.sort.bam
samtools index $outputFile.sort.bam



samtools mpileup -B -f $ref $outputFile.sort.bam | \
        java -jar /home/raymond/devel/varscan/varscan/VarScan.v2.4.2.jar mpileup2cns \
                --min-coverage 10 \
                --min-var-freq 0.1 \
                --variants \
                --output-vcf > $outputFile.vcf
