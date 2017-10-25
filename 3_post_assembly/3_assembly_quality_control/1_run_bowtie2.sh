#!/bin/bash


set -e 

#in1 and in2 are the 100x validation data R1, R2
in1=$1
in2=$2
outputFile=$3
#ref is the assembly fasta
ref=$4
threads=$5


bowtie2-build $ref $ref

bowtie2 \
    -q \
    -x $ref \
    -1 $in1 \
    -2 $in2 \
    -p $threads \
    -S $outputFile 

output=${outputFile/"sam"/"bam"}
outbamsort=${output/"bam"/"sort.bam"}
samtools sort -@ $threads  -o $outbamsort $outputFile
samtools index $outbamsort




