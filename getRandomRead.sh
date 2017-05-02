#!/bin/bash

#############################################################
#dir contains R1 and R2, in fq.gz format
inputDir='Sample_14710_PE_450bp_GBR_UNSW_H3NJJBCXX'
outputDir='separate/'
script='getRandomRead.py'
#0.1 means keep 90% reads in large
fraction=0.1
#############################################################

for R1 in $inputDir/*R1*fastq.gz
do
    echo $R1
    R2=${R1//"R1"/"R2"}
    outputR1=$outputDir/$(basename $R1 ".fastq.gz")
    outputR2=$outputDir/$(basename $R2 ".fastq.gz")
    ./$script $R1 $R2 $outputR1.large.fastq $outputR1.small.fastq $outputR2.large.fastq $outputR2.small.fastq $fraction
    gzip $outputR1.large.fastq
    gzip $outputR1.small.fastq
    gzip $outputR2.large.fastq
    gzip $outputR2.small.fastq

done

