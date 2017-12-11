#!/bin/bash

#assume the coverage of short read (R1, R2) is 1000x each
inputR1='R1.fastq.gz'
inputR2='R2.fastq.gz'

outputDir='output'
read_coverage=(5 8 10 20 40 60 80 100 200 300 400 500)
split_per=(995 992 990 980 960 940 920 900 800 700 600 500)

length=${#read_coverage[@]}

for ((i=0;i<$length;i++))
do
    outputR1=$outputDir/${read_coverage[$i]}.R1.fastq
    outputR2=$outputDir/${read_coverage[$i]}.R2.fastq
    splitReads_paired_not_test.sh \
        $inputR1 \
        $inputR2 \
        $outputR1 \
        $outputR2 \
        ${split_per[$i]}


done
