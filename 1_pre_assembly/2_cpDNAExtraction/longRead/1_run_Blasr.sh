#!/bin/bash


set -e 

source /home/raymond/devel/blasr/install/setup-env.sh


inputDir='query/'
outputDir='result'
threads=40
ref='ref/cp_double_up_combine.fasta'
minMatch=15
minAlnLength=5000


for in1 in $inputDir/*fastq
do
    echo "Processing " $in1
    
    # output file setup
    f1=$(basename $in1)
    id=${f1%%.fastq}
    outputFile=$outputDir/$id".out"
    start=$(date +%s)

    blasr \
          $in1 \
          $ref \
          --nproc $threads \
          --bestn 1 \
          -m 1 \
          --minMatch $minMatch \
          --minAlnLength $minAlnLength \
          --out $outputFile
         

    end=$(date +%s)
    DIFF=$(( $end - $start ))
    echo $outputDir/$id "Excution time: " $DIFF " seconds." >> runTime.log



done  


