#!/bin/bash


set -e 

source /home/raymond/devel/C/gcc/setenv.sh


inputDir='/home/raymond/work/Eucalyptus_pauciflora/genome/bin/raw_read_check/8_branches/short_read/3_correction/data/'
outputDir='result'
threads=40
ref='ref/cp_double_up_combine'




for in1 in $inputDir/*R1*fastq.gz
do
    in2=${in1/"R1"/"R2"}
    echo "mapping files: "
    echo $in1
    echo $in2
    
    # output file setup
    f1=$(basename $in1)
    id=${f1%%_R1*}
    outputFile=$outputDir/$id".sam"
    start=$(date +%s)
    bowtie2 \
        -q \
        -x $ref \
        -1 $in1 \
        -2 $in2 \
        -p $threads \
        -S $outputFile > $(basename $in1).log

    end=$(date +%s)
    DIFF=$(( $end - $start ))
    echo $outputDir/$id "Excution time: " $DIFF " seconds." >> runCPTime.log

    output=${outputFile/"sam"/"bam"}
    outbamsort=${output/"bam"/"sort.bam"}
    echo $outbamsort
    samtools sort -@ $threads  -o $outbamsort $outputFile
    samtools index $outbamsort

done



