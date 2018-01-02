#!/bin/bash

set -e

#the ref using for bowtie2 mapping
ref='/home/raymond/work/Eucalyptus_pauciflora/genome/bin/double_up/cp_double_up_combine.fasta'
#the dir of bowtie2 output results, bam file, sorted
inputDir='result'
#output dir for storing cp fastq read.
outputDir='map_to_chl'



for inputFile in $inputDir/*.sort.bam
do
    echo "Processing " $inputFile
    temp=temp_$(basename $inputFile)
    outputFileR1=$outputDir/$(basename ${inputFile/"sort.bam"/"R1.fastq"})
    outputFileR2=$outputDir/$(basename ${inputFile/"sort.bam"/"R2.fastq"})
    echo $outputFileR1
    echo $outputFileR2


    for chlSpecies in $(awk '{if (NR % 2 == 1) {print substr($1,2)}}' $ref)
    do
        samtools view -b $inputFile $chlSpecies > $temp$chlSpecies
        samtools sort -n -o $temp$chlSpecies.sort.bam $temp$chlSpecies
        wait
        bamToFastq -i $temp$chlSpecies.sort.bam -fq R1chl_$chlSpecies.temp -fq2 R2chl_$chlSpecies.temp
    done
    wait
    cat R1chl*.temp > $outputFileR1
    cat R2chl*.temp > $outputFileR2
    pigz $outputFileR1
    pigz $outputFileR2
    wait
    rm *temp*
done



