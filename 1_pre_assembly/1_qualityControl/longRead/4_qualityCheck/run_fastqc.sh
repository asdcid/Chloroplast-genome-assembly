#!/bin/bash


export PATH='fastqcPath':$PATH

inputFile='/home/raymond/work/Eucalyptus_pauciflora/genome/data/nanopore/RB7_C4.fastq.gz'
outputDir='result/'
threads=20



fastqc $inputFile -o $outputDir -t $threads

