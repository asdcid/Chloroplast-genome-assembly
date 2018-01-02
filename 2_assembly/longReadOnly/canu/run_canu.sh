#!/bin/bash

set -e

export PATH='/home/raymond/devel/canu/canu/Linux-amd64/bin/':$PATH


inputFile=$1
outputDir=$2
name='Epau'
genomeSize=160kb
#If the coverage is lower than 40x, the correctedErrorRate is set to 0.154
correctedErrorRate=0.144
threads=40
gnuplotPath=/home/raymond/devel/gnuplot/install/bin/gnuplot
echo $inputFile
echo $outputDir

canu \
    -p $name \
    -d $outputDir \
    correctedErrorRate=$correctedErrorRate \
    genomesize=$genomeSize \
    -nanopore-raw $inputFile\
    gnuplot=$gnuplotPath \
    useGrid=false \
    maxThreads=$threads




