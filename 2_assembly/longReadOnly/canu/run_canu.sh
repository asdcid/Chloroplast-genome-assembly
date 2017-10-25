#!/bin/bash

set -e

#export PATH='/home/raymond/devel/canu/canu/Linux-amd64/bin/':$PATH

export PATH='canu/path':$PATH

inputFile=$1
outputDir=$2
name=$3
genomeSize=$4
corOutCoverage=$5
correctedErrorRate=$6
threads=$7
#/home/raymond/devel/gnuplot/install/bin/gnuplot
gnuplotPath=$8
echo $inputFile
echo $outputDir

canu \
    -p $name \
    -d $outputDir \
    corOutCoverage=$corOutCoverage \
    correctedErrorRate=$correctedErrorRate \
    genomesize=$genomeSize \
    -nanopore-raw $inputFile\
    gnuplot=$gnuplotPath \
    useGrid=false \
    maxThreads=$threads




