#!/bin/bash

set -e

#export PATH="/home/raymond/devel/mummer3/MUMmer3.23/":$PATH
#export PATH="/home/raymond/devel/gnuplot/install/bin/":$PATH

export PATH="mummer/path":$PATH
export PATH="gnuplot/path":$PATH

#assembly result
inputFile=$1
#the genome you want to compare to
ref=$2
#output delta result
outputFile=$3
#output the dot-plot
outputPlotFile=$4


nucmer \
    --mum \
    -prefix=$outputFile  \
    $ref \
    $inputFile \


mummerplot \
    --png \
    -prefix=$outputPlotFile \
    $outputFile.delta

#get the coord for mummer_direction.py (optional)
show-coords -r $outputFile.delta > $outputFile.coord
