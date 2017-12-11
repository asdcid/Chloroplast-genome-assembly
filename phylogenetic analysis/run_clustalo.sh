#!/bin/bash



export PATH='path/of/clustalo':$PATH


inputFile=$1
outputFile=$2
threads=40

clustalo \
    -i $inputFile \
    -o $outputFile \
    --outfmt='clustal' \
    --thread $threads
