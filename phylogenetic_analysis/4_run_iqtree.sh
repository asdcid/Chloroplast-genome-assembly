#!/bin/bash


export PATH='path/of/iqtree':$PATH


inputFile='fixed_whole_genome_alignment'
nexus='nexus.file'
outputFile='all.aln'

iqtree\
    -s $inputFile \
    -spp $nexus \
    -pre $outputFile \
    -bb 1000

