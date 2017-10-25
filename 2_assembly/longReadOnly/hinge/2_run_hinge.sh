#!/bin/bash


set -e
#PPWD='/home/raymond/devel/hinge/HINGE/'
PPWD='hinge/path'
export PATH="$PATH:$PPWD/thirdparty:$PPWD/thirdparty/DALIGNER:$PPWD/thirdparty/DAZZ_DB:$PPWD/thirdparty/DEXTRACTOR/:$PPWD/thirdparty/DASCRUBBER"
export PATH="$PATH:$PPWD/inst/bin"
export MANPATH="$MANPATH:$PPWD/inst/share/man"

#export PATH='/home/raymond/devel/polish/racon/racon/bin/':$PATH
#export PATH='/home/raymond/devel/miniasm/minimap/':$PATH

inputFile=$1
outputDir=$2
coverage=$3
nominal=$4

cd $outputDir 
# reads.fasta should be in data/$coverage
fasta2DB $coverage $inputFile
DBsplit -x500 -s100 $coverage     
HPC.daligner -t5 $coverage | csh -v
DASqv -c100 $coverage $coverage.las
##echo 1
# Run filter

mkdir log
hinge filter --db $coverage --las $coverage.las -x $coverage --config $nominal
echo 'hinge maximal'
# Get maximal reads

hinge maximal --db $coverage --las $coverage.las -x $coverage --config $nominal
echo 'hinge layout'
# Run layout

hinge layout --db $coverage --las $coverage.las -x $coverage --config $nominal -o $coverage
echo 'clip-nanopore'
# Run postprocessing

hinge clip-nanopore $coverage.edges.hinges $coverage.hinge.list _temp_
echo 'draft assembly'

# get draft assembly 

hinge draft-path . $coverage ${coverage}_temp_.G2.graphml
hinge draft --db $coverage --las $coverage.las --prefix $coverage --config $nominal --out $coverage.draft
echo 'get consensus assembly'
# get consensus assembly

hinge correct-head $coverage.draft.fasta $coverage.draft.pb.fasta draft_map.txt
fasta2DB draft $coverage.draft.pb.fasta 
HPC.daligner $coverage draft | zsh -v  
hinge consensus draft $coverage draft.$coverage.las $coverage.consensus.fasta $nominal
hinge gfa . $coverage $coverage.consensus.fasta
echo 7
#results should be in $coverage_consensus.gfa


