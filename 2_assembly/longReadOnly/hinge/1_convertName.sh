#!/bin/bash


#from normal name to PabBio specific

#inputDir='/home/raymond/work/Eucalyptus_pauciflora/genome/bin/assembly/hinge/data/'
#outputDir='data/pacbioName'

##the file in inputDir should be fasta, fasta.gz format has not tested
inputDir=$1
outputDir=$2

cd $outputDir
for inputFile in $inputDir/*.fasta
do
    outputFile=$(basename ${inputFile%%.fasta}).pacbioName.fasta
    ./falcon_name_fasta.pl -i $inputFile -o $outputFile


done


