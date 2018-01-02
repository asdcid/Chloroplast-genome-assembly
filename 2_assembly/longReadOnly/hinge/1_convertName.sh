#!/bin/bash


#from normal name to PabBio specific
##the file in inputDir should be fasta, fasta.gz format has not tested

inputDir='/home/raymond/work/Eucalyptus_pauciflora/genome/bin/assembly/hinge/data/'
outputDir='data/pacbioName'



cd $outputDir
for inputFile in $inputDir/*.fasta
do
    outputFile=$(basename ${inputFile%%.fasta}).pacbioName.fasta
    ./falcon_name_fasta.pl -i $inputFile -o $outputFile


done


