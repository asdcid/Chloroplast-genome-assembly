#!/bin/bash



inputFile=$1
outputDir='result/'
threads=20



fastqc $inputFile -o $outputDir -t $threads

