#!/bin/bash

export PATH='/home/raymond/devel/python3/install/bin/':$PATH
export PATH='/home/raymond/devel/unicycler/Unicyclero/':$PATH
export PATH='/home/raymond/devel/unicycler/three_party/SPAdes/SPAdes-3.10.1-Linux/bin/':$PATH
export PATH='/home/raymond/devel/polish/pilon/':$PATH



inputFile_R1=$1
inputFile_R2=$2
outputDir=$3
threads=$4

start=$(date +%s)

unicycler-runner.py \
    -1 $inputFile_R1 \
    -2 $inputFile_R2 \
    -o $outputDir \
    --no_correct \
    -t $threads

end=$(date +%s)
DIFF=$(( $end - $start ))
echo $outputDir/$id "Excution time: " $DIFF " seconds." >> runTime.log
