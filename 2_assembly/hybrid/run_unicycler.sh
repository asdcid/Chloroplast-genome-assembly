#!/bin/bash


#export PATH='/home/raymond/devel/python3/install/bin/':$PATH
#export PATH='/home/raymond/devel/unicycler/Unicycler_old/':$PATH
#export PATH='/home/raymond/devel/unicycler/three_party/SPAdes/SPAdes-3.10.1-Linux/bin/':$PATH
#export PATH='/home/raymond/devel/polish/pilon/':$PATH

export PATH='unicycler/path':$PATH
export PATH='SPAdes/path':$PATH
export PATH='pilon/path':$PATH

inputFiles_longRead=$1
inputFile_R1=$2
inputFile_R2=$3
outputDir=$4
threads=$5
start=$(date +%s)

unicycler-runner.py \
    -1 $inputFile_R1 \
    -2 $inputFile_R2 \
    -l $inputFiles_longRead \
    -o $outputDir \
    -t $threads

end=$(date +%s)
DIFF=$(( $end - $start ))
echo $outputDir/$id "Excution time: " $DIFF " seconds." >> runTime.log




