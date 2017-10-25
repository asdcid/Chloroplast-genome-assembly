#!/bin/bash



#export PATH='/home/raymond/devel/nanopolish/nanopolish/':$PATH
export PATH='nanopolish/path':$PATH



polish()
{
    inputFile=$1
    outputDir=$2
    ref=$3
    threads=$4

    mkdir $outputDir

    echo $ref
    echo $inputFile
    refDir=$outputDir/ref
    bwaDir=$outputDir/bwa
    resultDir=$outputDir/result
    tempDir=$outputDir/temp

    mkdir $refDir
    mkdir $bwaDir
    mkdir $resultDir
    mkdir $tempDir

    coverage=$(basename ${inputFile%%.fastq})

    #index
    ln $ref $refDir
    ref=$refDir/*fa
    build_index $ref

    #run bwa
    outputBwa=$bwaDir/$coverage.sort.bam
    run_bwa $inputFile $bwaDir $ref $threads $coverage $outputBwa

    #run nanopolish
    outputNanopolish=$resultDir/$coverage.polished.nanopolish.fa
    tempFile=$tempDir/temp
    run_nanopolish $inputFile $outputNanopolish $outputBwa $ref $tempFile $threads
    echo $ref
    echo $outputBwa
}

build_index()
{
    indexFile=$1
    bwa index $indexFile
}

run_bwa()
{
    inputFile=$1
    outputDir=$2
    ref=$3
    threads=$4
    coverage=$5
    outputFile=$6

    bwa mem \
        -x ont2d \
        -t $threads \
        $ref \
        $inputFile > $outputDir/$coverage.sam

    samtools sort \
        -o $outputFile \
        $outputDir/$coverage.sam

    samtools index \
        $outputFile
}

run_nanopolish()
{
echo run_nanopolish
inputFile=$1
outputFile=$2
bwa_inputFile=$3
ref=$4
temp=$5
threads=$6

# path of nanopolish_makerange.py, should be change
python /home/raymond/devel/nanopolish/nanopolish/scripts/nanopolish_makerange.py \
    $ref \
    | \
    parallel \
    --results $temp.results \
    -P $threads \
    nanopolish variants \
    --consensus $temp.{1}.fa \
    -w {1} \
    -r $inputFile \
    -b $bwa_inputFile \
    -g $ref \
    -t $threads \
    --min-candidate-frequency 0.1 

# path of nanopolish_merge.py, should be change
python /home/raymond/devel/nanopolish/nanopolish/scripts/nanopolish_merge.py $temp.*.fa > $outputFile
}


inputFile=$1
outputDir=$2
ref=$3
threads=$4


cd $outputDir


n=1

outputDir=round-$n
echo polish
polish $inputFile $outputDir $ref $threads
echo "finish 1"

ref=round-$n/result/*fa
n=$(( $n + 1 ))
echo $n
outputDir=round-$n
polish $inputFile $outputDir $ref $threads
echo "finish 2"

LAST=round-$(( $n - 1 ))/result/*.fa 
echo $LAST
CURRENT=round-$n/result/*.fa
echo $CURRENT
DIFF=$(diff -q $LAST $CURRENT)
while [ "$DIFF" ]
do
    echo "begin round"
    ref=round-$n/result/*fa
    n=$(( $n + 1 ))
    echo round $n
    outputDir=round-$n
    polish $inputFile $outputDir $ref $threads

    LAST=round-$(( $n - 1 ))/result/*.fa 
    CURRENT=round-$n/result/*fa
    DIFF=$(diff -q $LAST $CURRENT)
    echo 'diff' $DIFF
done

cp round-$n/result/* ../final/
#rm -r round*
