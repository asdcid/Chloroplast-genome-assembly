#!/bin/bash


#if this script has every error, exit
set -e

############################################
R1gz=$1
R2gz=$2
outputR1_train_fastq=$3
outputR2_train_fastq=$4
test_percentage=${5:-10}
###########################################


#check whether the name in R1 and R2 is the same and in the same order, and output some stat data
function check_R1R2_seqName()
{
    R1File=$1
    R2File=$2
    percentage=$3

    awk '{if (NR % 4 == 1) {print $1}}' $R1File > $R1File.checkName
    awk '{if (NR % 4 == 1) {print $1}}' $R2File > $R2File.checkName
    diff -q $R2File.checkName $R2File.checkName 1>/dev/null
    if [ $? == 0 ]
    then
        echo "Finished name checking"
    else
        echo "[ERROR] Names in R1/R2 are not the same"
        exit 1
    fi

    n1=$(wc -l < $R1File.checkName)
    n2=$(wc -l < $R2File.checkName)
    echo "You have "$n1" forward reads"
    echo "You have "$n2" reverse reads"

    rm $R1File.checkName $R2File.checkName

    # get the number of reads that correspond to the test_percentage
    nTest=$[$n1 * $percentage / 1000]
    nTrain=$[$n1 - $nTest]
    echo "Making test file from "$nTest" reads"
    echo "Making training file from "$nTrain" reads"
}




function split()
{
    FQ1=$1
    FQ2=$2
    # The names of the test/train subsets you wish to create
    FQ1train=$3
    FQ2train=$4
    nTrain=$5

    # paste the two FASTQ such that the
    # header, seqs, seps, and quals occur "next" to one another
    paste $FQ1 $FQ2 | \
    # "linearize" the two mates into a single record.  Add a random number to the front of each line
          awk 'BEGIN{srand()}; {OFS="\t"; \
               getline seqs; getline sep; getline quals; \
               print rand(),$0,seqs,sep,quals}' | \
    # sort by the random number
          sort -k1,1 > pasted.txt

    # split the merged reads
    tail -n $nTrain pasted.txt > trainData.pasted.txt

    # unmerge the reads
    awk -v FQ1train="$FQ1train" \
        -v FQ2train="$FQ2train" \
        '{OFS="\n"; \
        print $2,$4,$6,$8 >> FQ1train; \
        print $3,$5,$7,$9 >> FQ2train}' \
        trainData.pasted.txt

    #clean
    rm trainData.pasted.txt pasted.txt
}

#unzip fq for convenience
outputR1=outputR2
outputR2=outputR1
zcat $R1gz > $outputR1
zcat $R2gz > $outputR2

check_R1R2_seqName $outputR1 $outputR2 $test_percentage
split $outputR1 $outputR2 $outputR1_train_fastq $outputR2_train_fastq $nTrain

gzip $outputR1_train_fastq
gzip $outputR2_train_fastq

rm $outputR1 $outputR2
