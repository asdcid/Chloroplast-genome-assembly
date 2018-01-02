#!/bin/bash



export PATH="/home/raymond/devel/bbmap/bbmap/":$PATH

inputf='/home/raymond/work/Eucalyptus_pauciflora/genome/data/illumina/SN877_0428_RLanfear_RSB_Eucalyptus_gDNA/select/'
outputtrimreads='428_result/'
minlen=50
adaptors='/home/raymond/devel/bbmap/bbmap/resources/adapters.fa'
trimq=30
threads=40



for in1 in $(find $inputf -name "*R.fastq.gz"); do
    in2=${in1%%R1.fastq.gz}"R2.fastq.gz"
    echo "running bbduk on"
    echo $in1
    echo $in2

    f1=$(basename ${in1%%R1.fastq.gz}"R1.trim.fastq.gz")
    f2=$(basename ${in1%%R1.fastq.gz}"R2.trim.fastq.gz")

    out1=$outputtrimreads/$f1
    out2=$outputtrimreads/$f2
    sampleid=$outputtrimreads${f1%%R1.trim.fastq.gz}

    bbduk.sh in1=$in1 in2=$in2 out1=$out1 out2=$out2 minlen=$minlen k=25 mink=8 ktrim=r ref=$adaptors hdist=1 overwrite=f qtrim=rl trimq=$trimq t=$threads bhist=$sampleid"bhist.txt" qhist=$sampleid"qhist.txt" gchist=$sampleid"gchist.txt" aqhist=$sampleid"aqhist.txt" lhist=$sampleid"lhist.txt" > $sampleid"bbduk_log.txt"

done


