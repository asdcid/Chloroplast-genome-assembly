export PATH="path/to/karect/":$PATH

threads=40
outputDir='result'
temp='temp'

mkdir $temp

karect \
    -correct \
     -inputfile='R1_short_read.fastq' \
     -inputfile='R2_short_read.fastq' \
    -resultdir=$outputDir \
    -tempdir=$temp \
    -threads=$threads \
    -celltype='haploid' \
    -matchtype='hamming'
