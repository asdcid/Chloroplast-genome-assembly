# Chroloplast-genome-assembly

## Requirements
- Fastqc
- Bbduk v37.31
- Porechop v0.2.1
- Nanofilt v1.2.0
- Bowtie2 v2.2.6
- Blasr v5.1
- Bwa v0.7.15
- samtools v1.5
- Hinge
- Canu v1.6
- Unicycler v0.3.1
- Mummer v3.23
- Racon
- Nanopolish v0.8.1
- Pilon v1.22
- Qualimap v2.2.1

## Procedure

1_pre_assembly: the input could be the raw long/short reads. Reads will be trimed the low quality region and adapters. Fastqc is used to visualize the read quality. Next, the processed read could be mapped to cp reference genomes (custom) to get the cp reads. Assuming the original reads contain non-cp genome, such at mtDNA or nuclear DNA.

2_assembly: 100x short reads are separated for assembly quality check in the 3_post_assembly. https://github.com/roblanf/splitreads provides scripts for randomly select short/long reads for assembly. The long read only assembl were using Hinge and Canu, whereas short read only and hybrid assembly were using Unicycler

3_post_assembly: the raw assembly was processed to create a single contig, which is no duplication and has same structure with cp reference genome (details see 3_post_assembly/README. Next, the 100x short reads kept from previous step was used to mapped to processed assembly to evaulate the assembly quality by checking the mapping rate, error rate and mismatch/insertion/deletion using Qualimap. 

# Running
A demo run for assembling cp genome is the following:
## Long read only assembly

#The inputFile is long.fastq.gz, in a dir ~/data/

### QualityControl

run fastqc to check the quality
```
cd 1\_pre\_assembly/1\_qualityControl/longRead/

mkdir 1\_qualityCheck/result

./1\_qualityCheck/run\_fastqc.sh ~/data/long.fastq.gz 1\_qualityCheck/result
```
trim adaptor
```
mkdir 2\_adapterTrim/result

./2\_adapterTrim/run\_porechop.sh ~/data/long.fastq.gz 2\_adapterTrim/result/long.trim.fastq.gz
```
trim low qualty region (<9) and read < 1kb
```
mkdir 3\_qualityTrim/result

./3\_quaityTrim\run\_nanoFilt.sh 2\_adapterTrim/result/long.trim.fastq.gz 3\_qualityTrim/result 9 1000
```
rerun fastqc to check the data again
```
mkdir 4\_qualityCheck/result

./4\_qualityCheck/run\_fastqc.sh 3\_qualityTrim/result/long.trim.fastq.gz 4\_qualityCheck/result
```

## cp\_DNA\_extraction

assume the ref.fa (other cp genomes, should be double-up, in case read maps to the 'cut-point') is in 1\_pre\_assembly/2\_cpDNAExtraction/longRead/ref/
```
cd ../../2\_cpDNAExtraction/longRead

mkdir result

```

#use 10 threads, minMatch is 15, minAlnLength is 1kb
```
./1\_run\_Blasr.sh ../../1\_qualityControl/longRead/3\_qualityTrim/result/ result 10 ref/ref.fa 15 1000  
```
