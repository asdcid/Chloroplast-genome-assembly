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
cd 1_pre_assembly/1_qualityControl/longRead/
mkdir 1_qualityCheck/result
./1_qualityCheck/run_fastqc.sh ~/data/long.fastq.gz 1_qualityCheck/result
```
trim adaptor
```
mkdir 2_adapterTrim/result
./2_adapterTrim/run_porechop.sh ~/data/long.fastq.gz 2_adapterTrim/result/long.trim.fastq.gz
```
trim low qualty region (<9) and read < 1kb
```
mkdir 3_qualityTrim/result
./3_quaityTrim/run_nanoFilt.sh 2_adapterTrim/result/long.trim.fastq.gz 3_qualityTrim/result 9 1000
```
rerun fastqc to check the data again
```
mkdir 4_qualityCheck/result
./\_qualityCheck/run_fastqc.sh 3_qualityTrim/result/long.trim.fastq.gz 4_qualityCheck/result
```

## cp\_DNA\_extraction

assume the ref.fa (other cp genomes, should be double-up, in case read maps to the 'cut-point') is in 1\_pre\_assembly/2\_cpDNAExtraction/longRead/ref/
```
cd ../../2_cpDNAExtraction/longRead
mkdir result

```
#use 10 threads, minMatch is 15, minAlnLength is 1kb (Blasr does not support the gz format, ungzip first, can be gzipped again after mapping)
```
pigz -d ../../1_qualityControl/longRead/3_qualityTrim/result/*gz
./1_run_Blasr.sh ../../1_qualityControl/longRead/3_qualityTrim/result/ result 10 ref/ref.fa 15 1000  
```
