# Chroloplast-genome-assembly

This pipeline is used to assemble cp genome with short read (Unicycler) or long read (Canu and Hinge) or hybrid (Unicycler). For long read only assembly, according to our analysis, Hinge is better than Canu, but the min coverage should be **>=20x**. In our analysis, long read only assembly result needs to be manually removed duplication and rearranged structure. For short read only assembly, we suggest that the short read coverage should be >= 20x. **When the short read coverage >=20x, Unicycler is able to correct distinguish short single copy (ssc), long singel copy (lsc) and inverst repeats (ir), which returns three individual, complete contigs (ssc, lsc and ir)**. In terms to hybrid assembly, we recommend using **>=20x coverage short read and >=20x coverage of long read**, which can provide **a single, no duplication and complete cp genome**. In general, the direction of lsc or ssc in hybrid assembly result could be different compared to other published cp genome. This could be cause by heteroplasmy. The structure of the hybrid assembly result could be easy changed according to the direction of gene exclusively located in lsc/ssc/ir (we provide an example script to do this (3\_post\_assembly/1\_same\_structure/direction.py), but this script cannot fix all conditions.) 

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
- IQ-TREE v1.5.5
- Cluster Omega v1.2.4

## Procedure

1\_pre\_assembly: the input could be the raw long/short reads. Reads will be trimed the low quality region and adapters. Fastqc is used to visualize the read quality. Next, the processed read could be mapped to cp reference genomes (custom) to get the cp reads. Assuming the original reads contain non-cp genome, such at mtDNA or nuclear DNA.

2\_assembly: 100x short reads are separated for assembly quality check in the 3_post_assembly. https://github.com/roblanf/splitreads provides scripts for randomly select short/long reads for assembly. The long read only assembl were using Hinge and Canu, whereas short read only and hybrid assembly were using Unicycler

3\_post_assembly: the raw assembly was processed to create a single contig, which is no duplication and has same structure with cp reference genome (details see 3_post_assembly/README. Next, the 100x short reads kept from previous step was used to mapped to processed assembly to evaulate the assembly quality by checking the mapping rate, error rate and mismatch/insertion/deletion using Qualimap. 

# Running
A demo run for assembling cp genome is the following:
## LONG READ ONLY ASSEMBLY

#The inputFile is long.fastq.gz, in a dir ~/data/, species is E.pau
### 1\_pre\_assembly
### QualityControl

run fastqc to check the quality
```
cd 1_pre_assembly/1_qualityControl/longRead/
mkdir 1_qualityCheck/result
./1_qualityCheck/run_fastqc.sh \
  ~/data/long.fastq.gz \
  1_qualityCheck/result
```
trim adaptor
```
mkdir 2_adapterTrim/result
./2_adapterTrim/run_porechop.sh \
   ~/data/long.fastq.gz \
   2_adapterTrim/result/long.trim.fastq.gz
```
trim low qualty region (<9) and read <= 1kb
```
mkdir 3_qualityTrim/result
./3_quaityTrim/run_nanoFilt.sh \
  2_adapterTrim/result/long.trim.fastq.gz \
  3_qualityTrim/result 9 1000
```
rerun fastqc to check the data again
```
mkdir 4_qualityCheck/result
./_qualityCheck/run_fastqc.sh \
  3_qualityTrim/result/long.trim.fastq.gz \
  4_qualityCheck/result
```

### cp\_DNA\_extraction

assume the ref.fa (other cp genomes, should be double-up, in case read maps to the 'cut-point') is in 1\_pre\_assembly/2\_cpDNAExtraction/longRead/ref/
```
cd ../../2_cpDNAExtraction/longRead
mkdir result
```
use 10 threads, minMatch is 15, minAlnLength is 1kb (Blasr does not support the gz format, ungzip first, can be gzipped again after mapping)
```
pigz -d ../../1_qualityControl/longRead/3_qualityTrim/result/*gz
./1_run_Blasr.sh \
   ../../1_qualityControl/longRead/3_qualityTrim/result/ \
   result \
   10 \
   ref/ref.fa \
   15 \
   1000  
```
get cp read from the Blasr output
```
./2_resultParse.py \
  result/long.trim.out \
  ~/data/long.fastq \
  result/long.fasta
```
the cp reads are in the long.fasta

### 2\_assembly
```
cd ../../../2_assembly/longReadOnly
```
canu assembly, assume cp genome size is 160kb, corOutCoverage is 40, correctedErrorRate is 0.154, the path of gnuplot is gunPlotPath, use 30 threads. The final assembly is Epau.contigs.fasta
```
cd canu
./run_canu.sh \
  ../../../1_pre_assembly/2_cpDNAExtraction/longRead/result/long.fasta \
  result \
  Epau \
  160kb \
  40 \
  0.154 \
  30 \
  gunPlotPath
```
hinge assembly. If MinION reads, run 1\_convertName.sh to get the correct header that hinge can recoginzed.
```
cd ../hinge
./1_convertName.sh ../../../1_pre_assembly/2_cpDNAExtraction/longRead/result .
```
the reads with new header are in long.trim.pacbioName.fasta. Assume the coverage is 20x nominal is nominal.ini (can be found in hinge install dir) in this dir. The final assembly is 40coverage.consensus.fasta
```
2_run_hinge.sh \
  long.trimm.pacbioName.fasta \
  result \
  20 \
  nominal.ini
```
### 3\_post\_assembly
```
cd ../../../3_post_assembly
```
### mummer
first, use mummer to check the contig alignment. Assume to cp genome which is used to compare is 3\_post\_assembly/1\_same\_structure/ref/ref.fa, the assembly result is 3\_post\_assembly/1\_same\_structure/assembly.contig.fasta (can be the canu/hinge result)
```
cd 1_same_structure
mkdir result
./mummer_plot.sh \
  assembly.contig.fasta \
  ref/ref.fa \
  result/assembly \
  result/assembly.png 
```
the alignment fig is result/assembly.png. According to the alignment suitation, choose direction.py or mummer\_direction.py to create a single contig which is no duplication and has the same sturcture of cp ref genome. The direction.py script is recommended for assembly which has clear lsc/ir/ssc contigs or a complete contig but the structure is different from cp ref genome. The direction.py changed the direction based on genes, the genes used to ensure direction can be changed in the script. The mummer\_direction.py is basing on mummer alignement result to merge the congits. **These two scripts are VERY VERY VERY VERY depending on the original assembly, NOT SUITABLE FOR ALL SITUATIONS.**

run direction.py
```
python direction.py \
  assembly.contig.fasta \
  assembly_one_contig.fa
```
run mummer\_plot.sh
```
./mummer_direction.sh \
  result/assembly.coord  \
  assembly_one_contig.fa \
  assembly.contig.fasta
```
### polish
```
cd ../2_polish
```
we use Racon+Nanopolish here. Racon runs 10 iterations, whereas Nanopolish runs until the result unchanged.
Run Racon first, and then use the Racon-polish result as input to run Nanopolish. The number of iteration can be changed if the while loop : _if [ $n -gt 10 ]_. Nanopolish is MinION specific, if data is from Pacbio, Nanopolish can be changed to another polisher or just skip.

run Racon, use 10 threads. The path of inputFile (the cp reads) and ref (assembly result) should be absolute path. The ref should end with 'fa', if end with 'fasta', change all 'fa' into 'fasta' in the code. In addition, the read used to assemble should be changed to fastq format, which can be obtained read from 1\_pre\_assembly/3\_qualityTrim/result/long.trim.fastq.gz according to the read name. 
```
mkdir result_racon
./run_racon.sh \
  1_pre_assembly/2_cpDNAExtraction/longRead/result/long.fastq \
  result_racon \
  3_post_assembly/1_same_structure/assembly_one_contig.fa \
  10
```
run Nanopolish with 10 threads. Nanopolish needs the index data (link to original Fast5 data, details see [nanopolish] (https://github.com/jts/nanopolish)). Assume the final polish result of Racon is  result\_racon/round10/result/assembly.polished.racon.fa
```
nanopolish index \
  -d /path/to/raw_fast5s/ \
  1_pre_assembly/2_cpDNAExtraction/longRead/result/long.fasta
mkdir result_nanopolish
./run_nanopolish.sh \
  1_pre_assembly/2_cpDNAExtraction/longRead/result/long.fasta \
  result_nanopolish \
  result\_racon/round10/result/assembly.polished.racon.fa 10
```
### assembly\_quality\_control
finally, we use the 100x short read (or other coverage of short read) to remap to the assembly to assess its quality. If no short read, this step can be ignored. Using 10 threads. Qualimap is used to grep the mapping information.
```
cd ../3_assembly_quality_control/
mkdir result_bowtie2
./1_run_bowtie2.sh \
  /path/to/shortReadR1 \
  /path/to/shortReadR2 \
  result_bowtie2 \
  /path/to/polished_assembly \
  10
mkdir result_qualimap
./2_run_qualimap.sh \
  result_bowtie2 \
  result_qualimap \
  10
```

## SHORT READ ONLY ASSEMBLY

#The inputFiles are R1.fastq.gz, R2.fastq.gz in a dir ~/data/
### 1\_pre\_assembly
### QualityControl

run fastqc to check the quality
```
cd 1_pre_assembly/1_qualityControl/shortRead/
mkdir 1_qualityCheck/result
./1_qualityCheck/run_fastqc.sh ~/data/R1.fastq.gz 1_qualityCheck/result
./1_qualityCheck/run_fastqc.sh ~/data/R2.fastq.gz 1_qualityCheck/result
```
trim adaptor and low quality region (<30). Read length <50 bp will be removed. Parameters can be changed in the script. Script will auto get the R2 read if they have the same name (R\*.fastq.gz). 10 threads will be used.
```
mkdir 2_adapterTrim/result
./2_adapterTrim/run_bbduk.sh \
  ~/data/R1.fastq.gz \
  2_adapterTrim/result/ \
  50 \
  /path/to/bbduk/adapterDB \
  10
```
rerun fastqc to check the data again
```
mkdir 3_qualityCheck/result
./_qualityCheck/run_fastqc.sh 3_qualityTrim/result/R1.trim.fastq.gz 3_qualityCheck/result
./_qualityCheck/run_fastqc.sh 3_qualityTrim/result/R2.trim.fastq.gz 3_qualityCheck/result
```

### cp\_DNA\_extraction

assume the ref.fa (other cp genomes, should be double-up, in case read maps to the 'cut-point') is in 1\_pre\_assembly/2\_cpDNAExtraction/shortRead/ref/
```
cd ../../2_cpDNAExtraction/shortRead
mkdir result
```
Bowtie2 is used, use 10 threads.R1 and R2 should have the same name, such as R\*.trim.fastq.gz.
```
mkdir result
./1_run_bowtie2.sh \
  ../../1_qualityControl/shortRead/2_adapterTrim/result/R1.trim.fastq.gz \
  result \
  10 \
  ref/ref.fa  
```
get cp read from the Bowtie2 output
```
mkdir cpRead
./2_getCPRead.py \
  ref/ref.fa \
  result \
  cpRead
```
the cp reads are in cpRead/R\*.trim.fastq.gz

### 2\_assembly
```
cd ../../../2_assembly/shortReadOnly
```
Unicycler is used to do the short read only assembly. In the default setting, the read is corrected by SPAdes, and assembly is polished by Pilon. The final assembly is result/assembly.fasta
```
#10 threads are used
./run_shortRead_unicycler.sh \
 ../../1_pre_assembly/2_cpDNAExtraction/shortRead/cpRead/R1.trim.fastq.gz \
 ../../1_pre_assembly/2_cpDNAExtraction/shortRead/cpRead/R2.trim.fastq.gz \
 result \
 10
```
### 3\_post\_assembly
```
cd ../../3_post_assembly
```
### mummer
as described above in long read only assembly part.
### polish
```
cd 2_polish
```
we use Pilon to polish the assembly, run until result unchanged. Using 10 threads. 
```
mkdir result_pilon
./run_pilon.sh \
  ../../1_pre_assembly/2_cpDNAExtraction/shortRead/cpRead/R1.trim.fastq.gz \
  ../../1_pre_assembly/2_cpDNAExtraction/shortRead/cpRead/R2.trim.fastq.gz \
  result_pilon \
  /path/to/final_assembly_single_contig.fa \
  10 \
  checkDiff.py
```
### assembly\_quality\_control
As described above in the long read only assembly part. NOTE: the short read used to remap should be **unuse** in the assembly. Read randomly separate can use script in https://github.com/roblanf/splitreads.

## HYBRID ASSEMBLY
how to get cp short read and cp long read are described above (1\_pre\_assembly). In general, every step in the hybrid assembly is the same as short read only assembly (including the 3\_post\_assembly). Only the assembly script is different.
```
cd 2_assembly/hybrid/
./run_unicycler.sh \
  ../../1_pre_assembly/2_cpDNAExtraction/longRead/result/long.fasta \
  ../../1_pre_assembly/2_cpDNAExtraction/shortRead/cpRead/R1.trim.fastq.gz \
  ../../1_pre_assembly/2_cpDNAExtraction/shortRead/cpRead/R2.trim.fastq.gz \
  result \
  10
```

## PHYLOGENETIC ANALYSIS
First, split the different chloroplast genomes by exon and non-exon, and then put the same region together (in fasta format), run run\_clustalo.sh first to get the alignment from different chloroplast genome, and then check and fix the alignment manually. After that, combine all the fragment fasta file together to create a corrected whole genome alignment (in fasta format), and create a nexus file (optional) according to the alignment. Finally, run run\_iqtree.sh with the corrected whole genome alignment file nexus to get the phylogenetic tree. 
