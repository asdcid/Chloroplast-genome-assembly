# Chroloplast-genome-assembly

This pipeline is used to assemble Eucalyptus pauciflora chloroplast genome with short read (Unicycler) or long read (Canu and Hinge) or hybrid (Unicycler). If anyone want to use, they should change the input/output file in scripts. We randomly selected 5x, 8x, 10x, 20x, 40x, 60x, 80x, 100x, 200x, 300x, 400x and 500x coverage of short/long read to do the assembly to compare the difference. For long read only assembly, according to our analysis, Hinge is better than Canu, but the min coverage should be **>=20x**. In our analysis, long read only assembly result needs to be manually removed duplication and rearranged structure. For short read only assembly, we found that **when the short read coverage >=20x, Unicycler is able to correct distinguish short single copy (ssc), long singel copy (lsc) and inverst repeats (ir), which returns three individual, complete contigs (ssc, lsc and ir)**. In terms to hybrid assembly, we found that **>=20x coverage short read and >=20x coverage of long read** can provide **a single, no duplication and complete chloroplast genome**.) 

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

1\_pre\_assembly: the input could be the raw long/short reads. Reads will be trimed the low quality region and adapters. Fastqc is used to visualize the read quality. Next, the processed read could be mapped to chloroplast reference genomes (custom) to get the chloroplast reads. Assuming the original reads contain non-chloroplast genome, such at mtDNA or nuclear DNA.

2\_assembly: 100x short reads are separated for assembly quality check in the 3\_post\_assembly. 5x, 8x, 10, 20x, 40x, 60x, 80x, 100x, 200x, 300x, 400x and 500x coverage of short/long read were randomly selected for assembly. The long read only assembly was using Hinge and Canu, whereas short read only and hybrid assembly were using Unicycler.

3\_post\_assembly: the raw assembly was processed to create a single contig, which is no duplication and has same structure with cp reference genome (details see 3_post_assembly/README. Next, the 100x short reads kept from previous step was used to mapped to processed assembly to evaulate the assembly quality by checking the mapping rate, error rate and mismatch/insertion/deletion using Qualimap. 

# Running
## LONG READ ONLY ASSEMBLY

### 1\_pre\_assembly
### QualityControl

run fastqc to check the quality
```
./1_pre_assembly/1_qualityControl/longRead/1_qualityCheck/run_fastqc.sh 
```
trim adaptor
```
./1_pre_assembly/1_qualityControl/longRead/2_adapterTrim/run_porechop.sh 
```
trim low qualty region (<9) and read <= 5kb
```
./1_pre_assembly/1_qualityControl/longRead/3_quaityTrim/run_nanoFilt.sh 
```
rerun fastqc to check the data again
```
./1_pre_assembly/1_qualityControl/longRead/4_qualityCheck/run_fastqc.sh 
```

### cp\_DNA\_extraction

Mapping all trimmed reads to refs (31 known Eucalyptus chloroplast genomes, all double-up) to get the chloroplast reads.
```
./1_pre_assembly/2_cpDNAExtraction/longRead/1_run_Blasr.sh 
```
get chloroplast reads from the Blasr output, from Blasr output to fasta
```
./1_pre_assembly/2_cpDNAExtraction/longRead/2_resultParse.py 
```

### 2\_assembly
Randomly selected 5x, 8x, 10x, 20x, 40x, 60x, 80x, 100x, 200x, 300x 400x and 500x coverage reads (assume chloroplast genome size is 160kb).
```
./2_assembly/randomSelection/split_longRead.py
```
Run assembly with different coverage of long-read with Canu/Hinge.

Canu assembly, the final assembly is Epau.contigs.fasta
```
./2_assembly/longReadOnly/canu/run_canu.sh 
```
Hinge assembly. Due to the data is MinION read, run 1\_convertName.sh to get the correct header that hinge can recoginzed.
```
./2_assembly/longReadOnly/hinge/1_convertName.sh 
```
Using the reads with new header to do the assembly. nominal.ini (can be found in hinge install dir) should be in this dir. The final assembly is xx.consensus.fasta
```
./2_assembly/longReadOnly/hinge/2_run_hinge.sh 
```
### 3\_post\_assembly
### mummer
first, use mummer to check the contig alignment. E.reg, the most close known species to E.pau, was used as the reference to compare. 
```
./3_post_assembly/1_same_structure/mummer_plot.sh 
```
 According to the alignment suitation (alignment figure), choose direction.py or mummer\_direction.py to create a single contig which is no duplication and has the same sturcture of cp ref genome. The direction.py script is recommended for assembly which has clear lsc/ir/ssc contigs or a complete contig but the structure is different from cp ref genome. The direction.py changed the direction based on genes, the genes used to ensure direction can be changed in the script. The mummer\_direction.py is basing on mummer alignement result to merge the congits, which may produce extra errors. In general, all long-read assembly used mummer\_plot.sh, because lots of duplicates were observed in the alignment. However, for the short-read only and hybrid assembly, direction.py can be used in the assembly with >= 20x short-read coverage input. Since those assemlies have clear lsc/ssc/ir contigs or has one single contig without overlap (according to mummer result)

run mummer\_plot.sh
```
./3_post_assembly/1_same_structure/mummer_direction.sh 
```
### polish

we polished the genome from last step useing Racon, Nanopolish and Racon+Nanopolish, respectively. Racon run 10 iterations (the result keeps change after 10 rounds), whereas Nanopolish run until the result unchanged. The aligner used here is bwa mem.

run Racon. The path of inputFile (the chloroplast long-reads) and ref (assembly result) should be absolute path. The ref should end with 'fa'. In addition, the read used to assemble should be changed to fastq format. 
```
./3_post_assembly/2_polish/run_racon.sh 
```
run Nanopolish. Nanopolish needs the index data (link to original Fast5 data, details see [nanopolish] (https://github.com/jts/nanopolish)) first.
```
./3_post_assembly/2_polish/
./3_post_assembly/2_polish/run_nanopolish.sh 
```
For Racon+Nanopolish, run Racon first, and then use the Racon-polish result as input to run Nanopolish.

### assembly\_quality\_control

finally, we use the 100x short read (randomly selected first, not used in assembly, method see below) to remap to the assembly to assess its quality. Qualimap is used to grep the mapping information.
```
./3_assembly_quality_control/1_run_bowtie2.sh 
./3_assembly_quality_control/2_run_qualimap.sh 
```

## SHORT READ ONLY and hybrid ASSEMBLY

#The inputFiles are R1.fastq.gz, R2.fastq.gz in a dir ~/data/
### 1\_pre\_assembly
### QualityControl

run fastqc to check the quality
```
./1_pre_assembly/1_qualityControl/shortRread/1_qualityCheck/run_fastqc.sh
```
trim adaptor and low quality region (<30). Read length <50 bp will be removed. Parameters can be changed in the script. Script will auto get the R2 read if they have the same name (R\*.fastq.gz). 10 threads will be used.
```
./1_pre_assembly/1_qualityControl/shortRead/2_adapterTrim/run_bbduk.sh 
```
rerun fastqc to check the data again
```
./1_pre_assembly/3_qualityCheck/run_fastqc.sh 
```

### cp\_DNA\_extraction

assume the ref.fa (other cp genomes, should be double-up, in case read maps to the 'cut-point') is in 1\_pre\_assembly/2\_cpDNAExtraction/shortRead/ref/
Bowtie2 is used, use 10 threads.R1 and R2 should have the same name, such as R\*.trim.fastq.gz.
```
./1_pre_assembly/2_cpDNAExtraction/shortRead/1_run_bowtie2.sh 
```
get cp read from the Bowtie2 output
```
./1_pre_assembly/2_cpDNAExtraction/shortRead/2_getCPRead.py 
```
the cp reads are in cpRead/R\*.trim.fastq.gz

### 2\_assembly
Unicycler is used to do the short read only assembly. In the default setting, the read is corrected by SPAdes, and assembly is polished by Pilon. The final assembly is result/assembly.fasta
```
./2_assembly/shortReadOnly/2_run_shortRead_unicycler.sh 
```
### 3\_post\_assembly
### mummer
as described above in long read only assembly part.

run direction.py
```
python ./3_post_assembly/1_same_structure/direction.py 
```
run mummer\_plot.sh
```
./3_post_assembly/1_same_structure/mummer_direction.sh 
```
### polish
we use Pilon to polish the assembly, run until result unchanged. Using 10 threads. 
```
./3_post_assembly/2_polish/run_pilon.sh 
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
