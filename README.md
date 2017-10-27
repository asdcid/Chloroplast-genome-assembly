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

1\_pre\_assembly: the input could be the raw long/short reads. Reads will be trimed the low quality region and adapters. Fastqc is used to visualize the read quality. Next, the processed read could be mapped to cp reference genomes (custom) to get the cp reads. Assuming the original reads contain non-cp genome, such at mtDNA or nuclear DNA.

2\_assembly: 100x short reads are separated for assembly quality check in the 3_post_assembly. https://github.com/roblanf/splitreads provides scripts for randomly select short/long reads for assembly. The long read only assembl were using Hinge and Canu, whereas short read only and hybrid assembly were using Unicycler

3\_post_assembly: the raw assembly was processed to create a single contig, which is no duplication and has same structure with cp reference genome (details see 3_post_assembly/README. Next, the 100x short reads kept from previous step was used to mapped to processed assembly to evaulate the assembly quality by checking the mapping rate, error rate and mismatch/insertion/deletion using Qualimap. 

# Running
A demo run for assembling cp genome is the following:
## Long read only assembly

#The inputFile is long.fastq.gz, in a dir ~/data/, species is E.pau

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
trim low qualty region (<9) and read <= 1kb
```
mkdir 3_qualityTrim/result
./3_quaityTrim/run_nanoFilt.sh 2_adapterTrim/result/long.trim.fastq.gz 3_qualityTrim/result 9 1000
```
rerun fastqc to check the data again
```
mkdir 4_qualityCheck/result
./_qualityCheck/run_fastqc.sh 3_qualityTrim/result/long.trim.fastq.gz 4_qualityCheck/result
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
./1_run_Blasr.sh ../../1_qualityControl/longRead/3_qualityTrim/result/ result 10 ref/ref.fa 15 1000  
```
get cp read from the Blasr output
```
./2_resultParse.py result/long.trim.out ~/data/long.fastq result/long.fasta
```
the cp reads are in the long.fasta

### 2\_assembly
```
cd ../../../2_assembly/longReadOnly
```
canu assembly, assume cp genome size is 160kb, corOutCoverage is 40, correctedErrorRate is 0.154, the path of gnuplot is gunPlotPath, use 30 threads. The final assembly is Epau.contigs.fasta
```
cd canu
./run_canu.sh ../../../1_pre_assembly/2_cpDNAExtraction/longRead/result/long.fasta result Epau 160kb 40 0.154 30 gunPlotPath
```
hinge assembly. If MinION reads, run 1\_convertName.sh to get the correct header that hinge can recoginzed.
```
cd ../hinge
./1_convertName.sh ../../../1_pre_assembly/2_cpDNAExtraction/longRead/result .
```
the reads with new header are in long.trim.pacbioName.fasta. Assume the coverage is 20x nominal is nominal.ini (can be found in hinge install dir) in this dir. The final assembly is 40coverage.consensus.fasta
```
2_run_hinge.sh long.trimm.pacbioName.fasta result 20 nominal.ini
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
./mummer_plot.sh assembly.contig.fasta ref/ref.fa result/assembly result/assembly.png 
```
the alignment fig is result/assembly.png. According to the alignment suitation, choose direction.py or mummer\_direction.py to create a single contig which is no duplication and has the same sturcture of cp ref genome. The direction.py script is recommended for assembly which has clear lsc/ir/ssc contigs or a complete contig but the structure is different from cp ref genome. The direction.py changed the direction based on genes, the genes used to ensure direction can be changed in the script. The mummer\_direction.py is basing on mummer alignement result to merge the congits. **These two scripts are VERY VERY VERY VERY depending on the original assembly, NOT SUITABLE FOR ALL SITUATIONS.**

run direction.py
```
python direction.py assembly.contig.fasta assembly_one_contig.fa
```
run mummer\_plot.sh
```
./mummer_direction.sh result/assembly.coord  assembly_one_contig.fa assembly.contig.fasta
```
### polish
```
cd ../2_polish
```
we use Racon+Nanopolish here. Racon runs 10 iterations, whereas Nanopolish runs until the result unchanged.
Run Racon first, and then use the Racon-polish result as input to run Nanopolish. The number of iteration can be changed if the while loop : _if [ $n -gt 10 ]_. Nanopolish is MinION specific, if data is from Pacbio, Nanopolish can be changed to another polisher or just skip.

run Racon, use 10 threads. The path of inputFile (the cp reads) and ref (assembly result) should be absolute path. The ref should end with 'fa', if end with 'fasta', change all 'fa' into 'fasta' in the code.
```
mkdir result_racon
./run_racon.sh 1_pre_assembly/2_cpDNAExtraction/longRead/result/long.fasta result_racon 3_post_assembly/1_same_structure/assembly_one_contig.fa 10
```
run Nanopolish with 10 threads. Nanopolish needs the index data (link to original Fast5 data, details see [nanopolish] (https://github.com/jts/nanopolish)). Assume the final polish result of Racon is  result\_racon/round10/result/assembly.polished.racon.fa
```
nanopolish index -d /path/to/raw_fast5s/ 1_pre_assembly/2_cpDNAExtraction/longRead/result/long.fasta
mkdir result_nanopolish
./run_nanopolish.sh 1_pre_assembly/2_cpDNAExtraction/longRead/result/long.fasta result_nanopolish result\_racon/round10/result/assembly.polished.racon.fa 10
```
### assembly\_quality\_control
finally, we use the 100x short read (or other coverage of short read) to remap to the assembly to assess its quality. If no short read, this step can be ignored. Using 10 threads. Qualimap is used to grep the mapping information.
```
cd ../3_assembly_quality_control/
mkdir result_bowtie2
./1_run_bowtie2.sh /path/to/shortReadR1 /path/to/shortReadR2 result_bowtie2 /path/to/polished_assembly 10
mkdir result_qualimap
./2_run_qualimap.sh result_bowtie2 result_qualimap 10
```

