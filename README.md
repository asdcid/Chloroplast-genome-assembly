# Chroloplast-genome-assembly

This pipeline is used to assemble Eucalyptus pauciflora chloroplast genome with short-read (Unicycler) or long-read (Canu and Hinge) or hybrid (Unicycler). If anyone want to use, they should change the input/output file in scripts. We randomly selected 5x, 8x, 10x, 20x, 40x, 60x, 80x, 100x, 200x, 300x, 400x and 500x coverage of long-/short- read to do the assembly to compare the difference. For long-read only assembly, according to our analysis, Hinge is better than Canu, but the min coverage should be **>=20x**. In our analysis, long-read only assembly result needs to be manually removed duplication and rearranged structure. For short-read only assembly, we found that **when the short-read coverage >=20x, Unicycler is able to correct distinguish short single copy (ssc), long singel copy (lsc) and inverst repeats (ir), which returns three individual, complete contigs (ssc, lsc and ir)**. In terms to hybrid assembly, we found that **>=20x coverage shor-read and >=20x coverage of long read** can provide **a single, no duplication and complete chloroplast genome**.) 

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

1\_pre\_assembly: the input could be the raw long-/short- reads. Reads will be trimed the low quality region and adapters. Fastqc is used to visualize the read quality. Next, the processed read could be mapped to chloroplast reference genomes (custom) to get the chloroplast reads, because the original reads contain non-chloroplast genome, such at mtDNA or nuclear DNA.

2\_assembly: 100x short-reads were separated for assembly quality check in the 3\_post\_assembly. 5x, 8x, 10, 20x, 40x, 60x, 80x, 100x, 200x, 300x, 400x and 500x coverage of long-/short- read were randomly selected for assembly. The long-read only assembly was using Hinge and Canu, whereas short-read only and hybrid assembly were using Unicycler.

3\_post\_assembly: the raw assembly was processed to create a single contig, which is no duplication and has same structure with cp reference genome (details see 3_post_assembly/README. Next, the 100x short-reads kept from previous step was used to mapped to processed assembly to evaulate the assembly quality by checking the mapping rate, error rate and mismatch/insertion/deletion using Qualimap. 

# Running
## LONG-READ ONLY ASSEMBLY

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

Mapping all trimmed reads to refs (31 known Eucalyptus chloroplast genomes, all double-up, in case read maps to the 'cut-point') to get the chloroplast reads, using Blasr.
```
./1_pre_assembly/2_cpDNAExtraction/longRead/1_run_Blasr.sh 
```
Get chloroplast reads from the Blasr output, from Blasr output to fasta.
```
./1_pre_assembly/2_cpDNAExtraction/longRead/2_resultParse.py 
```

### 2\_assembly
Randomly selected 5x, 8x, 10x, 20x, 40x, 60x, 80x, 100x, 200x, 300x 400x and 500x coverage reads (assume chloroplast genome size is 160kb).
```
./2_assembly/randomSelection/split_longRead.py
```
Run assembly with different coverage of long-read with Canu/Hinge.

Canu assembly, the final assembly is Epau.contigs.fasta.
```
./2_assembly/longReadOnly/canu/run_canu.sh 
```
Hinge assembly. Due to the data is MinION read, run 1\_convertName.sh to get the correct header that hinge can recoginze.
```
./2_assembly/longReadOnly/hinge/1_convertName.sh 
```
Using the reads with new header to do the assembly. nominal.ini (can be found in hinge install dir) should be in this dir. The final assembly is xx.consensus.fasta
```
./2_assembly/longReadOnly/hinge/2_run_hinge.sh 
```
### 3\_post\_assembly
### mummer
First, used mummer to check the contig alignment. E.reg, the most close known species to E.pau, was used as the reference to compare. 
```
./3_post_assembly/1_same_structure/mummer_plot.sh 
```
 According to the alignment suitation (alignment figure), choose direction.py or mummer\_direction.py to create a single contig which is no duplication and has the same sturcture of cp ref genome. The direction.py script is recommended for assembly which has clear lsc/ir/ssc contigs or a complete contig but the structure is different from cp ref genome. The direction.py changed the direction based on genes, the genes used to ensure direction can be changed in the script. The mummer\_direction.py is basing on mummer alignement result to merge the congits, which may produce extra errors. In general, all long-read assembly used mummer\_plot.sh, because lots of duplicates were observed in the alignment. However, for the short-read only and hybrid assembly, direction.py can be used in the assembly with >= 20x short-read coverage input. Since those assemlies have clear lsc/ssc/ir contigs or has one single contig without overlap (according to mummer result)

run mummer\_plot.sh
```
./3_post_assembly/1_same_structure/mummer_direction.sh 
```
### polish

We polished the genome from last step using Racon, Nanopolish and Racon+Nanopolish, respectively. Racon run 10 iterations (the result keeps change after 10 rounds), whereas Nanopolish run until the result unchanged. The aligner used here was bwa mem.

Run Racon. The path of inputFile (the chloroplast long-reads) and ref (assembly result) should be absolute path. The ref should end with 'fa'. In addition, the read used to assemble should be changed to fastq format. 
```
./3_post_assembly/2_polish/run_racon.sh 
```
Run Nanopolish. Nanopolish needs the index data (link to original Fast5 data, details see [Nanopolish] (https://github.com/jts/nanopolish)) .
```
./3_post_assembly/2_polish/
./3_post_assembly/2_polish/run_nanopolish.sh 
```
For Racon+Nanopolish, run Racon first, and then used the Racon-polish result as input to run Nanopolish.

### assembly\_quality\_control

Finally, we used the 100x short read (randomly selected first, not used in assembly, method see below) to remap to the assembly to assess its quality. Qualimap is used to grep the mapping information.
```
./3_post_assembly/3_assembly_quality_control/1_run_bowtie2.sh 
./3_post_assembly/3_assembly_quality_control/2_run_qualimap.sh 
```

## SHORT-READ ONLY ASSEMBLY

### 1\_pre\_assembly
### QualityControl

Run fastqc to check the quality
```
./1_pre_assembly/1_qualityControl/shortRread/1_qualityCheck/run_fastqc.sh
```
Trim adaptor and low quality region (<30). Read length <50 bp will be removed. 
```
./1_pre_assembly/1_qualityControl/shortRead/2_adapterTrim/run_bbduk.sh 
```
Rerun fastqc to check the data again
```
./1_pre_assembly/1_qualityControl/shortRead/3_qualityCheck/run_fastqc.sh 
```

### cp\_DNA\_extraction

Mapping all trimmed reads to refs (31 known Eucalyptus chloroplast genomes, all double-up, in case read maps to the 'cut-point') to get the chloroplast reads, using bowtie2.
```
./1_pre_assembly/2_cpDNAExtraction/shortRead/1_run_bowtie2.sh 
```
Get chloroplast read from the Bowtie2 output
```
./1_pre_assembly/2_cpDNAExtraction/shortRead/2_getCPRead.py 
```

### 2\_assembly
Short-reads were randomly selected (5x, 8x, 10x, 20x, 40x, 60x, 80x, 100x, 200x, 300x, 400x and 500x, assuming the genome size is 160kb). 100x coverage of short-read was seprated first as the validation data which did not use in assembly.
```
./2_assembly/randomSelection/split_pair_read.sh
```
Unicycler is used to do the short-read only assembly with different coverage. In the default setting, the read is corrected by SPAdes. In this study, we have tried three different read correction: SPAdes, Karect and SPAdes+Karect. However, the different error-correction pipelinne performed very similar.

get Karect-correct read
```
./2_assembly/run_Karect_correction.sh
```
For SPAdes-correct read assembly, using normal read as input:
```
./2_assembly/shortReadOnly/run_shortRead_unicycler.sh 
```
For Karect-correct read assembly, using the Karect-correct read as input:
```
./2_assembly/shortReadOnly/run_shortRead_unicycler_noSPAdes.sh 
```
For Karect-SPAdes-correct read assembly, using the Karect-correct read as input:
```
./2_assembly/shortReadOnly/run_shortRead_unicycler.sh 
```
### 3\_post\_assembly
### mummer
as described above in long-read only assembly part.

Run direction.py
```
./3_post_assembly/1_same_structure/direction.py 
```
Run mummer\_plot.sh
```
./3_post_assembly/1_same_structure/mummer_direction.sh 
```
### polish
we used Pilon to polish the assembly, run until result unchanged.
```
./3_post_assembly/2_polish/run_pilon.sh 
```
### assembly\_quality\_control
As described above in the long-read only assembly part，we used the 100x short-read (randomly selected first, not used in assembly, method see below) to remap to the assembly to assess its quality. Qualimap was used to grep the mapping information.
```
./3_assembly_quality_control/1_run_bowtie2.sh 
./3_assembly_quality_control/2_run_qualimap.sh 
```

## HYBRID ASSEMBLY
The method to get chloroplast long-/short- long read were described above (1\_pre\_assembly). In general, every step in the hybrid assembly is the same as short read only assembly (including the 3\_post\_assembly). Only the 2\_assembly script is different.

We used different combination of coverage (5X, 8X, 10X, 20X, 40X, 60X, 80X, 100X, 200X, 300X, 400X and 500x) of long-/short- read to do the assembly to find out the optimal combination, using Unicycler. Three different types of error-correction, SPAdes/Karect/Karect-SPAdes, were also compared in this hybrid assembly. 

For SPAdes-correct read assembly:
```
./2_assembly/hybrid/run_unicycler.sh 
```
For Karect-correct read assembly:
```
./2_assembly/hybrid/run_unicycler_noSPAdes.sh
```
For Karect-SPAdes-correct read assembly:
```
./2_assembly/hybrid/run_unicycler.sh 
```

## ANNOTATION
We use [GeSeq] (https://chlorobox.mpimp-golm.mpg.de/geseq.html) to do the annotation. The genes with very short exon were manually annotated according to other known Eucalyptus chloroplast genome annotations. 

## PHYLOGENETIC ANALYSIS
We drew the phylogenetic tree using E.pau chloroplast genome (this study generated) and other 31 known Eucalyptus chloroplast genomes. 

First, split the 32 chloroplast genomes by exon and non-exon according to our annotation (for E.pau) and genbank annotation (for other 31 known Eucalyptus):
```
./phylogenetic_analysis/1_gb_gff_fasta.py
```
Totally 312 fragments for each species were generated. We created 312 alignments using clustalo.
```
./phylogenetic_analysis/2_run_clustalo.sh
```
The 312 alignments (fasta format) were checked and fixed manually. After that, combined the 312 fasta files together to create a corrected whole genome alignment (in fasta format), and created a nexus file according to the alignment.
```
./phylogenetic_analysis/3_getFasta.py
```
 Finally, run IQ-TREE with the corrected whole genome alignment (fasta and nexus format) to get the phylogenetic tree. 
 ```
./phylogenetic_analysis/4_run_iqtree.sh
```

The final alignments are in phylogenetic\_analysis/iqtree\_result/InputFiles, the IQ-TREE results are in phylogenetic\_analysis/iqtree\_result/OutputFiles.
