# Chroloplast-genome-assembly

# Requirements
Fastqc;
Bbduk v37.31;
Porechop v0.2.1;
Nanofilt v1.2.0;
Bowtie2 v2.2.6;
Blasr v5.1;
samtools v1.5;
Hinge;
Canu v1.6;
Unicycler v0.3.1;
Mummer v3.23;
Racon;
Nanopolish v0.8.1;
Pilon v1.22;
Qualimap v2.2.1

# Procedure
1_pre_assembly: the input could be the raw long/short reads. Reads will be trimed the low quality region and adapters. Fastqc is used to visualize the read quality. Next, the processed read could be mapped to cp reference genomes (custom) to get the cp reads. Assuming the original reads contain non-cp genome, such at mtDNA or nuclear DNA.

2_assembly: 100x short reads are separated for assembly quality check in the 3_post_assembly. https://github.com/roblanf/splitreads provides scripts for randomly select short/long reads for assembly. The long read only assembl were using Hinge and Canu, whereas short read only and hybrid assembly were using Unicycler

3_post_assembly: the raw assembly was processed to create a single contig, which is no duplication and has same structure with cp reference genome (details see 3_post_assembly/README. Next, the 100x short reads kept from previous step was used to mapped to processed assembly to evaulate the assembly quality by checking the mapping rate, error rate and mismatch/insertion/deletion using Qualimap. 
