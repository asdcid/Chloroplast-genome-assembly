## Random reads
1. 100x illumina reads for assembly quality assessment were separte first using scripts in https://github.com/roblanf/splitreads
2. The random reads (can be set for different coverage) were generated using scripts in https://github.com/roblanf/splitreads (optional). Also, the scripts I used to random select different coverage of long reads and short reads can be found in randomSelection folder.

## long read assembly
1. Long read only assembly were using Hinge and Canu
2. Sometimes (in low coverage) Hinge failed to assemble a genome, change cut\_off in nominal.ini from 300 to 10 could help
3. The fasta2DB in Hinge can recoginze pacbio read header only, change the header (if not use pacbio data) using hinge/1\_convertName.sh

## short read and Hybrid assembly
1. Assembly is polished by Pilon in Unicycler by default setting
2. Unicycler is supported by python3
3. The hybrid assembly and short read only assembly was using Unicycler. The input read is corrected by SPAdes by default setting in Unicycler. It can be turned off using --no\_correct. An example script of using Karect correction pipeline is provided in run\_Karect\_correction.sh

