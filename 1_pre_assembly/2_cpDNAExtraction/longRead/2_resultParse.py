#!/usr/bin/env python
# encoding: utf-8
"""
getRead.py
convert blasr output to fasta
Created by www on 12:27 pm, Sep 08, 2017
"""
import numpy as np 
import sys
import os
from Bio import SeqIO
 
 
def main():
    inputFile   = 'result/RB7_C4.fasta.out'
    reads       = 'query/RB7_C4.fasta'
   
    #the readsFormat should be fasta or fastq
    readsFormat = 'fasta'
    outputFile  = 'result/blasr_RB7_C4.fasta'


    o = open(outputFile, 'w+')
    mappedReads   = {}
    with open(inputFile) as f:
        for line in f:
            info    = line.split()
            read    = info[0]
            mappedReads[read] = ''
    for r in SeqIO.parse(reads, readsFormat):
        if str(r.id) in mappedReads:
            o.write('>%s\n%s\n' % (r.id, r.seq))
    o.close()

 
if __name__ == '__main__':
    main()

