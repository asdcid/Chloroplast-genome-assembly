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
    #inputFile   = 'result/RB7_C4.fasta.out'
    #reference   = 'query/RB7_C4.fasta'
    #outputFile  = 'result/blasr_RB7_C4.fasta'

    #blasr output
    inputFile   = sys.argv[1]
    #blasr inputFile
    reference   = sys.argv[2]
    #output parse result
    outputFile  = sys.argv[3]

    o = open(outputFile, 'w+')
    reads   = {}
    with open(inputFile) as f:
        for line in f:
            info    = line.split()
            read    = info[0]
            reads[read] = ''
    for r in SeqIO.parse(reference, 'fasta'):
        if str(r.id) in reads:
            o.write('>%s\n%s\n' % (r.id, r.seq))
    o.close()

 
if __name__ == '__main__':
    main()

