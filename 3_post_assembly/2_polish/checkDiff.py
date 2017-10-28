#!/usr/bin/env python
# encoding: utf-8
"""
checkDiff.py

Created by www on  2:25 am, Oct 06, 2017
"""
import numpy as np 
import sys
import os
from Bio import SeqIO
 
 
def main():
    file1   = sys.argv[1]
    file2   = sys.argv[2]
    for r in SeqIO.parse(file1, 'fasta'):
        seq1    = r.seq
    for r in SeqIO.parse(file2, 'fasta'):
        seq2    = r.seq
    if seq1 == seq2:
        sys.stdout.write('same')
    else:
        sys.stdout.write('different')
if __name__ == '__main__':
    main()

