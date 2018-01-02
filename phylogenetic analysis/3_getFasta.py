#!/usr/bin/env python
# encoding: utf-8
"""
getFasta.py

Created by www on  2:38 pm, Aug 24, 2017
"""
import numpy as np 
import sys
import os
from Bio import SeqIO


def loadFile(inputFile, dict):

    for record in SeqIO.parse(inputFile, 'fasta'):
        speciesName     = str(record.id)
        if not speciesName:
            continue
        if not speciesName in dict:
            dict[speciesName] = ''
        dict[speciesName] += str(record.seq)
        length = len(record.seq)
    return dict, length
    
def output(outputDir, info, dict, name):
    o_fasta = os.path.join(outputDir, name + '.fasta')
    f   = open(o_fasta, 'w+')
    
    for species in dict:
        f.write('>%s\n%s\n' % (species, dict[species]))
    f.close()
    
    if name != 'nonCoding':
        d   = open(os.path.join(outputDir, name + '.info'), 'w+')
        d.write('#nexus\nbegin sets;\n')

        for i in info:
            gene, start, end = i
            gene = gene.replace('-', '_')
            d.write('\tcharset\t%s\t=\t%d-%d;\n' % (gene, start, end))  
   
        d.write('end;\n')
        d.close()


 
 
def main():
    inputDir    = '../data/fix_fasta/'
    outputDir   = '../data/final_data/'

    exon        = {}
    nonCoding   = {}
    all         = {}
    exon_info       = []
    nonCoding_info  = []
    all_info        = []

    f  = os.walk(inputDir)
    for root, dirs, names in f:
        for name in names:
            inputFile       = os.path.join(root, name)
            ID              = name.split('_')[1].split('.fas')[0]

            all, all_len    = loadFile(inputFile, all)
            if all_info:
                start   = all_info[-1][2] + 1
                end     = start + all_len - 1
                all_info.append([ID, start, end])
            else:
                all_info.append([ID, 1, all_len])

            if 'nocoding' in name:
                nonCoding, nonCoding_len = loadFile(inputFile, nonCoding)
            else:
                exon, exon_len  = loadFile(inputFile, exon)
                if exon_info:
                    start   = exon_info[-1][1]
                    end     = start + exon_len - 1
                    exon_info.append([ID, start, end])
                else:
                    exon_info.append([ID, 1, exon_len])

    output(outputDir, exon_info, exon, 'exon')
    output(outputDir, nonCoding_info, nonCoding, 'nonCoding')
    output(outputDir, all_info, all, 'all')


if __name__ == '__main__':
    main()

