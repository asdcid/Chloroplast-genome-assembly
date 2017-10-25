#!/usr/bin/env python
# encoding: utf-8
"""
getDirection.py

Created by www on  5:30 pm, Sep 13, 2017
"""
import numpy as np 
import sys
import os
from Bio import SeqIO
 
 
def getReference(reference):
    refs    = {}
    for r in SeqIO.parse(reference, 'fasta'):
        refs[str(r.id)] = r.seq
    return refs

def loadFile(inputFile, outputFile, refs):
    o   = open(outputFile, 'w+')
    coords  = {}
    for id in refs:
        coords[id] = np.zeros(len(refs[id]))

    with open(inputFile) as f:
        match = False
        lastRefStart = 0
        lastRefEnd   = 0
        seq = ''

        for line in f:
            if line[0] == '=':
                match = True
                continue
            if match:
                info    = line.split()
                refStart = int(info[0])
                refEnd   = int(info[1])
                length   = int(info[6])
                start   = int(info[3])
                end     = int(info[4])
                id      = info[-1]

                if end > start:
                    strand  = '+'
                else:
                    start, end = end, start
                    strand  = '-'
                if start < 50:
                    start = 1
                #first one
                if not lastRefStart:
                    lastRefStart = refStart
                    lastRefEnd   = refEnd
                    lastStart    = start
                    lastEnd      = end
                    coords[id][start - 1 : end] = 1
                    if strand == '+':
                        seq += str(refs[id][start - 1 : end])
                    else:
                        seq += str(refs[id][start - 1 : end].reverse_complement())
                    
                else:
                    if refStart >= lastRefEnd:
                        lastRefStart = refStart
                        lastRefEnd   = refEnd
                        length = end - start + 1

                        #check duplicate
                        #if sum(coords[id][start - 1 : end]) > length * 0.5:
                        #    continue                 
                        coords[id][start - 1 : end] = 1
                        if strand == '+':
                            seq += str(refs[id][start - 1 : end])
                        else:
                            seq += str(refs[id][start - 1 : end].reverse_complement())
                        
                    elif lastRefStart <= refStart < lastRefEnd:
                        if refEnd < lastRefEnd:
                            continue
                        else:
                            position = refEnd - lastRefEnd
                            if strand == '+':
                                start    = end - position + 1
                            else:
                                end      = start + position - 1

                            #check duplicate
                            #if sum(coords[id][start - 1 : end]) > length * 0.5:
                            #    continue                 
                            coords[id][start - 1 : end] = 1
                            if strand == '+':
                                seq += str(refs[id][start - 1 : end])
                            else:
                                seq += str(refs[id][start - 1 : end].reverse_complement())
                             
                            lastRefStart = refStart
                            lastRefEnd   = refEnd

                    else:
                        print inputFile 
                        print line
                        print 'ERROR'
                        sys.exit()  

        o.write('>1\t%d\n%s\n' % (len(seq), seq))
        o.close()    


def main():
    #the coord file, produce by mummer_plot.sh
    inputFile   = sys.argv[1]
    #the merged contig
    outputFile  = sys.argv[2] 
    #the original assembly, fasta
    reference   = sys.argv[3]
    refs    = getReference(reference) 

    loadFile(inputFile, outputFile, refs)

if __name__ == '__main__':
    main()

