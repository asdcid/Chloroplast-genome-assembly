#!/usr/bin/env python
# encoding: utf-8
"""
split_single.py

Created by www on  2:00 pm, Jun 07, 2017
"""
import numpy as np
import sys
import os
import gzip
import random


def loadFile(inputFile, outputDir, coverage, lengthCutoff):
    data = []
    bp = 0
    percent = 35
    for n in xrange(len(coverage)):
        with gzip.open(inputFile, 'rb') as f:
            while 1:
                name    = f.readline()
                if not name:
                    print 'ERROR, the end of file'
                    break
                seq     = f.readline()
                #info    = f.readline()
                #qual    = f.readline()
                #if len(seq) >= 10000:
                if random.randrange(1, 101) <= percent:
                    data.append([name, seq])
                    #data.append([name, seq, info, qual])
                    bp += len(seq.strip())
                    if bp >= lengthCutoff[n]:
                        print bp, coverage[n]
                        output(data, outputDir, coverage[n])
                        data = []
                        bp = 0
                        break

def output(data, outputDir, outputCoverage):
    outputFile = os.path.join(outputDir, str(outputCoverage)+'coverage.fasta')
    #outputFile = os.path.join(outputDir, str(outputCoverage)+'coverage.fastq')
    o   = open(outputFile, 'w+')
    for read in data:
        name, seq = read
        #name, seq, info, qual = read
        o.write(name)
        o.write(seq)
        #o.write(info)
        #o.write(qual)
    o.close()




def getCoverageLength(coverage, refLength):
    lengthCutoff    = []
    for i in coverage:
        lengthCutoff.append(i * refLength)
    return lengthCutoff


def main():
    inputFile   = 'longRead.fasta.gz'
    outputDir   = 'longRead/'
    coverage    = [5, 8, 10, 20, 40, 60, 80, 100, 200, 300, 400, 500]
    refLength   = 160000
    lengthCutoff=getCoverageLength(coverage, refLength)
    loadFile(inputFile, outputDir, coverage, lengthCutoff)

if __name__ == '__main__':
    main()
