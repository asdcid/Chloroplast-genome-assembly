#!/usr/bin/env python
# encoding: utf-8
"""
getRandomRead.py

Created by www on  1:52 pm, Apr 28, 2017
"""
import numpy as np 
import sys
import os
import gzip
import random
import copy
 
def loadFile(inputFile):
    seqs    = {}
    names   = {}
    with gzip.open(inputFile, 'rb') as f:
        while True: 
            name = f.readline()
            if not name:
                break
            ID   = name
            #the name is different between R1 and R2: @ASD 1:NXXXX   @ASD 2:Nxxxx'
            names[name.split()[0]] = ''
            seq  = f.readline()
            description     = f.readline()
            quality         = f.readline()
            if not name in seqs:
                seqs[name] = ''
            seqs[name.split()[0]] = [ID, seq, description, quality]
    return seqs, names

def separate(names, fraction):
    nRemove = int(len(names) * fraction)
    print nRemove
    smallItems     = random.sample(names, nRemove)
    print len(smallItems)
    largeItems     = copy.deepcopy(names)
    for i in smallItems:
        del largeItems[i]
    #map(largeItems.pop, filter(lambda k : k in smallItems, largeItems))
    print len(largeItems)
    return largeItems, smallItems


def output(seqs, large, small, outputLarge, outputSmall):
    om   = open(outputLarge, 'w+')
    ol   = open(outputSmall, 'w+')
    for name in large:
        om.write('%s%s%s%s' % (seqs[name][0], seqs[name][1], seqs[name][2], seqs[name][3]))
    for name in small:
        ol.write('%s%s%s%s' % (seqs[name][0], seqs[name][1], seqs[name][2], seqs[name][3]))
    om.close()
    ol.close()

def main():
    if 1: 
        if len(sys.argv) != 8 :
            print 'Usage: %s inputR1 inputR2 outputR1.large outputR1.small ouputR2.large outputR2.small fraction for separate (0-1.0)' % os.path.basename(sys.argv[0])
            sys.exit() 
    inputR1   = sys.argv[1]
    inputR2   = sys.argv[2]
    outputLargeR1   = sys.argv[3]
    outputSmallR1   = sys.argv[4]
    outputLargeR2   = sys.argv[5]
    outputSmallR2   = sys.argv[6]
    fraction  = float(sys.argv[7]) 
    seqsR1, namesR1    = loadFile(inputR1)
    print 'finished loading R1'
    seqsR2, namesR2    = loadFile(inputR2)
    print 'finished loading R2'
    large, small   = separate(namesR1, fraction)
    print 'finished separate'
    output(seqsR1, large, small, outputLargeR1, outputSmallR1)
    output(seqsR2, large, small, outputLargeR2, outputSmallR2)


if __name__ == '__main__':
    main()

