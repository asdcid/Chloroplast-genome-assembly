#!/usr/bin/env python
# encoding: utf-8
"""
checkOverlap.py
remove one rpl22, orf113
Created by www on 12:55 pm, Aug 21, 2017
"""
import numpy as np 
import sys
import os
from Bio  import SeqIO
 
def loadFile(inputFile):
    genes = [] 
    with open(inputFile) as f:
        for line in f:
            line    = line.strip()
            if not line:
                continue
            if line[0] == '#':
                continue
            info    = line.split()
            type    = info[2]
            start   = int(info[3])
            end     = int(info[4])
            strand  = info[6]
            feature = info[8]
            if type != 'exon':
                continue
            geneID  = feature.split('Parent=')[1].split(';')[0]
            if 'number' in feature:
                number = feature.split('number=')[1]
            else:
                number = ''
            genes.append([start, end, strand, geneID, number])
    return genes


def loadGb(inputFile):
    genes   = []
    gb      = SeqIO.read(inputFile, 'genbank')
    for feature in gb.features:
        if feature.type != 'CDS' and feature.type != 'tRNA' and feature.type != 'rRNA' and feature.type != 'gene':
            continue
        if feature.type == 'gene' and not 'pseudo' in feature.qualifiers:
            continue
        if 'gene' in feature.qualifiers:
            geneID      = feature.qualifiers['gene'][0]
            if geneID == 'ORF113':
                continue
        else:
            geneID      = feature.qualifiers['product'][0]
        strand = feature.strand
        if strand == -1:
            strand = '-'
        elif strand == 1:
            strand = '+'
        #remove one rpl22
        else:
            continue 
            
 
        location    = feature.location.parts
        if len(location) > 1:
            number = 0
            for i in location:
                start   = i.nofuzzy_start + 1
                end     = i.nofuzzy_end
                number += 1
                genes.append([start, end, strand, geneID, number])
        else:
            number  = ''
            start   = location[0].nofuzzy_start + 1
            end     = location[0].nofuzzy_end
            genes.append([start, end, strand, geneID, number])

    return genes, str(gb.seq)

def output(outputDir, genes, refs):
    n = 0
    for i in xrange(len(genes['E.pauciflora'])):
        n += 1
        gene    = genes['E.pauciflora'][i][3]
        geneID  = '%03d_%s' % (n, gene)
        outputFile  = os.path.join(outputDir, geneID + '.fasta')
        o           = open(outputFile, 'w+')
        for species in genes:
            start     = genes[species][i][0]
            end       = genes[species][i][1]
            seq       = refs[species][start - 1 : end]
            o.write('>%s\n%s\n' % (species, seq))
        o.close()




def orderGene(genes):
    for species in genes:
        genes[species].sort(key = lambda x : x[0])
        for i in xrange(len(genes[species]) - 1):
            currentGene  = genes[species][i]
            nextGene     = genes[species][i + 1]
            currentStart    = currentGene[0]
            currentEnd      = currentGene[1]
            nextStart       = nextGene[0]
            nextEnd         = nextGene[1]
        
            if currentEnd > nextEnd:
                print 'Outlier', currentGene
                print 'Outlier', nextGene
            if currentEnd > nextStart:
                #print currentGene
                #print nextGene
                genes[species][i + 1][0] = currentEnd + 1
    for species in genes:
        new_geneSet = []
        for i in xrange(len(genes[species]) - 1):
            currentGene  = genes[species][i]
            nextGene     = genes[species][i + 1]
            currentStart    = currentGene[0]
            currentEnd      = currentGene[1]
            nextStart       = nextGene[0]
            nextEnd         = nextGene[1]

            new_geneSet.append(currentGene)
            new_geneSet.append([currentEnd + 1, nextStart - 1, 'none', 'nocoding', ''])
        genes[species] = new_geneSet

    return genes


def getRef(EpauRef):
    seq = ''
    for r in SeqIO.parse(EpauRef, 'fasta'):
        seq = str(r.seq)
    return seq


def grandis_add_missGene(data):
    
    data.append([9718, 9740, '+', 'tRNA-Gly', 1])
    data.append([10485, 10533, '+', 'tRNA-Gly', 2])

    data.append([5075, 5268, '-', 'rps16', 2])
    data.append([6145, 6184, '-', 'rps16', 1])

    data.append([17880, 22057, '-', 'rpoC2', ''])

    data.append([22223, 23838, '-', 'rpoc1', 2])
    data.append([24572, 25024, '-', 'rpoc1', 1])
    
    data.append([53530, 54382, '-', 'ndhK', ''])

    data.append([81097, 81104, '+', 'petD', 1])
    data.append([81875, 82349, '+', 'petD', 2])

    data.append([61165, 62659, '+', 'accD', ''])

    data.append([129216, 134841, '-', 'ycf1', ''])

    data.append([88932, 89365, '-', 'rpl2', 2])
    data.append([90030, 90420, '-', 'rpl2', 1])

    data.append([125058, 125598, '-', 'ndhA', 2])
    data.append([126659, 127209, '-', 'ndhA', 1])

    data.append([158590, 158980, '+', 'rpl2', 1])
    data.append([159645, 160078, '+', 'rpl2', 2])


    return data 

def main():
    genes   = {}
    refs    = {}

    EpauInput   = '../data/Epau/Epau.gff3'
    inputDir    = '../data/other_cp/'
    outputDir   = '../data/fasta'
    EpauRef     = '../../unicycler/ref/Epau.fasta'

    genes['E.pauciflora']   = loadFile(EpauInput) 
    refs['E.pauciflora']    = getRef(EpauRef)

    f   = os.walk(inputDir)
    for root, dirs, names in f:
        for name in names:
            if not '.gb' in name:
                continue
            species         = name.split('.gb')[0]
            inputFile       = os.path.join(root, name)
            genes[species], refs[species]  = loadGb(inputFile)
            if species == 'E.grandis':
                genes[species] = grandis_add_missGene(genes[species])
    genes   = orderGene(genes) 
    #for i in genes:
    #    for j in xrange(len(genes[i])):
    #        pass
    #        print genes['E.pauciflora'][j], genes['E.grandis'][j]
    output(outputDir, genes, refs)   


 
if __name__ == '__main__':
    main()

