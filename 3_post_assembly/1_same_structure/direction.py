#!/usr/bin/env python
# encoding: utf-8
"""
direction.py

Created by www on  6:14 pm, Sep 04, 2017
"""
import numpy as np 
import sys
import os
from Bio import SeqIO
from Bio import Seq 


lsc_psbI = 'ATGCTTACTCTCAAACTCTTTGTTTACACAGTAGTGATATTTTTTGTTTCTCTCTTCATCTTTGGATTCCTATCTAATGATCCAGGACGTAATCCTGGACGTGAAGAATAA'
ssc_ccsA = 'ATGATATTTTCCACCTTAGAACATATATTAACTCATATATCCTTTTCGATCATTTCAATAGTAATTACTATTTTTTTGATAAGCTTATCGGTCGATGAAATCGTAGGACTATATGATTCGTCAGAAAAAGGCATGATAGCTACTTTTTTCTGTATAACAGGATTATTAGTCACTCGTTGGATTTATTCGGGACATTTCCCGTTAAGTGATTTATATGAATCATTAATTTTTCTTTCATGGAGTTTCTCCGTTATTCATATGGTTCCGTATTTTAAAAAACATAAAAATTATTTAAGCACAATAACCGCGCCAAGTACTATTTTTACCCAAGGCTTTGCTACTTCGGGTCTTTTAACTGAAATGCATCAATCCGAAATAGTAGTACCTGCTCTCCAATCCCAATGGTTAATGATGCATGTAAGTATGATGATATTGGGCTATGCAGCTCTTTTATGTGGATCATTATTATCAGTAGCTCTCTTAGTCATTACATTGCGAAAAGCCATAAGGGTTTTTAGTAAAAAAAACAATTTTTTAAATGAGTCATTTTCCTTTGTCGAGATCCAATACATGAATGAAAGAAGCAATGTTTTACTAACCACTTCTTTTTGTTCTTCTCGAAATTATTACAGGGCTCAACTGATTCAACAACTGGATCAGTGGAGTTATCGTATTATTAGTCTAGGGTTTATCTTTTTAACCATAGGTATTCTTTCGGGAGCAGTATGGGCTAATGAGGCATGGGGGTCATATTGGAATTGGGACCCAAAGGAAACTTGGGCATTTATTACTTGGACCCTATTTGCGATTTATTTACATACTCGAACAAATAAAAATTTGGAAAGTTTAAATTGCGCAATTGTGGCTTCTATAGGCTTTCTTATAATTTGGATATGCTATTTTGGGGTTAATTTATTAGGAATAGGATTACATAGTTATGGTTCATTTAATTTACATTAA'
ir_rpl23 = 'ATGGATGGAATCAAATATGCAGTATTTACAGACAAAAGTATTCGGTTATTGGGGAAAAATCAATATACTTTTAATGTCGAATCAGGATCAACTAGGACAGAAATAAAGCATTGGGTCGAACTCTTCTTTGGTGTCAAGGTAAAAGCTATGAATAGTCATCGACTCCCGGGAAAGGGTAGAAGAATGGGACCTATTCTGGGACATACAATGCATTACAGACGTATGATCATTACGCTTCAACCGGGTTATTCTATTCCACCTCTTAGAAAGAAAAGAACTTAA'

lsc_trnH = 'GCGGATGTAGCCAAGTGGATCAAGGCAGTGGATTGTGAATCCACCATGCGCGGGTTCAATTCCCGTCGTTCGCCC'
lsc_psbZ = 'ATGACTATTGCTTTTCAATTGGCTGTTTTTGCATTAATTGCTACTTCATTAATCTTACTGATTAGTGTACCCGTTGTATTTGCTTCTCCTGACGGTTGGTCGAGTAACAAAAATGTTGTATTTTCTGGTACATCATTATGGATTGGATTAGTCTTTCTGGTGGGTATCCTTAATTCTCTCATCTCTTGA'


def loadFile(inputFile):
    data = {}
    
    for r in SeqIO.parse(inputFile, 'fasta'):
        data[r.id] = r.seq
    #print len(data)
    if len(data) == 1:
        seq     = r.seq
        mark = True
    elif len(data) == 3:
        lsc = data['1']
        ir  = data['2']
        ssc = data['3']
        if lsc_psbI in lsc:
            lsc_direction = '+'
        elif Seq.Seq(lsc_psbI).reverse_complement() in lsc:
            lsc_direction = '-'
            lsc = lsc.reverse_complement()
        if ssc_ccsA in ssc:
            ssc_direction = '+'
        elif Seq.Seq(ssc_ccsA).reverse_complement() in ssc:
            ssc_direction = '-'
            ssc = ssc.reverse_complement()
        if ir_rpl23 in ir:
            ir_direction = '-'
            ir = ir.reverse_complement()
        elif Seq.Seq(ir_rpl23).reverse_complement() in ir:
            ir_direction = '+'

        #seq     = str(lsc) + str(ir) + str(ssc)
        seq     = str(lsc) + str(ir) + str(ssc) + str(ir.reverse_complement())
        mark    = False
        #seq     = Seq.Seq(str(data['1']) + str(data['2']) + str(data['3']) + str(data['2'].reverse_complement()))
    elif len(data) == 4:
        ir  = data['2']
        ssc = data['3']
        print inputFile
        print '444444444444'
        if lsc_psbZ in str(data['1']):
            pass
        elif str(Seq.Seq(lsc_psbZ).reverse_complement()) in str(data['1']):
            data['1'] = data['1'].reverse_complement()
        if lsc_trnH in str(data['4']):
            data['4'] = data['4'].reverse_complement()
        elif str(Seq.Seq(lsc_trnH).reverse_complement()) in str(data['4']):
            pass
        else:
            print 'ERROR, not lsc_trnH in >4'
        if ssc_ccsA in ssc:
            ssc_direction = '+'
        elif Seq.Seq(ssc_ccsA).reverse_complement() in ssc:
            ssc_direction = '-'
            ssc = ssc.reverse_complement()
        if ir_rpl23 in ir:
            ir_direction = '-'
            ir = ir.reverse_complement()
        elif Seq.Seq(ir_rpl23).reverse_complement() in ir:
            ir_direction = '+'

        mark    = False
        #seq     = Seq.Seq(str(data['4']) + str(data['1']) + str(data['2']) + str(data['3']))
        #seq     = Seq.Seq(str(data['4']) + str(data['1']) + str(data['2']) + str(data['3']) + str(data['2'].reverse_complement()))
        seq     = Seq.Seq(str(data['4']) + str(data['1']) + str(ir) + str(ssc))

    return seq, mark


def direction(seq, outputFile):
    position = str(seq).find(str(Seq.Seq(lsc_trnH).reverse_complement()))
    if position != -1:
        if position > 400:
            seq      = str(seq)
            seq = Seq.Seq(seq[position : ] + seq[ : position])

    
    reverse_seq     = str(seq.reverse_complement())
    begin           = 26000
    ir_frag         = reverse_seq[ : begin]
    sequence        = str(seq)
    match   = sequence.find(ir_frag) 
    n  = 0
    while 1:
        if sequence[n + match + begin] == reverse_seq[begin + n]:
            n += 1
        else:
            break
    irLength  = begin + n
    lsc = seq[ : match]
    ira = seq[match : match + begin + n]
    ssc = seq[match + begin + n : -irLength]
    irb = seq[-irLength : ]
    
    #debug
    if ira != irb.reverse_complement():
        print len(lsc), len(ira), len(ssc), len(irb)
        print 'ERROR, different ir'
        sys.exit()

    if lsc_psbI in lsc:
        lsc_direction = '+'
    elif Seq.Seq(lsc_psbI).reverse_complement() in lsc:
        lsc_direction = '-'

    if ssc_ccsA in ssc:
        ssc_direction = '+'
    elif Seq.Seq(ssc_ccsA).reverse_complement() in ssc:
        ssc_direction = '-'
   
    if ir_rpl23 in ira:
        ira_direction = '-'
    elif Seq.Seq(ir_rpl23).reverse_complement() in ira:
        ira_direction = '+' 


    if ir_rpl23 in irb:
        irb_direction = '+'
    elif Seq.Seq(ir_rpl23).reverse_complement() in irb:
        irb_direction = '-' 
    

    #print lsc_direction, ira_direction, ssc_direction, irb_direction

    #debug
    if ira_direction != irb_direction:
        print 'ERROR, the same direction for ir'
        sys.exit()
    
    #output
    with open(outputFile, 'w+') as o:
        final_seq = ''
        o.write('>1\t%d\n' % len(seq))

        if lsc_direction == '+':
            final_seq += str(lsc)
        else:
            final_seq += str(lsc.reverse_complement())

        if ira_direction == '+':
            final_seq += str(ira)
        else:
            final_seq += str(ira.reverse_complement())
 
        if ssc_direction == '+':
            final_seq += str(ssc)
        else:
            final_seq += str(ssc.reverse_complement())

        if irb_direction == '+':
            final_seq += str(irb)
        else:
            final_seq += str(irb.reverse_complement())
        
        o.write('%s\n' % final_seq)


def main():
    #the original assembly, fasta
    inputFile   = sys.argv[1]
    #the merged contig
    outputFile  = sys.argv[2]
    #print inputFile
    seq, mark     = loadFile(inputFile)
    if mark:
        direction(seq, outputFile) 
    else:
        o = open(outputFile, 'w+')
        o.write('>1\t%d\n%s\n' % (len(seq), seq))
        o.close()

if __name__ == '__main__':
    main()

