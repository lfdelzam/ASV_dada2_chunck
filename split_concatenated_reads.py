#!/usr/bin/env python3

import os
import argparse
import re

usage = 'python split_concatenated_reads.py -i -o'
description = 'This program splits concatenated ASVs into forward and revers reads'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i', dest='i', help='input file ASV_seqs.tsv', required=True)
parser.add_argument('-o', dest='o', help='output file',  required=True)

args = parser.parse_args()

def get_sequences(filein, fileout):
    #ASV     Sequence
    #ASV1    TGAGGAATATTGGACAATGGGCGAGAGCCTGATCCAGCCATGCCGCGTGCAGGAAGACTGCCCTATGGGTTGTAAACTGCTTTTATACAGGAAGAATAAGCCTTACGTGTAAGGTGATGACGGTACTGTAAGAATAAGGACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTCCGAGCGTTATCCGGAATTATTGGGTTTAAAGGGTCCGTAGGCGGATGNNNNNNNNNNCGGAATTATTGGGTTTAAAGGGTCCGTAGGCGGATGATTAAGTCAGGGGTGAAAGTTTGCAGCTCAACTGTAAAATTGCCTTTGATACTGGTCATCTTGAGTTGTATTGAAGTAGGCGGAATATGTAGTGTAGCGGTGAAATGCATAGATATTACATAGAACACCAATTGCGAAGGCAGCTTACTAAGTACTAACTGACGCTGATGGACGAAAGCGTGGGTAGCGAACA
    first_line=True
    with open(filein, "r") as fin, open(fileout, "w") as fout:
        for line in fin:
            line=line.rstrip()
            if first_line:
                first_line=False
            else:
                line=line.split("\t")
                asv=line[0]
                seq_split=line[1].split("NNNNNNNNNN")
                print(">{}.fwd\n{}\n>{}.rev\n{}".format(asv,seq_split[0],asv,seq_split[1]), file=fout)


get_sequences(args.i, args.o)
