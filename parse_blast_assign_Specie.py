#!/usr/bin/env python3

import os
import argparse
import re
import random

usage = 'parse_blast_assign_Specie.py -i -a -o'
description = 'This program parses a blast file of splited concatenated reads, from amplicon 16S, against a reference database and prints out hits for unique species where both reads matches 100% identity and alingment'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i', dest='i', help='input file, blast hits table', required=True)
parser.add_argument('-a', dest='a', help='input file ASV sequences table (from Dada2 output), .tsv', required=True)
parser.add_argument('-g', dest='g', help='path to reference database', required=True)
parser.add_argument('-b', dest='b', help='name of reference database used, options GTDB,silva,custom', required=True)
parser.add_argument('-o', dest='o', help='output prefix',  required=True)


args = parser.parse_args()

random.seed(123)

def get_sequences(file):

    asv_length={}
    with open(file, "r") as fin:
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                asv=line[1:]
            else:
                asv_length[asv]=len(line)
    return asv_length


def clean_name(name):
    name = re.sub(".chromosome.*$", "", name, count=1)
    name = re.sub(".genome.*$", "", name, count=1)
    name = re.sub(".plasmid.*$", "", name, count=1)
    name = re.sub(".DNA.*$", "", name, count=1)
    name = re.sub(".complete.*$", "", name, count=1)
    name = re.sub("sp._", "sp.", name, count=1)
    name=re.sub("contig.*$", "", name, count=1)
    return name

def get_selected_hits(file, asvlen):
# Fields: query acc., subject acc., % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 82 hits found
#1724    NZ_CP016320.1_Vibrio_vulnificus_SE_VV_18_11_contig00043 94.614  427     23      0       1       427     1224    798     0.0     666

    hits = {}
    hits_vulni = {}
    with open(file, "r") as fin:
        for line in fin:
            line = line.rstrip()
    # Fields: query acc., subject acc., % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
            if not line.startswith("#"):
                line = line.split()
                hit = line[1:]
                asv=line[0]
                #16S_rRNA::NC_002505.1_Vibrio_vulnificus_strain_DK_VV_18_01
                id = line[1].split("::")[1] #NC_002505.1_Vibrio_vulnificus_strain_DK_VV_18_01
                id = id.split("_")
                id = "_".join(id[2:]) #Vibrio_vulnificus_strain_DK_VV_18_01
                id=clean_name(id)
                bs=float(line[11])
                pi=float(line[2])
                al=(float(line[3])*100/asvlen[asv])

                if al >=99.99 and pi == 100 :
                    if asv in hits:
                        hits[asv].append((id, bs, hit))
                    else:
                        hits[asv]=[(id, bs, hit)]

    return hits

def get_both_matches(dict,fileo):
    with open(fileo+"_uniq_species", "w") as foutv, open(fileo+"_strain_level", "w") as fouts:
        real_asvs=set([ a.split(".")[0] for a in dict.keys()])
        for r in real_asvs:
            if  r+".fwd" in dict and r+".rev" in dict:
                fullnameF=set([ l[0] for l in dict[r+".fwd"] ])
                fullnameR=set([ l[0] for l in dict[r+".rev"]  ])
                fullname=fullnameF&fullnameR #both fwd and rev match
                if len(fullname) > 0:
                    commun_sps=set( [ " ".join(l.split("_")[:2]) for l in fullname ]  )
                    if len(commun_sps) == 1 : #only if hits corresponds to a specie
                        print("{}\t{}".format(r, ",".join(commun_sps)), file=foutv )
                        print("{}\t{}".format(r, ",".join(fullname)), file=fouts )




def get_selected_hits_gtdb(file, asvlen):
# ASV0.fwd        RS_GCF_001399455.2~NZ_LLEI02000013.1    100.000 229     0       0       1       229     366     594     4.29e-115       414

    hits = {}
    hits_vulni = {}
    with open(file, "r") as fin:
        for line in fin:
            line = line.rstrip()
    # Fields: query acc., subject acc., % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
            if not line.startswith("#"):
                line = line.split()
                hit = line[1:]
                asv=line[0]
                id = line[1]
                bs=float(line[11])
                pi=float(line[2])
                al=(100*float(line[3])/asvlen[asv])

                if al >=99.99 and pi == 100 :
                    if asv in hits:
                        hits[asv].append((id, bs, hit))
                    else:
                        hits[asv]=[(id, bs, hit)]

    return hits

def names_from_gtdb(filein):

    long_names={}
    with open(filein, "r") as fin:
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
        #>RS_GCF_002903645.1~NZ_PDGP01000109.1 Vibrio vulnificus | Silva >HG530070.1.1349 Trueperella pyogenes
                line=line.split()
                long_names[line[0][1:]]=" ".join(line[1:]) #new
    return long_names


def get_both_matches_gtdb(dict,fileo, name):

    with open(fileo+"_uniq_species", "w") as foutv, open(fileo+"_strain_level", "w") as fouts :
        real_asvs=set([ a.split(".")[0] for a in dict.keys()])
        for r in real_asvs:
            if  r+".fwd" in dict and r+".rev" in dict:
                speciesF=set( [ l[0] for l in dict[r+".fwd"] ] )
                speciesR=set( [ l[0] for l in dict[r+".rev"] ] )
                commun_sps=speciesF&speciesR
                if len(commun_sps)>0:
                    commun=set([ " ".join(name[i].split(" ")[:2]) for i in commun_sps ])
                    if len(commun) == 1 :
                         print("{}\t{}".format(r, ",".join(commun)), file=foutv )
                         Lcommun=set([ name[i] for i in commun_sps ])
                         print("{}\t{}".format(r, ",".join(Lcommun)), file=fouts )


asv_len=get_sequences(args.a)

if args.b == "GTDB" or args.b == "silva":
    print("Using "+args.g)
    selected_hits=get_selected_hits_gtdb(args.i, asv_len)
    long_species_names=names_from_gtdb(args.g)
    get_both_matches_gtdb(selected_hits, args.o,long_species_names)

elif args.b == "custom":
    print("Using custom database")
    selected_hits=get_selected_hits(args.i, asv_len)
    get_both_matches(selected_hits, args.o)
