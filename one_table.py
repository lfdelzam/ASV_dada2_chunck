#!/usr/bin/env python3
import os
import argparse
import re

usage = 'python one_table.py [options]'
description = 'This program creates a table combining counts and taxonomy table from dada2 output'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-t',dest='t',help='taxonomy table up to genus or species level', required=True)
parser.add_argument('-s',dest='s',help='species taxonomy file list')
parser.add_argument('-c', dest='c', help='count table', required=True)
parser.add_argument('-o', dest='o', help='prefix output file', required=True)
parser.add_argument('-d', dest='d', help='database - silva,gtdb,pr2', required=True)
parser.add_argument('-m', dest='m', help='metadata file',required=True)
args = parser.parse_args()

### new
def add_species(gen_table,spe_file):
    fname=str(gen_table).split(".")[0]
    with open(fname+"_species_level.tsv", "w") as fout:
        if args.d == "gtdb":
            print("ASV\tKingdom\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecie\tBoot_Kingdom\tBoot_Domain\tBoot_Phylum\tBoot_Class\tBoot_Order\tBoot_Family\tBoot_Genus", file=fout)

        if args.d == "silva":
            print("ASV\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecie\tBoot_Domain\tBoot_Phylum\tBoot_Class\tBoot_Order\tBoot_Family\tBoot_Genus", file=fout)

        if args.d == "pr2":
            print("ASV\tDomain\tSupergroup\tDivision\tSubdivision\tClass\tOrder\tFamily\tGenus\tSpecies\tBoot_Domain\tBoot_Supergroup\tBoot_Division\tBoot_Subdivision\tBoot_Class\tBoot_Order\tBoot_Family\tBoot_Genus\tBoot_Species", file=fout)    

        species={}
        with open(spe_file, "r") as fin:
            for line in fin:
                line=line.rstrip()
                line=line.split("\t")
                asv=line[0]
                g=set( [ i.split()[0] for i in line[1].split(",") ] )
                if len(g) == 1:
                    s=set( [ i.split()[1] for i in line[1].split(",") ] )
                    species[asv]={ list(g)[0]: "/".join(s) }
        #        else:
        #            print(line)

        with open(gen_table, "r") as fing:
            first_line=True
            for line in fing:
                line=line.rstrip()
                line=line.split("\t")
                if first_line:
                    first_line=False
                else:
                    ASV=line[0]
                    if args.d == "gtdb":
                        g=line[7]
                        if ASV in species:
                            ghit=[ k for k in species[ASV].keys() ][0]

                            if g == ghit:
                                print("{}\t{}\t{}".format("\t".join(line[:8]), species[ASV][g], "\t".join(line[9:]) ), file=fout)
                            else:
                                print("{}\t{}\t{}".format("\t".join(line[:8]), "NA", "\t".join(line[9:]) ), file=fout)

                        else:
                            print("{}\t{}\t{}".format("\t".join(line[:8]), "NA", "\t".join(line[9:]) ), file=fout)
                    if args.d == "silva":
                        g=line[6]
                        if ASV in species:
                            ghit=[ k for k in species[ASV].keys() ][0]

                            if g == ghit:
                                print("{}\t{}\t{}".format("\t".join(line[:7]), species[ASV][g], "\t".join(line[8:]) ), file=fout)
                            else:
                                print("{}\t{}\t{}".format("\t".join(line[:7]), "NA", "\t".join(line[8:]) ), file=fout)

                        else:
                            print("{}\t{}\t{}".format("\t".join(line[:7]), "NA", "\t".join(line[8:]) ), file=fout)                            

######

def sample_name(filem):
    first_line=True
    samples_true_names={}
    with open(filem, "r") as fin:
        for line in fin:
            line=line.rstrip()
            if first_line:
                first_line = False

            else:
                line=line.split("\t")
                if line[0] != "":
                    sm=re.sub("-",".",line[0])
                    samples_true_names[sm]=line[1]
    return(samples_true_names)

def get_asv_taxonomy(file, db):
    asv_tax={}
    first_line=True
    with open(file, "r") as fin:
        for line in fin:
            line=line.rstrip()
            if first_line:
                first_line = False

            else:
                line=line.split()
                #ASV;d_Domain;p_Phylum;c_Class;o_Order;f_Family;g_Genus;s_Species
                if db == "gtdb":
                    leves=["d","p","c","o","f","g","s"]
                    texto=[i+"_"+j for i,j in zip(leves,line[2:]) if j != "NA"]
                if db == "silva":
                    leves=["d","p","c","o","f","g","s"]
                    texto=[i+"_"+j for i,j in zip(leves,line[1:]) if j != "NA"]
                if db == "pr2":
                    leves=["D","S","D","d","C","O","F","G","S"]
                    texto=[i+"_"+j for i,j in zip(leves,line[1:]) if j != "NA"]
                asv_tax[line[0]]=";".join(texto)

    return(asv_tax)

if args.s:
    print("Using "+args.s)
    add_species(args.t, args.s)
    fname=str(args.t).split(".")[0]
    fileIn=fname+"_species_level.tsv"
    ASV_tax=get_asv_taxonomy(fileIn, args.d)
else:
    ASV_tax=get_asv_taxonomy(args.t, args.d)

Samples=sample_name(args.m)


unique_ASV_table={}
original_ASVs={}
fl=True
with open(args.c, "r") as fin, open(args.o, "w") as fout:
    for line in fin:
        line=line.rstrip()
        if fl:
            line=line.split("\t")
            header=[Samples[re.sub("_S.*","",s)] for s in line[1:]]
            print("{}\t{}".format(line[0], "\t".join(header)), file=fout)
            fl = False
        else:
            line=line.split()
            tax=ASV_tax[line[0]]
            print("{}\t{}".format(line[0]+";"+tax, "\t".join(line[1:])), file=fout)
            if tax in unique_ASV_table:
                original_ASVs[tax]+=[line[0]]
                for i in range(len(unique_ASV_table[tax])):
                    unique_ASV_table[tax][i]+=int(line[1+i])
            else:
                unique_ASV_table[tax]=[int(i) for i in line[1:] ]
                original_ASVs[tax]=[line[0]]

with open(args.o+"_merged_counts", "w") as foutm, open(args.o+"_ASV_merging_report", "w") as foutr:
    print("{}\t{}".format("ASV", "\t".join(header)), file=foutm)
    counter=1
    for k,v in unique_ASV_table.items():
        print("{}\t{}".format("ASV_"+str(counter)+";"+k, "\t".join(map(str,v))), file=foutm)
        print("{} groups {} ASVs including:\t{}".format("ASV_"+str(counter),len(original_ASVs[k]),",".join(original_ASVs[k]) ), file=foutr)
        counter +=1
