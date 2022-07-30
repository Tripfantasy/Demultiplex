#!/usr/bin/env python3.10
import gzip
import argparse
# This script points to target fastq file's sequence lines,
# reverse compliments them, and adds each revcomp to a list.
def get_args():
    parser = argparse.ArgumentParser(description = "Program to select filename")
    parser.add_argument("-f","--filename", help ="File name", type = str, required = True)
    return parser.parse_args()

def reversecomp(dna):
#This function will return reverse complemented DNA strings. 
    reversedna = dna[::-1]
    complement = reversedna.replace('A','t').replace('T','a').replace('G','c').replace('C','g').upper()
    return complement

args = get_args()
with gzip.open(args.filename,'rt') as myfile:
    i = 0 
    complements = []
    for line in (myfile):
        i += 1
        if i%4 == 2: #Qual lines .fq
            line = line.strip("\n")
            #print(line)
            line = reversecomp(line)
            complements.append(line)
            print(line)
print(complements)
