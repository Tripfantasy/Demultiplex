#!/usr/bin/env python3.10
import numpy as np
import matplotlib.pyplot as plt
import argparse
import bioinfo
import gzip
import itertools as it
# Interactive environment command: srun --account=bgmp --partition=bgmp --nodes=1 --ntasks-per-node=1 --time=2:00:00 --cpus-per-task=1 --pty bash
# Number of reads        | File Name
# R1 reads = 363,246,735 | 1294_S1_L008_R1_001.fastq.gz
# R2 reads = 363,246,735 | 1294_S1_L008_R2_001.fastq.gz
# R3 reads = 363,246,735 | 1294_S1_L008_R3_001.fastq.gz
# R4 reads = 363,246,735 | 1294_S1_L008_R4_001.fastq.gz
# Absolute path is /projects/bgmp/shared/2017_sequencing

#Argparse to get distribution for different read file. 
def get_args():
    parser = argparse.ArgumentParser(description = "Program to select filename and quantity of reads")
    parser.add_argument("-f","--filename", nargs = 4, help ="Input File name", type = str, required = True)
    parser.add_argument("-ind","--indfile",help = "input file containing indexes",type = str, required = True)
    return parser.parse_args()
args = get_args()

index_set=set()
with open(args.indfile, 'rt') as indexfile: #stores indexes from file in a set
    for line in indexfile:
        index = line.strip().split()[4]
        index_set.add(index)
#print("Index set:", index_set)
    #print(index_set, sep="\n")

indexes = {}
#Open the 48 barcode files
for index in index_set:
    knownR1 = open(index + '_R1.fq','wt')
    knownR2 = open(index + '_R2.fq','wt')
    indexes[index] = [knownR1, knownR2]
#print(indexes)

#Open 4 files R1 and R2 hopped and unknown
hopped_list = [open('Hopped_R1.fq','w'), open('Hopped_R2.fq','w')]
unknown_list = [open('Unknown_R1.fq','w'), open('Unknown_R2.fq','w')]

#Find all possible products using itertools
all_pairs = {}
products = it.product(index_set, index_set)
#print(products)
for pair in products:
    all_pairs[pair] = 0 #tuple with forward and reverse complimentary pairs
#print(all_pairs)

#Revcomp function defined
bases = {'A': 'T','T': 'A','C': 'G','G': 'C', 'N': 'N'}
def revcomp(seq: str) -> str:
    rev_com = ""
    for base in seq:
        rev_com += bases[base]
    return(rev_com[::-1])

#bioinfo.convert_phred safe check
# def bioinfo.convert_phred(letter: str) -> int:
#      """Converts a single character into a phred score"""
#     qual_score = (ord(letter) -33)
#     return qual_score

#Open input files R1,R2,R3,R4
R1 = gzip.open(args.filename[0],"rt")
I1 = gzip.open(args.filename[1],"rt")
I2 = gzip.open(args.filename[2],"rt")
R2 = gzip.open(args.filename[3],"rt")

#Counter definitions
count_hopped = {"Hopped_count:": 0}
count_unknown = {"Unknown_count:": 0}
count_matched = {"Match_count:": 0}

while True: #Accessing records in input files per-line.
    header_R1 = R1.readline().strip()
    header_I1 = I1.readline().strip()
    header_I2 = I2.readline().strip()
    header_R2 = R2.readline().strip()
    if header_R1 == "":
        break
    seq_R1 = R1.readline().strip()
    seq_I1 = I1.readline().strip()
    seq_I2 = I2.readline().strip()
    seq_R2 = R2.readline().strip()

    plus_R1 = R1.readline().strip()
    plus_I1 = I1.readline().strip()
    plus_I2 = I2.readline().strip()
    plus_R2 = R2.readline().strip()

    qs_R1 = R1.readline().strip()
    qs_I1 = I1.readline().strip()
    qs_I2 = I2.readline().strip()
    qs_R2 = R2.readline().strip()



    #Adds index + revcomp(index) to header line for records.
    ind_header_R1 = header_R1 + ' ' + seq_I1 + ' ' + '+' + revcomp(seq_I2)
    ind_header_R2 = header_R2 + ' ' + seq_I1 + ' ' + '+' + revcomp(seq_I2)

    #Conditionals for demux will become more selective as it goes down. Calling the bool "Funnel" to guide this.
    funnel = False
    quality_cutoff = 30 

    if seq_I1 not in index_set:  #Is index 1 in known indexes?
        unknown_list[0].write(ind_header_R1 + "\n" + seq_R1 + "\n" + plus_R1 + "\n" + qs_R1 + "\n")
        unknown_list[1].write(ind_header_R2 + "\n" + seq_R2 + "\n" + plus_R2 + "\n" + qs_R2 + "\n")
        count_unknown["Unknown_count:"] += 1
    elif revcomp(seq_I2) not in index_set: #Is the reverse comp of index 2 in known indexes?
        unknown_list[0].write(ind_header_R1 + "\n" + seq_R1 + "\n" + plus_R1 + "\n" + qs_R1 + "\n")
        unknown_list[1].write(ind_header_R2 + "\n" + seq_R2 + "\n" + plus_R2 + "\n" + qs_R2 + "\n")
        count_unknown["Unknown_count:"] += 1
    elif seq_I1 == revcomp(seq_I2): #If index 1 and revcomp(index 2) match, check quality score.
        for score in qs_I1:
            if bioinfo.convert_phred(score) < quality_cutoff:
                unknown_list[0].write(ind_header_R1 + "\n" + seq_R1 + "\n" + plus_R1 + "\n" + qs_R1 + "\n")
                unknown_list[1].write(ind_header_R2 + "\n" + seq_R2 + "\n" + plus_R2 + "\n" + qs_R2 + "\n")
                count_unknown["Unknown_count:"] += 1
                funnel = True
                break
        if funnel == False:
            for score in qs_I2:
                if bioinfo.convert_phred(score) < quality_cutoff:
                    unknown_list[0].write(ind_header_R1 + "\n" + seq_R1 + "\n" + plus_R1 + "\n" + qs_R1 + "\n")
                    unknown_list[1].write(ind_header_R2 + "\n" + seq_R2 + "\n" + plus_R2 + "\n" + qs_R2 + "\n")
                    count_unknown["Unknown_count:"] += 1
                    funnel = True
                    break
            if funnel == False: #If record passes quality score check and matches, add to matched.
                indexes[seq_I1][0].write(ind_header_R1 + "\n" + seq_R1 + "\n" + plus_R1 + "\n" + qs_R1 + "\n")
                indexes[seq_I1][1].write(ind_header_R2 + "\n" + seq_R2 + "\n" + plus_R2 + "\n" + qs_R2 + "\n")
                count_matched["Match_count:"] += 1
                all_pairs[(seq_I1, revcomp(seq_I2))] += 1 
    else: #If the record is not matched, or known, check for quality score.
        for qualscore in qs_I1:
            if bioinfo.convert_phred(qualscore) < quality_cutoff:
                unknown_list[0].write(ind_header_R1 + "\n" + seq_R1 + "\n" + plus_R1 + "\n" + qs_R1 + "\n")
                unknown_list[1].write(ind_header_R2 + "\n" + seq_R2 + "\n" + plus_R2 + "\n" + qs_R2 + "\n")
                count_unknown["Unknown_count:"] += 1
                funnel = True
                break
        if funnel == False: 
            for qualscore in qs_I2:
                if bioinfo.convert_phred(qualscore) < quality_cutoff:
                    unknown_list[0].write(ind_header_R1 + "\n" + seq_R1 + "\n" + plus_R1 + "\n" + qs_R1 + "\n")
                    unknown_list[1].write(ind_header_R2 + "\n" + seq_R2 + "\n" + plus_R2 + "\n" + qs_R2 + "\n")
                    count_unknown["Unknown_count:"] += 1
                    funnel = True
                    break
            if funnel == False: #If the quality score is not less than quality cutoff, add to hopped.
                hopped_list[0].write(ind_header_R1 + "\n" + seq_R1 + "\n" + plus_R1 + "\n" + qs_R1 + "\n")
                hopped_list[1].write(ind_header_R2 + "\n" + seq_R2 + "\n" + plus_R2 + "\n" + qs_R2 + "\n")
                count_hopped["Hopped_count:"] += 1 
                all_pairs[(seq_I1, revcomp(seq_I2))] += 1 
                funnel = True

#Close files: 2 unknown 2 hopped
unknown_list[0].close()
unknown_list[1].close()
hopped_list[0].close()
hopped_list[1].close()

#Close 48 index files
for index in indexes:
    indexes[index][0].close()
    indexes[index][1].close()


#Statistical analysis
matched_val = sum(count_matched.values())
print("Total matched indexes", matched_val)

hopped_val = sum(count_hopped.values())
print("Total hopped indexes", hopped_val)

unknown_val = sum(count_unknown.values())
print("Total unknown indexes", unknown_val)

total = matched_val + hopped_val + unknown_val
print("Total records", total)



with open ("stats.md", "wt") as stats:
    stats.write("\~SUMMARY~/")
    stats.write("\n")
    stats.write("Dual-matched: " + str(matched_val)+"\n")
    stats.write("Hopped: " + str(hopped_val)+"\n")
    stats.write("Unknown: " + str(unknown_val)+"\n")
    stats.write("Total: " + str(total)+"\n")
    stats.write("\n")
    stats.write("\~PERCENTAGES~/")
    stats.write("\n")
    stats.write("Dual-matched %: " + str((matched_val/total)*100) + "%" + "\n")
    stats.write("Hopped %: " + str((hopped_val/total)*100) + "%" + "\n")
    stats.write("Unknown %: " + str((unknown_val/total)*100) + "%" + "\n")
    stats.write("\n")
    stats.write("\~All Possible Index-Pairs~/")

    for pair in all_pairs: #Refers to itertools product pairs. 
        stats.write(pair[0]+"-"+pair[1]+": "+str(all_pairs[pair])+"\n")
    


# matched_percent = (matched_val / total) * 100 
# hopped_percent = (hopped_val / total) * 100
# unknown_percent = (unknown_val / total) * 100


                



    

        
    
                