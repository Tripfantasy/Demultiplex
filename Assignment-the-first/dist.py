#!/usr/bin/env python3.10
import numpy as np
import matplotlib.pyplot as plt
import argparse
import bioinfo
import gzip
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
    parser.add_argument("-f","--filename", help ="File name", type = str, required = False)
    parser.add_argument("-r","--reads", help = "Number of reads in fastq", type = int, required = True)
    parser.add_argument("-l","--length", help = "Length of fastq reads", type = int, required = True)
    return parser.parse_args()


#all_qscores_sum = np.empty([args.length, args.reads])

args = get_args()
with gzip.open(args.filename,'r') as myfile:
    all_qscores_sum = np.empty(args.length)
    mean = np.empty(args.length)
    i = 0 
    for line in (myfile):
        i += 1
        if i%4 == 0: #Qual lines .fq
            line = line.decode("ascii")
            line = line.strip("\n")
            x = 0
            for char in line:
                score = bioinfo.convert_phred(char)
                all_qscores_sum[x] += score
                x += 1
    #print(len(all_qscores_sum)) 
    #print(all_qscores_sum)
    for index in range(len(all_qscores_sum)):
        mean[index] = all_qscores_sum[index] / args.reads
print('All qscores list:')
print(all_qscores_sum)
print('Mean list:')
print(mean)
x = range(args.length)
y = mean
plt.bar(x,y, color = "tomato")
plt.ylabel("Quality Score")
plt.xlabel("Position")
plt.title("R4 Avg.Quality Scores per NT Position")
plt.savefig(f"R4.png")
    
    

