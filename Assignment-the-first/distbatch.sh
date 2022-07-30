#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --job-name=dist.py(R4)
#SBATCH --time=0-05:00:00

cd /projects/bgmp/dmarro/bioinfo/Bi622/Demultiplex/Assignment-the-first
conda activate bgmp_py310
/usr/bin/time -v ./dist.py \
 -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -r 363246735 -l 101