#!/bin/bash
#SBATCH --account=bgmp         ### SLURM account which will be charged for the job
#SBATCH --job-name=Demultiplex   ### Job Name
#SBATCH --output=Demultiplex_%j.out         ### File in which to store job output
#SBATCH --error=Demultiplex-%j.err          ### File in which to store job error messages
#SBATCH --time=0-24:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --cpus-per-task=4  ### Number of cpus (cores) per task
#SBATCH --partition=bgmp          ### partition to run things

cd /projects/bgmp/dmarro/bioinfo/Bi622/Demultiplex/Assignment-the-third
conda activate bgmp_py310
/usr/bin/time -v ./demux.py \
-f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -ind /projects/bgmp/shared/2017_sequencing/indexes.txt
