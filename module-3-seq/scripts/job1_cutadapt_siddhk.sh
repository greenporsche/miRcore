#!/bin/bash

#SBATCH --job-name=job1_cutadapt.sh
#SBATCH --mail-user=uniqname@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100m
#SBATCH --time=2:00:00
#SBATCH --account=mircore_project2
#SBATCH --partition=standard
#SBATCH --output=/home/%u/syg/%x-%j.log

module load Bioinformatics
module load cutadapt

cd /home/uniqname/syg

cutadapt -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o cuniqname.1.fq.gz -p cuniqname.2.fq.gz uniqname.1.fq.gz uniqname.2.fq.gz

