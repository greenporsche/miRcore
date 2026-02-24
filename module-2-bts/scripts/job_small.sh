#!/bin/bash
#SBATCH --job-name=job_small.sh
#SBATCH --mail-user=uniqname@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000m
#SBATCH --time=10:00
#SBATCH --account=mircore_project1
#SBATCH --partition=standard
#SBATCH --output=/home/%u/%x-%j.log

cd /home/uniqname/bts2025

module load Bioinformatics
module load bowtie2
module load samtools

bowtie2 -x sars2 -U small.fq -S sars2_small.sam
samtools view -Sb sars2_small.sam -o sars2_small.bam
samtools sort sars2_small.bam -o sars2_small.sorted.bam
samtools index sars2_small.sorted.bam
