#!/bin/bash
#SBATCH --job-name=job_hg38_small.sh
#SBATCH --mail-user=uniqname@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4000m
#SBATCH --time=10:00
#SBATCH --account=mircore_project1
#SBATCH --partition=standard
#SBATCH --output=/home/%u/%x-%j.log

cd /home/uniqname/bts2025

module load Bioinformatics
module load bowtie2
module load samtools

bowtie2 -x /scratch/mircore_project_root/mircore_project1/GRCh38_noalt_as/GRCh38_noalt_as -U small.fq -S hg38_small.sam > hg38_small.txt
samtools view -Sb hg38_small.sam -o hg38_small.bam
samtools sort hg38_small.bam -o hg38_small.sorted.bam
samtools index hg38_small.sorted.bam
