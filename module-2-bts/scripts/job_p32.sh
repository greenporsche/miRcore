#!/bin/bash
#SBATCH --job-name=job_p32.sh
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

bowtie2 -x /scratch/mircore_project_root/mircore_project1/GRCh38_noalt_as/GRCh38_noalt_as -U p32.fq -S p32_hg38.sam --un un_p32_hg38.fq
samtools view -Sb p32_hg38.sam -o p32_hg38.bam
samtools sort p32_hg38.bam -o p32_hg38.sorted.bam
samtools index p32_hg38.sorted.bam

bowtie2 -x sars2 -U un_p32_hg38.fq -S p32_sars2.sam
samtools view -Sb p32_sars2.sam -o p32_sars2.bam
samtools sort p32_sars2.bam -o p32_sars2.sorted.bam
samtools index p32_sars2.sorted.bam
