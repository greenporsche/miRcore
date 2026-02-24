#!/bin/bash
#SBATCH --job-name=job_p20.sh
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

bowtie2 -x /scratch/mircore_project_root/mircore_project1/GRCh38_noalt_as/GRCh38_noalt_as -U p20.fq -S p20_hg38.sam --un un_p20_hg38.fq
samtools view -Sb p20_hg38.sam -o p20_hg38.bam
samtools sort p20_hg38.bam -o p20_hg38.sorted.bam
samtools index p20_hg38.sorted.bam

bowtie2 -x sars2 -U un_p20_hg38.fq -S p20_sars2.sam
samtools view -Sb p20_sars2.sam -o p20_sars2.bam
samtools sort p20_sars2.bam -o p20_sars2.sorted.bam
samtools index p20_sars2.sorted.bam
