#!/bin/bash

#SBATCH --job-name=job3_samtools.sh
#SBATCH --mail-user=uniqname@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000m
#SBATCH --time=2:00:00
#SBATCH --account=mircore_project2
#SBATCH --partition=standard
#SBATCH --output=/home/%u/syg/%x-%j.log

module load Bioinformatics
module load samtools

cd /home/uniqname/syg/

samtools view -Sb hg38_uniqname.sam -o hg38_uniqname.bam
samtools sort hg38_uniqname.bam -o hg38_uniqname.sorted.bam
samtools index hg38_uniqname.sorted.bam
