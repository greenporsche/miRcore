#!/bin/bash

#SBATCH --job-name=job4_mpileup.sh
#SBATCH --mail-user=uniqname@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=600m 
#SBATCH --time=5:00:00
#SBATCH --account=mircore_project2
#SBATCH --partition=standard
#SBATCH --output=/home/%u/syg/%x-%j.log

module load Bioinformatics
module load bcftools

cd /home/uniqname/syg/
bcftools mpileup -f /scratch/mircore_project_root/mircore_project2/hg38/hg38.fa hg38_uniqname.sorted.bam -o hg38_uniqname.bcf
