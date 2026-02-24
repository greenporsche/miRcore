#!/bin/bash

#SBATCH --job-name=job2_bowtie2.sh
#SBATCH --mail-user=uniqname@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4000m 
#SBATCH --time=15:00:00
#SBATCH --account=mircore_project2
#SBATCH --partition=standard
#SBATCH --output=/home/%u/syg/%x-%j.log

module load Bioinformatics
module load bowtie2

cd /home/uniqname/syg/

bowtie2 --no-unal -x /scratch/mircore_project_root/mircore_project2/GRCh38_noalt_as/GRCh38_noalt_as -1 cuniqname.1.fq.gz -2 cuniqname.2.fq.gz -S hg38_uniqname.sam --un-conc-gz unhg38_uniqname.fq.gz 2> map_hg38_uniqname.txt

