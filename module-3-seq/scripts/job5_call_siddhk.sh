#!/bin/bash

#SBATCH --job-name=job5_call.sh
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
module load bcftools

cd /home/uniqname/syg/
bcftools index hg38_uniqname.bcf
bcftools call -mv -Oz hg38_uniqname.bcf -o hg38_uniqname.vcf.gz
bcftools index hg38_uniqname.vcf.gz
bcftools view -r chr9:133255176-133275214 hg38_uniqname.vcf.gz > ABO_uniqname.vcf
