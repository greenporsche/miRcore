#!/bin/bash

module load Bioinformatics
module load samtools

samtools view hg38_uniqname.sorted.bam chr9:133,255,176-133,275,214 -b > ABO_uniqname.sorted.bam
samtools index ABO_uniqname.sorted.bam

