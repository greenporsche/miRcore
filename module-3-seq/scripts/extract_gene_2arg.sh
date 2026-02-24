#!/bin/bash

module load Bioinformatics
module load samtools

samtools view hg38_uniqname.sorted.bam $1 -b > ${2}_uniqname.sorted.bam
samtools index ${2}_uniqname.sorted.bam

