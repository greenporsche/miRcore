# SEQ Day 2 Cheat Sheet: Exome Data Processing

## Commands & Tools
| Tool/Command | What It Does | Example Usage |
|---|---|---|
| `bwa mem` | Map/align reads to the human reference genome (hg38) | `bwa mem ref.fa read1.fq read2.fq > aligned.sam` |
| `samtools view` | Convert SAM to BAM (binary) format | `samtools view -bS aligned.sam > aligned.bam` |
| `samtools sort` | Sort BAM file by genomic coordinate | `samtools sort aligned.bam -o sorted.bam` |
| `samtools index` | Create index (.bai) for a sorted BAM file | `samtools index sorted.bam` |
| `samtools flagstat` | Summary statistics of alignment | `samtools flagstat sorted.bam` |
| `gatk HaplotypeCaller` | Call variants from aligned reads | `gatk HaplotypeCaller -R ref.fa -I sorted.bam -O output.vcf` |
| `sbatch` | Submit pipeline job to SLURM | `sbatch job2_mapping.sh` |
| `module load` | Load BWA, samtools, GATK modules | `module load BWA samtools GATK` |
| `ls -lh` | List files with human-readable sizes | `ls -lh *.bam` |

## Key Concepts
- **Read mapping/alignment**: Process of matching short sequencing reads to their location on the reference genome
- **Reference genome (hg38)**: The current human reference assembly (GRCh38) used as the coordinate backbone for alignment
- **SAM (Sequence Alignment/Map)**: Human-readable text format storing alignment information for each read
- **BAM**: Binary compressed version of SAM; much smaller file size, requires `samtools` to view
- **BAI (BAM index)**: Index file that enables fast random access to regions within a sorted BAM file
- **Mapping quality (MAPQ)**: Score indicating confidence that a read is mapped to the correct location
- **Variant calling**: Computational process of identifying positions where your DNA differs from the reference genome
- **Pipeline**: Series of bioinformatics steps (trim --> align --> sort --> index --> call variants) run sequentially

## File Formats
- **SAM (.sam)**: Tab-delimited text; each line is one read alignment with chromosome, position, CIGAR string, quality
- **BAM (.bam)**: Binary compressed SAM; primary working format for aligned data
- **BAI (.bai)**: Index for BAM enabling rapid regional queries; must match the BAM filename
- **VCF (.vcf)**: Variant Call Format; output of variant calling with position, ref/alt alleles, quality, genotype
- **BED (.bed)**: Tab-delimited genomic intervals (0-based start); used for exome target regions

## Databases & URLs
| Resource | URL | Used For |
|---|---|---|
| UCSC Genome Browser | `https://genome.ucsc.edu` | Visualizing genomic regions, gene structures, variant tracks |
| UCSC Table Browser | `https://genome.ucsc.edu/cgi-bin/hgTables` | Downloading annotation data for specific regions |
| Ensembl | `https://ensembl.org` | Alternative genome browser and gene annotation resource |
| SAM format specification | `https://samtools.github.io/hts-specs/SAMv1.pdf` | Detailed SAM/BAM format documentation |

## Common Pitfalls
- Forgetting to sort the BAM before indexing (samtools index requires a coordinate-sorted BAM)
- Running out of disk space on scratch when generating large BAM files -- check with `ls -lh`
- Not loading all required modules (BWA, samtools, GATK) in the job script
- Confusing 0-based BED coordinates with 1-based VCF/SAM coordinates when looking up positions
- Accidentally overwriting sorted BAM files by using the same filename for input and output
- Not verifying alignment quality with `samtools flagstat` before proceeding to variant calling

## Quick Check
1. What is the difference between SAM and BAM formats?
2. Why must a BAM file be sorted before it can be indexed?
3. What reference genome assembly does the SEQ camp pipeline use?
4. In the UCSC Genome Browser, how do you navigate to a specific gene?
5. What is the purpose of the variant calling step, and which tool performs it in our pipeline?
