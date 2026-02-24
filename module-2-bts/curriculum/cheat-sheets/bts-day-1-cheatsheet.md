# BTS Day 1 Cheat Sheet: Biotech Intro and Linux Basics

## Commands & Tools
| Tool/Command | What It Does | Example Usage |
|---|---|---|
| `ssh` | Connect to Great Lakes cluster remotely | `ssh uniqname@greatlakes.arc-ts.umich.edu` |
| `pwd` | Print working (current) directory | `pwd` |
| `ls` | List files and folders in a directory | `ls` or `ls -l` |
| `mkdir` | Create a new directory | `mkdir my_project` |
| `cd` | Change directory | `cd my_project` or `cd ..` |
| `cp` | Copy a file | `cp source.txt destination.txt` |
| `cat` | Display contents of a file | `cat sequence.fasta` |
| Duo Mobile | Two-factor authentication for UM login | Approve push notification on phone |

## Key Concepts
- **Biotechnology**: The use of biological systems or living organisms to develop products and technologies (e.g., insulin production, CRISPR gene editing, PCR testing, Illumina sequencing)
- **Wet Lab vs. Dry Lab**: Wet lab involves physical experiments with chemicals and biological materials; dry lab uses computers for data analysis (bioinformatics)
- **SARS-CoV-2**: The virus that causes COVID-19; contains a single-stranded RNA genome of ~30,000 nucleotides
- **Spike Protein (S-protein)**: The protein on the surface of SARS-CoV-2 that binds to human ACE2 receptors, enabling viral entry into cells; located at bases 21,563-25,384 in the genome
- **NCBI (National Center for Biotechnology Information)**: A public database for accessing genomic sequences, protein structures, and biological data
- **FASTA Format**: A text-based format for representing nucleotide or amino acid sequences, starting with a `>` header line
- **PDB (Protein Data Bank)**: A database of 3D structural data of proteins and nucleic acids
- **Great Lakes Cluster**: University of Michigan's high-performance computing (HPC) cluster used for bioinformatics analysis
- **Linux**: An open-source operating system commonly used in bioinformatics; the Great Lakes cluster runs Linux
- **Reference Genome**: A standardized representative genome sequence used as a baseline for comparison (e.g., SARS-CoV-2 reference sequence NC_045512.2)

## File Formats
- **FASTA (.fasta, .fa)**: Text file with a header line starting with `>` followed by sequence data; used for storing nucleotide or protein sequences
- **PDB (.pdb)**: Protein Data Bank file containing 3D atomic coordinates of a protein structure

## Databases & URLs
| Resource | URL | Used For |
|---|---|---|
| NCBI | https://www.ncbi.nlm.nih.gov/ | Searching genomes, genes, proteins |
| NCBI Nucleotide | https://www.ncbi.nlm.nih.gov/nuccore/ | Looking up SARS-CoV-2 reference genome |
| Protein Data Bank (PDB) | https://www.rcsb.org/ | Viewing 3D protein structures |
| Great Lakes HPC | https://greatlakes.arc-ts.umich.edu/ | Remote Linux computing cluster |

## Common Pitfalls
- Forgetting to authenticate with Duo Mobile when logging into the Great Lakes cluster
- Confusing `cd` (change directory) with `cp` (copy) --- both start with "c"
- Not realizing you need to type `ls` to see what is in your current directory before navigating
- Typing `cd` without a path returns you to your home directory, which may be confusing
- Mixing up FASTA format (sequences) with FASTQ format (sequences + quality scores, introduced later)
- Forgetting that Linux commands are case-sensitive (`LS` is not the same as `ls`)

## Quick Check
1. What is the difference between a wet lab and a dry lab?
2. What Linux command shows you which directory you are currently in?
3. What does the `>` character indicate at the beginning of a line in a FASTA file?
4. Where on the SARS-CoV-2 genome is the spike protein located (base pair range)?
5. How do you create a new folder called `results` in Linux?
