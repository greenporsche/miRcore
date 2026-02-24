# BTS Day 2 Cheat Sheet: NGS and FASTQ Analysis

## Commands & Tools
| Tool/Command | What It Does | Example Usage |
|---|---|---|
| `grep` | Search for patterns in a file | `grep "ATCG" reads.fastq` |
| `grep -c` | Count the number of matching lines | `grep -c "@" reads.fastq` |
| `wc -l` | Count the number of lines in a file | `wc -l reads.fastq` |
| `head` | Display the first lines of a file | `head -n 4 reads.fastq` |
| `tail` | Display the last lines of a file | `tail -n 4 reads.fastq` |
| `scp` | Securely copy files between local and remote machines | `scp user@greatlakes:~/file.txt .` |
| `module load` | Load a software module on Great Lakes | `module load Bioinformatics` |
| `less` | View a file one screen at a time (press `q` to quit) | `less reads.fastq` |

## Key Concepts
- **Next-Generation Sequencing (NGS)**: High-throughput DNA sequencing technology that generates millions of short reads simultaneously; used by platforms like Illumina
- **Reads**: Short DNA sequences (typically 50-300 bp) produced by NGS machines from fragmented DNA/RNA
- **FASTQ Format**: A text file format that stores both the nucleotide sequence and its per-base quality scores; each read occupies exactly 4 lines
- **Phred Quality Score**: A measure of the confidence that a base was called correctly; encoded as ASCII characters in FASTQ files (higher = better)
- **GC Content**: The percentage of guanine (G) and cytosine (C) bases in a sequence; important because GC pairs have 3 hydrogen bonds (stronger than AT's 2)
- **Nucleotide Hybridization**: The process by which complementary single-stranded DNA/RNA molecules bind together (A-T/U, G-C)
- **Melting Temperature (Tm)**: The temperature at which 50% of double-stranded DNA denatures into single strands; higher GC content = higher Tm
- **Spectrometry/Spectrophotometry**: A technique to measure the concentration and purity of DNA/RNA using light absorbance (e.g., 260nm wavelength)
- **Business Pitch**: A concise presentation to potential investors about a biotech company, product, or idea; includes problem, solution, market, and team

## File Formats
- **FASTQ (.fastq, .fq)**: Four lines per read: (1) `@` header, (2) sequence, (3) `+` separator, (4) quality scores as ASCII characters
- **FASTA (.fasta, .fa)**: Two lines per sequence: (1) `>` header, (2) sequence data (no quality info)

### FASTQ 4-Line Format
```
@SEQ_ID            <-- Line 1: Header (starts with @)
GATTTGGGGTTCAAAG   <-- Line 2: DNA Sequence
+                  <-- Line 3: Separator (starts with +)
!''*((((***+))%%   <-- Line 4: Quality Scores (ASCII encoded)
```

## Databases & URLs
| Resource | URL | Used For |
|---|---|---|
| NCBI | https://www.ncbi.nlm.nih.gov/ | Reference sequences and genome data |
| Illumina | https://www.illumina.com/ | NGS platform manufacturer |
| Great Lakes HPC | https://greatlakes.arc-ts.umich.edu/ | Running bioinformatics analyses |

## Common Pitfalls
- Forgetting that FASTQ files have **exactly 4 lines per read** --- dividing total line count by 4 gives you the number of reads
- Confusing FASTA (sequences only) with FASTQ (sequences + quality scores)
- Using `grep -c ">"` to count reads in a FASTQ file instead of `grep -c "@"` (FASTA uses `>`, FASTQ uses `@`)
- Misunderstanding quality scores: they are ASCII characters, not numbers; each character maps to a Phred score
- Forgetting to use `module load Bioinformatics` before trying to use bioinformatics tools on Great Lakes
- Mixing up GC content calculation --- count only G and C bases, divide by total bases, multiply by 100

## Quick Check
1. How many lines does a single read occupy in a FASTQ file?
2. If a FASTQ file has 2,000,000 lines, how many reads does it contain?
3. What Linux command would you use to count how many reads are in a FASTQ file?
4. Why does higher GC content lead to a higher melting temperature?
5. What is the difference between FASTA and FASTQ formats?
