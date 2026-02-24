# BTS Day 3 Cheat Sheet: SAM/BAM Files and Alignment

## Commands & Tools
| Tool/Command | What It Does | Example Usage |
|---|---|---|
| `bowtie2` | Align reads to a reference genome | `bowtie2 -x ref_index -U reads.fastq -S output.sam` |
| `samtools view` | Convert SAM to BAM format | `samtools view -bS output.sam > output.bam` |
| `samtools sort` | Sort a BAM file by genomic position | `samtools sort output.bam -o output.sorted.bam` |
| `samtools index` | Create an index for a sorted BAM file | `samtools index output.sorted.bam` |
| `nano` | Simple text editor in the terminal | `nano myfile.txt` |
| `cut` | Extract specific columns from a file | `cut -f1,3,4 output.sam` |
| `sort` | Sort lines of a text file | `sort -k3 output.sam` |
| `module load Bioinformatics` | Load bioinformatics software on Great Lakes | `module load Bioinformatics` |
| `module load Bowtie2` | Load Bowtie2 alignment tool | `module load Bowtie2` |
| `module load SAMtools` | Load SAMtools for SAM/BAM processing | `module load SAMtools` |

## Key Concepts
- **Sequence Alignment**: The process of matching short sequencing reads against a reference genome to determine where each read originated
- **Bowtie2**: A fast and memory-efficient tool for aligning sequencing reads to a reference genome; uses a Burrows-Wheeler index
- **Reference Index**: A pre-built data structure (from `bowtie2-build`) that allows Bowtie2 to quickly search the reference genome
- **SAM (Sequence Alignment/Map)**: A tab-delimited text format storing alignment information; human-readable
- **BAM (Binary Alignment/Map)**: The compressed binary version of SAM; not human-readable but much smaller and faster to process
- **Mapping Quality (MAPQ)**: A score indicating the confidence that a read is aligned to the correct position; higher is better
- **FLAG Field**: A bitwise integer in SAM files encoding properties of the alignment (e.g., whether the read is unmapped, reverse complemented)
- **CIGAR String**: A compact representation of how a read aligns to the reference (e.g., `50M` = 50 bases match, `3I` = 3 base insertion)
- **Business Pitch Presentations**: Students present pitches on real biotech companies (e.g., CRISPR Therapeutics, rapid COVID testing)

## File Formats
- **SAM (.sam)**: Tab-separated text file with header lines (`@`) and alignment records; each record has 11+ mandatory fields (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, etc.)
- **BAM (.bam)**: Binary compressed version of SAM; requires `samtools view` to read
- **BAI (.bai)**: Index file for sorted BAM; required for visualization in IGV

### SAM File Key Columns
| Column | Field | Description |
|---|---|---|
| 1 | QNAME | Read/query name |
| 2 | FLAG | Bitwise alignment flag |
| 3 | RNAME | Reference sequence name |
| 4 | POS | 1-based leftmost mapping position |
| 5 | MAPQ | Mapping quality score |
| 6 | CIGAR | Alignment description string |
| 10 | SEQ | Read sequence |
| 11 | QUAL | Base quality scores |

## Databases & URLs
| Resource | URL | Used For |
|---|---|---|
| Bowtie2 Manual | https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml | Reference for Bowtie2 options |
| SAM Format Spec | https://samtools.github.io/hts-specs/SAMv1.pdf | Official SAM format documentation |
| SAMtools | https://www.htslib.org/ | SAM/BAM processing toolkit |

## Common Pitfalls
- Forgetting to load both `Bioinformatics` and `Bowtie2` modules before running alignment
- Confusing the Bowtie2 index prefix with an actual file --- the index consists of multiple `.bt2` files, but you specify only the prefix name
- Forgetting the `-bS` flag in `samtools view` when converting SAM to BAM (without it, you get text output instead of binary)
- Not sorting the BAM file before indexing --- `samtools index` requires a position-sorted BAM
- Trying to open a BAM file with `cat` or `nano` --- BAM is binary and will show garbled text
- Mixing up the pipeline order: FASTQ -> SAM (via Bowtie2) -> BAM (via samtools view) -> sorted BAM (via samtools sort) -> index (via samtools index)

## Quick Check
1. What is the correct order of the alignment pipeline from raw reads to a viewable alignment?
2. What command converts a SAM file to BAM format?
3. Why must a BAM file be sorted before it can be indexed?
4. What does a MAPQ score of 0 mean for an aligned read?
5. In a SAM file, which column contains the reference genome name that a read mapped to?
