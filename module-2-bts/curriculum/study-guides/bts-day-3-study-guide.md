# BTS Day 3 Study Guide: SAM/BAM Files and Alignment

## Learning Objectives
By the end of this day, you should be able to:
- Explain the purpose of sequence alignment and why it is the central step in the bioinformatics pipeline
- Use Bowtie2 to align FASTQ reads against a reference genome
- Interpret Bowtie2 alignment summary statistics (alignment rate, concordance)
- Describe the SAM file format and identify key fields (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR)
- Convert SAM to BAM, sort, and index using SAMtools
- Use `cut`, `sort`, and `nano` for file manipulation on the command line
- Deliver a group business pitch presentation on a real biotech company

## Key Vocabulary
| Term | Definition |
|---|---|
| Sequence Alignment | The computational process of matching short reads to positions on a reference genome |
| Bowtie2 | A fast, memory-efficient short read aligner that uses a Burrows-Wheeler Transform index |
| Reference Index | A pre-computed data structure built from the reference genome that enables fast alignment searching |
| SAM | Sequence Alignment/Map format; a tab-delimited text file storing alignment information |
| BAM | Binary Alignment/Map format; the compressed binary equivalent of SAM |
| BAI | BAM Index file; enables random access to specific regions of a sorted BAM file |
| MAPQ | Mapping Quality; a Phred-scaled score indicating the confidence that a read is correctly placed |
| CIGAR | Compact Idiosyncratic Gapped Alignment Report; a string describing how a read aligns to the reference |
| FLAG | A bitwise integer in SAM encoding alignment properties (mapped/unmapped, strand, etc.) |
| SAMtools | A suite of command-line tools for manipulating SAM/BAM files |
| nano | A simple terminal-based text editor |
| cut | A Linux command for extracting columns from tab-delimited files |
| sort | A Linux command for sorting lines of text files |
| Pipeline | A series of computational steps where the output of one step feeds into the next |

## Concept Summaries

### Sequence Alignment with Bowtie2
Sequence alignment is the most critical step in the BTS bioinformatics pipeline. The goal is to take millions of short reads from a FASTQ file and determine where each one originated on a reference genome. This process answers a fundamental question: does this patient sample contain SARS-CoV-2?

Bowtie2 is the alignment tool used in this camp. It works by first loading a pre-built index of the reference genome (created using `bowtie2-build`), then systematically comparing each read against the indexed reference to find the best matching position. The index uses a data structure called the Burrows-Wheeler Transform, which compresses the reference genome in a way that makes searching extremely fast --- much faster than a simple text search.

The Bowtie2 command has three essential arguments: `-x` specifies the reference index prefix, `-U` specifies the input FASTQ file (for unpaired/single-end reads), and `-S` specifies the output SAM file. After alignment, Bowtie2 prints a summary to the screen showing the total number of reads processed, how many aligned zero times (unmapped), how many aligned exactly once (uniquely mapped), and how many aligned more than once (multi-mapped). The overall alignment rate is the most important number: a high rate to SARS-CoV-2 suggests the sample contains viral sequences.

### The SAM File Format
SAM (Sequence Alignment/Map) is the standard text format for storing alignment results. It is tab-delimited and human-readable, meaning you can view it with commands like `cat`, `head`, or `less`. A SAM file has two parts: header lines (beginning with `@`) that describe the reference sequences and alignment metadata, and alignment records (one per line) containing the actual mapping information.

Each alignment record has at least 11 mandatory fields. The most important for this camp are: QNAME (the read name), FLAG (a bitwise integer encoding alignment properties like whether the read is mapped or unmapped), RNAME (the reference sequence name the read aligned to), POS (the 1-based leftmost mapping position), MAPQ (the mapping quality score, with higher values indicating greater confidence), CIGAR (a compact string describing the alignment, such as `150M` for 150 matching bases), SEQ (the read sequence), and QUAL (the base quality scores).

Understanding SAM fields allows you to extract meaningful information. For example, you can use `cut -f3` to extract the reference name column and see which genome your reads mapped to, or `cut -f5` to examine the distribution of mapping quality scores.

### SAM to BAM Processing with SAMtools
While SAM files are human-readable, they are very large and slow to process. BAM is the binary compressed version of SAM --- it contains identical information but takes up 3-5 times less disk space and can be processed much faster. The conversion from SAM to BAM is performed using `samtools view -bS`.

After conversion, the BAM file must be sorted by genomic position using `samtools sort`. Sorting arranges all reads in order of their mapping position along the reference genome. This is a prerequisite for the final step: indexing with `samtools index`. The index (a `.bai` file) creates a lookup table that enables tools like IGV to quickly jump to any position in the genome without reading the entire file.

The complete pipeline for this step is always the same four commands in sequence:
1. `bowtie2 -x ref -U reads.fastq -S output.sam` (align)
2. `samtools view -bS output.sam > output.bam` (convert)
3. `samtools sort output.bam -o output.sorted.bam` (sort)
4. `samtools index output.sorted.bam` (index)

This pipeline is so fundamental that students will use it repeatedly for the rest of the camp, eventually automating it in shell scripts and SLURM job scripts.

### Business Pitch Presentations
On Day 3, student groups deliver their first business pitch presentations. Groups research real biotech companies or technologies --- examples from the camp include CRISPR Therapeutics (gene editing for sickle cell disease and cancer) and rapid COVID-19 testing companies. Each group presents to a panel of "investors" (instructors and peers), covering the company's mission, its core technology, target market, competitive advantages, and potential for growth.

The pitch component teaches communication skills that complement the technical bioinformatics work. Scientists and engineers must be able to explain complex technologies in accessible terms, justify the market need for their innovations, and make a compelling case for investment. Instructor feedback after each pitch helps students refine their presentation skills for the final Day 5 pitch.

### Linux Text Manipulation: cut, sort, and nano
Day 3 introduces additional Linux tools for working with structured text data. The `cut` command extracts specific columns from tab-delimited files, which is invaluable for parsing SAM files (e.g., `cut -f1,3,4,5` extracts the read name, reference name, position, and mapping quality). The `sort` command arranges lines alphabetically or numerically, and when combined with `uniq -c`, can count occurrences of each unique value.

The `nano` text editor is a beginner-friendly alternative to more powerful but complex editors like `vim`. Students use nano to create notes files, write simple scripts, and edit configuration files. Key nano shortcuts include `Ctrl+O` to save, `Ctrl+X` to exit, and `Ctrl+K` to cut a line.

## How It Connects
- **Previous Day**: Day 2 introduced the FASTQ format and NGS concepts. Students now take those FASTQ files and align them to a reference genome, transforming raw sequences into positional information.
- **Next Day**: Day 4 extends the alignment pipeline to the human genome (HG38), introduces IGV for visualization, and covers SLURM job submission for running alignments on compute nodes rather than interactively.

## Review Questions
1. What three arguments are essential when running a Bowtie2 alignment command?
2. Explain what an overall alignment rate of 95% means in the context of patient sample analysis.
3. In a SAM file, what do columns 3 and 4 represent?
4. What is the purpose of converting SAM to BAM format?
5. Why must BAM files be sorted before they can be indexed?
6. Describe the complete 4-step pipeline from FASTQ to an indexed BAM file.
7. What does a MAPQ score of 0 tell you about a read's alignment?
8. How could you use the `cut` command to extract read names and their mapping positions from a SAM file?

## Further Reading
- Bowtie2 Manual: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
- SAM Format Specification: https://samtools.github.io/hts-specs/SAMv1.pdf
- SAMtools Documentation: https://www.htslib.org/doc/samtools.html
- Understanding SAM Flags: https://broadinstitute.github.io/picard/explain-flags.html
- Burrows-Wheeler Transform Explained: https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform
