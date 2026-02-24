# BTS Day 2 Study Guide: NGS and FASTQ Analysis

## Learning Objectives
By the end of this day, you should be able to:
- Explain what next-generation sequencing (NGS) is and how it differs from Sanger sequencing
- Describe the FASTQ file format and its 4-line structure
- Interpret Phred quality scores and explain their significance
- Use `grep`, `wc`, and `head` to analyze FASTQ files on the command line
- Calculate GC content and explain its biological significance
- Explain nucleotide hybridization and the concept of melting temperature (Tm)
- Transfer files between Great Lakes and your local machine using `scp`
- Describe the basic structure and components of a biotech business pitch

## Key Vocabulary
| Term | Definition |
|---|---|
| Next-Generation Sequencing (NGS) | High-throughput DNA sequencing technology that reads millions of short DNA fragments simultaneously |
| Reads | Short DNA sequences (typically 50-300 base pairs) produced by an NGS instrument |
| FASTQ | A file format that stores both nucleotide sequences and their per-base quality scores in a 4-line-per-read structure |
| Phred Quality Score | A numerical score indicating the probability that a given base call is incorrect; Q30 means 1 in 1,000 chance of error |
| GC Content | The percentage of bases in a DNA sequence that are guanine (G) or cytosine (C) |
| Nucleotide Hybridization | The process by which complementary single-stranded nucleic acid molecules bind together via hydrogen bonds |
| Melting Temperature (Tm) | The temperature at which 50% of a double-stranded DNA population has separated into single strands |
| Complementary Base Pairing | A pairs with T (2 H-bonds); G pairs with C (3 H-bonds); in RNA, A pairs with U |
| Spectrophotometry | A technique to measure DNA/RNA concentration and purity by absorbance of light at specific wavelengths |
| scp | Secure Copy Protocol, a command for transferring files between local and remote machines via SSH |
| grep | A Linux command for searching text patterns in files |
| Business Pitch | A structured presentation designed to persuade investors to fund a biotech company or product |

## Concept Summaries

### Next-Generation Sequencing
Next-generation sequencing (NGS) represents a revolutionary advance over the original Sanger sequencing method. While Sanger sequencing reads one DNA fragment at a time, NGS platforms like Illumina can read millions of short DNA fragments in parallel. This massive parallelization makes it possible to sequence an entire human genome in days rather than years, at a fraction of the cost.

In the context of the BTS camp, NGS is the technology that generated the patient data students analyze. When a patient sample is collected, the DNA or RNA is extracted, fragmented into small pieces, and sequenced. The result is millions of short "reads" --- typically 50 to 300 base pairs long --- stored in FASTQ files. The challenge then becomes computational: using bioinformatics tools to align these millions of short reads back to a reference genome to determine what organism the DNA came from and whether any mutations exist.

The quality of sequencing data is not perfect. Each base in a read has an associated confidence score, which is why the FASTQ format includes quality information alongside the sequence. Understanding these quality scores is essential for interpreting sequencing results correctly.

### The FASTQ File Format
The FASTQ format is the standard output format for NGS instruments. Each read consists of exactly four lines, making the format highly structured and predictable:

Line 1 begins with `@` and contains the read identifier (a unique name for the read, often including instrument and run information). Line 2 contains the actual DNA sequence (A, T, G, C characters). Line 3 begins with `+` and serves as a separator (it may optionally repeat the read identifier). Line 4 contains the quality scores encoded as ASCII characters, with exactly one character per base.

A key principle: **total lines divided by 4 equals total reads**. This simple math is the fastest way to count reads in a FASTQ file. For example, a file with 4,000,000 lines contains exactly 1,000,000 reads.

The quality scores on line 4 are Phred scores encoded as ASCII characters. A Phred score of 20 (Q20) means a 1% chance the base call is wrong. A Phred score of 30 (Q30) means a 0.1% chance. Most modern Illumina data has average qualities of Q30 or above, meaning the vast majority of base calls are highly reliable. The ASCII encoding maps each quality value to a printable character, allowing the quality string to be stored compactly on a single line.

### GC Content and Melting Temperature
GC content is the percentage of bases in a DNA sequence that are guanine (G) or cytosine (C). It is calculated as: GC% = (G + C) / (A + T + G + C) * 100. GC content is biologically meaningful because G-C base pairs are held together by three hydrogen bonds, while A-T base pairs have only two. This makes GC-rich regions of DNA more thermodynamically stable.

The melting temperature (Tm) of a DNA molecule is the temperature at which half of the double-stranded molecules have separated into single strands (denatured). Higher GC content correlates with higher melting temperature because more energy (heat) is needed to break the three hydrogen bonds of GC pairs compared to the two bonds of AT pairs. This concept is important in laboratory techniques like PCR, where primer design must account for Tm to ensure proper hybridization.

SARS-CoV-2 has a relatively low GC content (~38%), which is a characteristic shared by many coronaviruses. When analyzing patient reads, GC content can serve as a quick sanity check --- reads from SARS-CoV-2 should have a GC content near this expected value.

### Nucleotide Hybridization and Spectrophotometry
Nucleotide hybridization is the fundamental biochemical process by which complementary single-stranded nucleic acids find each other and form stable double-stranded structures. This occurs through Watson-Crick base pairing: adenine (A) pairs with thymine (T) in DNA (or uracil (U) in RNA), and guanine (G) pairs with cytosine (C). The specificity of this pairing is what makes all sequencing and alignment work possible.

Spectrophotometry (specifically UV spectrophotometry at 260nm) is a laboratory technique for measuring DNA/RNA concentration. Nucleic acids absorb UV light at 260nm, and the absorbance is proportional to the concentration. This technique is used in wet labs to quantify how much DNA was extracted from a sample before sequencing. While students in the BTS camp focus on dry lab analysis, understanding how the data was generated in the wet lab provides important context for interpreting quality and reliability.

### Introduction to Business Pitches
The BTS camp has a dual focus: bioinformatics research skills and biotech startup/entrepreneurship. On Day 2, students are introduced to the business pitch component, where they work in groups to research a real biotech company or product and prepare a pitch presentation.

A strong business pitch has several key components: a clear problem statement (what unmet need exists), a proposed solution (the product or technology), the target market (who will buy it and how large is the market), the competitive advantage (why this solution is better than alternatives), and the ask (how much funding is needed and what it will accomplish). Students learn that communicating scientific innovation to non-scientists and potential investors is a critical skill in the biotech industry.

## How It Connects
- **Previous Day**: Day 1 introduced biotechnology concepts, SARS-CoV-2 biology, NCBI navigation, and basic Linux commands. Those foundational skills are directly used on Day 2 for working with FASTQ files on the command line.
- **Next Day**: Day 3 introduces sequence alignment with Bowtie2, SAM/BAM file formats, and SAMtools. Students will take the FASTQ files they explored on Day 2 and align them to a reference genome, producing alignment results for downstream analysis.

## Review Questions
1. How does NGS differ from Sanger sequencing in terms of throughput?
2. Draw or describe the 4-line structure of a single read in a FASTQ file.
3. What does a Phred quality score of 30 tell you about a base call?
4. If a FASTQ file has 8,000,000 lines, how many reads does it contain?
5. Why does a DNA molecule with 70% GC content have a higher melting temperature than one with 30% GC content?
6. Write the Linux command to count the number of reads in a file called `sample.fastq`.
7. What is the purpose of the `scp` command, and why can you not run it from within your Great Lakes session?
8. Name three essential components of a biotech business pitch.

## Further Reading
- Illumina Sequencing Technology Overview: https://www.illumina.com/science/technology/next-generation-sequencing.html
- FASTQ Format Specification: https://maq.sourceforge.net/fastq.shtml
- Understanding Phred Quality Scores: https://drive5.com/usearch/manual/quality_score.html
- DNA Melting Temperature Calculator: https://www.biophp.org/minitools/melting_temperature/
- Biotech Startup Pitch Guide: https://www.ycombinator.com/library/6q-how-to-pitch-your-startup
