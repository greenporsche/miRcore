# SEQ Day 2 Study Guide: Exome Data Processing

## Learning Objectives
By the end of this day, you should be able to:
- Explain the purpose of read alignment and why a reference genome is needed
- Describe what the BWA aligner does and how it maps short reads to genomic coordinates
- Distinguish between SAM, BAM, and BAI file formats and explain why each exists
- Understand the steps from alignment through sorting, indexing, and variant calling
- Navigate the UCSC Genome Browser to look up gene locations, exon structures, and genomic features
- Explain the difference between 0-based (BED) and 1-based (VCF/SAM) coordinate systems
- Interpret alignment quality metrics from samtools flagstat

## Key Vocabulary
| Term | Definition |
|---|---|
| Read alignment (mapping) | The process of determining where each short sequencing read originated on the reference genome by finding the best-matching position |
| Reference genome (hg38/GRCh38) | The current consensus human genome assembly (~3.2 billion base pairs) used as the coordinate framework for alignment |
| BWA (Burrows-Wheeler Aligner) | A fast and accurate tool for aligning short sequencing reads to a large reference genome |
| SAM (Sequence Alignment/Map) | A human-readable, tab-delimited text format where each line represents one read's alignment information |
| BAM | The binary, compressed version of SAM; much smaller file size and faster to process, but not human-readable |
| BAI (BAM index) | An index file for a sorted BAM that enables rapid random access to any genomic region without reading the entire file |
| CIGAR string | A compact representation of how a read aligns to the reference (e.g., 150M = 150 matching bases, 5M1I144M = 5 match, 1 insertion, 144 match) |
| Mapping quality (MAPQ) | A Phred-scaled confidence score for an alignment; MAPQ 60 means the aligner is highly confident the read maps to the correct position |
| Variant calling | The computational process of identifying positions where the sequenced individual differs from the reference genome |
| GATK HaplotypeCaller | A sophisticated variant caller from the Broad Institute that performs local re-assembly to accurately call SNPs and short indels |
| VCF (Variant Call Format) | The standard file format for recording genetic variants with position, alleles, quality, and genotype information |
| BED format | A tab-delimited format using 0-based, half-open coordinates to describe genomic intervals; commonly used for exome target regions |

## Concept Summaries

### Read Alignment: Mapping Reads to the Genome
After trimming adapters on Day 1, you now have clean FASTQ reads that need to be placed on the reference genome. This is the alignment step, performed by BWA (Burrows-Wheeler Aligner). BWA takes each read and searches the reference genome (all 3.2 billion bases) to find the position where that read matches best. For paired-end data, BWA uses both reads from a fragment simultaneously, leveraging the known insert-size range to improve mapping accuracy.

The output is a SAM file, where each line represents one read and includes: the read name, a bitwise FLAG encoding various properties (is it paired? did it map? is it first or second in the pair?), the chromosome and position where it mapped, the mapping quality score, and the CIGAR string describing the alignment pattern. Because SAM files are enormous (often 20-50 GB for WES), we immediately convert them to BAM (binary compressed), which reduces file size by 3-5x while retaining all information.

### Sorting, Indexing, and the Pipeline
After alignment, the BAM file must be sorted by genomic coordinate. Sorting is required because downstream tools (variant callers, genome browsers) need to efficiently access reads in a specific region. Without sorting, finding all reads at a particular position would require scanning the entire file. After sorting, we create an index (.bai) that acts like a book's table of contents, allowing tools to jump directly to any chromosomal region.

The complete pipeline on Day 2 is: BWA alignment (FASTQ --> SAM) --> SAM to BAM conversion (samtools view) --> coordinate sorting (samtools sort) --> indexing (samtools index) --> variant calling (GATK HaplotypeCaller, which produces the VCF file). Each step depends on the previous one completing successfully, which is why we chain them in a single job script.

### The UCSC Genome Browser
The UCSC Genome Browser is a web-based tool for visualizing the human genome and its annotations. You can search for any gene by name (e.g., TAS2R38) or navigate to specific coordinates (e.g., chr7:141,972,631-141,973,841). The browser displays multiple data tracks stacked vertically: gene models showing exon-intron structure, conservation across species, known SNPs from dbSNP, regulatory elements, and more.

Understanding coordinates is critical. The browser displays positions in 1-based coordinates (the first base of a chromosome is position 1). However, the underlying BED format used for many tracks is 0-based and half-open: the first base is position 0, and an interval [start, end) includes the start but not the end. This discrepancy is one of the most common sources of off-by-one errors in bioinformatics. VCF and SAM use 1-based coordinates; BED uses 0-based. When converting between formats, always add 1 to the BED start to get the VCF position.

## How It Connects
- **Previous Day (Day 1)**: You trimmed adapters from raw FASTQ reads using cutadapt. Those clean reads are the input for today's alignment. If adapter trimming was incomplete, unmapped reads or misalignments will increase.
- **Next Day (Day 3)**: Tomorrow you will open the VCF file produced today and begin identifying specific variants, starting with the TAS2R38 bitter taste gene. Today's pipeline generates the raw material (the VCF) that makes variant exploration possible.

## Review Questions
1. Why can't we simply compare raw sequencing reads directly to each other to find variants? Why do we need a reference genome?
2. Explain the relationship between SAM, BAM, and BAI files. Why do we need all three?
3. What does a CIGAR string of `100M1D50M` tell you about how a read aligns to the reference?
4. You check `samtools flagstat` and see that only 80% of reads mapped. What might be wrong, and what steps would you take to investigate?
5. Explain the difference between 0-based and 1-based coordinate systems. A variant is at BED position [99, 100) -- what is its position in VCF format?
6. In the UCSC Genome Browser, you search for TAS2R38. How many exons does it have, and why is this unusual for a human gene?
7. Why is sorting a BAM file necessary before creating an index?
8. Describe what GATK HaplotypeCaller does and why it is more accurate than simply counting mismatches.

## Further Reading
- BWA manual: https://bio-bwa.sourceforge.net/bwa.shtml
- SAM format specification: https://samtools.github.io/hts-specs/SAMv1.pdf
- UCSC Genome Browser user guide: https://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html
- GATK Best Practices: https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows
- Coordinate systems in bioinformatics (Genome Browser blog): https://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/
