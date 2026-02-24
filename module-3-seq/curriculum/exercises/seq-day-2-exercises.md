# SEQ Day 2 Exercises: Exome Data Processing

## Warm-Up Questions
1. After adapter trimming on Day 1, what format are your cleaned reads in? What is the next step in the bioinformatics pipeline?
2. Why do we need a reference genome to analyze exome sequencing data? What reference assembly does our pipeline use?
3. Explain the difference between a SAM file and a BAM file. Why do we convert from one to the other?

## Hands-On Exercises

### Exercise 1: Understanding the Pipeline
**Objective**: Trace the complete data processing pipeline from trimmed reads to variant calls.
**Instructions**:
Draw or write out the complete pipeline in order, listing:
1. Each step name
2. The tool used for that step
3. The input file format
4. The output file format

The steps are: Adapter trimming, Read alignment, SAM-to-BAM conversion, BAM sorting, BAM indexing, Variant calling.

**Expected Output**:
```
cutadapt:    FASTQ --> trimmed FASTQ
bwa mem:     trimmed FASTQ + reference --> SAM
samtools view: SAM --> BAM
samtools sort: unsorted BAM --> sorted BAM
samtools index: sorted BAM --> BAM + BAI
gatk HaplotypeCaller: sorted BAM + BAI + reference --> VCF
```

### Exercise 2: Submitting the Alignment Job
**Objective**: Run the BWA alignment and post-processing pipeline on your own data.
**Instructions**:
1. Open the mapping job script: `nano job2_mapping.sh`
2. Verify that:
   - Your uniqname replaces all `[UID]` placeholders
   - Input files point to your cutadapt output from Day 1
   - The reference genome path is correct: `/scratch/miRcore_project/reference/hg38.fa`
3. Submit the job: `sbatch job2_mapping.sh`
4. This job takes longer than cutadapt. Monitor with: `squeue -u [your_uniqname]`
5. After completion, check the email notification and log file.
6. List your output files with `ls -lh` and record the file sizes of your SAM, BAM, and sorted BAM files.

**Expected Output**: A sorted, indexed BAM file (and corresponding BAI file) in your output directory. BAM files are typically 2-5 GB for WES.

### Exercise 3: Checking Alignment Quality
**Objective**: Use samtools flagstat to assess how well your reads mapped to the reference genome.
**Instructions**:
1. Run: `samtools flagstat /scratch/miRcore_project/[your_uniqname]/output/sorted.bam`
2. Record the following from the output:
   - Total number of reads
   - Number (and percentage) of reads that mapped
   - Number of properly paired reads
   - Number of reads with a mapping quality of 0
3. Answer: Is your mapping rate acceptable? (Generally >95% is good for WES.)
4. What might cause reads to fail to map?

**Expected Output**: A flagstat report showing >95% mapping rate; understanding that unmapped reads may come from adapter remnants, contamination, or regions not in the reference.

### Exercise 4: Exploring the UCSC Genome Browser
**Objective**: Learn to navigate the UCSC Genome Browser and understand genomic coordinates.
**Instructions**:
1. Go to https://genome.ucsc.edu and select "Genome Browser."
2. Make sure the assembly is set to **GRCh38/hg38**.
3. In the search box, type `TAS2R38` and press Go.
4. Answer these questions:
   - On which chromosome is TAS2R38 located?
   - What are the start and end coordinates of the gene?
   - How many exons does TAS2R38 have?
   - What strand is the gene on (+ or -)?
5. Zoom in to exon 1 by clicking on the exon in the gene track. What is the approximate length of this exon?
6. Try searching for `ABO` gene. What chromosome is it on?

**Expected Output**: TAS2R38 is on chromosome 7 (approximately chr7:141,972,631-141,973,841 in hg38), has 1 coding exon, on the minus strand. ABO is on chromosome 9.

### Exercise 5: Coordinate System Comparison
**Objective**: Understand the critical difference between 0-based and 1-based coordinate systems.
**Instructions**:
Consider a short DNA sequence: `ATCGATCG` at the beginning of a chromosome.

1. In a **1-based coordinate system** (used by VCF, SAM, UCSC display):
   - What position is the first `A`?
   - What position is the `G` at position 4?
   - What are the coordinates for the substring `GATC` (middle four bases)?

2. In a **0-based, half-open coordinate system** (used by BED format):
   - What position is the first `A`?
   - What are the start and end coordinates for the substring `GATC`?
   - Explain what "half-open" means (start is inclusive, end is exclusive).

3. A VCF file lists a variant at position 100. What would the equivalent BED interval be?

**Expected Output**: 1-based: first A is at 1, G is at 4, GATC is positions 4-7. 0-based half-open: first A is at 0, GATC is [3, 7). VCF position 100 = BED interval [99, 100).

## Challenge Problems

### Challenge 1: SAM Format Deep Dive
**Objective**: Parse and interpret a real SAM alignment record.
**Instructions**:
Here is a SAM record (columns separated by tabs):
```
READ001  99  chr7  141972800  60  150M  =  141973000  350  ACGTACGT...  FFFFFFFF...
```
1. What is the read name?
2. What does FLAG 99 mean? (Hint: decode the bitwise flag -- use https://broadinstitute.github.io/picard/explain-flags.html)
3. What chromosome and position is this read aligned to?
4. What is the mapping quality, and what does it mean?
5. What does the CIGAR string `150M` mean?
6. What does `=` in the mate reference column indicate?
7. The insert size is 350. What does this represent?

### Challenge 2: File Size Analysis
**Objective**: Understand why different file formats have vastly different sizes.
**Instructions**:
After your pipeline runs, compare the file sizes:
1. Record the size of: raw FASTQ (R1+R2), SAM file, BAM file, sorted BAM, VCF file.
2. Calculate the compression ratio of BAM vs. SAM.
3. Why is the VCF file so much smaller than the BAM file?
4. If you had run whole genome sequencing instead of whole exome, estimate how much larger each file would be (hint: exome is ~1-2% of the genome).

---

## Answer Key

### Warm-Up Answers
1. Cleaned reads are in FASTQ format (trimmed FASTQ). The next step is aligning/mapping them to the human reference genome using BWA.
2. A reference genome provides the coordinate framework that tells us where each read belongs. Without it, we would have millions of unordered short sequences with no positional context. Our pipeline uses GRCh38/hg38.
3. SAM is a human-readable text format (very large). BAM is the binary compressed equivalent (much smaller, typically 3-5x compression). We convert to BAM to save disk space and enable faster processing by downstream tools.

### Exercise Answers
1. **Exercise 1**: See the pipeline diagram in Expected Output above. Key understanding: each step transforms data from one format to another, building toward variant calls.
2. **Exercise 2**: Students should produce sorted.bam and sorted.bam.bai files. File sizes vary by student but BAM is typically 2-5 GB for WES.
3. **Exercise 3**: Mapping rate should be >95%. Low mapping rates could indicate: adapter contamination remaining, low-quality reads, sample contamination from other organisms, or reads from regions not represented in hg38.
4. **Exercise 4**: TAS2R38 is on chr7, approximately chr7:141,972,631-141,973,841. It has 1 coding exon (unusual for a gene). It is on the minus (-) strand. ABO is on chr9.
5. **Exercise 5**: 1-based: A=1, G=4, GATC=4-7. 0-based half-open: A=0, GATC=[3,7). Half-open means the start position is included but the end position is not. VCF pos 100 = BED [99, 100). This is the most common source of coordinate errors in bioinformatics.

### Challenge Answers
1. **Challenge 1**: Read name = READ001. FLAG 99 = read paired (1) + read mapped in proper pair (2) + mate reverse strand (32) + first in pair (64) = paired, properly paired, mate on reverse strand, this is R1. Chromosome = chr7, position = 141,972,800. MAPQ = 60 (highest confidence, uniquely mapped). CIGAR `150M` = 150 bases all match/mismatch (no indels or clips). `=` means mate maps to the same chromosome. Insert size 350 = distance from leftmost mapped base of this read to rightmost mapped base of its mate.
2. **Challenge 2**: Typical compression BAM/SAM ratio is 3-5x. VCF is much smaller because it only records positions where the sample differs from the reference (typically 20,000-30,000 variants for WES), while BAM stores every individual read alignment. For WGS, files would be approximately 50-100x larger since the exome represents ~1-2% of the genome.
