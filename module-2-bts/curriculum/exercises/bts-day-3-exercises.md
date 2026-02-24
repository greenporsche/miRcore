# BTS Day 3 Exercises: SAM/BAM Files and Alignment

## Warm-Up Questions
1. What is the purpose of aligning sequencing reads to a reference genome?
2. What is the difference between a SAM file and a BAM file?
3. Why is it important to sort a BAM file before indexing it?

## Hands-On Exercises

### Exercise 1: Running Bowtie2 Alignment
**Objective**: Align sequencing reads to the SARS-CoV-2 reference genome using Bowtie2.
**Instructions**:
1. Log into Great Lakes and navigate to your working directory
2. Load the required modules:
   ```bash
   module load Bioinformatics
   module load Bowtie2
   ```
3. Run Bowtie2 to align reads against the SARS-CoV-2 reference index:
   ```bash
   bowtie2 -x /path/to/sarscov2_index -U small.fastq -S small.sam
   ```
4. Observe the alignment summary printed to the screen --- record the overall alignment rate
5. Examine the first few lines of the output SAM file: `head -n 20 small.sam`

**Expected Output**: Bowtie2 will print an alignment summary showing total reads processed and alignment rate (e.g., "95.5% overall alignment rate"). The SAM file will contain header lines starting with `@` followed by alignment records.

### Exercise 2: Understanding SAM Format
**Objective**: Parse and interpret the fields in a SAM file.
**Instructions**:
1. View the first non-header alignment line in your SAM file:
   ```bash
   grep -v "^@" small.sam | head -n 1
   ```
2. The line is tab-separated. Identify these fields:
   - Column 1 (QNAME): Read name
   - Column 2 (FLAG): Alignment flag number
   - Column 3 (RNAME): Reference name the read mapped to
   - Column 4 (POS): Position on the reference where the read aligned
   - Column 5 (MAPQ): Mapping quality score
   - Column 6 (CIGAR): How the read aligns (e.g., `150M` means 150 matching bases)
3. Use `cut` to extract specific columns:
   ```bash
   grep -v "^@" small.sam | cut -f1,3,4,5 | head -n 10
   ```
4. Record the reference name and position for 3 different reads

**Expected Output**: A table showing read name, reference name (e.g., NC_045512.2 for SARS-CoV-2), position, and mapping quality for several reads.

### Exercise 3: SAM to BAM Pipeline
**Objective**: Convert a SAM file to a sorted, indexed BAM file using SAMtools.
**Instructions**:
1. Load SAMtools: `module load SAMtools`
2. Convert SAM to BAM:
   ```bash
   samtools view -bS small.sam > small.bam
   ```
3. Sort the BAM file:
   ```bash
   samtools sort small.bam -o small.sorted.bam
   ```
4. Index the sorted BAM:
   ```bash
   samtools index small.sorted.bam
   ```
5. Verify all files were created:
   ```bash
   ls -lh small.*
   ```
6. Compare file sizes: How much smaller is the BAM compared to the SAM?

**Expected Output**: Four files: `small.sam`, `small.bam`, `small.sorted.bam`, `small.sorted.bam.bai`. The BAM file should be significantly smaller than the SAM file (typically 3-5x smaller).

### Exercise 4: Extracting Alignment Statistics
**Objective**: Use SAMtools and Linux commands to analyze alignment results.
**Instructions**:
1. Count total reads in the BAM file:
   ```bash
   samtools view -c small.sorted.bam
   ```
2. Count only mapped reads (exclude unmapped):
   ```bash
   samtools view -c -F 4 small.sorted.bam
   ```
3. Count unmapped reads:
   ```bash
   samtools view -c -f 4 small.sorted.bam
   ```
4. Calculate the mapping rate: mapped reads / total reads * 100
5. Extract reads that map to a specific region (e.g., the spike protein region):
   ```bash
   samtools view small.sorted.bam NC_045512.2:21563-25384 | wc -l
   ```

**Expected Output**: Total read count, mapped read count, unmapped read count, and the calculated mapping percentage. A count of reads mapping to the spike protein region.

### Exercise 5: Editing Files with nano
**Objective**: Practice using the nano text editor to create and modify files.
**Instructions**:
1. Create a new file: `nano my_notes.txt`
2. Type the following content:
   ```
   BTS Day 3 Notes
   ================
   Pipeline: FASTQ -> SAM -> BAM -> Sorted BAM -> Index
   Tool: Bowtie2 for alignment
   Tool: SAMtools for processing
   ```
3. Save the file: press `Ctrl+O`, then `Enter`
4. Exit nano: press `Ctrl+X`
5. Verify the contents: `cat my_notes.txt`

**Expected Output**: The file should display exactly what you typed.

## Challenge Problems

### Challenge 1: SAM Column Analysis
**Objective**: Investigate the distribution of mapping qualities in your alignment.
**Instructions**:
1. Extract only the MAPQ column (column 5) from the SAM file:
   ```bash
   grep -v "^@" small.sam | cut -f5 > mapq_scores.txt
   ```
2. Sort the scores and count unique values:
   ```bash
   sort mapq_scores.txt | uniq -c | sort -rn | head -n 10
   ```
3. Answer: What is the most common mapping quality score? What does MAPQ=0 mean? What does MAPQ=42 mean?

### Challenge 2: Multi-Step Pipeline Script
**Objective**: Write a complete alignment pipeline as a single script.
**Instructions**:
1. Open nano and create a file called `pipeline.sh`
2. Write the full pipeline:
   ```bash
   #!/bin/bash
   module load Bioinformatics
   module load Bowtie2
   module load SAMtools

   # Step 1: Align
   bowtie2 -x /path/to/sarscov2_index -U small.fastq -S small.sam

   # Step 2: Convert to BAM
   samtools view -bS small.sam > small.bam

   # Step 3: Sort
   samtools sort small.bam -o small.sorted.bam

   # Step 4: Index
   samtools index small.sorted.bam

   echo "Pipeline complete!"
   ```
3. Make it executable: `chmod +x pipeline.sh`
4. Run it: `bash pipeline.sh`

---

## Answer Key

### Warm-Up Answers
1. Aligning reads to a reference genome allows us to determine where each sequenced fragment came from in the genome. This is essential for detecting which organism the DNA belongs to (e.g., is it SARS-CoV-2 or human?) and for identifying mutations or variants.
2. SAM (Sequence Alignment/Map) is a human-readable text format --- you can open it with `cat` or `nano`. BAM (Binary Alignment/Map) is the compressed binary version that is much smaller but cannot be read directly as text. They contain the same information; BAM is just more efficient for storage and processing.
3. Sorting organizes the reads by their genomic position, which is required for indexing. The index allows tools like IGV to quickly jump to any genomic location without reading the entire file. Without sorting, the index cannot be built because it relies on reads being in positional order.

### Exercise Answers
1. **Bowtie2 Alignment**: The alignment summary shows the total reads, how many aligned exactly once, how many aligned more than once, and the overall alignment rate. A high alignment rate (>90%) to SARS-CoV-2 suggests the sample contains viral sequences.
2. **SAM Format**: Column 1 is the read name, column 3 is the reference genome name (NC_045512.2 for SARS-CoV-2), column 4 is the 1-based position, and column 5 is the MAPQ score. The CIGAR string (column 6) describes the alignment pattern.
3. **SAM to BAM Pipeline**: The correct order is: `samtools view -bS` (SAM to BAM) -> `samtools sort` (sort by position) -> `samtools index` (create .bai index). The BAM is typically 3-5x smaller than the SAM.
4. **Alignment Statistics**: Mapping rate = (mapped reads / total reads) * 100. The `-F 4` flag excludes unmapped reads (FLAG bit 4 = unmapped). The `-f 4` flag includes only unmapped reads.
5. **nano Editor**: Key shortcuts: `Ctrl+O` to save (write Out), `Ctrl+X` to exit, `Ctrl+K` to cut a line, `Ctrl+U` to paste.

### Challenge Answers
1. **MAPQ Analysis**: MAPQ=0 means the read maps equally well to multiple locations (ambiguous). MAPQ=42 is a high-confidence unique mapping. The most common score depends on the data, but for well-mapped reads, high MAPQ values should dominate.
2. **Pipeline Script**: The script should run all four steps sequentially. The `#!/bin/bash` line tells the system to use bash as the interpreter. Making the script executable with `chmod +x` allows you to run it directly. The `echo` at the end confirms completion.
