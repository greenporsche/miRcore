# BTS Day 2 Exercises: NGS and FASTQ Analysis

## Warm-Up Questions
1. What does NGS stand for, and why is it called "next-generation" compared to Sanger sequencing?
2. Explain in your own words why GC content affects the melting temperature of DNA.
3. How many lines does each individual read occupy in a FASTQ file, and what information does each line contain?

## Hands-On Exercises

### Exercise 1: Counting Reads in a FASTQ File
**Objective**: Determine the number of sequencing reads in a FASTQ file using Linux commands.
**Instructions**:
1. On Great Lakes, navigate to the directory containing the sample FASTQ file (e.g., `small.fastq`)
2. Count the total number of lines: `wc -l small.fastq`
3. Divide the total line count by 4 to get the number of reads
4. Verify your count using: `grep -c "^@" small.fastq`
5. Record both numbers --- do they match?

**Expected Output**: The line count divided by 4 should equal the `grep -c` count. For example, 40,000 lines = 10,000 reads.

### Exercise 2: Examining FASTQ Format
**Objective**: Understand the 4-line structure of a FASTQ file by manually inspecting reads.
**Instructions**:
1. Display the first 8 lines of the FASTQ file (2 reads): `head -n 8 small.fastq`
2. For each read, identify:
   - Line 1: The read identifier (starts with `@`)
   - Line 2: The DNA sequence
   - Line 3: The separator line (`+`)
   - Line 4: The quality scores (ASCII characters)
3. Count the number of bases in the sequence line
4. Verify that the quality score line has the same number of characters as the sequence line
5. Look up what the first quality character means using a Phred score table

**Expected Output**: A written breakdown of 2 reads showing all 4 components, with matching sequence and quality string lengths.

### Exercise 3: GC Content Calculation
**Objective**: Calculate the GC content of a DNA sequence using the command line.
**Instructions**:
1. Take the sequence from one read in your FASTQ file
2. Count the total number of bases (length of the sequence)
3. Count the number of G and C bases using:
   ```bash
   echo "ATCGGCTAGCTA" | grep -o "[GC]" | wc -l
   ```
   (Replace the sequence with your actual read sequence)
4. Calculate GC content: (number of G+C) / (total bases) * 100
5. Is your GC content above or below 50%? What might that suggest?

**Expected Output**: A GC percentage for the selected read. SARS-CoV-2 has a relatively low GC content (~38%), so reads from this virus should reflect that.

### Exercise 4: Pattern Searching with grep
**Objective**: Use `grep` to search for specific patterns in a FASTQ file.
**Instructions**:
1. Search for reads containing the motif "ATCGATCG": `grep "ATCGATCG" small.fastq`
2. Count how many sequence lines contain "AAAA" (poly-A stretch): `grep -c "AAAA" small.fastq`
3. Search for a specific read by its identifier: `grep "@READ_NAME" small.fastq`
4. Find all header lines: `grep "^@" small.fastq | head -n 5`
5. Record how many results each search returned

**Expected Output**: Varying counts depending on the FASTQ file. The `^@` search should return one result per read.

### Exercise 5: File Transfer with scp
**Objective**: Transfer a file from the Great Lakes cluster to your local computer.
**Instructions**:
1. On Great Lakes, note the full path to a small file you want to transfer (use `pwd` and `ls`)
2. Open a **new** terminal window on your **local** machine (do not use the Great Lakes terminal)
3. Run: `scp uniqname@greatlakes.arc-ts.umich.edu:~/path/to/file.txt ~/Desktop/`
4. Authenticate with your password and Duo
5. Verify the file arrived on your Desktop

**Expected Output**: The file should appear on your local Desktop with the same contents as on Great Lakes.

## Challenge Problems

### Challenge 1: FASTQ Quality Analysis
**Objective**: Write a series of commands to analyze the quality distribution of a FASTQ file.
**Instructions**:
1. Extract only the quality score lines (every 4th line, starting at line 4) from a FASTQ file:
   ```bash
   sed -n '4~4p' small.fastq > quality_lines.txt
   ```
2. Count the total number of quality lines: `wc -l quality_lines.txt`
3. Search for lines containing low-quality indicators (e.g., `!` which is Phred score 0, or `#` which is Phred score 2):
   ```bash
   grep -c '!' quality_lines.txt
   ```
4. What percentage of reads contain at least one very low-quality base?

### Challenge 2: Business Pitch Research
**Objective**: Prepare the foundation for your group's biotech business pitch.
**Instructions**:
1. Choose a biotech company or technology (e.g., Illumina sequencing, CRISPR therapeutics, rapid COVID testing)
2. Research and write answers to the following:
   - What problem does this company/technology solve?
   - Who is the target market?
   - What is the competitive advantage?
   - What would a 2-minute elevator pitch sound like?
3. Prepare 3 slides: Problem, Solution, and Market Opportunity

---

## Answer Key

### Warm-Up Answers
1. NGS stands for Next-Generation Sequencing. It is called "next-generation" because it can sequence millions of DNA fragments simultaneously (massively parallel), unlike Sanger sequencing which sequences one fragment at a time. This makes NGS much faster and cheaper for large-scale genomic studies.
2. GC base pairs are held together by 3 hydrogen bonds, while AT base pairs have only 2 hydrogen bonds. More hydrogen bonds means more energy (higher temperature) is required to separate the strands. Therefore, DNA with higher GC content has a higher melting temperature because it takes more heat to break apart those stronger bonds.
3. Each read occupies exactly 4 lines: (1) Header line starting with `@` containing the read ID, (2) DNA sequence, (3) Separator line with `+`, (4) Quality scores as ASCII characters (one character per base).

### Exercise Answers
1. **Counting Reads**: If `wc -l` returns 40,000 lines, there are 40,000 / 4 = 10,000 reads. The `grep -c "^@"` command should also return 10,000 (or close to it --- note that `@` can also appear in quality lines, so `grep -c "^@"` may slightly overcount in rare cases).
2. **FASTQ Format**: Each read should have exactly 4 lines. The sequence length and quality string length must be identical (one quality character per base). For example, a 150-base read will have 150 characters on both the sequence and quality lines.
3. **GC Content**: For SARS-CoV-2 reads, expect GC content around 38%. The calculation is: count of (G + C) bases divided by total bases, multiplied by 100.
4. **grep Searching**: Results vary by file. The `^@` pattern matches lines starting with `@`, which are header lines. Note: `grep "ATCGATCG"` will search all lines, not just sequence lines, so some matches could theoretically be in quality lines (unlikely but possible).
5. **scp Transfer**: The key is to run `scp` from your **local** terminal, not from the Great Lakes session. The syntax is `scp remote:path local_path`.

### Challenge Answers
1. **Quality Analysis**: The `sed -n '4~4p'` command extracts every 4th line starting from line 4 (the quality lines). The percentage of reads with low-quality bases = (lines matching `!`) / (total quality lines) * 100. A high percentage might indicate sequencing issues.
2. **Business Pitch**: Answers will vary. A strong pitch should clearly identify the problem (e.g., "genetic diseases affect millions"), the solution (e.g., "CRISPR-based gene therapy"), the market size, and a competitive advantage (e.g., "first FDA-approved CRISPR treatment").
