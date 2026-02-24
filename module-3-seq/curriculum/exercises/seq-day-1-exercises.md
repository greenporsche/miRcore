# SEQ Day 1 Exercises: Sequencing Fundamentals

## Warm-Up Questions
1. What is the difference between whole genome sequencing (WGS) and whole exome sequencing (WES)? Why might a researcher choose WES over WGS?
2. Explain in your own words why adapter sequences end up in sequencing reads, and why they need to be removed.
3. You have completed both the CB (Computational Biology) and BTS (Bench to Sequencer) camps. Name one thing you learned in each that directly connects to what we are doing in SEQ.

## Hands-On Exercises

### Exercise 1: Reading a FASTQ File
**Objective**: Understand the four-line structure of FASTQ format and interpret quality scores.
**Instructions**:
Given the following FASTQ record:
```
@SEQ_READ_001
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```
1. Identify what each of the four lines represents.
2. The quality character `F` corresponds to ASCII value 70. Using Phred+33 encoding, calculate the Phred quality score for each base. What is the error probability per base?
3. If a base had a quality character of `5` (ASCII 53), what would its Phred score and error probability be?

**Expected Output**: Line-by-line explanation; Phred score of 37 for `F` (error rate ~0.02%); Phred score of 20 for `5` (error rate 1%).

### Exercise 2: Connecting to Great Lakes
**Objective**: Practice logging into the HPC cluster and navigating the file system.
**Instructions**:
1. Connect to the U-M VPN.
2. Open a terminal and SSH into Great Lakes: `ssh uniqname@greatlakes.arc-ts.umich.edu`
3. Navigate to your project directory: `cd /scratch/miRcore_project/[your_uniqname]`
4. Run `pwd` to confirm your location.
5. Run `ls -la` to list all files, including hidden ones.
6. Record: How many files are in your directory? What are their names and sizes?

**Expected Output**: A list of files including your raw FASTQ files (R1 and R2 paired-end reads).

### Exercise 3: Examining Your Raw Data
**Objective**: Inspect your actual sequencing data and understand paired-end reads.
**Instructions**:
1. From your project directory, view the first 8 lines of your R1 FASTQ file:
   ```
   head -8 your_R1.fq.gz | zcat
   ```
   (If the file is gzipped, use `zcat` first: `zcat your_R1.fq.gz | head -8`)
2. Do the same for your R2 file.
3. Compare the read headers between R1 and R2. What do you notice about the read names?
4. Count the total number of reads in your R1 file:
   ```
   zcat your_R1.fq.gz | wc -l
   ```
   Divide by 4 to get the number of reads. Record this number.

**Expected Output**: R1 and R2 headers should match (same read name, different /1 and /2 suffix); total read count typically in the millions.

### Exercise 4: Running Cutadapt
**Objective**: Edit a job script, submit an adapter trimming job, and verify the output.
**Instructions**:
1. Open the cutadapt job script in nano: `nano job1_cutadapt.sh`
2. Verify the following are correctly set:
   - Your uniqname replaces ALL placeholder `[UID]` entries (including brackets)
   - Input file paths match your actual FASTQ filenames
   - The adapter sequence matches the Illumina universal adapter
3. Save the file (Ctrl+O, Enter, Ctrl+X).
4. Submit the job: `sbatch job1_cutadapt.sh`
5. Check your job status: `squeue -u [your_uniqname]`
6. After the job completes, check your email for the completion notification.
7. Examine the log file to confirm success.

**Expected Output**: A successful job completion email; trimmed output FASTQ files in your output directory.

### Exercise 5: Cutadapt Summary Analysis
**Objective**: Interpret the cutadapt summary report to assess data quality.
**Instructions**:
After cutadapt completes, examine the summary output (in the log file or stdout) and answer:
1. How many read pairs were processed?
2. What percentage of reads had adapters trimmed?
3. How many read pairs passed the quality filter?
4. What was the total number of base pairs before and after trimming?
5. Based on these numbers, would you say your data is high quality? Why or why not?

Fill out the Google Form provided by your group leader with your cutadapt results.

**Expected Output**: A completed Google Form entry; written answers showing understanding of trimming metrics.

## Challenge Problems

### Challenge 1: Quality Score Detective
**Objective**: Explore the relationship between sequencing technology and quality scores.
**Instructions**:
Research and answer the following:
1. Illumina sequencers typically produce reads of 150-300 bp. Why do quality scores tend to decrease toward the 3' end of reads?
2. A FASTQ file has average quality scores of Q15 across all reads. Would you proceed with analysis? What cutoff would you recommend and why?
3. Explain why paired-end sequencing gives better results than single-end for exome sequencing. Think about what happens when a read maps to a repetitive region.

### Challenge 2: Script Debugging
**Objective**: Practice troubleshooting common SLURM script errors.
**Instructions**:
The following job script has THREE errors. Find and fix all of them:
```bash
#!/bin/bash
#SBATCH --job-name=cutadapt_run
#SBATCH --output=/scratch/miRcore_project/[UID]/output/job1_%j.log
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

module load Bioinformatics
module load cutadapt

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
-o /scratch/miRcore_project/[UID]/output/trimmed_R1.fq
-p /scratch/miRcore_project/[UID]/output/trimmed_R2.fq
/scratch/miRcore_project/[UID]/raw/sample_R1.fq.gz /scratch/miRcore_project/[UID]/raw/sample_R2.fq.gz
```

---

## Answer Key

### Warm-Up Answers
1. WGS sequences the entire genome (~3.2 billion bases); WES sequences only exons (~1-2% of the genome, ~30 million bases). WES is chosen when the focus is on protein-coding variants because it is cheaper, faster, and produces smaller data files while capturing the most medically relevant regions.
2. During library preparation, short synthetic DNA sequences (adapters) are ligated to both ends of DNA fragments. When a fragment is shorter than the read length, the sequencer reads through the insert and into the adapter on the other end. These adapter sequences are not part of the person's genome and must be removed to avoid false alignments.
3. Answers will vary. Example: In CB, we learned about DNA structure, codons, and how mutations change amino acids. In BTS, we isolated DNA and learned about PCR amplification. Both connect to SEQ because we now analyze the digital output of sequencing our actual DNA.

### Exercise Answers
1. **Exercise 1**: Line 1 = read identifier (starts with @); Line 2 = nucleotide sequence; Line 3 = separator (+); Line 4 = quality scores (one character per base). Phred score for `F` = 70 - 33 = 37, error probability = 10^(-37/10) = ~0.02%. Phred score for `5` = 53 - 33 = 20, error probability = 10^(-20/10) = 1%.
2. **Exercise 2**: Students should successfully log in and list files. Typical contents include R1/R2 FASTQ files and possibly job scripts.
3. **Exercise 3**: R1 and R2 read headers share the same read name identifier (indicating they come from the same DNA fragment). Read count will vary per student, typically 20-50 million paired reads for WES.
4. **Exercise 4**: Students should verify: (a) all `[UID]` placeholders replaced, (b) cutadapt command is on a single line or uses proper line continuation (`\`), (c) job completes successfully per email and log.
5. **Exercise 5**: Answers depend on individual data. High quality is typically >90% reads passing filters and >80% of bases retained after trimming.

### Challenge Answers
1. **Challenge 1**: (1) Quality decreases at 3' ends because Illumina sequencing-by-synthesis chemistry has increasing error rates with each cycle due to signal decay and phasing effects. (2) Q15 is generally too low for reliable variant calling; Q20 (1% error) is a common minimum, and Q30 (0.1% error) is preferred. (3) Paired-end reads improve mapping in repetitive regions because even if one end maps ambiguously, the paired constraint (known insert size) helps the aligner determine the correct location.
2. **Challenge 2**: Three errors: (1) `[UID]` is never replaced with actual username -- appears four times. (2) The cutadapt command is broken across multiple lines without backslash (`\`) continuation characters, so the shell treats each line as a separate command. (3) The entire cutadapt command with all flags must be on one line or properly continued with `\` at the end of each line.
