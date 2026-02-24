# BTS Track (Biotech Startup) -- Detailed Analysis

> **Note:** The study guides, exercises, and cheat sheets in the `curriculum/` folder are my own self-generated study materials, not official miRcore content. See the [root README](../README.md) for details.

## Track Overview

The Biotechnology Startup (BTS) track is the second week of the miRcore summer camp curriculum. It targets high school students in grades 9 through 12, with most students having completed the CB track the previous week. BTS shifts from the database-driven analysis of CB to hands-on command-line bioinformatics and sequencing data processing. The track is co-led by Dr. Inhan Lee and Vina, with a team of group leaders providing breakout room support.

**Audience:** Students who have completed CB (recommended) or will take CB the following week. Some students enter BTS as their first miRcore experience, but the curriculum assumes familiarity with basic genomics concepts from CB.

**Goals:**
- Understand how DNA sequencing works (Illumina next-generation sequencing)
- Learn Linux command-line basics on a real high-performance computing cluster
- Execute a complete bioinformatics pipeline: FASTQ -> alignment -> SAM/BAM -> visualization
- Identify an unknown pathogen from patient sequencing data
- Develop and deliver a biotech startup business pitch

**Structure:** Like CB, each day has morning and afternoon sessions in hybrid format. BTS is distinctive in its dual-track nature: a research/technical track running alongside a business/entrepreneurship track. Students work in small groups throughout the week, culminating in both a research presentation and a business pitch on Day 5.

---

## Day-by-Day Breakdown

### Day 1: Biotech Introduction and Linux Basics

**Morning Session (Part 1)**

*Main Topics:*
- Welcome and camp context: Dr. Lee opens with genuine enthusiasm about teaching high school students to use the UM Linux cluster, setting an encouraging tone for the week
- Review of why students are here: learning to do real research, not school-style exercises
- Biotechnology overview: broad definition (using biological systems to develop products)
  - Insulin production: first biotech product, recombinant DNA technology in E. coli
  - CRISPR-Cas9: gene editing tool, molecular scissors, Nobel Prize 2020
  - PCR (Polymerase Chain Reaction): amplifying DNA, how it enabled COVID testing
  - Illumina sequencing: dominant NGS platform, sequencing by synthesis
- Wet lab vs. dry lab distinction: hands-on bench work vs. computational analysis
- COVID-19 biology deep dive:
  - Virus structure: lipid envelope, spike protein, RNA genome
  - How SARS-CoV-2 infects cells: spike protein binds ACE2 receptor
  - The SARS-CoV-2 genome: ~30,000 base pairs, positive-sense single-stranded RNA
  - Structural proteins: S (spike), E (envelope), M (membrane), N (nucleocapsid)
- NCBI exploration: finding the SARS-CoV-2 reference genome
  - Switching between text and graphics mode on NCBI
  - Locating the S-protein: genomic coordinates 21,563 to 25,384

*Key Concepts Introduced:*
- Next-generation sequencing technology
- Viral genome structure
- Reference sequences and genomic coordinates
- NCBI as a central genomics resource

*Hands-on Activities:*
- Exploring the SARS-CoV-2 genome on NCBI
- Finding the S-protein coordinates and first amino acid
- Group formation and initial team building

**Afternoon Session (Part 2)**

*Main Topics:*
- S-protein deep dive: amino acid sequence, protein structure
- PDB (Protein Data Bank): visualizing the 3D structure of the spike protein
- Introduction to the Great Lakes HPC cluster
  - Logging in with University of Michigan credentials
  - What is a computing cluster? Many computers working together
  - Difference between a login node and compute nodes
- First Linux commands:
  - `pwd` (print working directory): "Where am I?"
  - `ls` (list): "What files are here?"
  - `mkdir` (make directory): creating folders
  - `cd` (change directory): navigating the file system
  - `cp` (copy): duplicating files
- FASTA format: header line (starts with >) followed by sequence data
- File system concepts: home directory, paths, directory hierarchy

*Key Concepts Introduced:*
- High-performance computing and why bioinformatics needs it
- The Linux file system: directories, paths, navigation
- FASTA as the fundamental sequence file format
- Protein structure databases

*Tools/Databases Used:*
- **NCBI**: SARS-CoV-2 reference genome exploration
- **PDB**: 3D protein structure visualization
- **Great Lakes HPC**: First login and navigation
- **Linux CLI**: `pwd`, `ls`, `mkdir`, `cd`, `cp`

*Group Work:*
- Logging into the cluster for the first time (troubleshooting authentication issues)
- Practicing basic Linux commands in breakout rooms
- Creating project directories

---

### Day 2: NGS and FASTQ Analysis

**Morning Session (Part 1)**

*Main Topics:*
- Linux command review: practice with `pwd`, `ls`, `cd`, `mkdir`, `cp`
- Nucleotide hybridization: how complementary DNA/RNA strands bind
  - DNA-DNA duplex, DNA-RNA duplex, RNA-RNA duplex
  - Base pairing rules and hydrogen bonds
- Melting temperature (Tm): the temperature at which half the DNA duplexes are denatured
  - Factors affecting Tm: GC content (G-C has 3 hydrogen bonds vs. A-T's 2), length, salt concentration
- GC content: percentage of G and C bases in a sequence
  - Higher GC content = more stable duplex = higher melting temperature
  - Importance for sequencing and primer design
- Spectrometry: measuring DNA concentration using UV absorption (A260)
- Guest speaker: **Dr. Sally Choi**
  - Biotech industry career path
  - How academic research connects to biotech companies
  - Career advice for students interested in biotech

*Key Concepts Introduced:*
- Molecular basis of DNA hybridization (essential for understanding how sequencing works)
- GC content as a sequence quality metric
- Melting temperature and its practical significance
- Career pathways in biotechnology

**Afternoon Session (Part 2)**

*Main Topics:*
- Next-generation sequencing (NGS) technology in depth:
  - Library preparation: fragmenting DNA, adding adapters
  - Cluster generation: bridge amplification on a flow cell
  - Sequencing by synthesis: fluorescent nucleotides, imaging each cycle
  - How reads are produced: millions of short sequences (~150 bp each)
- FASTQ file format:
  - Four lines per read: header (@), sequence, separator (+), quality scores
  - Phred quality scores: encoding base-calling confidence as ASCII characters
  - Q30 = 1 in 1,000 chance of error; Q40 = 1 in 10,000
- File transfer on the cluster: moving and copying data files
- `grep` for pattern matching: searching for specific sequences or headers within FASTQ files
  - Counting reads: `grep -c "@"` to count headers
  - Finding specific patterns in sequence data
- GC content analysis of sequencing data: practical exercise
- Business pitch fundamentals:
  - Problem/solution framing: identify a real health problem, propose a biotech solution
  - Target market identification
  - Value proposition: what makes your solution unique
  - Introduction to the pitch format for Day 3 and Day 5 presentations

*Key Concepts Introduced:*
- The complete NGS workflow from sample to data
- FASTQ format as the raw data format of sequencing
- Quality scores and their meaning
- `grep` as a powerful pattern-matching tool
- Business pitch structure

*Tools/Databases Used:*
- **Great Lakes HPC**: File operations, `grep`, data exploration
- **Linux CLI**: `grep`, file handling, data analysis

*Group Work:*
- Exploring FASTQ files on the cluster
- Counting reads and analyzing GC content
- Beginning to formulate business pitch ideas in breakout rooms
- Cheat sheet updates with new commands and concepts

---

### Day 3: SAM/BAM Files and Alignment

**Morning Session (Part 1)**

*Main Topics:*
- Student business pitch presentations: groups present their first-draft biotech startup pitches
  - Example pitches observed in transcript: CRISPR Therapeutics (gene-editing therapies for sickle cell disease, beta-thalassemia, cancer immunotherapy, HIV), rapid diagnostic testing startup
  - Pitches include: company background, technology description, target market, competitive advantage
  - Peer feedback and discussion after each pitch
- Quiz review: reinforcing key bioinformatics concepts from Days 1--2
  - DNA structure, base pairing, FASTQ format, Linux commands
  - Addressing common misconceptions

*Key Concepts Reinforced:*
- Biotechnology as a bridge between science and business
- Presentation and public speaking skills
- Peer feedback as a learning tool

**Afternoon Session (Part 2)**

*Main Topics:*
- Introduction to sequence alignment: the core bioinformatics operation
  - Why align? To determine where sequencing reads originated in a reference genome
  - The alignment problem: finding the best match for millions of short reads
- Bowtie2: a fast short-read aligner
  - Reference genome indexing: why a genome must be indexed before alignment
  - Running Bowtie2: `bowtie2 -x <index> -U <reads.fastq> -S <output.sam>`
  - Alignment rates: what percentage of reads map to the reference
  - Interpreting alignment summary statistics
- SAM (Sequence Alignment/Map) format:
  - Header lines (start with @) and alignment lines
  - Key fields: QNAME (query name), FLAG (bitwise flag), RNAME (reference name), POS (position), MAPQ (mapping quality), CIGAR (alignment description)
  - Understanding the FLAG field: mapped/unmapped, strand, mate information
  - CIGAR strings: M (match), I (insertion), D (deletion), S (soft clip)
- SAMtools: the essential toolkit for alignment files
  - `samtools view`: converting SAM to BAM (binary format, smaller and faster)
  - `samtools sort`: sorting BAM files by genomic position (required for visualization)
  - `samtools index`: creating a BAI index file for fast random access
  - The full processing chain: SAM -> BAM -> sorted BAM -> indexed BAM
- Additional Linux tools:
  - `nano`: terminal text editor for editing files and scripts
  - `cut`: extracting specific columns from tabular data
  - `sort`: sorting data by specific fields

*Key Concepts Introduced:*
- Sequence alignment as the bridge between raw reads and biological meaning
- SAM/BAM format as the standard for alignment data
- The SAMtools processing pipeline
- Why sorting and indexing are necessary for downstream analysis

*Tools/Databases Used:*
- **Bowtie2**: Aligning reads to the SARS-CoV-2 reference genome
- **SAMtools**: `view`, `sort`, `index`
- **Linux CLI**: `nano`, `cut`, `sort`

*Group Work:*
- Running Bowtie2 alignment on example SARS-CoV-2 data
- Processing SAM files through the SAMtools pipeline
- Examining alignment fields and understanding the output
- Troubleshooting alignment issues with group leaders

---

### Day 4: Genome Mapping and IGV

**Morning Session (Part 1)**

*Main Topics:*
- Review of Bowtie2 and SAMtools workflow, building on the previous day's discovery of SARS-CoV-2 sequences
- Loading bioinformatics modules on the cluster: `module load` for Bowtie2 and SAMtools
- IGV (Integrative Genomics Viewer) introduction:
  - Desktop application for visualizing aligned sequencing data
  - Loading reference genomes and BAM files
  - Navigating to genomic regions: zoom, scroll, search by coordinate or gene name
  - Reading the visualization: coverage track (histogram of read depth), alignment track (individual reads), reference track
  - Color coding: gray = normal mapping, colored bases = mismatches to reference
  - Identifying mutations: positions where multiple reads show the same non-reference base
- Shell scripting for batch job submission:
  - Writing `.sh` scripts with SLURM directives (`#SBATCH` headers)
  - Job parameters: partition, time, memory, CPU cores, email notifications
  - `sbatch` command: submitting scripts to the job scheduler
  - Monitoring jobs: email notifications (begin, end, fail), resource usage reports
  - Converting interactive work to batch jobs -- once a command is confirmed to work interactively, it gets packaged into a SLURM script for submission
- Mapping reads to the human genome (HG38 reference):
  - Why HG38? The current human reference genome assembly
  - Much larger than viral genomes: requires more computational resources
  - Submitting as a batch job because it takes longer to run

*Key Concepts Introduced:*
- Genome visualization as a way to inspect alignment quality and find variants
- Batch computing: submitting long-running jobs to a scheduler
- Human reference genome and its scale compared to viral genomes
- The concept of read depth/coverage

**Afternoon Session (Part 2)**

*Main Topics:*
- Shell scripting deep dive: modifying scripts for different inputs
  - Changing file names in scripts to process different samples
  - Using scripts as templates for batch processing
- Project data analysis: the 8 patient samples
  - Context: 8 unknown patient samples from a mystery outbreak
  - Task: align each sample to the SARS-CoV-2 reference genome
  - Diagnostic approach: samples with high alignment rates contain SARS-CoV-2
  - Students discover which patients are infected by examining alignment percentages
- SARS-CoV-2 identification from alignment rates:
  - High alignment rate (e.g., >80%) = sample contains the virus
  - Low alignment rate = sample does not contain the virus (or contains a different pathogen)
- Mutation identification:
  - Using IGV to visualize aligned patient reads against the reference
  - Identifying positions where patient reads consistently differ from the reference
  - Connecting mutations to viral evolution and variants of concern
- Human genome mapping results:
  - Examining what happens when you align the same reads to the human genome
  - Reads from SARS-CoV-2 should NOT align well to human genome (confirming the pathogen is viral, not human)

*Key Concepts Introduced:*
- Using alignment rates as a diagnostic tool
- Pathogen identification through comparative genomics
- Viral mutation detection from sequencing data
- The logic of differential alignment (viral reads align to viral genome, not human genome)

*Tools/Databases Used:*
- **Bowtie2**: Alignment to both SARS-CoV-2 and HG38 reference genomes
- **SAMtools**: Processing all 8 patient samples
- **IGV**: Visualizing aligned data, identifying mutations
- **Shell scripting + sbatch**: Automating the pipeline for multiple samples

*Group Work:*
- Processing multiple patient samples through the pipeline
- Comparing alignment rates across samples
- Visualizing results in IGV
- Discussing findings in breakout rooms: which patients are infected?

---

### Day 5: Research Projects and Business Pitches

**Morning Session (Part 1)**

*Main Topics:*
- Camp recap: full pipeline review from Day 1 to Day 4
  - Day 1: Biotech intro, Linux basics, SARS-CoV-2 genome
  - Day 2: NGS concepts, FASTQ format, GC content
  - Day 3: Alignment with Bowtie2, SAM/BAM processing
  - Day 4: IGV visualization, patient sample analysis, pathogen identification
- Research project continuation:
  - Finalizing analysis of 8 patient samples
  - Preparing results for presentation
  - Groups compile alignment rates and mutation findings
- Business pitch preparation:
  - Target market and audience identification: who benefits from your biotech solution?
  - Revenue model: how does the company make money?
  - Refining the pitch based on Day 3 feedback
  - Practice delivery and timing

*Group Work:*
- Final research analysis and slide creation
- Business pitch rehearsal
- Cheat sheet completion

**Afternoon Session (Part 2)**

*Main Topics:*
- Final business pitch presentations with investment simulation:
  - Each group presents their biotech startup pitch (5--8 minutes)
  - Business pitches observed in sessions include: CRISPR-based therapeutics, rapid diagnostic testing, genomic data platforms
  - After all pitches: investment simulation where students allocate virtual investment dollars to the startups they find most compelling
  - Discussion of why certain pitches were more convincing (clear problem, viable solution, realistic market)
- Research presentations:
  - Groups present their pathogen identification findings
  - Presentation structure: project background, pipeline description, alignment results across 8 patients, mutation analysis, conclusions
  - Students explain the full FASTQ -> Bowtie2 -> SAM -> BAM -> IGV pipeline
  - Peer questions and feedback
- Parent presentation:
  - Dr. Lee presents the miRcore vision and program overview
  - Camp outcomes and student accomplishments
  - Volunteer program and continuing research opportunities
- Camp wrap-up:
  - Awards for cheat sheets and presentations
  - Information about miRcore volunteer program (MVP), GIDAS clubs, and future camps
  - Birthday celebrations and final community moments

---

## Learning Arc

```
Day 1: BIOLOGY      -> What is biotech? How do viruses work? Where is the data?
Day 2: SEQUENCING   -> How is DNA sequenced? What does raw data look like?
Day 3: ALIGNMENT    -> How do we map reads to a genome? What does the output mean?
Day 4: VISUALIZATION -> How do we see and interpret the results? Can we find mutations?
Day 5: COMMUNICATION -> How do we present research + pitch a business?
```

The BTS learning arc has a distinctive dual thread. The technical thread progresses from biology through computation to visualization. The entrepreneurship thread runs in parallel, progressing from understanding the biotech landscape (Day 1) through pitch fundamentals (Day 2) to practice pitches (Day 3) to final investment-ready pitches (Day 5).

---

## Technical Pipeline: FASTQ -> Bowtie2 -> SAM -> BAM -> IGV

The core technical achievement of the BTS track is that students execute a complete bioinformatics pipeline from raw sequencing data to visual interpretation. Here is each step in detail:

### Step 1: FASTQ (Raw Sequencing Data)
- **What it is:** The raw output of an Illumina sequencer
- **Format:** Four lines per read -- header, sequence, separator, quality scores
- **Example:**
  ```
  @SRR12345.1
  ATCGATCGATCGATCGATCG
  +
  IIIIIIIIIIIIIIIIIIIII
  ```
- **What students do:** Examine files with `head`, count reads with `grep -c "@"`, analyze GC content

### Step 2: Bowtie2 (Sequence Alignment)
- **What it does:** Maps each short read to the best-matching position in a reference genome
- **Command pattern:** `bowtie2 -x <genome_index> -U <reads.fastq> -S <output.sam>`
- **Key output:** Alignment rate percentage (critical for pathogen identification)
- **What students learn:** The concept of aligning millions of reads, interpreting alignment statistics

### Step 3: SAM (Sequence Alignment/Map)
- **What it is:** Human-readable text file describing where each read aligned
- **Key fields students learn:** QNAME, FLAG, RNAME, POS, MAPQ, CIGAR
- **What students do:** Examine SAM files with `head`, extract columns with `cut`, understand the flag field

### Step 4: BAM (Binary Alignment/Map)
- **What it is:** Compressed binary version of SAM (much smaller, faster to process)
- **Conversion:** `samtools view -bS input.sam > output.bam`
- **Sorting:** `samtools sort input.bam -o output.sorted.bam`
- **Indexing:** `samtools index output.sorted.bam` (creates .bai file)
- **Why sort and index:** Required for efficient random access and visualization

### Step 5: IGV (Visualization)
- **What it does:** Displays aligned reads against a reference genome in an interactive graphical interface
- **What students see:** Coverage histogram showing read depth, individual aligned reads with mismatches highlighted in color
- **What students find:** Mutations in SARS-CoV-2 sequences, regions of high/low coverage
- **Key insight:** Visual confirmation of what the statistics tell them -- they can literally see where mutations are

---

## Dual Track: Research Skills + Business Pitch

The BTS track's distinctive feature is its integration of scientific research with entrepreneurial thinking. This is not an afterthought -- the business component is woven throughout the week.

### Research Track Timeline
- **Day 1:** Learn the biology (SARS-CoV-2, viruses, biotech)
- **Day 2:** Learn the data format (FASTQ, sequencing)
- **Day 3:** Learn the analysis (Bowtie2, SAMtools)
- **Day 4:** Apply to real data (8 patient samples, pathogen identification)
- **Day 5:** Present findings (research presentation)

### Business Track Timeline
- **Day 1:** Introduction to biotech industry landscape
- **Day 2:** Pitch fundamentals: problem/solution framing, target market, value proposition
- **Day 3 Morning:** First pitch presentations (groups present draft pitches, receive peer feedback)
- **Day 5 Afternoon:** Final refined pitches with investment simulation

### Integration Points
The two tracks reinforce each other. Students' scientific understanding informs their business pitches (they can credibly discuss the technology underlying their startup ideas), while the business framing helps students think about real-world applications of bioinformatics (who needs this technology? why does it matter? how does it reach patients?).

---

## Pathogen Investigation: 8 Patient Samples

The central applied project of BTS is a pathogen identification exercise that mirrors real-world diagnostic bioinformatics.

### Scenario
Students are given FASTQ files from 8 unknown patient samples. Their task is to determine which patients are infected with SARS-CoV-2.

### Method
1. Align each patient's reads to the SARS-CoV-2 reference genome using Bowtie2
2. Examine the alignment rate for each sample
3. Samples with high alignment rates (reads successfully mapping to the viral genome) indicate the presence of SARS-CoV-2
4. Samples with low alignment rates do not contain the virus

### Additional Analysis
- Align the same reads to the human reference genome (HG38) as a control
- Viral reads should NOT align well to the human genome, confirming they are of viral origin
- Use IGV to visualize the aligned reads for infected samples
- Identify mutations: positions where patient reads consistently differ from the reference sequence
- Connect observed mutations to known SARS-CoV-2 variants

### Pedagogical Value
This exercise teaches students that bioinformatics is not abstract -- it solves real diagnostic problems. The scenario of identifying a pathogen from sequencing data is directly analogous to what public health laboratories do during outbreaks. Students experience the thrill of discovery when their alignment results reveal which patients are infected.

---

## Business Component: Startup Pitch Elements

### Problem Statement
- Identify a real, significant health problem
- Quantify the problem: how many people are affected, what is the economic burden
- Make it personal and compelling: why should investors care

### Solution
- Describe the biotech solution clearly
- Explain the underlying technology (students leverage their week's learning)
- Describe how the technology addresses the root cause, not just symptoms
- Examples from student pitches: CRISPR-based gene therapy for sickle cell disease, rapid point-of-care diagnostic testing

### Market Analysis
- Target market identification: who are the customers
- Market size estimation
- Competitive landscape: what existing solutions exist and why yours is better

### Team and Credibility
- Company background (for pitches based on real companies like CRISPR Therapeutics)
- Technical expertise and track record
- Key partnerships and locations

### Investment Simulation
On Day 5, after all pitches are presented, students participate in an investment simulation where they allocate virtual dollars to the startups they find most promising. This gamification element motivates students to create compelling pitches and think critically about what makes a biotech venture viable.

---

## Instructor Notes: Teaching Strategies and Pedagogical Approaches

### The Empowerment Moment
Dr. Lee explicitly aims for a transformative experience where students realize they can use a command line. The goal is for students to feel empowered and capable -- presenting Linux skills as exciting rather than intimidating. This is a deliberate strategy for engaging students who may have no prior programming experience.

### Incremental Pipeline Construction
The technical pipeline is built up incrementally over Days 2--4, with each day adding one or two new steps. Students do not see the full pipeline at once; they build it piece by piece, ensuring they understand each component before moving on. By Day 4, they can execute the complete workflow because they understand what each step does.

### Connecting to CB
For students who completed CB the previous week, BTS instructors frequently reference CB concepts. This reinforcement helps students see the connections between database-driven analysis (CB) and the data production pipeline (BTS).
