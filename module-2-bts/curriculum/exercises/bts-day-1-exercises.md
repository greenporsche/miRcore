# BTS Day 1 Exercises: Biotech Intro and Linux Basics

## Warm-Up Questions
1. What is the difference between a wet lab and a dry lab, and which type does bioinformatics fall under?
2. Name three real-world applications of biotechnology that were discussed in the lecture.
3. What role does the spike protein play in how SARS-CoV-2 infects human cells?

## Hands-On Exercises

### Exercise 1: Navigating NCBI
**Objective**: Use the NCBI database to locate the SARS-CoV-2 reference genome and identify key features.
**Instructions**:
1. Go to https://www.ncbi.nlm.nih.gov/
2. Search for the SARS-CoV-2 reference genome (accession: NC_045512.2)
3. Switch from FASTA view to Graphics view
4. Locate the S (spike) protein and record its genomic position (start and end base pairs)
5. Click on the S protein to find the first amino acid of the spike protein
6. Record the total length of the SARS-CoV-2 genome in base pairs

**Expected Output**: S protein location: 21,563-25,384 bp; first amino acid: M (methionine); genome length: ~29,903 bp.

### Exercise 2: Logging into Great Lakes
**Objective**: Successfully connect to the University of Michigan Great Lakes HPC cluster via SSH.
**Instructions**:
1. Open your terminal application
2. Type `ssh uniqname@greatlakes.arc-ts.umich.edu` (replace `uniqname` with your UM username)
3. Enter your password when prompted
4. Approve the Duo Mobile two-factor authentication push
5. Once logged in, type `pwd` and record the output
6. Type `ls` to see what files and directories exist in your home directory

**Expected Output**: After logging in, `pwd` should show something like `/home/uniqname`.

### Exercise 3: Basic Linux Navigation
**Objective**: Practice creating directories, navigating, and copying files using Linux commands.
**Instructions**:
1. Create a new directory called `bts_project` using `mkdir bts_project`
2. Navigate into the new directory: `cd bts_project`
3. Confirm your location with `pwd`
4. Create two subdirectories inside `bts_project`: `data` and `results`
5. Navigate into `data`: `cd data`
6. Go back up one directory: `cd ..`
7. List all contents to confirm both subdirectories exist: `ls`

**Expected Output**:
```
$ ls
data  results
```

### Exercise 4: Exploring FASTA Files
**Objective**: Understand the FASTA file format by examining a real sequence file.
**Instructions**:
1. On Great Lakes, navigate to the shared course data directory (as provided by instructors)
2. Use `cat` to display the contents of a FASTA file containing the SARS-CoV-2 genome
3. Identify the header line (starts with `>`) and note the accession number
4. Count how many lines the file has using `wc -l filename.fasta`
5. Use `head -n 5` to view just the first 5 lines

**Expected Output**: The file should start with a `>` header line containing the accession number, followed by lines of nucleotide sequence (A, T, G, C characters).

### Exercise 5: Viewing a 3D Protein Structure
**Objective**: Use the Protein Data Bank to visualize the SARS-CoV-2 spike protein structure.
**Instructions**:
1. Go to https://www.rcsb.org/
2. Search for SARS-CoV-2 spike protein (try PDB ID: 6VXX or 6VYB)
3. Use the 3D viewer to rotate and examine the protein structure
4. Identify where the receptor binding domain (RBD) is located on the spike
5. Write 2-3 sentences describing what you observe about the shape and structure

**Expected Output**: A written description noting the trimeric structure of the spike protein and the location of the RBD at the top of the protein.

## Challenge Problems

### Challenge 1: Linux File Organization
**Objective**: Create a complete project directory structure for the BTS camp using only the command line.
**Instructions**:
1. Starting from your home directory, create the following structure:
   ```
   bts_camp/
   ├── day1/
   │   ├── notes/
   │   └── data/
   ├── day2/
   │   ├── notes/
   │   └── data/
   ├── day3/
   ├── day4/
   └── day5/
   ```
2. Use `ls` at each level to verify your structure
3. Navigate from `day1/notes/` to `day5/` using a single `cd` command with `..` (relative path)

### Challenge 2: NCBI Deep Dive
**Objective**: Use NCBI to compare different coronaviruses.
**Instructions**:
1. On NCBI, look up the genome lengths of SARS-CoV-2, SARS-CoV (2003), and MERS-CoV
2. Record the genome length for each
3. Compare the spike protein lengths across all three viruses
4. Write a short paragraph about the similarities and differences you notice

---

## Answer Key

### Warm-Up Answers
1. A wet lab involves physical experiments with chemicals and biological materials (e.g., PCR, gel electrophoresis). A dry lab uses computers for data analysis. Bioinformatics is a dry lab discipline.
2. Examples include: insulin production using recombinant DNA, CRISPR gene editing for treating sickle cell disease, PCR (polymerase chain reaction) for COVID testing, and Illumina next-generation sequencing.
3. The spike protein (S-protein) on the surface of SARS-CoV-2 binds to the ACE2 receptor on human cells, allowing the virus to enter and infect the cell. This is why it is a key target for vaccines.

### Exercise Answers
1. **NCBI Exercise**: S protein is at positions 21,563-25,384 on the SARS-CoV-2 genome. The first amino acid is methionine (M), as is typical for most proteins. The total genome length is approximately 29,903 base pairs.
2. **Great Lakes Login**: `pwd` should output `/home/your_uniqname`. The `ls` command will show any existing files/directories in your home folder.
3. **Linux Navigation**: After running the commands, `ls` inside `bts_project` should show `data` and `results`. The `pwd` command from inside `bts_project` should show `/home/uniqname/bts_project`.
4. **FASTA Files**: The first line begins with `>` followed by an accession number and description. Subsequent lines contain nucleotide sequences using A, T, G, C characters.
5. **PDB Viewer**: The spike protein has a trimeric structure (three identical subunits). The receptor binding domain (RBD) is located at the top of each subunit and is the region that directly contacts the human ACE2 receptor.

### Challenge Answers
1. **Linux File Organization**: Commands needed:
   ```bash
   mkdir -p bts_camp/day1/notes bts_camp/day1/data bts_camp/day2/notes bts_camp/day2/data bts_camp/day3 bts_camp/day4 bts_camp/day5
   ```
   To navigate from `day1/notes/` to `day5/`: `cd ../../day5`

2. **NCBI Deep Dive**: SARS-CoV-2 genome is ~29,903 bp, SARS-CoV (2003) is ~29,751 bp, and MERS-CoV is ~30,119 bp. All three are positive-sense single-stranded RNA viruses with similar genome sizes. The spike proteins vary in length but share structural similarities as they all mediate host cell entry through different receptors.
