# miRcore 2025 Summer Camp -- Full Curriculum Map

> **Note:** The READMEs, study guides, exercises, and cheat sheets in this repository are **my own study materials**, created independently for personal review and practice after attending the miRcore summer camp. They are not official miRcore curriculum, publications, or instructional materials. All original course content remains the intellectual property of miRcore and Dr. Inhan Lee.

## Program Overview

miRcore is a nonprofit organization founded and led by Dr. Inhan Lee (PhD, Biomedical Engineering) that connects the public to real scientific research through hands-on bioinformatics education. Since 2011, miRcore has educated high school students through summer camps, volunteer programs, school outreach clubs (Genes in Disease and Symptoms / GIDAS), and a multi-year research pipeline that culminates in student presentations at professional research conferences alongside world-leading scientists.

The **2025 Summer Camp** runs three sequential one-week tracks over the summer. Each track meets for five days, with morning and afternoon sessions delivered in a hybrid format (in-person at the University of Michigan plus Zoom for remote participants). Students are organized into small groups led by trained group leaders -- typically miRcore alumni now attending top universities -- who facilitate breakout-room discussions, troubleshoot technical issues, and mentor students through their research projects. The camp culture emphasizes community: daily "brain reload" polls on Zulip, cheat-sheet competitions, and a welcoming environment where cameras stay on, questions are encouraged, and lunch menus are shared.

**Students** range from 9th through 12th graders. The CB and BTS tracks accept students across all high school grades, while the SEQ track is restricted to 11th and 12th graders and requires completion of the prior two camps. Many students participate in all three tracks across the summer, building progressively deeper skills.

The overarching mission is **precision medicine**: understanding how individual genomic differences influence health, disease, and treatment, and empowering the next generation to contribute to that future.

---

## Three-Track Comparison Table

| Dimension | CB (Computational Biology) | BTS (Biotechnology Startup) | SEQ (Sequencing Your Genome) |
|---|---|---|---|
| **Target Audience** | 9th--12th graders; no prerequisites | 9th--12th graders; CB recommended | 11th--12th graders only; CB + BTS required |
| **Duration** | 5 days (10 sessions) | 5 days (10 sessions) | 5 days (10 sessions) |
| **Core Question** | How do genes relate to disease, and how can we measure that? | How do we sequence and identify a pathogen from raw data? | What does my own genome tell me about myself? |
| **Learning Goals** | Use public databases to find disease-gene links, analyze gene expression with statistics, interpret biological pathways | Navigate Linux, run a bioinformatics pipeline from raw reads to visualization, develop a biotech business pitch | Process and interpret personal whole-exome sequencing data, understand variant significance, appreciate ethical dimensions |
| **Primary Tools** | OMIM, GEO, GEO2R, KEGG, microarray virtual lab | Great Lakes HPC cluster, Linux CLI, Bowtie2, SAMtools, IGV, NCBI | Great Lakes HPC cluster, Linux CLI, Bowtie2, SAMtools, cutadapt, UCSC Genome Browser, IGV, ClinVar, dbSNP, OMIM |
| **Key Datasets** | Publicly available GEO microarray datasets for chosen diseases | 8 patient FASTQ samples (pathogen identification exercise) | Students' own whole-exome sequencing (WES) FASTQ data |
| **Final Deliverable** | Group research presentation with slides on disease, gene expression analysis, and pathway interpretation | Research presentation on pathogen identification + biotech startup pitch with investment simulation | Individual presentation on personal genome variants and their phenotypic implications |
| **Assessment** | Cheat sheets, daily diagnostics, group presentations, peer feedback | Cheat sheets, business pitch judging, research presentations | Pre/post surveys, self-learning quizzes, individual presentations |

---

## Cross-Track Learning Progression

The three tracks form a deliberate **prerequisite chain** (CB then BTS then SEQ), each building on the knowledge and skills of the previous one.

### CB: Foundation Layer -- "What are genes, and how do we study them?"

CB establishes the conceptual and statistical foundations students need for all subsequent work. Students learn what genes are, how they relate to disease, and how gene expression is measured via microarrays. They are introduced to public databases (OMIM, GEO) and learn fundamental statistics (p-values, t-tests, fold change, volcano plots, FDR correction). They practice reading research papers, interpreting data visualizations, and communicating findings to peers. By the end of CB, students can navigate from a disease question to a dataset, run a GEO2R analysis, identify differentially expressed genes, and contextualize those genes in biological pathways using KEGG.

### BTS: Technical Layer -- "How does the data actually get produced and processed?"

BTS takes students behind the scenes of the sequencing data that CB treated as pre-processed. Students learn the biology of DNA sequencing (Illumina NGS, library prep, reads, quality scores), then gain hands-on Linux command-line skills on the University of Michigan's Great Lakes HPC cluster. They learn the full bioinformatics pipeline: FASTQ files, sequence alignment with Bowtie2, SAM/BAM format processing with SAMtools, and genome visualization with IGV. The applied project -- identifying a pathogen from 8 unknown patient samples -- gives students a realistic research scenario. A parallel entrepreneurship track has students develop biotech startup pitches, connecting science to real-world business applications.

### SEQ: Synthesis Layer -- "What does my own genome say about me?"

SEQ is the culminating experience. Students apply every skill learned in CB and BTS to their **own whole-exome sequencing data**. They process their personal FASTQ files through the same pipeline they learned in BTS (cutadapt for adapter removal, Bowtie2 for alignment, SAMtools for processing), then perform variant calling and interpretation. They look up their own variants in ClinVar, dbSNP, and OMIM, exploring genes like TAS2R38 (bitter taste), ABO (blood type), ABCC11 (earwax type), and MMP1 (sun sensitivity). The week has a strong ethics emphasis -- data security, penetrance, the need for professional confirmation, genetic discrimination (GINA Act), and the philosophical implications of knowing one's genome.

### The Full Arc

```
CB (Week 1)          BTS (Week 2)           SEQ (Week 3)
Question             Measurement            Personal Application
  |                      |                       |
  v                      v                       v
Databases --------> Sequencing Pipeline ---> Own Genome Analysis
Statistics -------> Data Processing -------> Variant Interpretation
Pathways ---------> Visualization ---------> Phenotype Prediction
Communication ----> Business Pitch --------> Ethical Reflection
```

---

## Shared Concepts

### Overlapping Across All Three Tracks

- **Genomics fundamentals**: DNA structure, base pairs, codons, amino acids, gene expression, the central dogma of molecular biology
- **Bioinformatics databases**: NCBI, OMIM (used in CB for disease-gene lookup; used in SEQ for variant annotation)
- **Scientific communication**: All three tracks require group or individual presentations with slides as final deliverables
- **Research methodology**: Asking questions, forming hypotheses, analyzing data, interpreting results
- **Precision medicine vision**: Dr. Lee's opening lecture on connecting genes to individualized treatment appears in all tracks
- **Ethics and data responsibility**: Genetic discrimination, GINA Act, data privacy discussed in CB and deeply in SEQ
- **Cheat sheets**: Students maintain running notes throughout each track, with awards for the most thorough ones
- **Zulip communication**: Daily brain reloads, polls, and group coordination through the Zulip messaging platform
- **Small-group mentorship**: Breakout rooms with group leaders are a consistent pedagogical pattern

### Unique to Each Track

| CB Only | BTS Only | SEQ Only |
|---|---|---|
| Microarray technology | Linux command line (pwd, ls, cd, mkdir, cp, nano, grep, cut, sort) | Personal genome data |
| GEO and GEO2R | Great Lakes HPC cluster and sbatch job submission | Adapter removal with cutadapt |
| Statistics deep dive (p-values, t-tests, Gaussian distribution, FDR) | FASTQ format and quality scores | VCF file format |
| Log2 fold change and volcano plots | Bowtie2 alignment (first encounter) | Variant calling and annotation |
| KEGG pathway analysis | SAM/BAM file processing (first encounter) | ClinVar and dbSNP databases |
| Microarray virtual lab | IGV genome browser (first encounter) | UCSC Genome Browser |
| Guest neurologist speaker | Biotech startup pitch and investment simulation | Personal phenotype prediction |
| -- | Pathogen identification (SARS-CoV-2) | Specific genes: TAS2R38, ABO, ABCC11, MMP1 |
| -- | Protein structure (PDB) | Penetrance and genetic counseling ethics |
| -- | Nucleotide hybridization and melting temperature | 0-based vs 1-based coordinate systems |

---

## Complete Tools and Databases Inventory

### Biological Databases

| Tool/Database | Description | Tracks Used |
|---|---|---|
| **OMIM** (Online Mendelian Inheritance in Man) | Catalog of human genes and genetic disorders; entry point for disease-gene relationships | CB, SEQ |
| **GEO** (Gene Expression Omnibus) | NCBI repository of gene expression datasets from microarray and sequencing experiments | CB |
| **GEO2R** | Interactive web tool for comparing groups of GEO samples; performs statistical analysis and generates volcano plots | CB |
| **KEGG** (Kyoto Encyclopedia of Genes and Genomes) | Pathway database showing how genes interact in biological processes; used for contextualizing differentially expressed genes | CB |
| **NCBI** (National Center for Biotechnology Information) | Central hub for genomic data; used for reference sequences, gene lookup, and SARS-CoV-2 genome exploration | BTS, SEQ |
| **ClinVar** | NCBI database of relationships between genetic variants and clinical conditions | SEQ |
| **dbSNP** | Database of single nucleotide polymorphisms and other short genetic variations | SEQ |
| **PDB** (Protein Data Bank) | 3D structural data for proteins; used to examine S-protein structure | BTS |
| **UCSC Genome Browser** | Web-based tool for visualizing genomic data in the context of a reference genome | SEQ |
| **UniProt** | Protein sequence and functional information database | SEQ |

### Bioinformatics Software

| Tool | Description | Tracks Used |
|---|---|---|
| **Bowtie2** | Fast short-read aligner; maps sequencing reads to reference genomes | BTS, SEQ |
| **SAMtools** | Suite for manipulating SAM/BAM alignment files (view, sort, index) | BTS, SEQ |
| **IGV** (Integrative Genomics Viewer) | Desktop application for interactive visualization of aligned genomic data | BTS, SEQ |
| **cutadapt** | Tool for removing adapter sequences and low-quality bases from sequencing reads | SEQ |
| **nano** | Terminal-based text editor used for editing scripts on the cluster | BTS, SEQ |
| **grep** | Command-line pattern matching tool; used to search within FASTQ and other files | BTS, SEQ |

### Computing Infrastructure

| Resource | Description | Tracks Used |
|---|---|---|
| **University of Michigan Great Lakes HPC Cluster** | High-performance computing cluster where students run bioinformatics pipelines | BTS, SEQ |
| **Linux command line** | Students learn essential Unix commands: `pwd`, `ls`, `cd`, `mkdir`, `cp`, `mv`, `cat`, `less`, `head`, `cut`, `sort`, `wc`, `nano`, `grep` | BTS, SEQ |
| **sbatch / SLURM** | Job scheduler for submitting batch jobs to the cluster | BTS, SEQ |
| **Shell scripting** | Students write `.sh` scripts for automating pipeline steps | BTS, SEQ |

### Communication and Collaboration

| Tool | Description | Tracks Used |
|---|---|---|
| **Zulip** | Team messaging platform used for polls, brain reloads, Q&A, and group coordination | CB, BTS, SEQ |
| **Zoom** | Video conferencing with breakout rooms for small-group discussions | CB, BTS, SEQ |
| **Google Slides** | Collaborative presentation tool for final research presentations | CB, BTS, SEQ |

### Educational Resources

| Resource | Description | Tracks Used |
|---|---|---|
| **Microarray virtual lab** | Interactive simulation of the microarray experiment process | CB |
| **P-value simulator** | Interactive tool for exploring the relationship between sample size, effect size, and statistical significance | CB |
| **Cheat sheets** | Student-maintained reference notes compiled throughout each track | CB, BTS, SEQ |

---

## Skill Progression: From Zero to Genome Analysis

The miRcore 2025 Summer Camp takes students from no bioinformatics background to analyzing their own genome in 15 days. Below is the full progression.

### Week 1 -- CB: Conceptual Foundations (Days 1--5)

**Day 1 -- Genes, Disease, and Precision Medicine**
- What is precision medicine and why does it matter
- How genes relate to diseases (one gene can affect many diseases; one disease involves many genes)
- First database: OMIM -- looking up diseases, finding associated genes
- Ethics: genetic discrimination, GINA Act, data ownership
- Group formation and disease selection for week-long research project

**Day 2 -- Gene Expression and Microarrays**
- DNA/RNA structure, nucleotides, base pairing, central dogma
- Transcription and translation: how genetic information becomes proteins
- Microarray technology: how gene expression is measured at scale
- Virtual microarray lab experience
- GEO database: finding and selecting appropriate datasets

**Day 3 -- Statistics and Volcano Plots**
- Null hypothesis and experimental design
- P-values, t-tests, and statistical significance
- Gaussian (normal) distribution and standard deviation
- Sample size effects on statistical power
- Log2 fold change, log10 p-values, and volcano plots
- Multiple testing problem and False Discovery Rate (FDR / Benjamini-Hochberg)
- GEO2R: running differential expression analysis on chosen datasets

**Day 4 -- GEO2R Analysis and Pathways**
- Guest speaker: Dr. Sami Barmada (neurologist, memory disorders research)
- Completing GEO2R analysis for group research projects
- Introduction to KEGG pathway analysis
- Understanding biological pathways: how genes work in networks, not isolation
- Research project work: interpreting results and building presentations

**Day 5 -- Research Presentations**
- Final diagnostic: synthesizing the week's learning
- Group presentation preparation and practice
- Student presentations to peers and parents
- Social impact discussion of disease research
- Awards for best cheat sheets
- Introduction to miRcore volunteer program (MVP) and future opportunities

### Week 2 -- BTS: Technical Pipeline (Days 6--10)

**Day 6 (BTS Day 1) -- Biotech Introduction and Linux Basics**
- Biotechnology overview: insulin production, CRISPR, PCR, Illumina sequencing
- Wet lab vs. dry lab distinction
- COVID-19 biology: virus structure, spike protein, SARS-CoV-2 genome on NCBI
- S-protein location, amino acid sequence, protein structure (PDB)
- First login to Great Lakes HPC cluster with UM credentials
- Linux fundamentals: `pwd`, `ls`, `mkdir`, `cd`, `cp`
- FASTA file format

**Day 7 (BTS Day 2) -- NGS and FASTQ Analysis**
- DNA hybridization, melting temperature, GC content
- Guest speaker: Dr. Sally Choi (biotech career path)
- Next-generation sequencing (NGS) technology: how reads are produced
- FASTQ file format: sequences + quality scores (Phred scores)
- File transfer on the cluster
- `grep` for pattern matching in sequence files
- GC content analysis of sequencing data
- Business pitch fundamentals: problem/solution framing

**Day 8 (BTS Day 3) -- SAM/BAM Files and Alignment**
- Student business pitch presentations (CRISPR Therapeutics, rapid testing startups)
- Quiz review on bioinformatics concepts
- Bowtie2: aligning reads to a reference genome (SARS-CoV-2)
- SAM file format: understanding alignment fields (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR)
- SAMtools: `view` (SAM to BAM conversion), `sort`, `index`
- `nano` text editor, `cut`, `sort` commands
- Processing aligned reads into analysis-ready BAM files

**Day 9 (BTS Day 4) -- Genome Mapping and IGV**
- Review of Bowtie2 and SAMtools workflow
- IGV (Integrative Genomics Viewer): loading and visualizing BAM files
- Shell scripting for batch job submission (`sbatch`)
- Mapping reads to human genome (HG38 reference)
- Project data analysis: processing 8 patient samples
- Identifying SARS-CoV-2 in patient samples through alignment rates
- Mutation identification in viral sequences

**Day 10 (BTS Day 5) -- Research Projects and Business Pitches**
- Camp recap and review of full pipeline
- Research project completion: finalizing pathogen identification analysis
- Target market and audience identification for business pitches
- Final business pitch presentations with investment simulation (students allocate virtual dollars)
- Research presentations to parents and peers
- Introduction to volunteer program and continuing research opportunities

### Week 3 -- SEQ: Personal Genome (Days 11--15)

**Day 11 (SEQ Day 1) -- Sequencing Fundamentals**
- Self-learning format introduction (more independent than CB/BTS)
- Pre-assessment quiz for baseline measurement
- Whole-exome sequencing (WES) overview: what is sequenced and why
- Receiving personal exome sequencing data files
- Data quality assessment
- Adapter removal with cutadapt
- Submitting first processing jobs on the cluster

**Day 12 (SEQ Day 2) -- Exome Data Processing**
- Processing personal WES data: mapping to human reference genome (HG38) with Bowtie2
- Understanding alignment results and mapping statistics
- SAMtools processing: converting, sorting, and indexing BAM files
- Introduction to UCSC Genome Browser
- BAM/BAI file formats and their relationship
- Exploring genomic regions in the browser context

**Day 13 (SEQ Day 3) -- Variant Identification**
- TAS2R38 gene: bitter taste receptor, haplotypes, and phenotype prediction
- Introduction to variant identification from aligned data
- OMIM and UCSC Genome Browser integration for variant lookup
- Understanding variant coordinates: 0-based vs. 1-based coordinate systems
- Phenotype-genotype connections: how a single nucleotide change can alter taste perception
- Variant ranking by phenotype impact
- Breakout room discussions on interpreting variant significance

**Day 14 (SEQ Day 4) -- Personal Genome Exploration**
- ABO blood type gene: determining personal blood type from variant data
- ABCC11 gene: earwax type (wet vs. dry) and its population genetics
- MMP1 gene: sun sensitivity and skin response
- IGV visualization of personal genome data
- Comparing personal variants to known database entries in ClinVar and dbSNP
- Navigating variant databases: understanding pathogenicity classifications
- Ethics deep dive: penetrance (having a variant does not guarantee the phenotype), the importance of professional genetic counseling, data security

**Day 15 (SEQ Day 5) -- Final Presentations**
- Post-camp survey (compared to pre-camp baseline)
- Final variant hunting and analysis
- Slide creation for individual presentations
- Practice presentations in breakout rooms
- Final presentations to parents and peers: each student presents findings from their own genome
- Camp conclusion: reflection on the journey from zero to genome analyst
- Introduction to Citizen Scientist Sequencing Initiative (CSSI)
- Pathways forward: volunteer program, GIDAS clubs, multi-year research, conference presentations

---

## Summary

The miRcore 2025 Summer Camp represents a carefully designed 15-day curriculum that transforms high school students with no bioinformatics background into confident genomic data analysts. The progression from database exploration (CB) through pipeline construction (BTS) to personal genome interpretation (SEQ) mirrors the actual workflow of modern genomics research. By the final day, students have not only mastered technical skills but have also grappled with the ethical, social, and personal dimensions of genomic data -- an experience that positions them to be thoughtful contributors to the precision medicine future that miRcore envisions.
