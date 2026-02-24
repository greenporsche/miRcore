# SEQ Track (Sequencing Your Genome) -- Detailed Analysis

## Track Overview

The Sequencing Your Genome (SEQ) track is the culminating experience of the miRcore summer camp curriculum. Restricted to 11th and 12th graders who have completed both the CB and BTS tracks, SEQ is where students apply everything they have learned to the most personal dataset imaginable: their own genome. Students process their own whole-exome sequencing (WES) data through a bioinformatics pipeline, identify variants in genes with known phenotypic effects, and grapple with the ethical implications of personal genomic knowledge.

**Audience:** 11th--12th graders only. CB and BTS completion required. This is the most advanced and selective of the three tracks.

**Goals:**
- Process personal whole-exome sequencing data through a complete bioinformatics pipeline
- Identify and interpret genetic variants in their own genome
- Understand the relationship between genotype and phenotype for specific well-characterized genes
- Develop a deep appreciation for the ethical dimensions of personal genomic data
- Present findings about their own genome to peers and parents

**Structure:** SEQ differs from CB and BTS in its emphasis on self-directed learning. As Dr. Lee explains on Day 1: "This genome camp is slightly different... It will be a little bit more your part." The first 30 minutes of each day feature a self-learning quiz where students work independently, without group leader guidance. The main camp sessions are more collaborative but still expect greater independence. The overall tone is more mature, reflecting the older student population and their accumulated experience.

---

## Day-by-Day Breakdown

### Day 1: Sequencing Fundamentals

**Morning Session (Part 1)**

*Main Topics:*
- Camp context: "You guys are doing our miRcore summer camp for the past two weeks... this is August, so it seems like almost like you are being on entire summer with miRcore."
- How SEQ differs from CB and BTS: older students only, more self-directed, more independent work
- Self-learning format introduction: first 30 minutes of each day are independent self-assessment
- Pre-assessment quiz: baseline measurement of student knowledge before the week begins
  - Covers: sequencing concepts, genomics vocabulary, bioinformatics tools
  - Purpose: not graded, but used to identify areas needing reinforcement
- Group formation: smaller groups given the more advanced level
- Whole-exome sequencing (WES) overview:
  - What is the exome? The ~1.5% of the genome that codes for proteins (~20,000 genes)
  - Why sequence only the exome? Most known disease-causing variants are in coding regions; more cost-effective than whole-genome sequencing
  - How WES works: capture probes selectively bind to exonic DNA, which is then sequenced
- Receiving personal data: students receive their own FASTQ files
  - Data was generated from actual saliva or blood samples provided by the students prior to camp
  - Each student has unique data representing their individual genome
  - Emphasis on data privacy: "sequencing your genome without sharing with anybody else"

*Key Concepts Introduced:*
- Whole-exome sequencing as a targeted approach
- The exome as the protein-coding portion of the genome
- Personal genomic data ownership and privacy

**Afternoon Session (Part 2)**

*Main Topics:*
- Review of sequencing fundamentals from self-learning quiz
- Data quality assessment: examining raw FASTQ files
  - Checking read counts, read lengths, quality score distributions
  - Identifying potential quality issues before processing
- Adapter removal with cutadapt:
  - What are adapters? Short synthetic sequences added during library prep that are not part of the genome
  - Why remove them? Adapters interfere with alignment if left in the reads
  - Running cutadapt: specifying adapter sequences, quality trimming thresholds
  - Examining cutadapt output: how many reads were trimmed, how many were discarded
- Submitting first processing jobs on the cluster:
  - Review of SLURM job submission from BTS
  - Writing shell scripts for cutadapt jobs
  - Email notifications for job start/completion
  - Examining log files: job resource usage (time, memory, CPU)

*Key Concepts Introduced:*
- Data quality control as the essential first step (reinforcing BTS: "cleaning data is probably 90% of your job")
- Adapter contamination and its effects on downstream analysis
- cutadapt as a preprocessing tool

*Tools/Databases Used:*
- **Great Lakes HPC cluster**: Job submission, file management
- **cutadapt**: Adapter removal and quality trimming
- **Linux CLI**: File examination, script writing

*Group Work:*
- Running cutadapt on personal data
- Comparing adapter removal statistics across group members
- Troubleshooting job submission issues

---

### Day 2: Exome Data Processing

**Morning Session (Part 1)**

*Main Topics:*
- Review of Day 1: adapter removal results, job submission mechanics
- Understanding job completion emails: begin notification, end notification, resource usage stats
- Log files: where output goes, how to interpret them
- Mapping personal exome data to the human reference genome (HG38) with Bowtie2:
  - Review of Bowtie2 from BTS (students already know this tool from the pathogen identification exercise)
  - Scale difference: human exome data is much larger than the SARS-CoV-2 data from BTS
  - Running Bowtie2 as a batch job: writing SLURM scripts with appropriate resource requests
  - Expected alignment rates for exome data (should be high since reads come from human DNA)
- Understanding alignment results and mapping statistics:
  - Overall alignment rate
  - Concordantly mapped reads vs. discordantly mapped
  - Unmapped reads and what they might represent
- File management: organizing output files, moving between scratch and home directories

*Key Concepts Introduced:*
- Scaling up: from viral genome to human genome alignment
- Alignment statistics for quality assessment
- File organization on an HPC system

**Afternoon Session (Part 2)**

*Main Topics:*
- SAMtools processing of personal data:
  - `samtools view`: SAM to BAM conversion
  - `samtools sort`: Position-based sorting
  - `samtools index`: Creating BAI index files
  - Review from BTS: students are familiar with these steps but now apply them to their own data
- BAM and BAI file formats:
  - BAM: compressed binary alignment file
  - BAI: index file that enables fast lookups of specific genomic regions
  - Why both files are needed for visualization and downstream analysis
- Introduction to UCSC Genome Browser:
  - Web-based genome visualization tool (different from IGV used in BTS)
  - Navigating the browser: searching by gene name or coordinates
  - Understanding genome tracks: reference sequence, known genes, conservation, variant databases
  - Zooming in and out of genomic regions
  - Customizing track visibility
- Exploring genomic regions in the browser context:
  - How to find a gene of interest
  - Understanding the visual representation of exons, introns, UTRs
  - Connecting what they see in the browser to the data in their BAM files

*Key Concepts Introduced:*
- UCSC Genome Browser as a rich genomic context tool
- The visual language of genome browsers: tracks, coordinates, annotations
- Relationship between BAM files and genome browser visualization

*Tools/Databases Used:*
- **SAMtools**: Full processing pipeline on personal data
- **UCSC Genome Browser**: First exploration
- **Great Lakes HPC**: Continued batch job management

---

### Day 3: Variant Identification

**Morning Session (Part 1)**

*Main Topics:*
- Recap of processing pipeline completed so far: FASTQ -> cutadapt -> Bowtie2 -> SAMtools -> BAM/BAI
- The importance of clean data: "If you have a clean nice data, the rest is really good"
- Discussion of data bias: lack of healthy/normal genome data vs. abundant patient data
- **TAS2R38 gene deep dive: the bitter taste receptor**
  - TAS2R38: a gene encoding a bitter taste receptor on chromosome 7
  - The PTC (phenylthiocarbamide) taste test: a classic genetics demonstration
  - Three key SNPs in TAS2R38 that determine bitter taste sensitivity:
    - Position 145: Proline (P) vs. Alanine (A) -- amino acid 49
    - Position 785: Alanine (A) vs. Valine (V) -- amino acid 262
    - Position 886: Isoleucine (I) vs. Valine (V) -- amino acid 296
  - Haplotypes: PAV (taster) vs. AVI (non-taster)
    - PAV/PAV: strong taster (homozygous taster)
    - PAV/AVI: taster (heterozygous, dominant)
    - AVI/AVI: non-taster (homozygous non-taster)
  - Students predict their own taster status from their variant data
- Introduction to variant identification from aligned data:
  - What is a variant? A position where the student's sequence differs from the reference genome
  - Types of variants: SNP (single nucleotide polymorphism), insertion, deletion
  - How to find variants: examining BAM files for positions with non-reference bases
- OMIM and UCSC Genome Browser integration:
  - Looking up TAS2R38 in OMIM to understand the gene's phenotypic associations
  - Navigating to TAS2R38 in the UCSC Genome Browser to see the gene structure
  - Finding the specific variant positions in the browser

*Key Concepts Introduced:*
- Haplotypes: combinations of variants that travel together
- Dominant vs. recessive inheritance in a concrete example
- Variant identification as the bridge from raw data to biological meaning
- Phenotype prediction from genotype

*Tools/Databases Used:*
- **OMIM**: Gene-phenotype lookup for TAS2R38
- **UCSC Genome Browser**: Navigating to gene regions, viewing variant positions
- **IGV or BAM file examination**: Finding variants in personal data

**Afternoon Session (Part 2)**

*Main Topics:*
- Understanding variant coordinates: **0-based vs. 1-based coordinate systems**
  - Different databases and tools use different conventions
  - BED format uses 0-based coordinates (start inclusive, end exclusive)
  - VCF format and SAM format use 1-based coordinates
  - This off-by-one difference is a common source of errors in bioinformatics
  - Students practice converting between coordinate systems
- Variant ranking by phenotype impact:
  - Not all variants are equal: some are benign, some are pathogenic
  - How to assess variant significance: population frequency, known associations, predicted protein impact
- Phenotype-genotype connections:
  - Moving from observing a variant to understanding what it means for the individual
  - The concept of effect size: how much does this variant change the phenotype?
- Breakout room discussions on variant interpretation:
  - Students share their TAS2R38 results with group members
  - Discussion: Did your genotype match your actual taste experience?
  - Understanding that genetics is probabilistic, not deterministic

*Key Concepts Introduced:*
- Coordinate system differences across bioinformatics tools
- Variant pathogenicity classification
- The probabilistic nature of genotype-phenotype relationships

*Group Work:*
- Identifying personal TAS2R38 variants
- Comparing results across the group
- Discussing the experience of discovering something genetic about themselves

---

### Day 4: Personal Genome Exploration

**Morning Session (Part 1)**

*Main Topics:*
- Review of TAS2R38 variant identification process
- Students share results from their Day 3 practice: "How was looking at your own data set? Isn't it exciting?"
- **ABO blood type gene:**
  - ABO gene on chromosome 9: encodes a glycosyltransferase enzyme
  - Three main alleles: A, B, O
  - Blood type determination: A/A or A/O = Type A; B/B or B/O = Type B; A/B = Type AB; O/O = Type O
  - Key variant positions that distinguish A, B, and O alleles
  - Students determine their own blood type from their variant data
  - Cross-referencing with known blood type (if students know theirs)
- **ABCC11 gene: earwax type**
  - ABCC11 on chromosome 16: encodes an ABC transporter
  - A single SNP (rs17822931, G>A) determines earwax type
  - G/G: wet earwax (dominant worldwide except East Asia)
  - A/A: dry earwax (common in East Asian populations)
  - G/A: wet earwax (wet is dominant)
  - Population genetics connection: geographic distribution of the alleles
  - Students identify their own genotype and predict earwax type
- **MMP1 gene: sun sensitivity**
  - MMP1 (matrix metalloproteinase 1) on chromosome 11
  - A promoter variant affects MMP1 expression levels
  - Higher expression associated with increased skin sensitivity to UV radiation
  - Connection to melanoma risk and sun protection behavior
  - Students identify their variant and discuss implications
- IGV visualization of personal genome data:
  - Loading personal BAM and BAI files into IGV
  - Navigating to each gene of interest (TAS2R38, ABO, ABCC11, MMP1)
  - Visually confirming variants: seeing the colored bases that indicate non-reference positions
  - Comparing heterozygous vs. homozygous variants in the read pileup

*Key Concepts Introduced:*
- Multiple genes with clear, well-understood phenotypic effects
- Dominant vs. recessive patterns across different genes
- Population genetics and geographic variation in allele frequencies
- Visual variant confirmation in IGV

*Tools/Databases Used:*
- **IGV**: Loading and visualizing personal BAM files at gene-specific locations
- **OMIM**: Looking up gene-phenotype relationships
- **ClinVar**: Checking variant clinical significance
- **dbSNP**: Looking up rs numbers and population frequencies

**Afternoon Session (Part 2)**

*Main Topics:*
- Comparing personal variants to known database entries:
  - **ClinVar**: NCBI database of variant-disease relationships
    - Pathogenicity classifications: benign, likely benign, uncertain significance (VUS), likely pathogenic, pathogenic
    - Understanding what each classification means
    - The challenge of Variants of Uncertain Significance (VUS)
  - **dbSNP**: database of known short genetic variations
    - rs numbers as unique variant identifiers
    - Population allele frequencies: how common is your variant in different populations?
    - Minor Allele Frequency (MAF): the frequency of the less common allele
- Navigating variant databases:
  - How to search by rs number, gene name, or genomic coordinate
  - Interpreting the evidence behind pathogenicity classifications
  - Understanding that classifications can change as more evidence accumulates
- **Ethics deep dive:**
  - **Penetrance:** Having a variant does not guarantee the phenotype will manifest
    - Incomplete penetrance: some carriers of pathogenic variants never develop the condition
    - Variable expressivity: same variant can cause different severity in different individuals
    - Why penetrance matters: a "pathogenic" variant is a risk factor, not a diagnosis
  - **Professional confirmation:** Personal genome analysis is educational, not clinical
    - Students' results are NOT medical diagnoses
    - Clinical genetic testing requires CLIA-certified labs and genetic counselor interpretation
    - Emphasis: "If you find something concerning, talk to a healthcare provider, do not self-diagnose"
  - **Data security:**
    - Who has access to your genomic data?
    - miRcore's approach: students own their data, miRcore does not retain copies
    - The Citizen Scientist Sequencing Initiative (CSSI): a vision for individual genome ownership
  - **Genetic discrimination:**
    - Review of GINA Act protections (introduced in CB, deepened here)
    - Limitations: GINA does not cover life insurance, disability insurance, or long-term care insurance
    - Military and government employment considerations
    - The tension between genetic knowledge and potential discrimination

*Key Concepts Introduced:*
- Variant pathogenicity classification systems
- The critical concept of penetrance
- The distinction between educational and clinical genetic analysis
- Data sovereignty and genomic privacy

*Tools/Databases Used:*
- **ClinVar**: Variant clinical significance lookup
- **dbSNP**: Population frequency data
- **OMIM**: Gene-disease relationship confirmation

*Group Work:*
- Looking up personal variants in ClinVar and dbSNP
- Discussing findings with group members (at each student's comfort level -- sharing is voluntary)
- Beginning to prepare presentation slides
- Ethical discussions in breakout rooms: what would you want to know? what would you not want to know?

---

### Day 5: Final Presentations

**Morning Session (Part 1)**

*Main Topics:*
- Post-camp survey: mirroring the pre-assessment from Day 1
  - Same or similar questions to measure growth over the week
  - Students can see how much their understanding has deepened
- Framing: "You are actually a pioneer. This is really, really just new things we are trying to build together."
- Final variant hunting and analysis:
  - Students complete their exploration of personal variants
  - Looking for additional interesting variants beyond the four core genes
  - Searching ClinVar for any variants of note in their exome data
- Slide creation for individual presentations:
  - Unlike CB (group presentations) and BTS (group research + group pitches), SEQ presentations have an individual component -- each student presents findings from their own genome
  - Presentation structure: introduce yourself, describe the pipeline, show your variants for each gene, explain the phenotype predictions, discuss what you learned about your own biology
  - Guidance on what to share and what to keep private: "You choose what to present; you do not have to share anything you are not comfortable with"

*Hands-on Activities:*
- Final data analysis
- Slide creation in Google Slides
- Writing presentation scripts or notes
- Practice presentations in breakout rooms

**Afternoon Session (Part 2)**

*Main Topics:*
- Final presentations to parents and peers:
  - Each student or group presents findings from their genome analysis
  - Presentations cover the full pipeline journey: raw data -> processing -> variant identification -> phenotype interpretation
  - Students show IGV screenshots of their own variants
  - Students explain what they learned about TAS2R38 (can they taste bitter compounds?), ABO (blood type), ABCC11 (earwax type), and MMP1 (sun sensitivity)
  - Q&A from peers and parents after each presentation
- Camp conclusion:
  - Reflection on the journey: from the first CB day knowing nothing about bioinformatics to analyzing their own genome
  - Discussion of what this experience means for their future
  - The significance of personal genome ownership in the precision medicine era
- Introduction to continuing opportunities:
  - Citizen Scientist Sequencing Initiative (CSSI): Dr. Lee's vision for a future where individuals own and understand their own genomic data
  - miRcore Volunteer Program (MVP): continuing research and leadership development
  - GIDAS clubs: bringing bioinformatics to their own schools
  - Research conference: presenting more developed research at the annual miRcore conference
  - Multi-year research program: some students continue for years, developing deep expertise

---

## Learning Arc

```
Day 1: ETHICS + PREPROCESSING -> Receive data, understand privacy, clean raw reads
Day 2: PROCESSING             -> Map to human genome, process alignments, explore browser
Day 3: VARIANT ANALYSIS       -> Identify first variants (TAS2R38), learn coordinates
Day 4: PERSONAL GENOME        -> Explore multiple genes (ABO, ABCC11, MMP1), ethics deep dive
Day 5: COMMUNICATION          -> Present personal findings, reflect on the journey
```

The SEQ arc is unique because the **personal stakes** increase throughout the week. On Day 1, students are processing data that happens to be theirs, but it feels abstract (FASTQ files look the same regardless of source). By Day 3, they are discovering whether they are bitter tasters. By Day 4, they know their blood type, earwax type, and sun sensitivity from their DNA. The week builds toward an increasingly personal and emotionally resonant experience, which is why the ethics discussions are strategically placed on Days 1 and 4 -- bookending the most personal discoveries.

---

## Personal Genome Analysis: Students Analyze Their Own Exome Sequencing Data

This is the defining feature of the SEQ track. Prior to camp, students provide biological samples (saliva or blood) that are sent for whole-exome sequencing. The sequencing is performed by a commercial service, and the resulting FASTQ files are made available to students on the Great Lakes cluster on Day 1.

### What Makes This Special
- **It is real.** These are not simulated datasets or classroom exercises. Students are working with their actual DNA sequence.
- **It is personal.** Every variant they find tells them something about their own biology.
- **It is private.** miRcore does not retain copies of student data. Students own their sequence. The camp philosophy aligns with the Citizen Scientist Sequencing Initiative: individuals should control their own genomic information.
- **It is empowering.** By the end of the week, students have done something that most adults have never done -- analyzed their own genome at the variant level.

### The Processing Pipeline (Applied to Personal Data)
1. **Receive FASTQ files** -- personal whole-exome sequencing data in their home directory
2. **Quality control** -- examine read counts and quality scores
3. **Adapter removal** -- cutadapt to remove library preparation artifacts
4. **Alignment** -- Bowtie2 to map reads to the human reference genome (HG38)
5. **Post-processing** -- SAMtools to convert, sort, and index BAM files
6. **Visualization** -- IGV and UCSC Genome Browser to view aligned data
7. **Variant identification** -- finding positions where personal reads differ from the reference
8. **Variant interpretation** -- looking up variants in ClinVar, dbSNP, and OMIM

---

## Key Genes Analyzed

### TAS2R38 -- Bitter Taste Receptor (Day 3)

**Gene:** TAS2R38, located on chromosome 7p34

**What it does:** Encodes a G-protein coupled receptor in taste receptor cells on the tongue. This receptor detects bitter compounds, specifically PTC (phenylthiocarbamide) and PROP (6-n-propylthiouracil), as well as naturally occurring bitter compounds in vegetables like broccoli and Brussels sprouts.

**The variants:** Three SNPs define the two main haplotypes:
- **PAV haplotype** (Proline-Alanine-Valine at positions 49, 262, 296): TASTER -- can strongly taste bitter compounds
- **AVI haplotype** (Alanine-Valine-Isoleucine at the same positions): NON-TASTER -- cannot taste PTC/PROP

**Inheritance pattern:** Tasting is dominant. PAV/PAV and PAV/AVI individuals are tasters; only AVI/AVI individuals are non-tasters.

**Why this gene is taught first:** TAS2R38 is the ideal introductory variant because (1) the phenotype is immediately verifiable (students can try PTC test strips), (2) it has no medical consequences, so there is no anxiety about the result, (3) it demonstrates the concept of haplotypes clearly, and (4) it is engaging -- students are excited to learn why they do or do not like broccoli.

**Population genetics:** Approximately 70% of people worldwide are tasters (PAV carriers). The frequency varies by population, providing an entry point for discussing human genetic diversity.

---

### ABO -- Blood Type (Day 4)

**Gene:** ABO, located on chromosome 9q34.2

**What it does:** Encodes a glycosyltransferase enzyme that adds sugar molecules to the H antigen on the surface of red blood cells. The specific sugars added (or not added) determine blood type.

**The variants:**
- **A allele:** Enzyme adds N-acetylgalactosamine -> Type A
- **B allele:** Enzyme adds galactose -> Type B (differs from A by key nucleotide changes)
- **O allele:** A frameshift deletion inactivates the enzyme -> no additional sugars added -> Type O

**Inheritance:** Codominant (A and B are codominant with each other; both are dominant over O)
- A/A or A/O = Type A
- B/B or B/O = Type B
- A/B = Type AB
- O/O = Type O

**Why this gene matters:** Blood type is one of the most familiar genetic concepts, and most students already know their blood type (or can find out). Being able to determine their blood type from their DNA data provides a powerful validation experience: "My sequencing data correctly predicted something I already know about myself."

---

### ABCC11 -- Earwax Type (Day 4)

**Gene:** ABCC11, located on chromosome 16q12.1

**What it does:** Encodes an ATP-binding cassette (ABC) transporter protein. The protein is expressed in ceruminous glands (earwax-producing glands) and apocrine sweat glands.

**The variant:** A single SNP, rs17822931 (c.538G>A), determines earwax type:
- **G/G (Gly180Gly):** Wet, sticky, honey-colored earwax
- **G/A:** Wet earwax (wet is dominant)
- **A/A (Gly180Arg):** Dry, flaky, gray earwax

**Population genetics:** This variant has one of the starkest geographic distributions of any human polymorphism. The dry earwax allele (A) reaches frequencies above 80--95% in East Asian populations but is rare (<5%) in European and African populations. This provides an excellent case study in population genetics, natural selection, and human migration.

**Additional association:** The same variant affects body odor and sweat composition, and has been associated with breast cancer risk through its effect on apocrine gland function.

**Why this gene is taught:** It is a fun, non-threatening phenotype that (1) most students have never thought about, (2) demonstrates pleiotropy (one gene affecting multiple traits), and (3) vividly illustrates population-level genetic variation.

---

### MMP1 -- Sun Sensitivity (Day 4)

**Gene:** MMP1 (Matrix Metalloproteinase 1), located on chromosome 11q22.2

**What it does:** Encodes an enzyme (collagenase-1) that breaks down collagen in the extracellular matrix. MMP1 is involved in tissue remodeling, wound healing, and the skin's response to UV radiation.

**The variant:** A promoter polymorphism affects the level of MMP1 gene expression. Higher expression leads to more collagen degradation in response to UV exposure, which is associated with:
- Greater sun sensitivity
- Faster photoaging (wrinkling)
- Potentially increased melanoma risk

**Why this gene is taught:** MMP1 introduces a different type of genetic variant -- a regulatory variant that affects gene expression levels rather than protein structure. This contrasts with the coding variants in TAS2R38 and ABO. It also has practical, relatable implications: students can use their result to inform their sunscreen habits. The connection to skin cancer risk also introduces the concept of risk factors vs. deterministic outcomes, reinforcing the penetrance discussion.

---

## Ethics Emphasis

Ethics is not a one-time lecture in the SEQ track -- it is woven throughout the entire week. The ethical framework has several key components:

### Data Security and Ownership
- Students process their data on the University of Michigan's secure HPC cluster
- Data is stored in the student's personal directory with standard HPC access controls
- miRcore does not retain copies of student genomic data after the camp
- The Citizen Scientist Sequencing Initiative (CSSI) represents miRcore's vision for a future where individuals own their genome data without it passing through commercial databases
- Dr. Lee: "We hope that there is some possible way individual can own their own genome without sharing with anybody else including us"

### Penetrance and Probabilistic Genetics
- A recurring theme: having a variant does not mean you will develop the condition
- Incomplete penetrance is discussed explicitly: some people carry "pathogenic" variants and never become sick
- Variable expressivity: the same variant can manifest differently in different people
- The importance of not jumping to conclusions from a single variant

### Professional Confirmation
- Repeatedly emphasized: camp results are educational, NOT clinical diagnoses
- Clinical genetic testing requires CLIA-certified laboratories and interpretation by genetic counselors
- Students are told: if you find something concerning, discuss it with a healthcare provider
- The camp intentionally focuses on non-medical, fun phenotypes (taste, earwax, blood type) to minimize anxiety while still teaching real variant interpretation skills

### Genetic Discrimination
- The GINA Act (introduced in CB, deepened in SEQ): protects against discrimination in health insurance and employment
- Limitations of GINA: does not cover life insurance, disability insurance, long-term care, or military service
- Discussion of the tension between genetic knowledge and its potential misuse
- Connection to the data ownership philosophy: if you control your data, you control who sees it

### The Right Not to Know
- Implicit in the presentation guidelines: "You choose what to share; you do not have to share anything you are not comfortable with"
- Recognition that some students may discover variants they did not expect
- Group leaders are trained to handle emotional reactions sensitively
- The voluntary nature of sharing reinforces the principle that genomic information is deeply personal

---

## Technical Depth

### VCF Format
While the primary focus of SEQ is on visual variant identification through IGV, students encounter Variant Call Format (VCF) concepts:
- VCF as the standard format for storing variant information
- Fields: chromosome, position, ID (rs number), reference allele, alternate allele, quality, filter, info
- How variant callers produce VCF files from BAM files
- The relationship between what they see in IGV and what would appear in a VCF file

### Variant Calling
- Students perform a simplified version of variant calling: visual identification in IGV
- They look at positions where their reads show non-reference bases
- Heterozygous variants: approximately 50% of reads show the alternate allele
- Homozygous variants: approximately 100% of reads show the alternate allele
- Understanding read depth: more reads at a position = more confidence in the variant call

### Coordinate Systems
The 0-based vs. 1-based coordinate distinction is a genuine stumbling block in bioinformatics, and SEQ addresses it directly:
- **1-based systems** (VCF, SAM, UCSC Genome Browser display): Position 1 is the first base. Used by most human-readable interfaces.
- **0-based systems** (BED format, BAM internal representation, programming languages): Position 0 is the first base. Used by many computational tools.
- Off-by-one errors are one of the most common mistakes in bioinformatics; teaching this early prevents future confusion.

### ClinVar Classification System
Students learn the five-tier classification:
1. **Benign:** Variant is not associated with disease
2. **Likely benign:** Evidence suggests variant is not disease-causing
3. **Uncertain significance (VUS):** Not enough evidence to classify
4. **Likely pathogenic:** Evidence suggests variant causes disease
5. **Pathogenic:** Variant is known to cause disease

The challenge of VUS is emphasized: most variants in any individual's genome are classified as VUS, which means the field's knowledge is still evolving.

### dbSNP and Population Frequencies
- rs numbers as universal variant identifiers
- Minor Allele Frequency (MAF): the frequency of the less common allele in a population
- Populations in gnomAD/dbSNP: African, European, East Asian, South Asian, Latino, etc.
- Why population frequency matters: common variants are usually benign; rare variants are more likely to be pathogenic (but not always)
- Students can see where their variants fall in the frequency spectrum and compare to global populations

---

## Instructor Notes: Teaching Strategies and Pedagogical Approaches

### The Self-Learning Model
SEQ's self-learning quiz format (first 30 minutes of each day, no group leader guidance) is a deliberate pedagogical choice for this advanced track. It develops independent learning skills, gives students time to process complex material at their own pace, and provides a natural assessment mechanism without the pressure of a traditional exam.

### Emotional Management
The instructors are clearly aware that analyzing one's own genome can be emotionally charged. The curriculum strategically begins with non-threatening phenotypes (bitter taste, earwax, blood type) before progressing to anything with health implications. The emphasis on penetrance, professional confirmation, and the right not to share serves as emotional scaffolding, ensuring students feel safe and supported.

### Bridging from BTS
Many BTS skills (Linux, Bowtie2, SAMtools, IGV) are directly reused in SEQ, but applied to personal data. This repetition-with-variation reinforces technical skills while keeping the material fresh. Students who struggled with commands in BTS get a second opportunity to master them in SEQ, and students who found BTS easy are challenged by the larger dataset sizes and more nuanced interpretation required.

### Pioneer Framing
Dr. Lee repeatedly frames SEQ students as pioneers: "You are actually a pioneer. This is really, really just new things we are trying to build together." This framing gives students a sense of significance and ownership over their experience, motivating deeper engagement with the material and the broader CSSI vision.

### Community of Practice
By the third week, students have formed genuine relationships with group leaders and peers. The transcripts show warmer, more relaxed interactions compared to CB Day 1 -- inside jokes, shared memories from previous weeks, and a deeper level of trust. This social capital enables more authentic discussions about the ethical and personal dimensions of genomic data.
