# CB Track (Computational Biology) -- Detailed Analysis

> **Note:** The study guides, exercises, and cheat sheets in the `curriculum/` folder are my own self-generated study materials, not official miRcore content. See the [root README](../README.md) for details.

## Track Overview

The Computational Biology (CB) track is the entry point to the miRcore summer camp curriculum. It targets high school students in grades 9 through 12 with no prerequisite bioinformatics experience. Over five days (10 sessions), students progress from foundational biology concepts to independent gene expression analysis using real public datasets. The track is led by Dr. Inhan Lee, with Vina (a Johns Hopkins alumna) serving as camp coordinator and facilitator, supported by a team of undergraduate and graduate student group leaders.

**Goals:**
- Understand how genes relate to diseases and the concept of precision medicine
- Learn to use bioinformatics databases (OMIM, GEO, KEGG) for research
- Gain practical statistics skills for interpreting gene expression data
- Perform a complete GEO2R differential expression analysis on a real dataset
- Communicate scientific findings through group presentations

**Structure:** Each day consists of a morning session (lecture + breakout room discussion) and an afternoon session (lecture + hands-on activity + breakout room work). Students work in small groups of 3--5 with a dedicated group leader throughout the week. The week culminates in group research presentations to peers and parents.

---

## Day-by-Day Breakdown

### Day 1: Genes, Disease, and Precision Medicine

**Morning Session (Part 1)**

*Main Topics:*
- miRcore's mission: connecting the public to scientists, bringing precision medicine closer to reality
- Camp structure: hybrid format (in-person at University of Michigan + Zoom), small groups with group leaders, Zulip for communication
- What is precision medicine? Current medicine gives the same treatment to everyone (some improve, some don't); precision medicine uses biomarkers to match individuals to the right treatment
- Genes as the molecular basis of individual differences
- Chromosomes and base pairs: 23 pairs, approximately 3 billion base pairs, 20,000--25,000 genes
- How genes relate to disease: one gene can be involved in many diseases, one disease can involve many genes
- Ethics of genomic data: genetic discrimination, insurance implications, employment concerns
- The GINA Act (Genetic Information Nondiscrimination Act): protections and limitations
- Data ownership: who controls your genomic data and why it matters

*Key Concepts Introduced:*
- Precision medicine vs. current medicine
- Gene-disease relationships (many-to-many)
- Biomarkers and molecular diagnostics
- Ethical frameworks for genetic information

*Hands-on Activities:*
- Icebreaker activities and group formation
- Discussion questions about precision medicine scenarios in breakout rooms
- Introduction to Zulip messaging platform

**Afternoon Session (Part 2)**

*Main Topics:*
- Introduction to OMIM (Online Mendelian Inheritance in Man)
- How to search OMIM for diseases and related genes
- Navigating OMIM entries: gene symbols, inheritance patterns, clinical descriptions
- Disease selection for group research projects
- Setting up group research infrastructure

*Key Concepts Introduced:*
- OMIM as a curated catalog of genetic disorders
- MIM numbers and how to navigate entries
- Phenotype-genotype relationships in OMIM

*Tools/Databases Used:*
- **OMIM** (https://omim.org): Students search for diseases, browse gene associations, and select a disease for their week-long research project

*Group Work:*
- Each group selects a disease to study for the week
- Groups explore OMIM entries for their chosen disease
- Begin populating cheat sheets with key terms and concepts

---

### Day 2: Gene Expression and Microarrays

**Morning Session (Part 1)**

*Main Topics:*
- Review of Day 1 concepts: genes, diseases, OMIM
- DNA structure: double helix, nucleotides (A, T, G, C), base pairing rules (A-T, G-C)
- RNA structure: single-stranded, uracil replaces thymine (A-U, G-C)
- The central dogma of molecular biology: DNA -> RNA -> Protein
- Transcription: DNA to mRNA in the nucleus
- Translation: mRNA to protein at the ribosome
- Codons: three-nucleotide sequences that code for amino acids
- Gene expression: not all genes are active in all cells at all times
- Why gene expression matters: same genome, different cell types, different expression patterns
- Introduction to microarray technology: measuring the expression of thousands of genes simultaneously

*Key Concepts Introduced:*
- Central dogma (transcription + translation)
- Gene expression as a quantitative, measurable phenomenon
- Microarrays as a high-throughput measurement tool
- The concept of probes hybridizing to complementary mRNA sequences

*Hands-on Activities:*
- Microarray virtual lab: students walk through the steps of a microarray experiment (RNA extraction, labeling, hybridization, scanning, data analysis)
- Breakout room discussions on how gene expression differences could relate to disease

**Afternoon Session (Part 2)**

*Main Topics:*
- Review of microarray virtual lab concepts
- Introduction to GEO (Gene Expression Omnibus)
- How to search GEO for relevant datasets
- Dataset selection criteria: organism (human), experiment type (expression profiling by array), relevant disease/condition, adequate sample sizes, control vs. disease groups
- Understanding GEO accession numbers (GSE series)
- Reading dataset descriptions and associated publications
- Data filtering: selecting appropriate samples within a dataset

*Key Concepts Introduced:*
- GEO as a public repository of gene expression data
- Experimental design: control group vs. treatment/disease group
- Sample selection and data quality considerations

*Tools/Databases Used:*
- **GEO** (https://www.ncbi.nlm.nih.gov/geo/): Browsing and selecting datasets
- Students learn to read GEO dataset pages, understand sample annotations, and identify appropriate comparisons

*Group Work:*
- Each group finds a GEO dataset related to their chosen disease
- Groups read the abstract of the associated research paper
- Group leaders help translate technical jargon into accessible language
- Students begin creating simplified summaries of the research context

---

### Day 3: Statistics and Volcano Plots

**Morning Session (Part 1)**

*Main Topics:*
- Why statistics matter: we need rigorous methods to distinguish real biological signal from random noise
- The null hypothesis: the assumption that there is no difference between groups
- P-value: the probability of observing results at least as extreme as what was measured, assuming the null hypothesis is true
- Common significance thresholds: p < 0.05 (95% confidence), p < 0.01 (99% confidence)
- T-tests: comparing means between two groups
- Gaussian (normal) distribution: the bell curve, mean, and standard deviation
- Standard deviation: measuring the spread of data
- Sample size effects: larger samples give more reliable results (demonstrated with interactive simulator)
- Interactive p-value simulator: students explore how changing sample size, effect size, and variability affects the p-value

*Key Concepts Introduced:*
- Null hypothesis and alternative hypothesis
- P-value interpretation (not "the probability the result is real" but "the probability of seeing this result by chance")
- Statistical significance vs. biological significance
- The relationship between sample size and statistical power

*Hands-on Activities:*
- Interactive p-value simulator exploration
- Breakout room: discussion of statistics concepts, working through example problems
- Practice interpreting p-values from example datasets

**Afternoon Session (Part 2)**

*Main Topics:*
- Log transformation: why we use logarithmic scales for gene expression data
- Log2 fold change: measuring the magnitude of expression differences between groups
  - Log2 fold change of 1 = 2x increase; log2 fold change of -1 = 2x decrease
  - Why log2 (not log10): intuitive doubling/halving interpretation
- Log10 p-value: compressing the p-value scale for visualization
  - Negative log10 transformation: smaller p-values become larger numbers (more significant = higher on the plot)
- Volcano plots: visualizing both fold change (x-axis) and statistical significance (y-axis) simultaneously
  - Reading volcano plots: upper-left = significantly downregulated, upper-right = significantly upregulated
  - Thresholds: horizontal line for p-value cutoff, vertical lines for fold-change cutoff
- Multiple testing problem: when testing thousands of genes, some will appear significant by chance alone
  - With 20,000 genes at p < 0.05, expect ~1,000 false positives
- False Discovery Rate (FDR): the Benjamini-Hochberg correction
  - Adjusted p-value vs. raw p-value
  - Why FDR correction is essential for genome-wide studies
- GEO2R walkthrough: step-by-step analysis
  - Defining groups (control vs. disease)
  - Running the analysis
  - Interpreting the results table and volcano plot

*Key Concepts Introduced:*
- Logarithmic transformations for biological data
- Volcano plots as the standard visualization for differential expression
- The multiple testing problem and FDR correction
- GEO2R as an accessible tool for differential expression analysis

*Tools/Databases Used:*
- **GEO2R**: Students begin using GEO2R on their group's selected dataset, defining comparison groups and running analysis

*Group Work:*
- Continue reading and translating the research paper abstract
- Begin setting up GEO2R analysis for the group dataset
- Discuss initial GEO2R results in breakout rooms

---

### Day 4: GEO2R Analysis and Pathways

**Morning Session (Part 1)**

*Main Topics:*
- Guest speaker: **Dr. Sami Barmada**, neurologist at the University of Michigan Memory Disorders Clinic
  - His research on ALS (amyotrophic lateral sclerosis) and FTD (frontotemporal dementia)
  - TDP-43 protein: found in 95%+ of ALS patients and ~50% of FTD patients
  - The 2006 breakthrough paper that connected ALS and FTD through TDP-43 antibody studies
  - Discussion of causation vs. correlation in disease biology (TDP-43 buildup: cause or consequence?)
  - Genetic therapies and antisense oligonucleotides (ASOs) as emerging treatments
  - Career path discussion: how basic science informs clinical medicine
  - Student Q&A about medical research careers, the research-to-clinic pipeline, and neurological diseases

*Key Concepts Introduced:*
- Translational research: from bench to bedside
- Protein aggregation in neurodegenerative disease
- Genetic therapies (antisense oligonucleotides)
- The interplay between clinical observation and molecular biology

*Continuing GEO2R Analysis:*
- Students complete their group's GEO2R analysis
- Interpreting the top differentially expressed genes table
- Understanding adjusted p-values (FDR) in the results
- Selecting genes of interest from the results

**Afternoon Session (Part 2)**

*Main Topics:*
- Introduction to KEGG (Kyoto Encyclopedia of Genes and Genomes)
- Pathway analysis: understanding how genes work together in biological networks
- Reading KEGG pathway diagrams: nodes (genes/proteins), edges (interactions), colored boxes for different functions
- Connecting GEO2R results to pathways: which pathways are affected by differentially expressed genes
- How to search KEGG for pathways related to the group's disease
- Understanding pathway enrichment: are certain pathways overrepresented among differentially expressed genes
- Research project work: integrating all findings into a coherent narrative

*Key Concepts Introduced:*
- Genes function in networks and pathways, not in isolation
- Pathway analysis as a method for biological interpretation of gene lists
- KEGG as a resource for understanding biological systems

*Tools/Databases Used:*
- **KEGG** (https://www.genome.jp/kegg/): Students explore pathways related to their disease and identify which of their differentially expressed genes appear in relevant pathways

*Group Work:*
- Research project work: combining OMIM background, GEO2R results, and KEGG pathway analysis
- Beginning slide creation for final presentations
- Group leaders guide students in building a narrative: disease background -> gene expression analysis -> pathway interpretation -> implications

*Guest Speaker: Guidance Club Leaders*
- miRcore alumni discuss how to start a Genes in Disease and Symptoms (GIDAS) club at their schools
- Scheduling meetings, teaching content, preparing for the Genes and Health Contest
- Connecting school clubs to the broader miRcore research community

---

### Day 5: Research Presentations

**Morning Session (Part 1)**

*Main Topics:*
- Final daily diagnostic: synthesis question covering the entire week
- Week in review: genes and diseases (Day 1) -> gene expression and microarrays (Day 2) -> statistics and volcano plots (Day 3) -> GEO2R analysis and pathways (Day 4) -> presentations (Day 5)
- Social impact of disease research: each group discusses the broader impact of their disease of focus
- Presentation expectations: clear slides, practice delivery, each group member presents a section
- Presentation structure guidance: introduce the disease, explain the dataset, present GEO2R results (volcano plot, top genes), show KEGG pathway connections, discuss implications

*Hands-on Activities:*
- Final breakout room sessions for presentation preparation
- Adding animations and visuals to slides
- Practice runs within groups

**Afternoon Session (Part 2)**

*Main Topics:*
- Brain reload activity
- Final presentation practice
- Student group presentations to peers and parents
  - Each group presents their disease, dataset, analysis, findings, and implications
  - Peer questions and feedback after each presentation
- Parent showcase: Dr. Lee presents miRcore's vision, history, and outcomes to parents
  - Alumni success stories (Harvard, Stanford, MIT, Johns Hopkins, Yale, Princeton, Columbia, Cornell, UC Berkeley)
  - Explanation of the volunteer program pipeline
  - Genes and Health Contest and research conference opportunities
- Award ceremony: best cheat sheets and presentation recognition
- Introduction to miRcore Volunteer Program (MVP): how students can continue their research journey
- Citizen Scientist Sequencing Initiative (CSSI): Dr. Lee's vision for individuals owning their own genome data

---

## Learning Arc

The CB track follows a deliberate five-day arc:

```
Day 1: QUESTION    -> What is the relationship between genes and disease?
Day 2: MEASUREMENT -> How do we measure gene expression? (microarrays, GEO)
Day 3: STATISTICS  -> How do we know if differences are real? (p-values, volcano plots)
Day 4: INTERPRETATION -> What do the results mean biologically? (pathways, KEGG)
Day 5: COMMUNICATION -> How do we share findings with others? (presentations)
```

This arc mirrors the actual scientific research process: start with a question, identify appropriate measurements, apply statistical rigor, interpret results in biological context, and communicate findings. Students experience each stage of this process within a single week, gaining an integrated understanding of how bioinformatics research works in practice.

---

## Key Databases: How Each Is Used and When

### OMIM -- Day 1 (entry point)
Students use OMIM as their first research tool. They search for diseases they are interested in, browse associated genes, read clinical descriptions, and select a disease for their week-long project. OMIM provides the biological context: what genes are involved, how they are inherited, and what is known about the molecular basis of the disease. This gives students a foundation before they begin looking at expression data.

### GEO -- Days 2--3 (data source)
GEO is introduced on Day 2 as the source of the gene expression data students will analyze. Students learn to search for datasets relevant to their disease, evaluate datasets based on experimental design (control vs. disease samples, adequate sample sizes, human subjects), and read the associated publication abstracts. The dataset selected from GEO becomes the raw material for the rest of the week's analysis.

### GEO2R -- Days 3--4 (analysis engine)
GEO2R is the primary analytical tool. Introduced late on Day 3 after students have learned the underlying statistics, GEO2R allows students to define comparison groups within their chosen GEO dataset and run a differential expression analysis. The tool produces volcano plots, ranked gene lists with fold changes and adjusted p-values, and downloadable results tables. Students spend parts of Days 3 and 4 working with GEO2R, learning to interpret its outputs and select genes of interest.

### KEGG -- Day 4 (biological interpretation)
KEGG is introduced on Day 4 as the final analytical layer. After students have identified differentially expressed genes through GEO2R, KEGG helps them understand what those genes do in the context of biological pathways. Students search for pathways related to their disease and identify which of their significant genes appear in those pathways, enabling them to build a more complete biological narrative for their presentations.

---

## Assessment Strategy

The CB track uses multiple assessment approaches, none of which are traditional exams.

### Cheat Sheets
Students maintain running reference documents throughout the week, recording key concepts, definitions, tools, and procedures. These cheat sheets serve both as a learning tool (the act of writing reinforces understanding) and as an assessment artifact. At the end of the week, awards are given for the most thorough, well-organized, and aesthetically pleasing cheat sheets. The instructors explicitly encourage creative, comprehensive cheat sheets from Day 1.

### Daily Diagnostics
Each morning begins with a brief diagnostic exercise -- typically a discussion question or reflection prompt related to the previous day's content. These serve as formative assessments, helping instructors gauge understanding and identify concepts that need reinforcement. They are low-stakes and discussion-oriented rather than graded.

### Group Presentations (Day 5)
The primary summative assessment is the group research presentation. Each group presents a 10--15 minute talk covering:
1. Disease background (from OMIM research)
2. Dataset description (from GEO)
3. Analysis methodology (GEO2R setup and execution)
4. Results (volcano plot, top genes, statistical significance)
5. Pathway connections (from KEGG)
6. Social and medical implications

Every group member must present a section, requiring all students to understand and articulate the material. Peers ask questions after each presentation, adding an element of real-time scientific discourse.

### Peer Feedback (from Final Presentations Transcript)
During the parent showcase, students hear each other's presentations and ask questions. The Q&A component tests students' depth of understanding beyond their prepared slides. The supportive group environment -- cultivated all week through brain reloads, lunch polls, and breakout room bonding -- ensures this feels more like a celebration of learning than a high-stakes evaluation.

---

## Instructor Notes: Teaching Strategies and Pedagogical Approaches

### Scaffolded Complexity
The curriculum is carefully scaffolded. Day 1 concepts (genes, diseases) are intuitive and accessible. Statistics on Day 3 is acknowledged as potentially challenging, with instructors noting the range of math experience in the room and encouraging students to ask questions. Technical vocabulary is introduced gradually and reinforced through the cheat sheet system.

### Translating Academic Language
A notable pedagogical practice observed in the sessions: group leaders help students translate research paper abstracts from technical jargon into "regular words." This explicit exercise in scientific literacy -- taking dense academic language and making it comprehensible -- is a valuable skill that most students do not encounter until much later in their education.

### Affective Support
Instructors invest significant time in community building. Daily brain reloads (Zulip polls about lunch, fun questions) serve a dual purpose: they get students active on the communication platform and they create a warm, personal atmosphere. The instructors' genuine concern for students' wellbeing -- including checking in on those who skip meals -- exemplifies the caring culture. This affective support reduces the anxiety that can accompany challenging technical content.

### Small-Group Mentorship Model
The breakout room structure -- 3--5 students per group leader -- ensures personalized attention. Group leaders are miRcore alumni, typically undergraduates at research universities, who understand the material deeply but are close enough in age to relate to the students. They function as peer mentors rather than traditional instructors, creating a less intimidating learning environment.

### Interactive Engagement
Dr. Lee and the facilitators use frequent interactive elements: polls, chat responses, thumbs-up checks, "roller coaster" activities to re-energize the room, and real-time Q&A during lectures. The hybrid format (in-person + Zoom) presents challenges, and the instructors visibly adapt -- asking online students to use chat, monitoring breakout rooms, and adjusting pacing based on visual feedback from camera feeds.

### Real-World Connections
The guest speaker on Day 4 (Dr. Barmada) connects classroom concepts to active clinical and research work, showing students that the skills they are learning have direct applications in medical research. The parent showcase highlights alumni who have gone on to Harvard, Stanford, MIT, and other institutions, creating aspirational role models. The emphasis on precision medicine throughout the week grounds abstract concepts in a compelling future vision.

### Growth Mindset Framing
The camp consistently frames learning as a journey. The Day 5 recap explicitly lists "all the things we did in just the past four days" to help students appreciate their own growth. The pre/post diagnostic structure (though more prominent in SEQ) reinforces the idea that learning is measurable and worth celebrating.
