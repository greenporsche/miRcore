# BTS Day 5 Exercises: Research Projects and Pitches

## Warm-Up Questions
1. What are the key steps in the bioinformatics pipeline you have used throughout this camp (from raw data to visualization)?
2. How can you determine from your alignment data whether a patient sample is positive for SARS-CoV-2?
3. What are the essential components of a successful biotech business pitch?

## Hands-On Exercises

### Exercise 1: Compiling Patient Results
**Objective**: Create a summary table of all 8 patient alignment results.
**Instructions**:
1. Review the Bowtie2 output (`.out` files) for each patient's SARS-CoV-2 alignment
2. For each patient, record:
   - Total reads processed
   - Number of reads aligned to SARS-CoV-2
   - SARS-CoV-2 alignment rate (%)
   - Number of reads aligned to HG38 (if available)
   - HG38 alignment rate (%)
3. Create a table in a text file using nano:
   ```
   Patient | Total Reads | SCV2 Rate | HG38 Rate | Diagnosis
   --------|-------------|-----------|-----------|----------
   1       | ...         | ...%      | ...%      | +/-
   2       | ...         | ...%      | ...%      | +/-
   ...
   ```
4. Based on the alignment rates, classify each patient as positive (+) or negative (-) for SARS-CoV-2

**Expected Output**: A completed table with all 8 patients and their diagnoses. Patients with high SARS-CoV-2 alignment rates should be classified as positive.

### Exercise 2: Mutation Identification in IGV
**Objective**: Use IGV to identify specific mutations in SARS-CoV-2-positive patient samples.
**Instructions**:
1. Open IGV and load the SARS-CoV-2 reference genome
2. Load a sorted BAM file from a SARS-CoV-2-positive patient
3. Navigate to the spike protein region: `NC_045512.2:21,563-25,384`
4. Zoom in until you can see individual bases in the reads
5. Look for colored bases (mismatches from the reference) --- these are potential mutations
6. For each mutation found, record:
   - Genomic position
   - Reference base
   - Variant base
   - Approximate percentage of reads showing the variant
7. Compare mutations across multiple positive patients --- do they share the same mutations?

**Expected Output**: A list of identified mutations with their positions and alleles. Shared mutations across patients may indicate a common viral lineage or variant of concern.

### Exercise 3: Research Findings Summary
**Objective**: Write a concise research summary of your bioinformatics analysis.
**Instructions**:
1. Create a file called `research_summary.txt` using nano
2. Include the following sections:
   - **Background**: Why are we analyzing these patient samples? (2-3 sentences)
   - **Methods**: What tools and pipeline did we use? (List the steps)
   - **Results**: Which patients tested positive? What mutations were found? (Reference your data)
   - **Conclusion**: What do the results suggest? (2-3 sentences)
3. Keep the entire summary under one page (approximately 250 words)

**Expected Output**: A structured research summary that could be presented to parents or peers.

### Exercise 4: Business Pitch Preparation
**Objective**: Finalize your group's biotech business pitch for the investment simulation.
**Instructions**:
1. Review and refine your pitch slides to include:
   - **Slide 1**: Company name, logo, tagline
   - **Slide 2**: The problem you are solving
   - **Slide 3**: Your solution / product
   - **Slide 4**: Target market and audience
   - **Slide 5**: Competitive advantage
   - **Slide 6**: Team (your group members and their roles)
   - **Slide 7**: The Ask (how much investment you need and what you will use it for)
2. Practice delivering the pitch in under 7 minutes
3. Prepare for Q&A --- anticipate 2-3 questions investors might ask

**Expected Output**: A polished slide deck and a rehearsed 5-7 minute presentation.

### Exercise 5: Investment Evaluation
**Objective**: Critically evaluate other groups' business pitches as a mock investor.
**Instructions**:
1. As each group presents, take notes on:
   - Clarity of the problem statement
   - Viability of the proposed solution
   - Size and accessibility of the target market
   - Strength of the competitive advantage
   - Quality of the presentation delivery
2. Rate each pitch on a scale of 1-5 for each criterion
3. Decide which group(s) you would invest your play money in
4. Write a brief justification (2-3 sentences) for your investment decision

**Expected Output**: A completed evaluation sheet with ratings and a written investment justification.

## Challenge Problems

### Challenge 1: Variant Significance Research
**Objective**: Investigate the biological significance of mutations you identified.
**Instructions**:
1. Take one or more mutations you found in Exercise 2
2. Use NCBI or published literature to determine:
   - Is this mutation found in known SARS-CoV-2 variants of concern (Alpha, Delta, Omicron)?
   - Does this mutation occur in the spike protein's receptor binding domain?
   - Could this mutation affect how the virus binds to human ACE2 receptors?
3. Write a short paragraph (100-150 words) explaining the potential significance of the mutation

### Challenge 2: Complete Pipeline Documentation
**Objective**: Create a comprehensive documentation of the entire BTS pipeline that a new student could follow.
**Instructions**:
1. Write a step-by-step guide covering:
   - Logging into Great Lakes
   - Navigating to the data directory
   - Running the alignment pipeline (all commands)
   - Submitting jobs with SLURM
   - Downloading results
   - Visualizing in IGV
   - Interpreting results
2. Include example commands, expected outputs, and common troubleshooting tips
3. The guide should be detailed enough that someone with no prior experience could follow it

---

## Answer Key

### Warm-Up Answers
1. The full pipeline: (1) Start with raw FASTQ reads from NGS, (2) Align reads to a reference genome using Bowtie2, (3) Convert the resulting SAM file to BAM using `samtools view`, (4) Sort the BAM file by position using `samtools sort`, (5) Index the sorted BAM using `samtools index`, (6) Visualize alignments in IGV to examine coverage and identify mutations.
2. Compare alignment rates: map the patient's reads to the SARS-CoV-2 reference genome and to the human genome (HG38). If a high percentage of reads align to SARS-CoV-2 (e.g., >80%), the sample likely contains viral RNA. If nearly all reads align to HG38 and very few to SARS-CoV-2, the patient is likely negative.
3. Essential pitch components: (1) Clear problem statement, (2) Proposed solution with scientific backing, (3) Target market identification and size, (4) Competitive advantage or unique value proposition, (5) Team qualifications, (6) Specific investment ask and how funds will be used.

### Exercise Answers
1. **Patient Results**: The summary table should show clear differentiation between positive and negative patients based on SARS-CoV-2 alignment rates. A threshold of ~50% alignment rate is a reasonable cutoff, though actual values will depend on the dataset.
2. **Mutation Identification**: Mutations appear as colored bases in IGV reads. A true mutation should be present in a high percentage of reads (>50%) at that position. Low-frequency variants could be sequencing errors. Common SARS-CoV-2 mutations include D614G in the spike protein.
3. **Research Summary**: The summary should follow standard scientific structure: Background (why), Methods (how), Results (what), Conclusion (so what). Key findings include which patients are positive and any notable mutations.
4. **Business Pitch**: A strong pitch stays within the time limit, clearly articulates the problem-solution fit, quantifies the market opportunity, and makes a specific, justified investment ask. Practice is essential for smooth delivery.
5. **Investment Evaluation**: Good investors evaluate both the business fundamentals and the presentation quality. A strong evaluation considers whether the market is large enough, the solution is feasible, and the team can execute.

### Challenge Answers
1. **Variant Significance**: Students should connect their observed mutations to known variants. For example, the D614G mutation (aspartate to glycine at position 614 of the spike protein) is found in virtually all circulating SARS-CoV-2 lineages and is associated with increased transmissibility.
2. **Pipeline Documentation**: The guide should be step-by-step and reproducible. Common troubleshooting tips include: checking that modules are loaded, verifying file paths, ensuring BAM is sorted before indexing, and confirming both BAM and BAI files are transferred for IGV.
