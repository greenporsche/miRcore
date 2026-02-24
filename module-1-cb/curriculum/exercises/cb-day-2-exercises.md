# CB Day 2 Exercises: Gene Expression and Microarrays

## Warm-Up Questions
1. What is the Central Dogma of molecular biology, and what are its three steps?
2. If a gene is "highly expressed" in liver cells but "not expressed" in skin cells, what does that mean at the molecular level?
3. Why would scientists want to measure the expression of thousands of genes at once instead of one at a time?

## Hands-On Exercises

### Exercise 1: Base Pairing Practice
**Objective**: Reinforce understanding of complementary base pairing in DNA and RNA.
**Instructions**:
1. Given the following DNA template strand, write the complementary DNA strand and the mRNA transcript:
   - DNA template: 3'- T A C G G A T T C A A G -5'
2. For the mRNA you wrote, identify which bases changed compared to the complementary DNA strand.
3. If the mRNA had 900 nucleotides (not counting the stop codon), how many amino acids would the resulting protein have?

**Expected Output**: The complementary DNA strand, the mRNA sequence, identification that U replaces T in RNA, and the calculation that 900 / 3 = 300 amino acids.

### Exercise 2: Exploring NCBI GEO
**Objective**: Learn to search for and evaluate gene expression datasets.
**Instructions**:
1. Go to [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/).
2. Search for a disease your group is studying (e.g., "Alzheimer's disease microarray").
3. Find a dataset (GEO Series, starting with "GSE") and record:
   - The GSE accession number
   - The organism studied
   - The number of samples
   - What comparison is being made (e.g., Alzheimer's brain vs. healthy brain)
   - The platform used (e.g., Affymetrix Human Genome U133)
4. Evaluate: Does this dataset have enough samples to give reliable results? (Aim for at least 3 samples per group.)

**Expected Output**: A completed information sheet for one GEO dataset with the GSE number, sample details, and a brief evaluation of its quality.

### Exercise 3: Virtual Microarray Interpretation
**Objective**: Practice reading microarray results.
**Instructions**:
Look at the following simplified microarray data table:

| Gene | Control Signal | Disease Signal | Ratio (Disease/Control) |
|---|---|---|---|
| Gene A | 100 | 800 | 8.0 |
| Gene B | 500 | 510 | 1.02 |
| Gene C | 700 | 175 | 0.25 |
| Gene D | 50 | 48 | 0.96 |
| Gene E | 200 | 1600 | 8.0 |

1. Which genes are up-regulated in the disease group? By how much?
2. Which genes are down-regulated? By how much?
3. Which genes show essentially no change?
4. If you could only study two genes further, which would you pick and why?

**Expected Output**: Identification of Genes A and E as up-regulated (~8x), Gene C as down-regulated (~4x), and Genes B and D as unchanged. Students should pick Genes A/E or Gene C for further study because they show the largest changes.

### Exercise 4: Finding Your Dataset
**Objective**: Select a GEO dataset for your group's research project.
**Instructions**:
1. As a group, agree on your disease of focus (assigned on Day 1).
2. Search GEO for datasets related to your disease.
3. Evaluate at least 3 candidate datasets using these criteria:
   - Human samples (Homo sapiens)
   - Microarray platform (not RNA-seq for this exercise)
   - At least 3 control and 3 disease samples
   - Clear comparison between healthy and diseased tissue
4. Select your best dataset and record the GSE number in your group research document.

**Expected Output**: A comparison of 3 datasets with the selected GSE number and justification for the choice.

## Challenge Problems

### Challenge 1: From Gene to Protein — Full Pathway
**Objective**: Trace the complete path from DNA to a functional outcome.
**Instructions**:
1. Use OMIM to look up the gene **TP53** (tumor protein p53).
2. Describe what protein this gene encodes and what its function is.
3. Explain what happens when TP53 is under-expressed (down-regulated) in cells.
4. If you saw TP53 as significantly down-regulated in a microarray comparing cancer tissue to healthy tissue, what would that tell you about the cancer?
5. Connect this back to the Central Dogma: where in the DNA -> RNA -> Protein pathway could something go wrong to cause TP53 under-expression?

### Challenge 2: Designing a Microarray Experiment
**Objective**: Think like a researcher and design a proper microarray study.
**Instructions**:
You want to study how a new drug affects gene expression in breast cancer cells. Design your experiment:
1. What is your control group? What is your treatment group?
2. How many samples would you use in each group, and why?
3. What would you expect to see on the microarray if the drug is working?
4. What potential problems could arise, and how would you address them?

---

## Answer Key

### Warm-Up Answers
1. The Central Dogma is: DNA -> RNA -> Protein. The three steps are: (1) Replication (DNA copies itself), (2) Transcription (DNA is copied into mRNA), and (3) Translation (mRNA is read by ribosomes to make protein).
2. At the molecular level, the gene's DNA is present in both cell types, but in liver cells the gene is actively transcribed into mRNA and translated into protein, while in skin cells the gene is not being transcribed — the gene is "silent."
3. Diseases are rarely caused by a single gene acting alone. Measuring thousands of genes simultaneously with microarrays reveals patterns, pathways, and networks of genes that are disrupted, providing a much more complete picture of the disease.

### Exercise Answers
1. **Base Pairing**:
   - DNA template: 3'- T A C G G A T T C A A G -5'
   - Complementary DNA: 5'- A T G C C T A A G T T C -3'
   - mRNA transcript: 5'- A U G C C U A A G U U C -3'
   - The difference is that T (thymine) in DNA is replaced by U (uracil) in RNA.
   - 900 nucleotides / 3 nucleotides per codon = 300 amino acids.

2. **GEO Exploration**: Answers will vary. Example: GSE5281 — Alzheimer's disease study, Homo sapiens, 161 samples (87 AD, 74 control), Affymetrix Human Genome U133 Plus 2.0 platform. With 87 disease and 74 control samples, this is a well-powered dataset.

3. **Microarray Interpretation**:
   - Up-regulated: Gene A (8x higher in disease) and Gene E (8x higher in disease)
   - Down-regulated: Gene C (4x lower in disease, ratio 0.25 = 1/4)
   - No change: Gene B (ratio ~1.02) and Gene D (ratio ~0.96)
   - Best genes to study: Gene A or E (strongly up-regulated) and Gene C (strongly down-regulated), because they show the largest expression changes and are most likely to be biologically meaningful.

4. **Dataset Selection**: Answers will vary by group. Key criteria: human samples, microarray platform, sufficient sample size (3+ per group), clear disease vs. control comparison.

### Challenge Answers
1. **TP53 Pathway**: TP53 encodes the p53 protein, known as the "guardian of the genome." It functions as a tumor suppressor by activating DNA repair, halting cell division, or triggering apoptosis (programmed cell death) when DNA damage is detected. When TP53 is down-regulated, cells with damaged DNA are not stopped from dividing, leading to uncontrolled growth (cancer). Seeing TP53 down-regulated in cancer tissue confirms that this tumor suppressor mechanism is impaired. Problems could occur at transcription (mutations in the promoter region preventing mRNA production), at the mRNA level (degradation of mRNA), or at translation (problems producing the protein).

2. **Experiment Design**: Control group = breast cancer cells treated with a placebo/vehicle; Treatment group = breast cancer cells treated with the drug. Use at least 3 biological replicates per group (ideally 5+) to achieve statistical power. If the drug works, you would expect to see changes in gene expression: tumor-promoting genes should be down-regulated and tumor-suppressing genes should be up-regulated. Potential problems include batch effects (run all samples at the same time), RNA quality degradation (use proper handling), and biological variability (use multiple replicates).
