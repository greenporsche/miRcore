# CB Day 2 Study Guide: Gene Expression and Microarrays

## Learning Objectives
By the end of this day, you should be able to:
- Describe the structure of DNA and RNA, including their differences
- Explain the Central Dogma of molecular biology (DNA -> RNA -> Protein)
- Define gene expression and explain why different cells express different genes
- Describe how microarray technology works and what it measures
- Navigate the NCBI Gene Expression Omnibus (GEO) database
- Evaluate and select an appropriate GEO dataset for your group's research
- Explain what a control group and experimental (disease) group are and why both are needed

## Key Vocabulary
| Term | Definition |
|---|---|
| DNA (Deoxyribonucleic acid) | A double-stranded molecule that stores genetic information; uses the bases adenine (A), thymine (T), cytosine (C), and guanine (G) |
| RNA (Ribonucleic acid) | A single-stranded molecule involved in protein synthesis; uses the bases A, uracil (U), C, and G — note U replaces T |
| mRNA (Messenger RNA) | The type of RNA that carries a copy of a gene's instructions from the nucleus to the ribosome |
| Transcription | The process of copying a gene's DNA sequence into an mRNA molecule; occurs in the nucleus |
| Translation | The process of reading the mRNA sequence at the ribosome to assemble a chain of amino acids (protein) |
| Central Dogma | The fundamental principle that genetic information flows from DNA to RNA to Protein |
| Gene Expression | The process by which a gene's information is used to produce a functional product (usually a protein); a gene that is being used is said to be "expressed" |
| Microarray | A laboratory tool consisting of a small chip containing thousands of DNA probes, used to measure the expression levels of many genes simultaneously |
| Probe | A short, known DNA sequence attached to a fixed position on the microarray chip; each probe corresponds to a specific gene |
| Hybridization | The process by which complementary DNA or RNA strands bind together; on a microarray, fluorescently labeled sample RNA binds to matching probes |
| Fluorescence | Light emitted by labeled molecules; on a microarray, brighter fluorescence indicates higher gene expression |
| GEO (Gene Expression Omnibus) | An NCBI public database that archives gene expression data from microarray and sequencing experiments worldwide |
| GSE (GEO Series) | An accession number for a complete experiment/study in GEO (e.g., GSE5281) |
| Control Group | The baseline group in an experiment (e.g., healthy tissue) used for comparison |
| Experimental Group | The group being tested (e.g., diseased tissue) that is compared against the control |
| Nucleotide | The basic building block of DNA and RNA, consisting of a sugar, a phosphate group, and a nitrogenous base |
| Codon | A sequence of three nucleotides in mRNA that codes for a specific amino acid during translation |
| Amino Acid | The building block of proteins; there are 20 standard amino acids |

## Concept Summaries

### DNA and RNA Structure
DNA and RNA are both nucleic acids made up of nucleotides, but they have important structural differences. DNA is double-stranded (the famous double helix), uses the sugar deoxyribose, and contains the bases A, T, C, and G. RNA is typically single-stranded, uses the sugar ribose, and replaces thymine (T) with uracil (U). These structural differences reflect their different roles: DNA serves as the long-term storage of genetic information, while RNA acts as a temporary working copy used to make proteins.

Base pairing is the key to how genetic information is stored and copied. In DNA, adenine (A) always pairs with thymine (T), and cytosine (C) always pairs with guanine (G). During transcription, these pairing rules are used to create an mRNA copy of a gene, except that uracil (U) pairs with adenine (A) instead of thymine. Understanding complementary base pairing is essential because it is also the principle that makes microarray technology work.

### The Central Dogma: DNA -> RNA -> Protein
The Central Dogma describes the flow of genetic information in cells. First, during transcription, one strand of DNA is used as a template to synthesize an mRNA molecule in the nucleus. The mRNA then travels out of the nucleus to the ribosome. During translation, the ribosome reads the mRNA in groups of three nucleotides (codons), and each codon specifies a particular amino acid. The chain of amino acids folds into a functional protein.

This pathway matters for our research because diseases can result from disruptions at any point along this process. A mutation in the DNA can produce a defective mRNA, which produces a defective protein (or no protein at all). By measuring how much mRNA is present for each gene (gene expression), we can identify which genes are over-active or under-active in diseased tissue compared to healthy tissue.

### Gene Expression — Why It Matters
Every cell in your body contains the same DNA, yet a liver cell behaves very differently from a brain cell. The difference is gene expression — which genes are turned "on" (being transcribed into mRNA and translated into protein) and which are turned "off." Gene expression is tightly regulated, and disruptions in this regulation are a hallmark of many diseases, especially cancer.

When we say a gene is "up-regulated" in a disease, we mean more mRNA is being produced from that gene in diseased cells compared to healthy cells. "Down-regulated" means less mRNA is being produced. By comparing expression levels across the entire genome between healthy and diseased samples, we can identify the genes (and biological pathways) most affected by the disease.

### How Microarrays Work
A microarray (also called a "gene chip") is a small glass or silicon chip containing thousands of tiny spots, each with a different known DNA probe. To measure gene expression, researchers extract mRNA from a biological sample, convert it to complementary DNA (cDNA), label it with fluorescent molecules, and allow it to hybridize (bind) to the probes on the chip. Each probe will capture its matching gene's cDNA through complementary base pairing.

After washing away unbound material, a scanner measures the fluorescence intensity at each spot. Spots that glow brightly indicate that the corresponding gene is highly expressed (lots of mRNA was present in the sample). Dim spots indicate low expression. By running microarrays on both a control sample and a disease sample, researchers can compare fluorescence intensities to determine which genes are differentially expressed. This technology allows scientists to measure the expression of 20,000+ genes in a single experiment.

### Navigating NCBI GEO
The Gene Expression Omnibus (GEO) is a public repository hosted by the National Center for Biotechnology Information (NCBI) where researchers worldwide deposit their gene expression data. When you search GEO, you can find thousands of experiments comparing gene expression in different conditions: healthy vs. diseased, treated vs. untreated, young vs. old, and more.

Each experiment is assigned a GEO Series accession number (GSE followed by digits). When evaluating a dataset for your research, you should check: (1) Is it from the correct organism (Homo sapiens for human studies)? (2) Does it use microarray technology? (3) Does it compare disease tissue to healthy control tissue? (4) Are there enough samples in each group (at least 3, ideally more) to provide statistical reliability? Selecting a good dataset is crucial because all of your subsequent analyses depend on the quality of the underlying data.

## How It Connects
- **Previous Day**: On Day 1, you learned about genes, diseases, and used OMIM to find which genes are associated with your disease. Today builds on that foundation by explaining HOW genes produce their effects (through expression) and introducing the technology (microarrays) that measures expression.
- **Next Day**: On Day 3, you will learn the statistical methods needed to determine which gene expression differences are real (significant) versus random noise. You will use GEO2R to analyze the dataset you selected today and produce a list of differentially expressed genes.

## Review Questions
1. List three structural differences between DNA and RNA.
2. What are the complementary base pairs in DNA? How does base pairing change when making mRNA from a DNA template?
3. Explain why a liver cell and a brain cell can have the same DNA but different functions.
4. Describe, step by step, how a microarray experiment works (from sample preparation to reading results).
5. What does it mean if a spot on a microarray chip is very bright? What if it is dim?
6. When searching for a dataset on GEO, what four criteria should you check before selecting it for your research?
7. Why is it essential to have both a control group and a disease group in a microarray experiment?
8. A microarray experiment finds that Gene X produces 10 times more mRNA in cancer tissue than in healthy tissue. In the context of the Central Dogma, what might this mean for the protein encoded by Gene X?

## Further Reading
- NCBI Gene Expression Omnibus: https://www.ncbi.nlm.nih.gov/geo/
- Nature Education — Gene Expression: https://www.nature.com/scitable/topicpage/gene-expression-14121669/
- Khan Academy — Central Dogma: https://www.khanacademy.org/science/biology/gene-expression-central-dogma
- How Microarrays Work (NHGRI): https://www.genome.gov/genetics-glossary/Microarray-Technology
