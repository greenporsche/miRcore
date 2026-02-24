# CB Day 2 Cheat Sheet: Gene Expression and Microarrays

## Commands & Tools
| Tool/Command | What It Does | Example Usage |
|---|---|---|
| GEO (Gene Expression Omnibus) | NCBI database of gene expression datasets from microarray and sequencing experiments | Search for a disease to find related expression datasets |
| Virtual Microarray Lab | Interactive simulation to learn how microarrays work step by step | Walk through hybridization, scanning, and data output |
| Google Sheets | Spreadsheet tool for organizing and sorting gene expression data | Sort genes by p-value or log fold change |

## Key Concepts
- **DNA**: Deoxyribonucleic acid — the double-stranded molecule that stores genetic information; uses bases A, T, C, G
- **RNA**: Ribonucleic acid — single-stranded molecule transcribed from DNA; uses bases A, U, C, G (uracil replaces thymine)
- **mRNA (messenger RNA)**: The RNA copy of a gene that carries instructions from DNA to the ribosome for protein synthesis
- **Transcription**: The process of copying DNA into mRNA in the nucleus
- **Translation**: The process of reading mRNA to build a protein at the ribosome
- **Central Dogma**: DNA -> RNA -> Protein (the flow of genetic information)
- **Gene Expression**: The process by which information from a gene is used to make a functional product (protein); genes can be "turned on" (expressed) or "turned off"
- **Microarray**: A technology that measures expression levels of thousands of genes simultaneously using a chip with DNA probes
- **Hybridization**: The process where complementary DNA/RNA strands bind together on the microarray chip
- **Probe**: A short DNA sequence fixed on the microarray chip that matches a known gene
- **Fluorescence**: The light signal emitted when labeled RNA binds to probes; brighter = more expression
- **Control vs. Disease Sample**: Healthy tissue (control) compared to diseased tissue to find differentially expressed genes
- **Nucleotide**: The building block of DNA and RNA, consisting of a sugar, phosphate group, and nitrogenous base

## File Formats
- **GEO Dataset files**: Expression data tables available for download from NCBI GEO
- **Spreadsheet (.xlsx/.csv)**: Used to organize and sort gene expression results

## Databases & URLs
| Resource | URL | Used For |
|---|---|---|
| NCBI GEO | https://www.ncbi.nlm.nih.gov/geo/ | Searching for gene expression datasets |
| NCBI GEO Datasets | https://www.ncbi.nlm.nih.gov/gds/ | Browsing curated datasets |
| OMIM | https://omim.org | Cross-referencing genes found in expression data |

## Common Pitfalls
- Confusing DNA and RNA bases — remember RNA uses Uracil (U) instead of Thymine (T)
- Thinking microarrays sequence DNA — they measure expression levels, not sequences
- Not understanding that brighter fluorescence = higher gene expression
- Choosing a GEO dataset without checking sample size — small sample sizes give less reliable results
- Forgetting that complementary base pairing is what makes microarrays work (A-T/U, C-G)
- Not filtering GEO search results properly — make sure to select datasets comparing disease vs. control

## Quick Check
1. What are the three steps of the Central Dogma?
2. What base replaces Thymine in RNA?
3. How does a microarray measure gene expression?
4. Why do we need both a control group and a disease group in a microarray experiment?
5. What does it mean if a spot on a microarray chip is very bright?
