# SEQ Day 4 Study Guide: Personal Genome Exploration

## Learning Objectives
By the end of this day, you should be able to:
- Determine your predicted ABO blood type from exome variant data
- Identify the ABCC11 earwax variant and predict wet vs. dry earwax phenotype
- Understand the MMP1 promoter variant and its association with sun sensitivity
- Load and navigate BAM files in IGV (Integrative Genomics Viewer)
- Visually inspect variant positions at the read alignment level
- Assess variant reliability by examining read depth and allele balance in IGV
- Compile a personal variant summary table for multiple genes

## Key Vocabulary
| Term | Definition |
|---|---|
| ABO gene | Gene on chromosome 9 encoding the glycosyltransferase that determines ABO blood type through modifications to red blood cell surface antigens |
| Codominance | A pattern of inheritance where both alleles (A and B) are fully expressed in heterozygotes, producing blood type AB |
| Frameshift deletion | A deletion of bases (not a multiple of 3) that shifts the reading frame, usually destroying protein function; rs8176719 in ABO creates the O allele |
| ABCC11 | Gene on chromosome 16 encoding an ATP-binding cassette transporter; the rs17822931 variant determines earwax type and correlates with body odor |
| MMP1 (Matrix Metalloproteinase 1) | Gene on chromosome 11 encoding a collagen-degrading enzyme; a promoter variant (1G/2G polymorphism) affects expression level and sun sensitivity |
| IGV (Integrative Genomics Viewer) | A desktop application from the Broad Institute for interactively visualizing aligned reads, variants, and genomic annotations |
| Read depth (coverage) | The number of sequencing reads overlapping a specific genomic position; higher depth provides more statistical confidence in variant calls |
| Allele balance | The ratio of reads supporting the reference allele versus the alternate allele; for heterozygous variants, expect approximately 50:50 |
| Exon | A segment of a gene that is retained in the mature mRNA after splicing; these are the regions captured by exome sequencing |
| Intron | A non-coding segment between exons that is spliced out during mRNA processing; generally not captured by exome sequencing |
| Promoter | A regulatory DNA region upstream of a gene that controls when and how much the gene is transcribed; usually NOT captured by exome sequencing |
| Population stratification | Differences in allele frequencies between human populations due to evolutionary history, drift, and selection |

## Concept Summaries

### ABO Blood Type from Exome Data
The ABO blood type system is one of the most medically important genetic traits, determining blood transfusion compatibility. The ABO gene on chromosome 9 encodes a glycosyltransferase enzyme that adds sugar molecules to a precursor molecule (H antigen) on red blood cells. The A allele adds N-acetylgalactosamine, the B allele adds galactose, and the O allele makes no functional enzyme due to a single-base deletion (rs8176719, also called 261delG) that shifts the reading frame.

Determining blood type from exome data requires checking multiple variants. First, look for rs8176719: if both copies of your ABO gene carry this deletion, you are type O. If at least one copy does not have the deletion, check rs8176746 and rs8176747, which distinguish the A allele from the B allele. Having one A and one B allele produces type AB (codominance). This exercise is more complex than TAS2R38 because ABO involves multiple variant types (deletion, missense) and a three-allele system (A, B, O) rather than a simple two-haplotype system.

### ABCC11 and Population Genetics
The ABCC11 gene provides a striking example of population-specific genetic variation. A single SNP (rs17822931, C>T) determines whether you produce wet or dry earwax. The C allele (ancestral) produces wet earwax, and the T allele produces dry earwax. Because C is dominant, only TT homozygotes have dry earwax. What makes this variant remarkable is its extreme variation across populations: the T allele frequency is 80-95% in East Asian populations but only 5-20% in European and African populations. This is one of the strongest signals of population-specific natural selection in the human genome, possibly related to thermoregulation in cold environments (reduced apocrine gland secretion).

This example teaches an important lesson about interpreting genetic variants: a variant that appears "rare" in one population may be very common in another. Variant databases that are biased toward European populations can lead to misclassification of variants in underrepresented groups.

### IGV: Seeing Your Reads
Until now, you have been working with summary data (VCF tells you a variant exists, but you cannot see the underlying evidence). IGV changes that. By loading your sorted BAM and BAI files, you can navigate to any position in the genome and see every individual sequencing read. Variant bases appear as colored letters against the gray reference background. For a heterozygous variant, you should see roughly half the reads with the reference allele and half with the alternate.

IGV is the gold standard for visually confirming variant calls. A strong variant call will show: high read depth (20+ reads), balanced allele ratio (close to 50:50 for het, close to 100% for hom alt), variants present on reads from both strands, and no systematic mapping artifacts. A suspicious call might show: low coverage, strand bias (variant only on forward or reverse reads), or the variant appearing only at read ends (possible alignment artifact). Learning to read IGV is an essential skill in genomics.

### MMP1 and the Limits of Exome Sequencing
MMP1 (matrix metalloproteinase 1) is included in today's exploration to teach an important limitation of exome sequencing. The key variant (rs1799750, the 1G/2G polymorphism) is in the gene's PROMOTER region -- the DNA that controls when and how much MMP1 protein is made. Because exome sequencing is designed to capture only coding exons, promoter variants are typically missed. Many students will NOT find this variant in their VCF. This is not a data quality problem -- it is a fundamental limitation of the WES approach. Whole genome sequencing would capture this variant, but at much higher cost. Understanding what your data CAN and CANNOT tell you is a crucial lesson in genomics.

## How It Connects
- **Previous Day (Day 3)**: Yesterday you learned to navigate VCF files and explored TAS2R38. Today you apply the same skills to new genes (ABO, ABCC11, MMP1) and add visual confirmation through IGV. The workflow is the same: find variant in VCF --> look up in databases --> confirm in IGV.
- **Next Day (Day 5)**: Tomorrow you will compile your findings into a group presentation for parents and peers. Today's variant summary table becomes the content backbone of your presentation. You will also discuss the ethical implications of the personal genomic analysis you have performed this week.

## Review Questions
1. Explain how the ABO blood type system works at the molecular level. What is the role of the glycosyltransferase, and how do the A, B, and O alleles differ?
2. Why is determining ABO blood type from exome data more complex than determining TAS2R38 haplotype?
3. The ABCC11 dry earwax allele (T) is at ~90% frequency in East Asian populations but ~10% in European populations. What evolutionary forces could explain this dramatic difference?
4. You load your BAM into IGV and navigate to a variant position. You see 30 reads total, but only 3 show the alternate allele. The VCF says this position is heterozygous (0/1). Should you trust the VCF call? Why or why not?
5. Why might the MMP1 promoter variant (rs1799750) be absent from your VCF? Is this a problem with your data?
6. What is the difference between read depth and mapping quality? Why do both matter for variant interpretation?
7. Explain why allele frequency databases need to include diverse populations. Give an example of how population bias could lead to a misinterpretation.
8. You have now explored four genes (TAS2R38, ABO, ABCC11, MMP1). For which gene did your predicted phenotype match your known trait? For which did it not? What might explain discrepancies?

## Further Reading
- ABO blood group system (NCBI GeneReviews): https://www.ncbi.nlm.nih.gov/books/NBK2267/
- ABCC11 and human earwax: Yoshiura et al., 2006, Nature Genetics (DOI: 10.1038/ng1733)
- IGV User Guide: https://software.broadinstitute.org/software/igv/UserGuide
- Population genetics and variant interpretation: Popejoy & Fullerton, 2016, Nature (DOI: 10.1038/nature19785)
- Limitations of exome sequencing: Meienberg et al., 2016, Human Genetics
