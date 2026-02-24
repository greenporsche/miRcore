# SEQ Day 3 Study Guide: Variant Identification

## Learning Objectives
By the end of this day, you should be able to:
- Read and interpret VCF file entries, including genotype fields
- Explain what an rs number is and use it to look up variants in ClinVar and dbSNP
- Describe the TAS2R38 gene, its three key SNPs, and the PAV/AVI haplotype system
- Predict PTC bitter taste phenotype from TAS2R38 genotype data
- Use ClinVar to assess clinical significance of a variant (Benign through Pathogenic)
- Explain the concepts of haplotype, penetrance, and genotype-phenotype relationships
- Navigate OMIM to understand gene-disease associations

## Key Vocabulary
| Term | Definition |
|---|---|
| VCF (Variant Call Format) | The standard format for reporting genomic variants; columns include CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, and sample genotypes |
| SNP (Single Nucleotide Polymorphism) | A single base-pair variation at a specific position, the most common type of human genetic variation |
| rs number | A Reference SNP identifier assigned by dbSNP (e.g., rs713598); a stable, unique ID for each cataloged variant |
| Genotype | The allele combination at a locus, encoded in VCF as 0/0 (homozygous reference), 0/1 (heterozygous), or 1/1 (homozygous alternate) |
| Haplotype | A set of alleles at linked loci inherited together on a single chromosome; TAS2R38 has major haplotypes PAV and AVI |
| TAS2R38 | The taste receptor gene on chromosome 7 that determines sensitivity to PTC (phenylthiocarbamide) and similar bitter compounds |
| PTC (phenylthiocarbamide) | A synthetic compound that tastes intensely bitter to people with the PAV haplotype but is tasteless to AVI/AVI non-tasters |
| ClinVar | An NCBI database that aggregates reports on the clinical significance of genomic variants |
| dbSNP | The NCBI database of short genetic variations, assigning rs numbers and cataloging allele frequencies |
| OMIM (Online Mendelian Inheritance in Man) | A comprehensive database of human genes and genetic disorders, focusing on gene-phenotype relationships |
| Pathogenicity classification | The five-tier system for variant clinical significance: Benign, Likely Benign, Variant of Uncertain Significance (VUS), Likely Pathogenic, Pathogenic |
| Penetrance | The proportion of individuals with a particular genotype who express the expected phenotype; incomplete penetrance means not everyone with the variant shows the trait |
| Missense variant | A single nucleotide change that alters the amino acid encoded at that position in the protein |
| Linkage disequilibrium | The non-random association of alleles at different loci, meaning certain allele combinations are inherited together more often than expected by chance |

## Concept Summaries

### Reading Your VCF File
The VCF file is the culmination of the computational pipeline you built over the first two days. It contains every position in your exome where your DNA differs from the reference genome. Each variant line has ten or more columns. The first five (CHROM, POS, ID, REF, ALT) tell you WHERE the variant is and WHAT it is. The QUAL and FILTER columns tell you how confident the variant caller is. The INFO field contains additional annotations. The FORMAT and sample columns tell you YOUR genotype at that position.

Understanding the genotype field is critical. A genotype of `0/1:45:99` with FORMAT `GT:DP:GQ` means: heterozygous (0/1), covered by 45 reads (DP=45), with a genotype quality of 99 (very high confidence). The read depth is particularly important -- a variant call supported by only 2-3 reads is much less reliable than one supported by 45 reads.

### TAS2R38: Your First Gene Exploration
TAS2R38 is the perfect first gene to explore because it has a clear, testable phenotype. This gene encodes a bitter taste receptor on chromosome 7. Three SNPs define two major haplotypes: the taster haplotype PAV (Proline-Alanine-Valine at positions 49, 262, 296) and the non-taster haplotype AVI (Alanine-Valine-Isoleucine). The PAV haplotype produces a functional receptor that strongly detects PTC and related bitter compounds. The AVI haplotype produces a receptor with dramatically reduced function.

If you are PAV/PAV (homozygous taster), PTC tastes extremely bitter to you. If you are PAV/AVI (heterozygous), you can still taste PTC but perhaps less intensely. If you are AVI/AVI (homozygous non-taster), PTC paper tastes like blank paper. This is one of the cleanest examples of a genotype-to-phenotype relationship in human genetics, which is why it has been a teaching staple for decades.

Importantly, these three SNPs are in strong linkage disequilibrium: they are almost always inherited together as a PAV or AVI block. However, rare recombinant haplotypes do exist (e.g., AAV, PVI), which is why examining all three SNPs rather than just one gives the most accurate prediction.

### Using Genomic Databases
This camp introduces you to the major public databases used by geneticists and clinicians worldwide. **ClinVar** is your first stop for understanding clinical significance: is a variant known to cause disease, or is it benign? Variants are classified on the five-tier ACMG scale from Benign to Pathogenic. **dbSNP** catalogs known variants with population allele frequencies, allowing you to see how common your variant is globally and in specific populations. **OMIM** provides deep context on gene function and associated genetic disorders. These databases work together: you find a variant's rs number in dbSNP, check its clinical significance in ClinVar, and understand the gene's role in disease via OMIM.

A critical lesson: just because a variant changes a protein does not mean it is pathogenic. The TAS2R38 variants are missense (they change amino acids) and functional (they alter taste perception), but they are classified as Benign because taste variation is a normal human trait, not a disease.

## How It Connects
- **Previous Day (Day 2)**: You ran the alignment and variant calling pipeline, producing the VCF file that you are now exploring. Understanding VCF quality metrics (QUAL, DP, GQ) depends on understanding the pipeline that generated them.
- **Next Day (Day 4)**: Tomorrow you will expand your exploration to additional genes (ABO blood type, ABCC11 earwax, MMP1 sun sensitivity) and use IGV to visually confirm variants. Today's TAS2R38 exercise establishes the workflow you will repeat with new genes.

## Review Questions
1. What are the ten standard columns in a VCF file? Which columns tell you the variant location, and which tell you your genotype?
2. Explain the difference between 0/0, 0/1, and 1/1 genotypes. What does 1/2 mean?
3. What are the three rs numbers for the TAS2R38 haplotype-defining SNPs? What amino acid change does each cause?
4. A student is heterozygous (0/1) at all three TAS2R38 SNPs. What is the most likely haplotype combination and why? What are the alternative possibilities?
5. You find a variant in your VCF with rs123456. ClinVar classifies it as "Variant of Uncertain Significance (VUS)." What does this mean, and should you be worried?
6. Why is allele frequency important for variant interpretation? If a variant is found in 30% of the population, is it likely to cause a rare disease?
7. Explain penetrance and give an example of a variant with incomplete penetrance.
8. What is linkage disequilibrium, and why are the three TAS2R38 SNPs in strong LD?

## Further Reading
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
- dbSNP: https://www.ncbi.nlm.nih.gov/snp/
- OMIM: https://www.omim.org
- TAS2R38 and PTC tasting (Nature Education): https://www.nature.com/scitable/topicpage/the-genetics-of-bitter-taste-525/
- ACMG/AMP variant interpretation standards: Richards et al., 2015, Genetics in Medicine
- gnomAD population database: https://gnomad.broadinstitute.org
