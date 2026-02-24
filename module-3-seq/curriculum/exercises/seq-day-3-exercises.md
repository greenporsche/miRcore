# SEQ Day 3 Exercises: Variant Identification

## Warm-Up Questions
1. What is a VCF file, and what are the key columns it contains?
2. Explain the difference between a genotype of 0/0, 0/1, and 1/1 in VCF format.
3. What is a haplotype? How is it different from a single SNP genotype?

## Hands-On Exercises

### Exercise 1: Exploring Your VCF File
**Objective**: Navigate and interpret your personal VCF file from the variant calling pipeline.
**Instructions**:
1. View the header of your VCF file:
   ```
   grep "^##" output.vcf | head -20
   ```
2. View the column headers:
   ```
   grep "^#CHROM" output.vcf
   ```
3. Count the total number of variants called:
   ```
   grep -v "^#" output.vcf | wc -l
   ```
4. View the first 5 variant lines:
   ```
   grep -v "^#" output.vcf | head -5
   ```
5. For each of those 5 variants, identify: chromosome, position, reference allele, alternate allele, and your genotype.

**Expected Output**: A count of total variants (typically 20,000-40,000 for WES) and a parsed table of 5 variants with their key fields.

### Exercise 2: Finding TAS2R38 Variants
**Objective**: Locate the three key TAS2R38 SNPs in your VCF and determine your bitter taste haplotype.
**Instructions**:
1. Search for TAS2R38 variants in your VCF:
   ```
   grep "chr7" output.vcf | grep -E "141972[0-9]{3}|141973[0-9]{3}"
   ```
   Or search by rs numbers:
   ```
   grep -E "rs713598|rs1726866|rs10246939" output.vcf
   ```
2. For each of the three SNPs, record:

| SNP | rs Number | Your Genotype (0/0, 0/1, or 1/1) | Your Alleles |
|---|---|---|---|
| SNP 1 (P49A) | rs713598 | | |
| SNP 2 (A262V) | rs1726866 | | |
| SNP 3 (V296I) | rs10246939 | | |

3. Based on your three genotypes, determine your TAS2R38 haplotype combination:
   - PAV = taster allele (Pro-Ala-Val at positions 49, 262, 296)
   - AVI = non-taster allele (Ala-Val-Ile)
4. Predict your PTC bitter taste phenotype: Can you taste PTC paper?
5. If you have done the PTC taste test strip in a previous camp, does your prediction match your actual experience?

**Expected Output**: A completed table with genotypes for all three SNPs, haplotype determination (PAV/PAV, PAV/AVI, or AVI/AVI), and phenotype prediction.

### Exercise 3: ClinVar Variant Lookup
**Objective**: Use ClinVar to investigate the clinical significance of variants in your VCF.
**Instructions**:
1. Go to https://www.ncbi.nlm.nih.gov/clinvar/
2. Search for `rs713598` (the first TAS2R38 SNP).
3. Record:
   - Clinical significance classification
   - Associated condition(s)
   - Review status (how many stars?)
   - Allele frequency in different populations
4. Now search for `rs1726866` and `rs10246939`. Record the same information.
5. Are any of these TAS2R38 variants classified as pathogenic? Why or why not?
6. Pick one OTHER variant from your VCF (choose one with an rs number in the ID column). Look it up in ClinVar. What did you find?

**Expected Output**: ClinVar entries showing TAS2R38 variants as "Benign" or "drug response" (related to PTC taste perception, not a disease). Understanding that common functional variants are not necessarily pathogenic.

### Exercise 4: Genotype-to-Phenotype Mapping
**Objective**: Connect VCF genotype data to observable traits using the TAS2R38 example.
**Instructions**:
Complete the following table based on your understanding of TAS2R38:

| Student (Example) | rs713598 Genotype | rs1726866 Genotype | rs10246939 Genotype | Haplotype(s) | Predicted Phenotype |
|---|---|---|---|---|---|
| Alice | 0/1 (C/G) | 0/1 (T/C) | 0/1 (C/T) | PAV/AVI | Taster |
| Bob | 1/1 (G/G) | 1/1 (C/C) | 1/1 (T/T) | AVI/AVI | Non-taster |
| You | ? | ? | ? | ? | ? |

1. Fill in your own row.
2. In your breakout room group, share your haplotypes (if comfortable). How many tasters vs. non-tasters does your group have?
3. The frequency of PAV in European populations is about 45%. If we assume Hardy-Weinberg equilibrium, what percentage of people would be: PAV/PAV, PAV/AVI, AVI/AVI?

**Expected Output**: Completed personal row; group tally; Hardy-Weinberg calculation: p=0.45, q=0.55, so PAV/PAV = 20.25%, PAV/AVI = 49.50%, AVI/AVI = 30.25%.

### Exercise 5: Variant Filtering and Prioritization
**Objective**: Learn to filter VCF files to find variants of interest.
**Instructions**:
1. Count how many variants are on each chromosome:
   ```
   grep -v "^#" output.vcf | awk '{print $1}' | sort | uniq -c | sort -rn
   ```
2. Which chromosome has the most variants? Does this correlate with chromosome size?
3. Filter for only heterozygous variants:
   ```
   grep -v "^#" output.vcf | grep "0/1" | wc -l
   ```
4. Filter for only homozygous alternate variants:
   ```
   grep -v "^#" output.vcf | grep "1/1" | wc -l
   ```
5. What is the ratio of heterozygous to homozygous alternate variants? For exome data, a ratio of approximately 1.5-2.0 is expected. How does your ratio compare?

**Expected Output**: Variant counts by chromosome (chr1 typically has the most); het/hom ratio close to 1.5-2.0, indicating good data quality.

## Challenge Problems

### Challenge 1: Haplotype Phase Ambiguity
**Objective**: Understand why determining haplotype phase from genotype data is not always straightforward.
**Instructions**:
Consider a person who is heterozygous (0/1) at all three TAS2R38 SNPs.
1. List ALL possible haplotype combinations that could produce this genotype pattern.
   (Hint: There are more possibilities than just PAV/AVI.)
2. In practice, why is PAV/AVI the most likely interpretation?
3. What technique could you use to definitively determine the phase? (Hint: Think about family members or long-read sequencing.)
4. Explain what "linkage disequilibrium" means and why TAS2R38 SNPs are in strong LD.

### Challenge 2: Ethical Scenario - Unexpected Findings
**Objective**: Grapple with the ethical dimensions of personal genome analysis.
**Instructions**:
Read the following scenario and answer the questions:

*While exploring your VCF file, you find a variant in the BRCA2 gene that ClinVar classifies as "Likely Pathogenic" for hereditary breast and ovarian cancer. You are 16 years old.*

1. Should you tell your parents about this finding? Why or why not?
2. Is this finding from your camp exome data reliable enough to make medical decisions? Why not?
3. What is the appropriate next step if someone finds a potentially serious variant in their research-grade data?
4. What does "likely pathogenic" mean, and how is it different from "pathogenic"?
5. Discuss the concept of penetrance: If a variant is pathogenic, does that mean you will definitely develop the disease?
6. How does GINA (Genetic Information Nondiscrimination Act) protect you, and what are its limitations?

---

## Answer Key

### Warm-Up Answers
1. VCF (Variant Call Format) records positions where the sample differs from the reference genome. Key columns: CHROM (chromosome), POS (position, 1-based), ID (rs number if known), REF (reference allele), ALT (alternate allele), QUAL (variant quality score), FILTER (pass/fail), INFO (additional annotations), FORMAT (genotype format), and the SAMPLE column (actual genotype data).
2. 0/0 = homozygous reference (both alleles match the reference). 0/1 = heterozygous (one reference allele, one alternate allele). 1/1 = homozygous alternate (both alleles are the alternate).
3. A haplotype is a specific combination of alleles at multiple linked loci inherited together on one chromosome. A single SNP genotype only tells you about one position. The TAS2R38 haplotype combines three SNPs (PAV or AVI) because they are inherited as a block.

### Exercise Answers
1. **Exercise 1**: Variant counts typically range from 20,000-40,000 for WES. Students should correctly identify CHROM, POS, REF, ALT, and genotype fields from VCF lines.
2. **Exercise 2**: Answers are personal to each student. The key learning is connecting three genotypes to a haplotype prediction. PAV/PAV or PAV/AVI = taster; AVI/AVI = non-taster.
3. **Exercise 3**: TAS2R38 variants in ClinVar are typically classified as "Benign" for disease but "drug response" for taste perception. They are not pathogenic because bitter taste variation is a normal human trait, not a disease.
4. **Exercise 4**: Hardy-Weinberg with p(PAV)=0.45: PAV/PAV = (0.45)^2 = 20.25%, PAV/AVI = 2(0.45)(0.55) = 49.50%, AVI/AVI = (0.55)^2 = 30.25%. In a class of 25 students, expect roughly 5 PAV/PAV, 12 PAV/AVI, 8 AVI/AVI (with variation).
5. **Exercise 5**: Chromosome 1 typically has the most variants because it is the largest chromosome. Expected het/hom ratio is 1.5-2.0 for exome data. Ratios significantly outside this range may indicate data quality issues, consanguinity (very low ratio), or contamination (very high ratio).

### Challenge Answers
1. **Challenge 1**: Possible haplotype pairs for triple heterozygous: PAV/AVI (most common), PAI/AVC, PVC/AAI, etc. PAV/AVI is most likely because these three SNPs are in very strong linkage disequilibrium -- they are inherited together >95% of the time. Definitively determining phase requires either parental genotyping (trio analysis), long-read sequencing that spans all three positions, or statistical phasing with a reference panel. Linkage disequilibrium means that alleles at different positions are inherited together more often than expected by chance, due to physical proximity on the same chromosome and limited historical recombination between them.
2. **Challenge 2**: (1) You should discuss with parents, but recognize this requires sensitivity -- the variant has implications for them too. (2) No -- camp data is research-grade, not clinical-grade. It has not been processed in a CLIA-certified lab and has not been validated. (3) The appropriate step is to consult a genetic counselor and potentially get confirmatory clinical testing. (4) "Likely pathogenic" means there is strong evidence (>90% certainty) that the variant causes disease, but it has not met the full threshold for "pathogenic." (5) Penetrance for BRCA2 is incomplete -- not everyone with a pathogenic BRCA2 variant will develop cancer. Lifetime risk is elevated but not 100%. (6) GINA prohibits health insurance companies and employers from discriminating based on genetic information, but it does NOT cover life insurance, disability insurance, or long-term care insurance.
