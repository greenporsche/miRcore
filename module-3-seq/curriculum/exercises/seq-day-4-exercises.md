# SEQ Day 4 Exercises: Personal Genome Exploration

## Warm-Up Questions
1. Name three genes (other than TAS2R38) whose variants produce observable phenotypes that you can check against your own traits.
2. What is the difference between looking at variants in a VCF file versus visualizing alignments in IGV? What does each approach tell you?
3. Why is read depth (coverage) important when interpreting whether a variant call is reliable?

## Hands-On Exercises

### Exercise 1: Determining Your ABO Blood Type
**Objective**: Use your VCF data to predict your ABO blood type from exome variants.
**Instructions**:
1. The ABO gene is on chromosome 9. Search your VCF for ABO region variants:
   ```
   grep "chr9" output.vcf | grep -E "133255[0-9]{3}|133256[0-9]{3}|133257[0-9]{3}"
   ```
2. Look for these key variants:
   - **rs8176719** (261delG): A single-base deletion that creates the O allele. If you are homozygous for this deletion, you are blood type O.
   - **rs8176746** and **rs8176747**: These distinguish the A allele from the B allele.
3. Record your ABO genotype findings in this table:

| Variant | rs Number | Your Genotype | Interpretation |
|---|---|---|---|
| 261delG | rs8176719 | | O allele indicator |
| Gly266Ala | rs8176746 | | B vs. A allele |
| Leu266Met | rs8176747 | | B vs. A allele |

4. Based on your variants, predict your blood type: A, B, AB, or O.
5. Do you know your actual blood type? If so, does your prediction match?

**Expected Output**: A completed table with personal genotypes and a blood type prediction. Note: ABO determination from exome data can be complex because some relevant variants may not be captured by exome sequencing.

### Exercise 2: ABCC11 Earwax Type
**Objective**: Find the earwax-determining variant in your data and predict your phenotype.
**Instructions**:
1. Search for the ABCC11 variant rs17822931:
   ```
   grep "rs17822931" output.vcf
   ```
   Or search by position on chromosome 16:
   ```
   grep "chr16" output.vcf | grep "48258198"
   ```
2. Record your genotype:
   - **CC** (homozygous reference): Wet earwax
   - **CT** (heterozygous): Wet earwax (C is dominant)
   - **TT** (homozygous alternate): Dry earwax
3. Look up the population frequency of the T allele in dbSNP (https://www.ncbi.nlm.nih.gov/snp/rs17822931). Record:
   - Global frequency of T allele
   - Frequency in East Asian populations
   - Frequency in European populations
4. Why is the allele frequency so different across populations? What does this tell us about human evolutionary history?
5. Check your own earwax type. Does it match your genotype prediction?

**Expected Output**: Genotype and phenotype prediction. T allele frequency is ~80-95% in East Asian populations but only ~5-20% in European and African populations, suggesting positive selection in East Asian ancestral populations.

### Exercise 3: MMP1 Sun Sensitivity
**Objective**: Investigate the MMP1 promoter variant associated with sun sensitivity and skin aging.
**Instructions**:
1. The MMP1 gene is on chromosome 11. Search for MMP1 variants:
   ```
   grep "chr11" output.vcf | grep -E "rs1799750|10234[0-9]{4}"
   ```
2. The key variant is **rs1799750** (commonly called the 1G/2G polymorphism) in the MMP1 promoter:
   - **1G allele**: Lower MMP1 expression, less collagen breakdown
   - **2G allele**: Higher MMP1 expression, more collagen breakdown, greater sun sensitivity
3. Record your genotype and predict your relative sun sensitivity:
   - 1G/1G: Lower sun sensitivity
   - 1G/2G: Intermediate
   - 2G/2G: Higher sun sensitivity
4. Note: This variant may not be in your VCF if it falls outside the exome capture region (it is in the promoter). If you cannot find it, explain why exome sequencing might miss promoter variants.

**Expected Output**: Genotype if found, or an explanation of why promoter variants can be missed by WES (exome capture targets coding exons, not regulatory regions like promoters).

### Exercise 4: IGV Visualization
**Objective**: Load your personal BAM file into IGV and visually inspect variant positions.
**Instructions**:
1. Download your sorted BAM and BAI files from Great Lakes to your local computer.
2. Open IGV (download from https://software.broadinstitute.org/software/igv/ if needed).
3. Set the reference genome to **hg38** in the dropdown.
4. Load your BAM file: File --> Load from File --> select your sorted.bam (BAI must be in the same folder).
5. Navigate to one of your TAS2R38 variant positions. In the search bar, type the coordinates for rs713598:
   ```
   chr7:141,972,804
   ```
6. Zoom in until you can see individual reads and base colors.
7. Answer:
   - Can you see the variant? What color are the mismatched bases?
   - How many reads cover this position (read depth)?
   - What fraction of reads show the reference allele vs. the alternate allele?
   - Does this visual match the genotype in your VCF?
8. Take a screenshot for your presentation.
9. Navigate to your ABO gene variant position and take another screenshot.

**Expected Output**: IGV screenshots showing variant positions with colored mismatches. For a heterozygous variant, approximately 50% of reads should show the alternate allele.

### Exercise 5: Variant Summary Table
**Objective**: Compile all personal findings into a comprehensive variant summary for presentation preparation.
**Instructions**:
Create a summary table of all the genes and variants you have explored:

| Gene | Chromosome | Variant(s) | Your Genotype | Predicted Phenotype | Confirmed? | Database Evidence |
|---|---|---|---|---|---|---|
| TAS2R38 | chr7 | rs713598, rs1726866, rs10246939 | (your data) | Taster/Non-taster | Y/N/Unknown | ClinVar: Benign |
| ABO | chr9 | rs8176719, rs8176746, rs8176747 | (your data) | Blood type: ? | Y/N/Unknown | dbSNP: common variant |
| ABCC11 | chr16 | rs17822931 | (your data) | Wet/Dry earwax | Y/N/Unknown | ClinVar: Benign |
| MMP1 | chr11 | rs1799750 | (your data) | Sun sensitivity: ? | Y/N/Unknown | dbSNP: common variant |

1. Fill in all personal data rows.
2. Star or highlight any discrepancies between predicted and confirmed phenotypes.
3. For any discrepancies, brainstorm possible reasons (incomplete penetrance, environmental factors, technical errors, complex inheritance).
4. Select 1-2 variants that you are comfortable sharing in your group presentation.

**Expected Output**: A completed summary table ready for use in the Day 5 group presentation.

## Challenge Problems

### Challenge 1: Exploring a Novel Variant
**Objective**: Investigate an unknown variant from your VCF using multiple databases.
**Instructions**:
1. Find a variant in your VCF that has no rs number (the ID column shows `.`).
2. Note its chromosome and position.
3. Look it up in these databases in order:
   a. **UCSC Genome Browser**: Is it in a gene? Which gene? In an exon or intron?
   b. **gnomAD** (https://gnomad.broadinstitute.org): Is it seen in any population? At what frequency?
   c. **ClinVar**: Is it in ClinVar? What is its classification?
4. Based on your research, classify this variant as: likely benign, uncertain, or potentially interesting.
5. Write a one-paragraph summary of your findings, as if you were writing a brief lab report.

### Challenge 2: Population Genetics of Personal Variants
**Objective**: Explore how your personal variant profile compares across global populations.
**Instructions**:
1. For each of the four genes explored today (TAS2R38, ABO, ABCC11, MMP1), look up the allele frequency of your variant(s) in different populations using gnomAD or dbSNP.
2. Create a comparison table:

| Gene | Your Allele | African Freq | European Freq | East Asian Freq | South Asian Freq |
|---|---|---|---|---|---|
| TAS2R38 | | | | | |
| ABO | | | | | |
| ABCC11 | | | | | |
| MMP1 | | | | | |

3. Which of your variants show the greatest frequency differences across populations?
4. Hypothesize why certain alleles might be at different frequencies in different populations (consider: genetic drift, natural selection, founder effects).
5. Discuss: Why is population diversity in genetic databases important for variant interpretation?

---

## Answer Key

### Warm-Up Answers
1. ABO (blood type, chr9), ABCC11 (earwax type, chr16), MMP1 (sun sensitivity, chr11). Other possibilities include MC1R (red hair/fair skin), HERC2/OCA2 (eye color), SLC24A5 (skin pigmentation).
2. VCF shows you a summary of variant calls (position, alleles, genotype, quality). IGV shows you the raw read alignments, letting you visually confirm: how many reads support the variant, whether the variant is on one strand or both, whether the region has good coverage, and whether there are any mapping artifacts. IGV is the "ground truth" check for VCF calls.
3. Read depth matters because variant calls based on few reads are less reliable. If only 3 reads cover a position and 1 shows an alternate allele, it could be a true heterozygous variant or a sequencing error. At 30x coverage with 15 reads showing the alternate, the call is much more confident.

### Exercise Answers
1. **Exercise 1**: ABO blood type determination requires checking multiple variants. rs8176719 (261delG) is the most important -- homozygous deletion = type O. If no deletion, rs8176746 and rs8176747 distinguish A from B alleles. Note: some students may find their exome does not fully cover all ABO variants.
2. **Exercise 2**: CC or CT = wet earwax, TT = dry earwax. The T allele (dry earwax) is at ~80-95% frequency in East Asian populations but only ~5-20% in European/African populations. This extreme frequency difference suggests positive selection, possibly related to reduced apocrine gland secretion in cold, dry climates.
3. **Exercise 3**: Many students will NOT find rs1799750 in their VCF because it is a promoter variant outside the exome capture region. This is an important lesson: WES only captures ~1-2% of the genome (coding regions) and misses regulatory variants, deep intronic variants, and structural variants.
4. **Exercise 4**: In IGV, variant bases appear as colored letters (different color than the gray reference). For heterozygous sites, roughly 50% of reads show the alternate allele. Read depth at exome targets is typically 30-100x.
5. **Exercise 5**: Individual results vary. Common discrepancies arise from: incomplete penetrance (genotype predicts one thing but phenotype is different), environmental modifiers, complex genetics (multiple genes contributing to a trait), or student uncertainty about their actual phenotype.

### Challenge Answers
1. **Challenge 1**: Answers vary by student. The key learning is the workflow: identify variant position --> check gene context (UCSC) --> check population frequency (gnomAD) --> check clinical significance (ClinVar). Most novel variants in coding regions will be rare missense or synonymous changes of uncertain significance.
2. **Challenge 2**: ABCC11 rs17822931 typically shows the greatest population frequency difference (East Asian T allele ~90% vs. African ~5%). This extreme difference is one of the most striking examples of population-specific selection in the human genome. Population diversity in databases matters because a variant that is "rare" in one population may be common in another -- interpreting pathogenicity without population context can lead to misdiagnosis, particularly for underrepresented populations.
