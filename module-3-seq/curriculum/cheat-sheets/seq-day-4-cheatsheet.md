# SEQ Day 4 Cheat Sheet: Personal Genome Exploration

## Commands & Tools
| Tool/Command | What It Does | Example Usage |
|---|---|---|
| IGV (Integrative Genomics Viewer) | Desktop tool to visualize aligned reads at specific positions | Load BAM + BAI, navigate to gene of interest |
| `grep` | Search your VCF for specific gene variants | `grep "chr9:133" output.vcf` (for ABO region) |
| ClinVar search | Look up pathogenicity of personal variants | Search rs number from your VCF |
| dbSNP search | Look up allele frequency of your variants | Search rs number to see global/population frequencies |
| UCSC Genome Browser | Visualize gene structure and variant location | Navigate to ABO, ABCC11, MMP1 gene regions |
| `samtools tview` | Terminal-based alignment viewer | `samtools tview sorted.bam ref.fa` |

## Key Concepts
- **ABO blood type**: Determined by variants in the ABO gene on chromosome 9; three main alleles (A, B, O) with codominance
- **ABCC11 gene**: Determines earwax type (wet vs. dry) and body odor; rs17822931 -- CC/CT = wet earwax, TT = dry earwax
- **MMP1 gene**: Matrix metalloproteinase 1; a promoter variant affects collagen breakdown and sun sensitivity/skin aging
- **Exon**: Protein-coding region of a gene that is retained in mature mRNA after splicing
- **Intron**: Non-coding region between exons that is spliced out; variants here usually have less impact
- **IGV (Integrative Genomics Viewer)**: Desktop application for visualizing BAM files and seeing individual read alignments
- **Read depth / Coverage**: Number of sequencing reads covering a specific position; higher depth = more confidence
- **Allele frequency**: How common a variant is in a population; common variants (>1%) are usually benign
- **Zygosity**: Whether you carry 0, 1, or 2 copies of an alternate allele at a position
- **Compound heterozygote**: Having two different variants on the same gene, one on each chromosome

## Key Genes Explored on Day 4
| Gene | Chromosome | Key Variant(s) | Trait | How to Check |
|---|---|---|---|---|
| ABO | chr9 | rs8176719, rs8176746, rs8176747 | Blood type (A, B, AB, O) | Look up ABO exon 6-7 variants in VCF |
| ABCC11 | chr16 | rs17822931 (C>T) | Earwax type (wet/dry) | Search VCF for chr16:48258198 |
| MMP1 | chr11 | rs1799750 (1G/2G) | Sun sensitivity / skin aging | Search VCF for MMP1 promoter region |
| TAS2R38 | chr7 | rs713598, rs1726866, rs10246939 | Bitter taste (PTC) | Review haplotype from Day 3 |

## ABO Blood Type Determination
| Genotype | Blood Type | Key Variant Pattern |
|---|---|---|
| AA or AO | Type A | Has A-specific variants, no B-specific variants |
| BB or BO | Type B | Has B-specific variants, no A-specific variants |
| AB | Type AB | Has both A and B specific variants |
| OO | Type O | Has frameshift deletion (rs8176719) on both alleles |

## File Formats
- **BAM + BAI**: Required pair for IGV visualization; BAI must be in the same directory as BAM
- **VCF**: Continue using your variant file from Day 2/3 pipeline to look up personal variants

## Databases & URLs
| Resource | URL | Used For |
|---|---|---|
| IGV Download | `https://software.broadinstitute.org/software/igv/` | Download IGV desktop application |
| ClinVar | `https://www.ncbi.nlm.nih.gov/clinvar/` | Variant clinical significance lookup |
| dbSNP | `https://www.ncbi.nlm.nih.gov/snp/` | Allele frequency and population data |
| OMIM | `https://www.omim.org` | Gene-disease associations |
| SNPedia | `https://www.snpedia.com` | Community-curated SNP phenotype information |

## Common Pitfalls
- Not having the BAI index file in the same directory as the BAM file when loading into IGV
- Confusing ABO blood type alleles with simple dominant/recessive -- ABO has codominance (A and B) and recessiveness (O)
- Assuming a single variant determines blood type -- ABO requires checking multiple positions across exons 6-7
- Forgetting that ABCC11 dry earwax (TT genotype) is most common in East Asian populations -- allele frequency varies dramatically by population
- Over-interpreting variants without checking read depth -- low coverage positions may have unreliable genotype calls
- Making medical conclusions from camp data -- these results are for educational purposes only, not clinical grade

## Quick Check
1. What three pieces of information do you need to load your data into IGV?
2. How do you determine your ABO blood type from exome data?
3. What does the rs17822931 variant in ABCC11 determine, and what are the possible phenotypes?
4. Why is read depth important when interpreting a variant you find in your VCF?
5. Why should personal genomic findings from this camp NOT be used for medical decisions?
