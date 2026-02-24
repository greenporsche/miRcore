# SEQ Day 3 Cheat Sheet: Variant Identification

## Commands & Tools
| Tool/Command | What It Does | Example Usage |
|---|---|---|
| `grep` | Search for specific variants or lines in a VCF file | `grep "rs713598" output.vcf` |
| `less` | Page through large VCF files | `less output.vcf` |
| `awk` | Extract specific columns from VCF | `awk '{print $1, $2, $4, $5}' output.vcf` |
| `wc -l` | Count total variant lines in a VCF | `grep -v "^#" output.vcf \| wc -l` |
| UCSC Genome Browser search | Look up gene location by name or coordinates | Search: `TAS2R38` or `chr7:141,972,631-141,973,841` |
| ClinVar web search | Look up clinical significance of a variant | Search by rs number: `rs713598` |
| dbSNP web search | Look up known SNP details | Search by rs number: `rs1726866` |
| OMIM web search | Look up gene-disease associations | Search by gene name: `TAS2R38` |

## Key Concepts
- **SNP (Single Nucleotide Polymorphism)**: A single base-pair change at a specific position; the most common type of genetic variation
- **rs number**: Reference SNP identifier from dbSNP (e.g., rs713598); a unique ID for each known variant
- **Genotype**: The combination of alleles at a locus; reported as 0/0 (homozygous reference), 0/1 (heterozygous), or 1/1 (homozygous alternate)
- **Haplotype**: A set of alleles inherited together on one chromosome; TAS2R38 has two major haplotypes: PAV (taster) and AVI (non-taster)
- **Phenotype**: The observable trait resulting from genotype (e.g., ability to taste PTC bitterness)
- **Pathogenicity**: Classification of a variant's clinical impact: Benign, Likely Benign, Uncertain Significance (VUS), Likely Pathogenic, Pathogenic
- **Penetrance**: The proportion of individuals with a given genotype who actually express the associated phenotype
- **Heterozygous advantage**: When carrying one copy of a variant confers a selective benefit (e.g., sickle cell trait and malaria resistance)
- **Missense variant**: A nucleotide change that results in a different amino acid in the protein
- **cDNA notation**: Naming convention for variants at the coding DNA level (e.g., c.145C>T means cytosine to thymine at position 145)

## File Formats
- **VCF (.vcf)**: Variant Call Format with header lines (##) and data columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE
- **VCF genotype field**: Format like `GT:DP:GQ` with values like `0/1:30:99` meaning heterozygous, 30 reads deep, genotype quality 99

## Databases & URLs
| Resource | URL | Used For |
|---|---|---|
| ClinVar | `https://www.ncbi.nlm.nih.gov/clinvar/` | Looking up clinical significance of variants |
| dbSNP | `https://www.ncbi.nlm.nih.gov/snp/` | Finding rs numbers, allele frequencies, population data |
| OMIM | `https://www.omim.org` | Gene-disease relationships and inheritance patterns |
| UCSC Genome Browser | `https://genome.ucsc.edu` | Visualizing variant positions in genomic context |
| GnomAD | `https://gnomad.broadinstitute.org` | Population allele frequencies across diverse groups |

## TAS2R38 Quick Reference
| SNP | rs Number | Taster (PAV) Allele | Non-Taster (AVI) Allele | Amino Acid Change |
|---|---|---|---|---|
| SNP 1 | rs713598 | C (Pro) | G (Ala) | P49A |
| SNP 2 | rs1726866 | T (Val) | C (Ala) | A262V |
| SNP 3 | rs10246939 | C (Val) | T (Ile) | V296I |

- **PAV/PAV**: Strong taster (homozygous taster haplotype)
- **PAV/AVI**: Taster (heterozygous)
- **AVI/AVI**: Non-taster (homozygous non-taster haplotype)

## Common Pitfalls
- Confusing VCF positions (1-based) with BED positions (0-based) -- off-by-one errors are very common
- Forgetting that VCF header lines start with `#` and must be skipped when counting variants
- Misreading the genotype field: `0/1` is heterozygous, NOT "has one copy of the reference and one alternate"
- Assuming all three TAS2R38 SNPs segregate independently -- they are inherited as a haplotype block
- Not checking allele frequency before assuming a variant is rare or disease-causing
- Confusing gene name searches (TAS2R38) with variant searches (rs713598) in databases

## Quick Check
1. What are the three SNPs that define the TAS2R38 haplotype, and what are their rs numbers?
2. How do you distinguish a heterozygous genotype from a homozygous alternate in VCF format?
3. What does it mean when ClinVar classifies a variant as "Variant of Uncertain Significance"?
4. If someone has the genotype PAV/AVI for TAS2R38, can they taste PTC? Why?
5. Why is it important to check population allele frequency when interpreting a variant?
