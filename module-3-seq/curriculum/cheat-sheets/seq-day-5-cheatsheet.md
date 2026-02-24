# SEQ Day 5 Cheat Sheet: Final Presentations & Ethics

## Commands & Tools
| Tool/Command | What It Does | Example Usage |
|---|---|---|
| Google Slides | Create and edit group presentation slides | Collaborative editing with breakout room team |
| IGV | Capture screenshots of personal variant alignments for slides | Navigate to variant, take screenshot |
| ClinVar / dbSNP | Final variant lookups for presentation evidence | Confirm pathogenicity and allele frequency |
| UCSC Genome Browser | Generate gene region images for slides | Screenshot of gene structure at locus of interest |

## Key Concepts
- **Genetic discrimination**: Using genetic information to discriminate in employment, insurance, or other areas
- **GINA (Genetic Information Nondiscrimination Act)**: U.S. federal law (2008) prohibiting genetic discrimination in health insurance and employment
- **GINA limitations**: Does not cover life insurance, disability insurance, long-term care insurance, or military
- **Direct-to-consumer (DTC) genetic testing**: Companies like 23andMe that provide genetic results directly to individuals without a physician intermediary
- **Clinical-grade sequencing**: Sequencing performed in a CLIA-certified lab with validated protocols; required for medical decisions
- **Research-grade sequencing**: Sequencing for educational/research purposes (like camp data); NOT suitable for clinical decisions
- **Informed consent**: The process of understanding what genetic testing involves, including risks of learning unexpected information
- **Incidental findings**: Unexpected discoveries in genomic data unrelated to the original purpose of testing
- **Penetrance (revisited)**: Not everyone with a pathogenic variant will develop the associated condition; incomplete penetrance is common
- **Genetic counseling**: Professional guidance for interpreting genetic test results and making informed health decisions

## Presentation Structure
| Slide | Content |
|---|---|
| Title Slide | Group name, member names, date |
| Pipeline Overview | Steps from raw FASTQ to variant calls (visual flowchart) |
| Gene/Trait of Interest | Background on chosen gene (TAS2R38, ABO, ABCC11, MMP1, or other) |
| Personal Findings | Summary of group members' variants and predicted phenotypes (only what each member consents to share) |
| Database Evidence | ClinVar classification, allele frequencies, population data |
| IGV Visualization | Screenshot showing read alignments at the variant position |
| Ethical Considerations | Discussion of privacy, GINA, limitations of research-grade data |
| Conclusion | What was learned, reflections on personal genomics |

## Ethics Discussion Points
| Topic | Key Consideration |
|---|---|
| Privacy of genetic data | Who should have access? Can you "un-know" results? |
| Consent for data sharing | Each student chooses what to share publicly in presentations |
| Insurance implications | GINA covers health insurance but not life/disability insurance |
| Family implications | Your genetic data reveals information about relatives who did not consent |
| Population-specific databases | Variant databases are biased toward European-descent populations |
| Actionability | Difference between interesting and medically actionable findings |

## Databases & URLs
| Resource | URL | Used For |
|---|---|---|
| GINA full text | `https://www.genome.gov/about-genomics/policy-issues/Genetic-Discrimination` | Understanding genetic nondiscrimination protections |
| ClinGen | `https://clinicalgenome.org` | Curated gene-disease validity and clinical actionability |
| ACMG Secondary Findings | `https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/` | List of genes recommended for return of incidental findings |
| NHGRI Genomics & Medicine | `https://www.genome.gov` | General genomics education and policy resources |

## Common Pitfalls
- Sharing personal genetic results in the presentation without explicit consent from the individual
- Making definitive health claims based on research-grade exome data (always caveat with "for educational purposes")
- Forgetting to cite databases and sources when presenting variant evidence
- Not explaining technical terms (VCF, haplotype, etc.) for a general audience (parents attend final presentations)
- Spending too much time on pipeline technical details and not enough on the biological story
- Overlooking the ethical dimensions -- presentations should include reflection on responsible use of genetic data

## Quick Check
1. What does GINA protect against, and what are its key limitations?
2. Why is it important to get consent before sharing someone else's genetic data in a presentation?
3. What is the difference between research-grade and clinical-grade sequencing?
4. Name two ethical concerns that arise when sequencing your own genome.
5. If you found a "pathogenic" variant in your exome data, what should your next step be (and what should it NOT be)?
