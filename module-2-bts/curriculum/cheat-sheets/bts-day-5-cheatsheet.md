# BTS Day 5 Cheat Sheet: Research Projects and Pitches

## Commands & Tools
| Tool/Command | What It Does | Example Usage |
|---|---|---|
| `bowtie2` | Align reads (review from pipeline) | `bowtie2 -x sarscov2 -U patient8.fastq -S patient8.sam` |
| `samtools flagstat` | Generate alignment statistics from a BAM file | `samtools flagstat patient.sorted.bam` |
| `samtools view -c -F 4` | Count only mapped reads in a BAM file | `samtools view -c -F 4 patient.sorted.bam` |
| `scp` | Transfer result files to local machine for IGV | `scp user@greatlakes:~/results/*.sorted.bam .` |
| IGV | Visualize alignment results to identify variants | Open sorted BAM against reference genome |

## Key Concepts
- **Research Project Completion**: Students finalize analysis of all 8 patient samples, comparing SARS-CoV-2 vs. human alignment rates to identify infected patients
- **Alignment Rate Comparison**: If a patient's reads align at a high rate to the SARS-CoV-2 genome, the sample likely contains viral RNA
- **Mutation Identification**: Using IGV to visually identify single nucleotide variations (SNVs) between patient sequences and the SARS-CoV-2 reference genome
- **Business Pitch (Final)**: A polished 5-7 minute presentation to a panel of "investors" (instructors and peers) covering problem, solution, market, team, and financials
- **Target Market**: The specific group of customers a biotech product or service is designed for; must be clearly defined in the pitch
- **Investment Simulation**: Students role-play as investors with play money, deciding which group's biotech pitch to invest in
- **Variant**: A difference between a patient's sequence and the reference genome; can be a SNV (single base change), insertion, or deletion
- **Parent Presentation**: Students present their research findings and bioinformatics pipeline results to parents at camp conclusion

## File Formats
- **Sorted BAM + BAI**: Final output files from the pipeline, used for IGV visualization and variant identification
- **Pitch Slide Deck**: Presentation slides covering biotech company overview, product, market, and investment ask

## Databases & URLs
| Resource | URL | Used For |
|---|---|---|
| NCBI BLAST | https://blast.ncbi.nlm.nih.gov/ | Verifying sequence identity |
| IGV | https://igv.org/ | Visualizing patient alignment results |
| GISAID | https://www.gisaid.org/ | SARS-CoV-2 variant tracking (reference) |

## Common Pitfalls
- Drawing conclusions from alignment data without comparing both SARS-CoV-2 and human genome mapping rates
- Confusing a low-quality alignment artifact in IGV with an actual mutation
- Forgetting to check mapping quality scores when identifying variants --- low MAPQ reads are unreliable
- Not preparing a clear "ask" in the business pitch (what specific amount of investment, and what it will fund)
- Assuming all 8 patient samples are positive for SARS-CoV-2 --- the analysis should reveal which are positive and which are negative
- Spending too much time on bioinformatics and neglecting pitch preparation, or vice versa

## Quick Check
1. How do you determine from alignment data whether a patient sample contains SARS-CoV-2?
2. What visual indicator in IGV suggests a mutation relative to the reference genome?
3. Name three essential components of a biotech business pitch.
4. Why is it important to map patient reads to both a viral genome and the human genome?
5. What is the full bioinformatics pipeline you learned this week, from raw data to visualization?
