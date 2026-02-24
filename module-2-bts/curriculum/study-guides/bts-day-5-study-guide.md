# BTS Day 5 Study Guide: Research Projects and Pitches

## Learning Objectives
By the end of this day, you should be able to:
- Compile and interpret alignment results for all 8 patient samples
- Classify patient samples as SARS-CoV-2 positive or negative based on alignment rate data
- Use IGV to identify mutations in SARS-CoV-2 positive samples
- Explain the significance of identified mutations in the context of viral biology
- Write a concise research summary covering background, methods, results, and conclusions
- Deliver a polished biotech business pitch to a panel of mock investors
- Critically evaluate business pitches from other groups as a mock investor
- Describe the complete bioinformatics pipeline from raw FASTQ data to biological conclusions

## Key Vocabulary
| Term | Definition |
|---|---|
| Variant | A difference between a sample's sequence and the reference genome (e.g., a single nucleotide change) |
| SNV (Single Nucleotide Variant) | A change in a single base at a specific position relative to the reference genome |
| Mutation | A heritable change in the nucleotide sequence of an organism's genome |
| Coverage Depth | The number of reads overlapping a specific genomic position; higher depth increases confidence in variant calls |
| Positive/Negative | Classification of a patient sample as containing (positive) or not containing (negative) SARS-CoV-2 |
| Alignment Rate | The percentage of total reads that successfully mapped to a given reference genome |
| Variant of Concern | A SARS-CoV-2 lineage with mutations that may increase transmissibility, immune evasion, or disease severity |
| Target Market | The specific group of customers a biotech product or service is designed for |
| Investment Ask | The specific amount of funding requested in a business pitch and how it will be used |
| Research Summary | A structured document describing the background, methods, results, and conclusions of a study |

## Concept Summaries

### Compiling and Interpreting Patient Results
Day 5 is when all the pieces of the BTS project come together. Students have aligned reads from 8 patient samples to both the SARS-CoV-2 genome and the human genome (HG38). The critical analytical step is comparing these alignment rates to determine which patients are infected.

The interpretation framework is straightforward: if a patient's reads align at a high rate to SARS-CoV-2 (e.g., >50-80%), the sample contains significant amounts of viral RNA, indicating infection. If the reads align primarily to HG38 with minimal SARS-CoV-2 alignment (<5%), the patient is likely uninfected and the low viral alignment is noise or contamination. Ambiguous cases (moderate alignment to both) require deeper investigation --- examining mapping quality scores, coverage uniformity, and visual inspection in IGV.

Students create a comprehensive results table listing each patient's alignment rates to both genomes and their diagnosis (positive or negative). This table is the centerpiece of the research summary and parent presentation.

### Mutation Identification with IGV
Beyond simply detecting the presence of SARS-CoV-2, students use IGV to look for specific mutations in positive samples. When viewing aligned reads in IGV, bases that differ from the reference genome appear as colored letters (instead of gray). These mismatches can represent real biological mutations (variants) or sequencing errors.

To distinguish real mutations from errors, students consider several factors: Is the variant present in a high percentage of reads at that position (e.g., >80%)? Is the position covered by many reads (high depth)? Do the reads supporting the variant have high mapping quality (MAPQ) scores? A true mutation will be consistently present across many high-quality reads, while a sequencing error will appear sporadically in only a few reads.

Of particular interest are mutations in the spike protein region (positions 21,563-25,384), because these can affect viral transmissibility and vaccine effectiveness. Students may identify mutations that correspond to known variants of concern, such as the D614G mutation that became dominant early in the pandemic. Connecting bioinformatics observations to real-world epidemiological significance is a powerful learning moment.

### Writing a Research Summary
The research summary consolidates the week's work into a structured scientific document. It follows a standard format: Background (why this analysis matters --- the COVID-19 pandemic and the need for rapid pathogen identification), Methods (the bioinformatics pipeline: FASTQ reads aligned with Bowtie2, processed with SAMtools, visualized in IGV), Results (which patients are positive, what mutations were found), and Conclusions (what the results mean in context).

This exercise introduces students to scientific writing, emphasizing clarity, precision, and evidence-based reasoning. The summary should reference specific numbers (alignment rates, mutation positions) rather than vague statements. It should also acknowledge limitations, such as the small sample size and the simplified analysis pipeline compared to clinical diagnostic workflows.

### Final Business Pitches and Investment Simulation
The business pitch presentations on Day 5 are the culmination of the entrepreneurship track. Groups deliver polished 5-7 minute presentations to the entire camp, with instructors and peers acting as an investment panel. The pitch structure typically includes: the problem (an unmet medical or scientific need), the solution (a biotech product or technology), the target market (who will use it and how large is the market), the competitive landscape (existing alternatives and the proposed advantage), the team (group members and their roles), and the ask (specific funding amount and allocation plan).

After all pitches are delivered, students participate in an investment simulation. Each student receives a fixed amount of play money and must decide which group(s) to invest in. This exercise teaches critical evaluation --- students must assess not just the presentation quality but the underlying business viability: Is the market real? Is the technology feasible? Does the team inspire confidence? The simulation creates an engaging, competitive atmosphere that reinforces the connection between scientific innovation and real-world commercialization.

### The Complete BTS Pipeline: A Week in Review
The full bioinformatics pipeline learned across all five days is:

1. **Data acquisition**: Patient samples are sequenced using NGS (Illumina), producing FASTQ files with millions of short reads
2. **Quality assessment**: Examining FASTQ files to understand read counts, quality scores, and GC content
3. **Alignment**: Using Bowtie2 to map reads against reference genomes (SARS-CoV-2 and HG38)
4. **File conversion**: SAM to BAM conversion using `samtools view -bS`
5. **Sorting**: Organizing BAM records by genomic position using `samtools sort`
6. **Indexing**: Creating a BAI index for fast random access using `samtools index`
7. **Visualization**: Loading sorted BAM + BAI into IGV for visual exploration
8. **Analysis**: Comparing alignment rates, identifying mutations, and classifying samples
9. **Interpretation**: Connecting computational results to biological and clinical conclusions

This pipeline represents a simplified but authentic version of what bioinformaticians do in research labs, clinical diagnostics facilities, and public health agencies around the world.

## How It Connects
- **Previous Day**: Day 4 introduced IGV visualization, SLURM job submission, and the dual mapping strategy. Day 5 brings all patient analyses to completion and adds the interpretation layer --- moving from raw results to biological conclusions.
- **Looking Ahead**: Students who continue in bioinformatics research will encounter more sophisticated pipelines involving variant calling software (e.g., GATK, bcftools), phylogenetic analysis, machine learning for genomic prediction, and multi-omics integration. The fundamentals learned in this camp --- Linux, alignment, SAM/BAM, visualization --- form the foundation for all of these advanced topics.

## Review Questions
1. Describe how you would determine whether a patient sample is positive for SARS-CoV-2 using only alignment rate data.
2. In IGV, how can you distinguish a real mutation from a sequencing error?
3. Why are mutations in the spike protein region of particular interest?
4. What are the four sections of a standard research summary?
5. List the complete bioinformatics pipeline steps from FASTQ to IGV visualization.
6. What makes a biotech business pitch persuasive? Name at least four components.
7. If all 8 patient samples showed 0% alignment to SARS-CoV-2, what would you conclude? What if all showed 100%?
8. How does the dual mapping approach (SARS-CoV-2 + HG38) provide stronger evidence than mapping to just one reference?

## Further Reading
- IGV Variant Visualization Guide: https://igv.org/doc/desktop/
- SARS-CoV-2 Variants of Concern (WHO): https://www.who.int/activities/tracking-SARS-CoV-2-variants
- Writing a Scientific Report: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3178846/
- GISAID SARS-CoV-2 Database: https://www.gisaid.org/
- Biotech Venture Capital Overview: https://www.nature.com/articles/nbt0908-978
- miRcore Organization: https://www.mircore.org/
