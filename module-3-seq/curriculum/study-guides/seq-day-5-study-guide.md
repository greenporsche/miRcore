# SEQ Day 5 Study Guide: Final Presentations & Ethics

## Learning Objectives
By the end of this day, you should be able to:
- Synthesize a week of bioinformatics work into a clear, audience-appropriate scientific presentation
- Communicate genomic concepts (pipeline, variants, genotype-phenotype relationships) to a non-technical audience
- Articulate the ethical dimensions of personal genomics including privacy, consent, and genetic discrimination
- Explain what GINA (Genetic Information Nondiscrimination Act) protects and its limitations
- Distinguish between research-grade and clinical-grade genetic testing results
- Reflect on the progression from CB through BTS to SEQ and your growth as a genomics student

## Key Vocabulary
| Term | Definition |
|---|---|
| GINA (Genetic Information Nondiscrimination Act) | U.S. federal law enacted in 2008 that prohibits discrimination based on genetic information in health insurance coverage and employment decisions |
| Genetic discrimination | The unfair treatment of individuals based on their genetic information, including decisions about insurance, employment, education, or social relationships |
| Informed consent | The process by which a participant is educated about the purpose, procedures, risks, and benefits of genetic testing before agreeing to participate |
| Incidental finding | An unexpected genetic result unrelated to the original purpose of testing, such as discovering a cancer predisposition variant while investigating bitter taste |
| Clinical-grade testing | Genetic testing performed in a CLIA-certified laboratory with validated protocols, quality assurance, proficiency testing, and regulatory oversight; results are suitable for medical decision-making |
| Research-grade testing | Genetic analysis performed for educational or research purposes without the rigorous quality controls required for clinical use; NOT suitable for medical decisions |
| Genetic counselor | A healthcare professional with specialized training in medical genetics and counseling who helps individuals and families understand and adapt to genetic test results |
| Direct-to-consumer (DTC) testing | Genetic tests marketed directly to consumers (e.g., 23andMe, AncestryDNA) without requiring a physician order, often with limited clinical validation |
| Penetrance | The probability that a person carrying a disease-associated variant will actually develop the disease; many pathogenic variants have incomplete penetrance |
| Actionability | Whether knowing about a genetic variant can lead to a specific medical intervention that improves health outcomes |

## Concept Summaries

### Science Communication: Presenting Your Findings
Today you synthesize everything from this week into a group presentation delivered to parents, peers, and camp instructors. Effective science communication requires translating technical work into accessible language without sacrificing accuracy. Your audience includes people who may have no background in genetics or computing.

A strong presentation tells a story with a clear arc: We started with saliva (something familiar), processed it through a pipeline (show the flowchart), found specific variants in our DNA (share the gene stories), and reflected on what it means (ethics). Use visuals: IGV screenshots, database entries, and the pipeline flowchart communicate more effectively than text-heavy slides. Define technical terms the first time you use them. Every group member should present at least one section, and practice at least once before the final session. Remember that the parents in the audience are hearing about what their child's DNA revealed -- be sensitive and thoughtful in how you frame findings.

### The Ethics of Personal Genomics
The SEQ camp places ethics at the center, not as an afterthought. When you analyze your own exome, you confront questions that society is still grappling with.

**Privacy and consent**: Your genetic data is uniquely identifying -- more so than a fingerprint. It also contains information about your biological relatives who did not consent to testing. When you find that you carry a specific variant, that information has implications for your parents and siblings. Deciding what to share, with whom, and in what context requires careful thought.

**GINA and its gaps**: The Genetic Information Nondiscrimination Act protects Americans from discrimination in health insurance and employment based on genetic information. However, GINA does NOT cover life insurance, disability insurance, long-term care insurance, or the military. This means that a person who discovers they carry a pathogenic variant could face higher life insurance premiums or denial of coverage. This gap in protection is why many experts recommend genetic counseling before any testing, and why the camp emphasizes that research-grade results should not be treated as clinical diagnoses.

**Research vs. clinical grade**: Perhaps the most important ethical message of the week is this: the exome data you analyzed is research-grade. It was not processed in a CLIA-certified clinical laboratory. The variant calls have not been independently validated. If you found something that concerns you, the appropriate response is to speak with a genetic counselor and pursue confirmatory clinical testing -- NOT to make medical decisions based on camp data. This distinction between research-grade and clinical-grade data is fundamental to responsible genomics.

### Reflection: The CB to BTS to SEQ Journey
The three-camp arc of the miRcore program represents a deliberate pedagogical progression. In Computational Biology (CB), you learned the fundamentals: DNA structure, the genetic code, how mutations change proteins, and basic computational tools. In Bench to Sequencer (BTS), you performed hands-on laboratory work: extracting DNA, running PCR, analyzing results on gels, and understanding how the wet lab generates data. In SEQ, you closed the loop: taking real sequencing data from your own genome through a complete bioinformatics pipeline to discover personal genetic variants.

This progression mirrors the real path of a geneticist or bioinformatician: you need to understand the biology (CB), the experimental methods (BTS), and the computational analysis (SEQ) to do meaningful work. The fact that you analyzed your own DNA makes this experience uniquely personal and motivating. You are among the relatively few high school students in the world who have gone from saliva to variant-level genomic analysis of their own data. That is genuinely pioneering.

## How It Connects
- **Previous Day (Day 4)**: Yesterday you completed your personal variant exploration across four genes and compiled a summary table. That table becomes the core content of today's presentation. The IGV screenshots you captured are your visual evidence.
- **The Full Week**: Day 1 (data processing fundamentals) --> Day 2 (alignment pipeline) --> Day 3 (variant identification) --> Day 4 (personal exploration) --> Day 5 (synthesis and communication). Each day built on the previous one, and today you demonstrate mastery by communicating the entire journey to others.
- **Beyond Camp**: The skills you learned this week -- Linux command line, HPC job submission, bioinformatics pipeline construction, database navigation, variant interpretation, IGV visualization, and ethical reasoning -- are directly applicable to careers in genomics, bioinformatics, genetic counseling, medicine, and bioethics.

## Review Questions
1. Why is it important to distinguish between research-grade and clinical-grade genetic testing? What could go wrong if someone made medical decisions based on research-grade data?
2. Explain GINA in your own words. What does it protect, and what are its key limitations?
3. Your friend says, "I found a pathogenic BRCA1 variant in my camp data. I need to get tested for breast cancer immediately." How would you respond?
4. You are presenting to parents who know nothing about genomics. How would you explain what a VCF file is in one sentence?
5. Why might a student choose NOT to share certain genetic findings in the group presentation? Is this decision ethically appropriate?
6. Explain the concept of incidental findings. Why are they a particular concern in exome and genome sequencing?
7. How does the miRcore three-camp progression (CB --> BTS --> SEQ) prepare you for understanding genomics more deeply than any single camp could?
8. If you could change one thing about genetic privacy laws in the United States, what would it be and why?

## Further Reading
- GINA and genetic discrimination: https://www.genome.gov/about-genomics/policy-issues/Genetic-Discrimination
- NHGRI: Genomics and Medicine: https://www.genome.gov/health/Genomics-and-Medicine
- Genetic counseling resources: https://www.nsgc.org
- ACMG recommendations on incidental findings: Green et al., 2013, Genetics in Medicine (DOI: 10.1038/gim.2013.73)
- Responsible use of genomic data (Presidential Commission for Bioethics): https://bioethicsarchive.georgetown.edu/pcsbi/
- Direct-to-consumer genetic testing: FDA overview at https://www.fda.gov/medical-devices/in-vitro-diagnostics/direct-consumer-tests
- The Immortal Life of Henrietta Lacks by Rebecca Skloot (book on consent and ethics in biomedical research)
