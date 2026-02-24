# SEQ Day 5 Exercises: Final Presentations & Ethics

## Warm-Up Questions
1. What is GINA, and why was it enacted? What does it protect, and what does it NOT protect?
2. Explain the difference between research-grade and clinical-grade genetic testing. Why does this distinction matter?
3. You have analyzed four genes this week (TAS2R38, ABO, ABCC11, MMP1). Which of these has the most direct medical relevance, and why?

## Hands-On Exercises

### Exercise 1: Presentation Preparation - Pipeline Slide
**Objective**: Create a clear, audience-friendly visual of the sequencing analysis pipeline.
**Instructions**:
1. In your group's Google Slides presentation, create a "Pipeline Overview" slide.
2. The slide should show the complete pipeline as a flowchart:
   ```
   Saliva Sample --> DNA Extraction --> Library Prep --> Illumina Sequencing
        --> Raw FASTQ --> Cutadapt (adapter trim) --> BWA (alignment)
        --> SAM --> BAM --> Sorted BAM --> GATK (variant calling) --> VCF
   ```
3. For EACH step, add a one-sentence annotation that a parent or non-scientist could understand.
4. Use color coding: green for wet lab steps, blue for computational steps.
5. Include approximate file sizes at key stages (FASTQ: ~5-10 GB, BAM: ~2-5 GB, VCF: ~5-10 MB).

**Expected Output**: A visually clear slide that communicates the pipeline to a non-technical audience.

### Exercise 2: Variant Highlight Slide
**Objective**: Present your most interesting personal variant finding with supporting evidence.
**Instructions**:
1. Choose ONE gene/variant from your analysis that you want to highlight (you must have consent from any group members whose data is shown).
2. Create a slide that includes:
   - Gene name and what it does (one sentence, plain language)
   - The variant(s) found and what they predict about the phenotype
   - An IGV screenshot showing the variant in your alignment data
   - Database evidence (ClinVar classification, allele frequency from dbSNP)
   - Whether the predicted phenotype matches reality
3. Add speaker notes explaining what you will say for this slide.

**Expected Output**: A polished slide with visual evidence (IGV screenshot), database references, and a clear genotype-to-phenotype story.

### Exercise 3: Ethics Discussion Preparation
**Objective**: Prepare thoughtful talking points about the ethical dimensions of personal genomics.
**Instructions**:
Your group presentation must include a section on ethics. Prepare responses to these prompts:

1. **Privacy**: If you found a variant associated with a serious disease in your own data, would you want to know? Would you tell your family? Why or why not?

2. **Consent**: Your genetic data contains information about your family members (parents, siblings) who did not consent to testing. How should this be handled?

3. **Data Security**: Your exome data is stored on a university server. Who should have access? What would happen if it were breached?

4. **Insurance**: Under GINA, health insurers cannot use genetic information against you. But life insurance companies can. Is this fair? Should GINA be expanded?

5. **Equity**: Genetic databases like ClinVar and gnomAD have historically over-represented people of European descent. How might this affect variant interpretation for people of other backgrounds?

Each group member should prepare to speak on at least one of these topics for 30-60 seconds.

**Expected Output**: Written bullet points for each ethical topic, with at least one assigned to each group member.

### Exercise 4: Peer Review Practice
**Objective**: Practice reviewing and giving constructive feedback on scientific presentations.
**Instructions**:
Before the final presentation, practice with your breakout room group:
1. Each group member presents their assigned slide(s) -- time yourself (aim for 1-2 minutes per slide).
2. After each person presents, other group members provide feedback using this rubric:

| Category | Score (1-5) | Feedback |
|---|---|---|
| Clarity: Could a non-scientist follow the explanation? | | |
| Accuracy: Are the scientific details correct? | | |
| Visuals: Are the slides clean, readable, with helpful images? | | |
| Engagement: Did the presenter make eye contact, speak clearly? | | |
| Ethics: Were responsible caveats included? | | |

3. Revise your slides and delivery based on feedback.
4. Do one final run-through as a group. Time the full presentation (target: 8-12 minutes total).

**Expected Output**: Feedback rubrics for each member; a timed, polished group presentation ready for the final session.

### Exercise 5: Reflection Journal
**Objective**: Reflect on the entire SEQ camp experience and what you have learned.
**Instructions**:
Write a short reflection (1-2 paragraphs) addressing the following:
1. What was the most surprising thing you learned about your own genome this week?
2. How has your understanding of genetics changed from CB --> BTS --> SEQ?
3. If you could sequence any additional part of your genome or analyze a different gene, what would it be and why?
4. How will you explain what you did this week to a friend or family member who knows nothing about genomics?
5. What is one ethical concern about personal genomics that you think society needs to address?

**Expected Output**: A thoughtful written reflection demonstrating growth across the three-camp progression.

## Challenge Problems

### Challenge 1: Designing a Genetic Testing Policy
**Objective**: Apply ethical reasoning to a real-world policy question.
**Instructions**:
You have been asked to advise a school board on whether to offer optional genetic testing to high school students as part of a biology curriculum. Write a one-page policy recommendation that addresses:
1. What genes/traits should be tested (and which should NOT)?
2. How should informed consent be handled for minors?
3. Who should have access to the results?
4. How should students with concerning findings be supported?
5. What safeguards are needed to prevent discrimination or stigmatization?
6. Should results be shared with parents? What if the student does not want to share?

### Challenge 2: The Future of Personal Genomics
**Objective**: Think critically about where personal genomics is heading.
**Instructions**:
Research and write a brief essay (300-500 words) on ONE of these topics:
1. **Pharmacogenomics**: How genetic variants affect drug metabolism. Find an example of a gene (e.g., CYP2D6, CYP2C19) where knowing your genotype could change which medication a doctor prescribes.
2. **Polygenic risk scores**: Unlike single-gene traits (TAS2R38), most diseases are influenced by hundreds of variants. How are polygenic risk scores calculated, and what are their current limitations?
3. **Gene editing and CRISPR**: If you could edit one variant in your genome, would you? Should this technology be used on embryos? Where do you draw the line?
4. **Genetic privacy in the age of forensic genealogy**: How was the Golden State Killer identified using genealogy databases? What does this mean for genetic privacy?

---

## Answer Key

### Warm-Up Answers
1. GINA (Genetic Information Nondiscrimination Act, 2008) prohibits genetic discrimination in health insurance and employment. It was enacted because advances in genetic testing raised fears that people would avoid testing if results could be used against them. GINA does NOT cover life insurance, disability insurance, long-term care insurance, the military, or companies with fewer than 15 employees.
2. Research-grade sequencing (like camp data) is performed for educational/research purposes without the rigorous quality controls required for medical decisions. Clinical-grade testing is done in CLIA-certified laboratories with validated protocols, quality assurance, and a chain of custody. Only clinical-grade results should be used for medical decisions. This matters because a variant call from research-grade data could be a false positive or false negative.
3. ABO blood type has the most direct medical relevance because blood type is critical for blood transfusions. Receiving incompatible blood can cause a fatal hemolytic reaction. TAS2R38 (taste) and ABCC11 (earwax) are interesting but not medically critical. MMP1 has some relevance to dermatology and cancer risk but is not directly actionable.

### Exercise Answers
1. **Exercise 1**: Pipeline slides should be visually clear with a logical flow. Common mistakes: too much text, missing the wet-lab-to-computational transition, not explaining file formats in plain language.
2. **Exercise 2**: Strong slides include: gene name with plain-language function, clear genotype notation, IGV screenshot with the variant highlighted, ClinVar/dbSNP evidence, and a statement about whether the prediction matches reality.
3. **Exercise 3**: There are no single right answers for ethics discussions. The key is demonstrating that students have thought critically about the issues. Strong responses acknowledge complexity and avoid black-and-white thinking.
4. **Exercise 4**: Feedback should be specific and constructive. Common presentation issues: speaking too fast, using jargon without defining it, slides with too much text, not mentioning ethical caveats, and not being mindful of consent when showing others' data.
5. **Exercise 5**: Reflections should show progression from CB (learning about DNA and genes) through BTS (hands-on lab work) to SEQ (computational analysis of personal data). Strong reflections connect the technical to the personal and ethical.

### Challenge Answers
1. **Challenge 1**: Strong policies will include: limiting testing to benign/educational traits (like TAS2R38, not BRCA), requiring parental consent for minors with student assent, results accessible only to the student (and parents, if the student agrees), referral pathways for concerning findings, clear statements that results are not clinical grade, and anonymous aggregate data only for classroom discussion.
2. **Challenge 2**: Answers vary by topic. Key elements for each: (1) Pharmacogenomics -- CYP2D6 poor metabolizers may need different doses of codeine, antidepressants; FDA now recommends pharmacogenomic testing for some drugs. (2) Polygenic scores combine effects of hundreds/thousands of variants; current limitations include poor transferability across populations, modest predictive power, and ethical concerns about use in embryo selection. (3) CRISPR raises questions about somatic vs. germline editing, consent of future generations, and the distinction between treating disease and enhancement. (4) Forensic genealogy used DNA uploaded to GEDmatch by relatives of the Golden State Killer; raises concerns about consent, genetic surveillance, and the erosion of genetic privacy for everyone in a DNA database.
