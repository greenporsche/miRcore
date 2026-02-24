# SEQ Day 1 Study Guide: Sequencing Fundamentals

## Learning Objectives
By the end of this day, you should be able to:
- Explain what whole exome sequencing (WES) is and why it is used instead of whole genome sequencing for this camp
- Describe the structure of a FASTQ file and interpret Phred quality scores
- Understand why adapter sequences appear in sequencing data and must be removed
- Navigate the Great Lakes HPC cluster via SSH and basic Linux commands
- Edit and submit a SLURM batch job script to run cutadapt for adapter trimming
- Check job status using squeue, email notifications, and log files

## Key Vocabulary
| Term | Definition |
|---|---|
| Whole Exome Sequencing (WES) | A targeted sequencing approach that captures only the protein-coding exons (~1-2% of the genome, ~30 Mb), providing deep coverage of medically relevant regions at lower cost than whole genome sequencing |
| FASTQ | A text-based file format that stores both nucleotide sequences and their per-base quality scores; each read occupies exactly four lines |
| Phred quality score | A logarithmic measure of base-calling accuracy: Q = -10 log10(P_error). Q20 = 1% error, Q30 = 0.1% error, Q40 = 0.01% error |
| Adapter sequence | Short synthetic oligonucleotides ligated to DNA fragments during library preparation; they are sequenced when the read length exceeds the insert size |
| Cutadapt | A bioinformatics tool that identifies and removes adapter sequences and low-quality bases from sequencing reads |
| Great Lakes | The University of Michigan's high-performance computing (HPC) cluster used for running computationally intensive bioinformatics analyses |
| SLURM | Simple Linux Utility for Resource Management; the job scheduling system on Great Lakes that allocates computing resources |
| sbatch | The SLURM command to submit a batch job script to the scheduler queue |
| Paired-end reads | Sequencing both ends of a DNA fragment, producing two files (R1/forward and R2/reverse) that together provide more accurate alignment |
| Next-Generation Sequencing (NGS) | High-throughput, massively parallel sequencing technologies (e.g., Illumina) that generate millions of reads simultaneously |
| Sanger sequencing | First-generation sequencing method that reads one DNA fragment at a time; considered the gold standard for variant validation |

## Concept Summaries

### The SEQ Camp Context
The SEQ (Sequencing Your Genome) camp is the culminating experience in the miRcore summer program. Students who have completed both the Computational Biology (CB) and Bench to Sequencer (BTS) camps are now ready to work with their own exome sequencing data. This is a unique educational experience: each student has had their saliva collected, DNA extracted, and whole exome sequenced by a clinical laboratory. Today you receive your raw data and begin processing it.

This camp is restricted to 11th and 12th graders because of the maturity needed to handle personal genetic information responsibly. You will learn not only the bioinformatics pipeline but also the ethical dimensions of personal genomics -- what it means to look at your own DNA, what the limitations are, and why you should never make medical decisions based on research-grade data.

### From Saliva to Data: The Sequencing Workflow
Before you ever touch a computer, a complex wet-lab process has already occurred. Your saliva sample was collected, DNA was extracted, and that DNA was fragmented into small pieces (typically 200-400 base pairs). Adapters were ligated to both ends of each fragment, creating a "library" that the sequencer can read. The Illumina sequencer then performs sequencing-by-synthesis: fluorescently labeled nucleotides are incorporated one at a time, and a camera records which base was added at each cycle. This produces millions of short reads (typically 150 bp each), stored as FASTQ files.

Each FASTQ read is four lines: the header (starting with @), the nucleotide sequence, a separator (+), and the quality string. The quality string encodes a Phred score for each base using ASCII characters. Higher scores mean higher confidence. Because the sequencer sometimes reads into the adapter at the end of short inserts, we must trim these adapter sequences before aligning reads to the reference genome. This is what cutadapt does.

### Working on Great Lakes HPC
Genomic data processing requires far more computing power than a laptop can provide. A single exome generates 5-10 GB of raw FASTQ data, and processing it involves reading and writing hundreds of millions of records. The Great Lakes cluster provides this computing power. You connect via VPN and SSH, navigate using Linux commands (cd, ls, pwd, cat), edit files with nano, and submit jobs using SLURM. Each job script specifies the resources needed (CPU, memory, time) and the commands to run. The scheduler manages the queue so that everyone's jobs get fair access to the cluster's resources.

Key practical skills: always replace [UID] placeholders with your actual username (including removing the brackets), make sure cutadapt commands stay on a single line in nano, and check job results via email notifications and log files.

## How It Connects
- **Previous Camps (CB + BTS)**: In CB, you learned about DNA, genes, mutations, and basic computational tools. In BTS, you performed hands-on laboratory work including DNA extraction, PCR, and gel electrophoresis. SEQ builds on both: you now analyze the digital output of the sequencing process.
- **Next Day (Day 2)**: Tomorrow you will take your trimmed FASTQ files and align them to the human reference genome using BWA, creating BAM files. Today's adapter trimming is a prerequisite -- poor trimming leads to poor alignment.

## Review Questions
1. Why does miRcore use whole exome sequencing rather than whole genome sequencing for this camp? List at least two practical reasons.
2. A FASTQ quality character `I` has ASCII value 73. Using Phred+33 encoding, what is the quality score? What is the probability of a base-calling error?
3. Describe the four-line structure of a FASTQ record. What information does each line contain?
4. What would happen if you forgot to trim adapters before aligning reads to the reference genome?
5. You submitted a job on Great Lakes and received an email saying "FAILED." Describe the steps you would take to troubleshoot.
6. Why must the cutadapt command in your job script be on a single line (or use proper line continuation)?
7. Explain why paired-end sequencing is preferred over single-end for exome sequencing.
8. What is the difference between submitting a job with `sbatch` and running a command directly on the login node?

## Further Reading
- Cutadapt documentation: https://cutadapt.readthedocs.io
- Illumina Sequencing Technology overview: https://www.illumina.com/science/technology/next-generation-sequencing.html
- FASTQ format specification: https://en.wikipedia.org/wiki/FASTQ_format
- Great Lakes HPC User Guide: https://arc.umich.edu/greatlakes/
- Phred quality scores explained: https://en.wikipedia.org/wiki/Phred_quality_score
