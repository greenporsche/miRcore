# BTS Day 4 Study Guide: Genome Mapping and IGV

## Learning Objectives
By the end of this day, you should be able to:
- Map patient reads to both the SARS-CoV-2 genome and the human genome (HG38) and compare results
- Write shell scripts to automate the alignment pipeline
- Convert a shell script into a SLURM job script with appropriate `#SBATCH` directives
- Submit batch jobs to the Great Lakes cluster using `sbatch` and monitor them with `squeue`
- Transfer sorted BAM and BAI files to your local machine using `scp`
- Load and navigate alignment data in IGV (Integrative Genomics Viewer)
- Begin analyzing all 8 patient samples to identify SARS-CoV-2 positive cases

## Key Vocabulary
| Term | Definition |
|---|---|
| HG38 (GRCh38) | The current human reference genome assembly, approximately 3.2 billion base pairs |
| IGV | Integrative Genomics Viewer, a desktop application for visually exploring genomic data |
| SLURM | Simple Linux Utility for Resource Management, the job scheduler on the Great Lakes cluster |
| sbatch | A SLURM command to submit a batch job script for execution on a compute node |
| squeue | A SLURM command to check the status of submitted jobs (PENDING, RUNNING, COMPLETED) |
| Shell Script | A text file containing a sequence of Linux commands that are executed in order |
| Job Script | A shell script with SLURM `#SBATCH` directives specifying resource requirements |
| Compute Node | A dedicated server in the HPC cluster where submitted jobs run (versus the shared login node) |
| Login Node | The server you connect to via SSH; used for file management, not heavy computation |
| Alignment Rate | The percentage of reads that successfully map to a reference genome |
| Coverage | The number of reads that overlap a given position in the genome; higher coverage = more confidence |
| Dual Mapping | Aligning the same reads to two different genomes to determine the source organism |
| BAI | BAM index file, required alongside a sorted BAM for IGV visualization |

## Concept Summaries

### Shell Scripts and Automation
By Day 4, students have run the alignment pipeline (Bowtie2 -> samtools view -> samtools sort -> samtools index) multiple times manually. Shell scripting automates this process by writing all commands into a single `.sh` file that can be executed with one command. A shell script begins with `#!/bin/bash` (called a "shebang line"), which tells the operating system to use the Bash shell to interpret the commands.

The key advantage of shell scripts is reusability and reproducibility. Instead of typing four or five commands each time, you run the script once. When you need to process a different patient sample, you only need to change the filename in the script. More advanced scripts can use variables and loops to process all patient samples automatically, as demonstrated in the for-loop approach: `for i in 1 2 3 4 5 6 7 8; do ... done`.

Students first test their pipeline interactively (running the script on the login node with a small file) to verify the commands are correct. Once confirmed, the script is converted to a SLURM job script for submission to compute nodes with full resources.

### SLURM Job Submission
The Great Lakes cluster uses SLURM to manage computational resources. When you run commands interactively on the login node, you are sharing resources with all other logged-in users, which means compute-heavy tasks like genome alignment can be slow and may impact others. SLURM solves this by allocating dedicated compute nodes for your jobs.

A SLURM job script is a shell script with special `#SBATCH` directives at the top. These directives specify: `--job-name` (a label for your job), `--account` (the billing account, e.g., the class account), `--partition` (which pool of nodes to use), `--nodes`, `--ntasks`, `--cpus-per-task` (how many compute resources to request), `--mem` (RAM allocation), `--time` (maximum wall clock time), and `--output` (where to save the job's console output). The `%x` and `%j` placeholders in the output filename expand to the job name and job ID, respectively.

Important: the `#SBATCH` lines look like comments to Bash (because they start with `#`), but SLURM reads them as configuration directives. Also, modules loaded in your interactive session are NOT inherited by submitted jobs --- you must include `module load` commands inside the job script itself.

### Dual Mapping Strategy: SARS-CoV-2 vs. HG38
The core analytical approach of the BTS project is dual mapping. Patient samples contain a mixture of nucleic acids: human DNA/RNA from the patient's own cells, and potentially viral RNA if the patient is infected with SARS-CoV-2. By aligning the same reads to both the SARS-CoV-2 reference genome and the human reference genome (HG38), students can determine the composition of each sample.

If a patient is infected, their reads will show a high alignment rate when mapped to the SARS-CoV-2 genome (because the sample contains viral RNA). The same reads will show a relatively lower alignment rate to HG38 (because many of the reads are viral, not human). Conversely, an uninfected patient's reads will show very low alignment to SARS-CoV-2 (essentially noise) and high alignment to HG38 (because the sample is predominantly human DNA).

Note that mapping to HG38 requires significantly more computational resources than mapping to SARS-CoV-2. The human genome is approximately 3.2 billion base pairs, while SARS-CoV-2 is only about 30,000 base pairs --- a factor of more than 100,000. This is why the HG38 job script requests more memory (8G vs. 4G) and more time (1 hour vs. 30 minutes).

### IGV (Integrative Genomics Viewer)
IGV is a powerful desktop application developed by the Broad Institute that enables visual exploration of genomic data. After running the alignment pipeline and producing sorted BAM files with indices, students transfer these files to their local computers and load them into IGV.

To use IGV, you must first select the correct reference genome. For SARS-CoV-2 analysis, load the virus genome; for human analysis, load HG38. Then load the sorted BAM file (IGV automatically looks for the `.bai` index in the same directory). The IGV display shows a coverage track at the top (a histogram of how many reads overlap each position) and individual read tracks below (horizontal bars representing each aligned read).

By zooming in to specific regions, students can see individual bases in each read. Colored bases indicate mismatches from the reference genome --- these are potential mutations or variants. The spike protein region (21,563-25,384) is of particular interest because mutations here can affect how the virus infects cells and how well vaccines work. IGV makes these abstract data points visible and intuitive, transforming numbers in a BAM file into a biological story.

### Processing All 8 Patient Samples
The BTS camp project involves analyzing 8 patient FASTQ files. By Day 4, students have the skills to process all of them. The workflow for each patient is identical: align to SARS-CoV-2, align to HG38, compare alignment rates, and visualize in IGV. Efficient students create scripts that can be quickly modified by changing the patient number, or use a for-loop to process all samples in a single job.

The goal is to build a comprehensive results table showing each patient's alignment rates to both genomes, allowing a clear classification of positive versus negative SARS-CoV-2 cases. This table becomes the foundation for the research summary and parent presentation on Day 5.

## How It Connects
- **Previous Day**: Day 3 introduced the alignment pipeline (Bowtie2 + SAMtools) and the SAM/BAM formats. Day 4 applies these same tools at scale, adding automation (scripts), resource management (SLURM), and visualization (IGV).
- **Next Day**: Day 5 is the final day. Students complete their analysis of all 8 patients, identify mutations using IGV, prepare their research summaries, and deliver final business pitch presentations to a panel of mock investors.

## Review Questions
1. What is the purpose of the `#!/bin/bash` line at the top of a shell script?
2. Name three `#SBATCH` directives and explain what each controls.
3. Why must `module load` commands be included inside a SLURM job script?
4. If a patient's reads show 92% alignment to SARS-CoV-2 and 15% alignment to HG38, what does this suggest?
5. Why does mapping to HG38 require more memory and time than mapping to SARS-CoV-2?
6. What two files must be in the same directory for IGV to open an alignment?
7. In IGV, what do colored bases in the read tracks indicate?
8. How would you check whether your submitted SLURM job is still running or has completed?

## Further Reading
- SLURM Documentation: https://slurm.schedmd.com/documentation.html
- IGV User Guide: https://igv.org/doc/desktop/
- Bash Scripting Tutorial: https://www.gnu.org/software/bash/manual/bash.html
- Human Reference Genome (GRCh38): https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/
- Great Lakes SLURM Guide: https://arc.umich.edu/greatlakes/slurm-user-guide/
