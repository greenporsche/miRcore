# SEQ Day 1 Cheat Sheet: Sequencing Fundamentals

## Commands & Tools
| Tool/Command | What It Does | Example Usage |
|---|---|---|
| `ssh` | Connect to Great Lakes HPC cluster | `ssh uniqname@greatlakes.arc-ts.umich.edu` |
| `cd` | Change directory | `cd /scratch/miRcore_project/uniqname` |
| `pwd` | Print current working directory | `pwd` |
| `ls` | List files in a directory | `ls -la` |
| `nano` | Text editor on Great Lakes | `nano job1_cutadapt.sh` |
| `cutadapt` | Remove adapter sequences from FASTQ reads | `cutadapt -a ADAPTER_SEQ -o output.fq input.fq` |
| `sbatch` | Submit a job to the SLURM scheduler | `sbatch job1_cutadapt.sh` |
| `squeue` | Check status of submitted jobs | `squeue -u uniqname` |
| `scancel` | Cancel a running job | `scancel JOB_ID` |
| `module load` | Load software on Great Lakes | `module load Bioinformatics cutadapt` |
| `cat` | View contents of a file | `cat output.log` |

## Key Concepts
- **Whole Exome Sequencing (WES)**: Sequencing only the protein-coding regions (exons) of the genome, roughly 1-2% of total DNA
- **FASTQ format**: Text-based format storing both nucleotide sequences and per-base quality scores from the sequencer
- **Phred quality score**: Logarithmic score representing base-call accuracy; Q30 = 99.9% accuracy, encoded as ASCII characters
- **Adapter sequences**: Short synthetic DNA sequences ligated during library prep that must be trimmed before analysis
- **Next-Generation Sequencing (NGS)**: Massively parallel sequencing technology (e.g., Illumina) producing millions of short reads simultaneously
- **Sanger sequencing**: First-generation method reading one fragment at a time; still the gold standard for validation
- **Paired-end reads**: Sequencing both ends of a DNA fragment (R1 = forward, R2 = reverse), improving mapping accuracy
- **HPC (High-Performance Computing)**: Remote computing cluster (Great Lakes at U-M) needed for processing large genomic datasets
- **SLURM**: Job scheduler on Great Lakes that manages resource allocation for submitted compute jobs

## File Formats
- **FASTQ (.fq, .fastq, .fq.gz)**: Raw sequencing reads with quality scores; four lines per read (header, sequence, +, quality)
- **Shell script (.sh)**: Batch script containing SLURM directives and commands for HPC job submission
- **Log file (.log)**: Output file from a completed job containing runtime info, errors, and resource usage

## Databases & URLs
| Resource | URL | Used For |
|---|---|---|
| Great Lakes HPC | `https://greatlakes.arc-ts.umich.edu` | Accessing the compute cluster via web portal |
| U-M VPN | `https://its.umich.edu/enterprise/wifi-networks/vpn` | Required to connect to Great Lakes remotely |
| Cutadapt documentation | `https://cutadapt.readthedocs.io` | Reference for adapter trimming parameters |
| Zulip chat | (Camp-specific link) | Daily communication and poll/attendance tool |

## Common Pitfalls
- Forgetting to connect to VPN before accessing Great Lakes
- Breaking a single-line command across two lines in `nano`, causing the job script to fail
- Forgetting to replace `[UID]` or `<uniqname>` placeholders (including the brackets) with your actual username
- Not checking your U-M email for job start/end notifications to verify whether jobs succeeded or failed
- Confusing the SLURM log file location with the directory where you submitted the job
- Running commands directly on the login node instead of submitting via `sbatch`

## Quick Check
1. What are the four lines of a FASTQ record, and what does each line contain?
2. Why do we need to remove adapter sequences before downstream analysis?
3. What does a Phred quality score of Q20 mean in terms of error probability?
4. How do you check if your submitted SLURM job succeeded or failed?
5. What is the difference between whole genome sequencing and whole exome sequencing?
