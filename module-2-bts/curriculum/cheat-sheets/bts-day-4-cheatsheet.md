# BTS Day 4 Cheat Sheet: Genome Mapping and IGV

## Commands & Tools
| Tool/Command | What It Does | Example Usage |
|---|---|---|
| `bowtie2 -x -U -S` | Align single-end reads to a reference | `bowtie2 -x hg38_index -U patient.fastq -S patient.sam` |
| `samtools view -bS` | Convert SAM to BAM | `samtools view -bS patient.sam > patient.bam` |
| `samtools sort` | Sort BAM by coordinate | `samtools sort patient.bam -o patient.sorted.bam` |
| `samtools index` | Index a sorted BAM | `samtools index patient.sorted.bam` |
| `sbatch` | Submit a job script to the SLURM scheduler | `sbatch my_job.sh` |
| `squeue -u` | Check the status of your submitted jobs | `squeue -u uniqname` |
| `cp` | Copy a file (used to create job script from shell script) | `cp small-shell.sh small-job.sh` |
| `nano` | Edit files in the terminal (e.g., add SLURM headers) | `nano small-job.sh` |
| `scp` | Transfer files between Great Lakes and local machine | `scp user@greatlakes:~/file.bam .` |
| IGV (desktop app) | Visualize alignments against a reference genome | Load genome, then load sorted BAM + BAI |

## Key Concepts
- **IGV (Integrative Genomics Viewer)**: A desktop application for visually exploring genomic data; used to view BAM alignments against a reference genome
- **HG38**: The current human reference genome assembly (GRCh38); used to determine if patient reads map to human DNA
- **SLURM**: The job scheduler on the Great Lakes cluster; manages computational resource allocation
- **sbatch**: A SLURM command to submit batch job scripts for execution on compute nodes
- **Shell Script (.sh)**: A text file containing a series of Linux commands that run sequentially; used to automate the alignment pipeline
- **Job Script**: A shell script with SLURM directives (`#SBATCH`) at the top specifying resources (time, memory, CPUs, account)
- **Patient Samples**: The BTS project analyzes 8 patient FASTQ files to determine which contain SARS-CoV-2
- **Dual Mapping Strategy**: Map reads first to the SARS-CoV-2 genome, then to the human genome (HG38) to compare alignment rates
- **Alignment Rate**: The percentage of reads that successfully map to a reference genome; a high rate to SARS-CoV-2 suggests infection

### SLURM Job Script Header
```bash
#!/bin/bash
#SBATCH --job-name=bowtie2_align
#SBATCH --account=bioinf545w24_class
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=%x-%j.out

module load Bioinformatics
module load Bowtie2
module load SAMtools

bowtie2 -x /path/to/index -U /path/to/reads.fastq -S output.sam
samtools view -bS output.sam > output.bam
samtools sort output.bam -o output.sorted.bam
samtools index output.sorted.bam
```

## File Formats
- **Shell Script (.sh)**: Plain text file with Linux commands; first line is `#!/bin/bash`
- **SLURM Output (.out)**: Text file capturing stdout/stderr from a submitted job
- **BAM (.bam) + BAI (.bai)**: Sorted BAM and its index are needed together to load into IGV

## Databases & URLs
| Resource | URL | Used For |
|---|---|---|
| IGV Download | https://igv.org/ | Download IGV desktop application |
| UCSC Genome Browser | https://genome.ucsc.edu/ | Alternative genome browser |
| Great Lakes Documentation | https://arc.umich.edu/greatlakes/ | SLURM usage and cluster info |

## Common Pitfalls
- Forgetting to add SLURM `#SBATCH` headers when converting a shell script to a job script
- Submitting a job without loading the required modules inside the job script (modules are not inherited from your interactive session)
- Using the wrong account name in `#SBATCH --account` --- this must match your class or lab allocation
- Not downloading BOTH the sorted BAM and BAI files to your local machine before opening in IGV
- Loading the wrong reference genome in IGV (e.g., loading hg19 instead of hg38, or forgetting to load the SARS-CoV-2 genome)
- Forgetting that `squeue` shows only running/pending jobs --- a missing job means it already completed (check the `.out` file)
- Trying to run a compute-heavy alignment interactively on the login node instead of submitting via `sbatch`

## Quick Check
1. What is the difference between running a shell script interactively and submitting it with `sbatch`?
2. What two files must you download to visualize an alignment in IGV?
3. Why do we map patient reads to both the SARS-CoV-2 genome and the human genome?
4. What command checks the status of your submitted SLURM jobs?
5. In IGV, what does it mean when you see many reads stacked up at a particular genomic location?
