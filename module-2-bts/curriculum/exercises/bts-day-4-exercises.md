# BTS Day 4 Exercises: Genome Mapping and IGV

## Warm-Up Questions
1. What is the difference between running commands interactively on Great Lakes versus submitting them as a batch job with `sbatch`?
2. Why do we map patient reads to both the SARS-CoV-2 genome and the human genome (HG38)?
3. What is IGV, and why is it useful for analyzing alignment data?

## Hands-On Exercises

### Exercise 1: Creating a Shell Script
**Objective**: Write a shell script that automates the alignment pipeline for a single patient sample.
**Instructions**:
1. Open nano to create a new script: `nano patient_align.sh`
2. Write the following content:
   ```bash
   #!/bin/bash
   module load Bioinformatics
   module load Bowtie2
   module load SAMtools

   # Align to SARS-CoV-2
   bowtie2 -x /path/to/sarscov2_index -U patient1.fastq -S patient1_sarscov2.sam

   # Convert, sort, index
   samtools view -bS patient1_sarscov2.sam > patient1_sarscov2.bam
   samtools sort patient1_sarscov2.bam -o patient1_sarscov2.sorted.bam
   samtools index patient1_sarscov2.sorted.bam

   echo "Alignment complete for patient 1"
   ```
3. Save and exit nano (`Ctrl+O`, `Enter`, `Ctrl+X`)
4. Test the script interactively first: `bash patient_align.sh`
5. Verify the output files exist: `ls -lh patient1_sarscov2.*`

**Expected Output**: Four output files (SAM, BAM, sorted BAM, BAI). The terminal should print "Alignment complete for patient 1."

### Exercise 2: Converting to a SLURM Job Script
**Objective**: Add SLURM headers to your shell script so it can be submitted as a batch job.
**Instructions**:
1. Copy your shell script to a new job script: `cp patient_align.sh patient_align_job.sh`
2. Open the job script: `nano patient_align_job.sh`
3. Add SLURM directives right after the `#!/bin/bash` line:
   ```bash
   #SBATCH --job-name=patient1_align
   #SBATCH --account=bioinf545w24_class
   #SBATCH --partition=standard
   #SBATCH --nodes=1
   #SBATCH --ntasks=1
   #SBATCH --cpus-per-task=1
   #SBATCH --mem=4G
   #SBATCH --time=00:30:00
   #SBATCH --output=%x-%j.out
   ```
4. Save and exit
5. Submit the job: `sbatch patient_align_job.sh`
6. Check job status: `squeue -u your_uniqname`

**Expected Output**: `sbatch` will return a job ID (e.g., "Submitted batch job 12345678"). `squeue` will show your job as PENDING or RUNNING.

### Exercise 3: Mapping to the Human Genome (HG38)
**Objective**: Align the same patient reads to the human reference genome to compare alignment rates.
**Instructions**:
1. Create a new script: `nano patient1_hg38.sh`
2. Write the alignment pipeline using the HG38 index instead:
   ```bash
   #!/bin/bash
   #SBATCH --job-name=patient1_hg38
   #SBATCH --account=bioinf545w24_class
   #SBATCH --partition=standard
   #SBATCH --nodes=1
   #SBATCH --ntasks=1
   #SBATCH --cpus-per-task=1
   #SBATCH --mem=8G
   #SBATCH --time=01:00:00
   #SBATCH --output=%x-%j.out

   module load Bioinformatics
   module load Bowtie2
   module load SAMtools

   bowtie2 -x /path/to/hg38_index -U patient1.fastq -S patient1_hg38.sam
   samtools view -bS patient1_hg38.sam > patient1_hg38.bam
   samtools sort patient1_hg38.bam -o patient1_hg38.sorted.bam
   samtools index patient1_hg38.sorted.bam
   ```
3. Submit: `sbatch patient1_hg38.sh`
4. After the job completes, compare alignment rates from the `.out` file with the SARS-CoV-2 alignment rate
5. Note: HG38 alignment needs more memory (8G) and time (1 hour) because the human genome is much larger

**Expected Output**: The `.out` file will contain Bowtie2's alignment summary. Compare the alignment rate to SARS-CoV-2: if the patient is infected, expect high alignment to SARS-CoV-2 and relatively lower alignment to HG38.

### Exercise 4: Visualizing Alignments in IGV
**Objective**: Load alignment data into IGV to visually inspect where reads map to the SARS-CoV-2 genome.
**Instructions**:
1. Transfer sorted BAM and BAI files to your local machine:
   ```bash
   scp uniqname@greatlakes.arc-ts.umich.edu:~/path/patient1_sarscov2.sorted.bam ~/Desktop/
   scp uniqname@greatlakes.arc-ts.umich.edu:~/path/patient1_sarscov2.sorted.bam.bai ~/Desktop/
   ```
2. Open IGV on your local computer
3. In IGV, load the SARS-CoV-2 reference genome:
   - Go to Genomes > Load Genome from Server (or load a local FASTA)
   - Select SARS-CoV-2 (NC_045512.2) if available
4. Load your BAM file: File > Load from File > select `patient1_sarscov2.sorted.bam`
5. Navigate to the spike protein region: type `NC_045512.2:21,563-25,384` in the search bar
6. Zoom in to see individual reads and look for mismatches (colored bases)

**Expected Output**: You should see reads stacked up along the genome (coverage track at top, individual reads below). Any colored bases indicate differences from the reference (potential mutations).

### Exercise 5: Analyzing Multiple Patient Samples
**Objective**: Set up alignment jobs for all 8 patient samples using shell script modification.
**Instructions**:
1. Start with your working patient1 job script
2. For each patient (2 through 8), create a copy and modify the filenames:
   ```bash
   cp patient_align_job.sh patient2_align_job.sh
   nano patient2_align_job.sh
   ```
3. In nano, change every occurrence of `patient1` to `patient2`
4. Submit each modified script: `sbatch patient2_align_job.sh`
5. Track all jobs: `squeue -u your_uniqname`
6. After all jobs complete, compare the alignment rates from each `.out` file

**Expected Output**: Eight sets of output files. By comparing the SARS-CoV-2 alignment rates across all 8 patients, you can identify which patients are likely infected.

## Challenge Problems

### Challenge 1: Automated Multi-Sample Script
**Objective**: Write a single shell script that processes all 8 patient samples using a for loop.
**Instructions**:
1. Create `all_patients.sh` with a loop:
   ```bash
   #!/bin/bash
   #SBATCH --job-name=all_patients
   #SBATCH --account=bioinf545w24_class
   #SBATCH --partition=standard
   #SBATCH --nodes=1
   #SBATCH --ntasks=1
   #SBATCH --cpus-per-task=1
   #SBATCH --mem=4G
   #SBATCH --time=02:00:00
   #SBATCH --output=%x-%j.out

   module load Bioinformatics
   module load Bowtie2
   module load SAMtools

   for i in 1 2 3 4 5 6 7 8; do
       echo "Processing patient $i..."
       bowtie2 -x /path/to/sarscov2_index -U patient${i}.fastq -S patient${i}.sam
       samtools view -bS patient${i}.sam > patient${i}.bam
       samtools sort patient${i}.bam -o patient${i}.sorted.bam
       samtools index patient${i}.sorted.bam
       echo "Patient $i complete."
   done
   ```
2. Submit the job and review all results

### Challenge 2: Alignment Rate Summary Table
**Objective**: Create a summary of alignment rates for all patient samples.
**Instructions**:
1. After all patient alignment jobs complete, examine each `.out` file
2. Create a table with columns: Patient Number, Reads, SARS-CoV-2 Rate, HG38 Rate
3. Based on the data, write a short paragraph identifying which patients are likely positive for SARS-CoV-2 and explain your reasoning

---

## Answer Key

### Warm-Up Answers
1. Running commands interactively executes them on the login node, which has limited resources and is shared with all users. Submitting with `sbatch` sends the job to a compute node with dedicated resources (memory, CPUs). Always use `sbatch` for compute-heavy tasks like genome alignment; the login node is only for file management and small tasks.
2. Mapping to both genomes helps differentiate the source of the reads. Patient samples contain both human DNA and potentially viral RNA. A high alignment rate to SARS-CoV-2 with the same reads showing low alignment to HG38 strongly suggests viral infection. If reads align well to HG38 but not SARS-CoV-2, the patient is likely not infected.
3. IGV (Integrative Genomics Viewer) is a desktop application that displays alignment data visually. It shows reads stacked along the genome, making it easy to see coverage depth, identify mutations (colored mismatches), and explore specific genomic regions. It is much more intuitive than reading SAM/BAM files as text.

### Exercise Answers
1. **Shell Script**: The script runs four commands sequentially: Bowtie2 alignment, SAM-to-BAM conversion, sorting, and indexing. Testing interactively first confirms the commands work before submitting as a job.
2. **SLURM Job Script**: The `#SBATCH` lines are not comments --- SLURM reads them as directives. Key parameters: `--account` (billing), `--mem` (RAM), `--time` (wall clock limit), `--output` (where stdout goes, `%x` = job name, `%j` = job ID).
3. **HG38 Mapping**: The human genome is ~3.2 billion bp versus ~30,000 bp for SARS-CoV-2, so alignment requires more memory and time. Comparing alignment rates between the two references is the core analytical approach of the BTS project.
4. **IGV Visualization**: Both the `.sorted.bam` and `.sorted.bam.bai` files must be in the same directory for IGV to work. Colored bases in the read display indicate mismatches from the reference, which may represent real mutations or sequencing errors.
5. **Multiple Patients**: Creating separate scripts per patient is a simple approach. The challenge problem shows how a for-loop can automate this. The key finding is comparing alignment rates across all 8 patients to identify infected samples.

### Challenge Answers
1. **Loop Script**: The `for` loop iterates over patient numbers 1-8, using `${i}` variable substitution to change filenames. This processes all samples in a single job, though each alignment runs sequentially within the job.
2. **Summary Table**: Students should create a table comparing rates. Patients with >80% SARS-CoV-2 alignment rate are likely positive. Patients with <5% SARS-CoV-2 rate but high HG38 rate are likely negative. Ambiguous cases may need further analysis.
