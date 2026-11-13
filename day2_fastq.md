<!--
INSTRUCTOR NOTES (not rendered on the page)
Data: RNA-seq of healthy dogs from the preclinical body atlas, ENA PRJNA906635.
6 dogs (1001, 1003, 1005 male; 1002, 1004, 1006 female) x 2 tissues (cervical spinal cord, dorsal root ganglia).
Students use the full data in /data/courses/courseB/RNAseq/fastq_full/<dog>_<tissue>_R1/R2.fastq.gz
(4 lanes merged per sample; prep/rnaseq_download_full.slurm with run table prep/rnaseq_runs.tsv).
Each student analyses one sample: assign the 12 samples to the students (e.g. list on the board).
Per sample: 45-75 million read pairs, ~2-3 GB per file. Estimated run times (not measured on the cluster):
FastQC ~10-20 min per file, fastp ~10-20 min per sample with 4 threads. fastp output ~4-5 GB per student (home quota 20 GB).
Checked 2026-11-12 on SRR22458006: reads 73-76 bp (already lightly trimmed), adapters ~0.01%, poly-G tails ~0.4% in R2,
read names start with NB501995 (NextSeq 500) although ENA lists "Illumina HiSeq 2000"; quality characters E A / < 6 #.
Answer to the dog 1002 questions (FastQC of all 24 files, checked 2026-11-13): dog 1002 (both tissues) differs from the other five dogs:
  - Adapter Content: ~2.6% Nextera transposase sequence at the read ends (other dogs: ~0.01%); --detect_adapter_for_pe removes it
  - Sequence Length Distribution: 57/60-76 bp (other dogs: 35-76 bp), so the reads were trimmed differently
  - Per base sequence quality: mean Q ~21 at base 1 (other dogs ~31), which is why its 4 files FAIL this module; the rest of the read is fine
  - fewer duplicates (~75-80% of reads remain after deduplication vs 47-71%)
  Likely explanation: these libraries were prepared and/or sequenced in a different batch (different kit or run) than the other dogs.
  Other FastQC findings, all samples: Per base sequence content FAIL (random priming, normal), poly-G reads up to ~3% (all-G, removed by --trim_poly_g and -l 36),
  overrepresented sequences look like abundant small RNAs (e.g. 7SL/SRP RNA); fastp settings need no change.
Module names below are from the 2025 course; check them with module avail before class:
  FastQC/0.11.9-Java-11, fastp/0.23.4-GCC-10.3.0, MultiQC (only used in "For the brave").
-->

# FASTQ Quality Control: Can We Trust These Reads?

What makes the spinal cord different from the nerve tissue right next to it? Using RNA-seq from healthy dogs, we will compare the **spinal cord** with the **dorsal root ganglia** (DRG), which contain the sensory nerve cells, and find the genes that are active in one but not the other.

We use RNA-seq data from a public atlas of healthy dog tissues ([ENA PRJNA906635](https://www.ebi.ac.uk/ena/browser/view/PRJNA906635)): spinal cord and DRG from the same **6 dogs** (3 males, 3 females). That makes 12 samples, and **each of us analyses one of them**. Before analysing the data, we have to know whether it is good. Today we look inside FASTQ files, check their quality with **FastQC**, and clean them with **fastp**.

## 1. Set up our working directory and choose a sample

Log in to the Bioinformatics server with VS Code and open a terminal.

```shell
cd ~/course
mkdir -p dataPreProcess/scripts dataPreProcess/raw_fastqc dataPreProcess/clean_fastq
cd dataPreProcess
```

The reads are in a shared course directory. They are large, so we don't copy them: a **link** makes them appear in our directory without using extra space.

```shell
ln -s /data/courses/courseB/RNAseq/fastq_full raw_data
ls -lh raw_data/
```

These are the 12 samples. Each sample has two files, `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz`:

| Dog | Sex | Spinal cord | Dorsal root ganglia |
|---|---|---|---|
| 1001 | male | `1001_SpinalCord` | `1001_DRG` |
| 1002 | female | `1002_SpinalCord` | `1002_DRG` |
| 1003 | male | `1003_SpinalCord` | `1003_DRG` |
| 1004 | female | `1004_SpinalCord` | `1004_DRG` |
| 1005 | male | `1005_SpinalCord` | `1005_DRG` |
| 1006 | female | `1006_SpinalCord` | `1006_DRG` |

**Each of us uses the sample the instructor assigns.** Write its name down; we need it all day. In the examples below, replace `1003_DRG` with your own sample.

At the end of the day our directory will look like this:

```text
course
└── dataPreProcess
    ├── raw_data -> /data/courses/courseB/RNAseq/fastq_full
    ├── raw_fastqc      FastQC reports of the raw reads
    ├── clean_fastq     the cleaned reads and the fastp report
    └── scripts
```

Questions:

1. How big are the two files of our sample?
2. Why does each sample have two files, R1 and R2?

## 2. Look inside a FASTQ file

A **variable** saves us from typing the sample name again and again:

```shell
SAMPLE=1003_DRG
echo $SAMPLE
zcat raw_data/${SAMPLE}_R1.fastq.gz | head -n 8
```

Each read takes **4 lines**:

| Line | Content |
|---|---|
| 1 | `@` and the read name: run, sequencing machine, flow cell, lane and position on the flow cell |
| 2 | the sequence |
| 3 | `+` |
| 4 | the **quality** of each base, one character per base |

Now look at the first read of R2:

```shell
zcat raw_data/${SAMPLE}_R2.fastq.gz | head -n 4
```

Questions:

1. Do the first reads of R1 and R2 have the same name? Why?
2. How long is the sequence? How long is the quality line?
3. The database says these samples were sequenced on an "Illumina HiSeq 2000". The machine name in the read name starts with `NB5`, which is the code of an **Illumina NextSeq**. What does this tell us about trusting metadata?

### What do the quality characters mean?

Each character stands for a number, the **Phred quality score** Q. The computer stores characters as numbers (ASCII codes), and Q is the code minus 33. We can look up a character's code with `printf`:

```shell
printf '%d\n' "'E"
printf '%d\n' "'A"
printf '%d\n' "'/"
printf '%d\n' "'#"
```

So `E` is 69 − 33 = **Q36**. Q tells us the chance that the base is wrong:

| Q | Chance the base is wrong | Accuracy |
|---|---|---|
| 2 | about 6 in 10 | basically unknown |
| 10 | 1 in 10 | 90% |
| 20 | 1 in 100 | 99% |
| 30 | 1 in 1,000 | 99.9% |
| 36 | about 1 in 4,000 | 99.97% |

Questions:

1. Which quality characters do we see? Why are there only a few different ones?

   **Answer:** We see only a handful of characters: `E` (Q36), `A` (Q32), `<` (Q27), `6` (Q21), `/` (Q14) and `#` (Q2). The NextSeq does not report every possible quality value. It groups the values into a few **bins** and gives all bases in a bin the same score. This makes the files smaller and is precise enough for the analysis.
2. What quality does the base `N` usually have?

## 3. FastQC, interactively

FastQC checks many aspects of quality and writes an HTML report for each file. Documentation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/

Our files have tens of millions of reads, which takes FastQC a while. To try it out, we first make a small test file with the first 100,000 reads, on a compute node. Ask SLURM for an interactive session:

```shell
srun --partition=pibu_el8 --account=class-26980-ws-2026 --reservation=class-26980-ws-2026 --cpus-per-task=1 --mem=4G --time=01:00:00 --pty /bin/bash
```

Look at the prompt: which node are we on now? We are in a new shell, so we set the sample again:

```shell
SAMPLE=1003_DRG
cd ~/course/dataPreProcess
zcat raw_data/${SAMPLE}_R1.fastq.gz | head -n 400000 > test_R1.fastq
module load FastQC/0.11.9-Java-11
fastqc test_R1.fastq
ls
exit
```

Questions:

1. Why `head -n 400000` for 100,000 reads?
2. Which files did FastQC write?
3. After `exit`, which node are we on?

Delete the test files when done: `rm test_R1*`

## 4. FastQC as a batch job

For the full files, write a job script: then the exact commands and settings are saved, and we can put them in the methods section of a paper.

Save this as `run_fastqc.slurm` in `scripts`, with **your own** sample name:

```shell
#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --partition=pibu_el8
#SBATCH --account=class-26980-ws-2026
#SBATCH --reservation=class-26980-ws-2026
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=fastqc.out

SAMPLE=1003_DRG

module load FastQC/0.11.9-Java-11

cd ~/course/dataPreProcess
fastqc --threads 2 -o raw_fastqc raw_data/${SAMPLE}_R1.fastq.gz raw_data/${SAMPLE}_R2.fastq.gz
```

```shell
cd ~/course/dataPreProcess/scripts
sbatch run_fastqc.slurm
squeue -u $USER
```

With 2 CPUs, FastQC works on R1 and R2 at the same time. It takes about 10 to 20 minutes: a good moment for a coffee, or to read the FastQC help pages.

## 5. Read the FastQC reports

The reports are HTML files in `raw_fastqc`. To open them in VS Code:

1. Install the **Live Preview** extension by Microsoft, if it isn't installed yet: open the Extensions panel (Ctrl+Shift+X or Cmd+Shift+X), search for "Live Preview" and click Install.
2. In the VS Code Explorer on the left, go to `course/dataPreProcess/raw_fastqc`.
3. Right-click `<your sample>_R1_fastqc.html` and choose **Show Preview**.

Start with **Basic Statistics** and **Per base sequence quality**.

Questions:

1. How many reads are in each file? How long are they? Why are the reads not all exactly the same length?
2. Look at **Per base sequence quality**. Where in the read is the quality lowest? Compare R1 and R2: which one is worse?
3. Which modules have a warning (orange) or fail (red)? Several of them are **normal for RNA-seq**. Read the FastQC help for them:
   - **Per base sequence content**: the first ~10 bases are uneven. RNA-seq libraries are made with random primers; why could that cause it?
   - **Sequence Duplication Levels**: many reads are duplicated. In RNA-seq, what does it mean when a gene is very highly expressed?
   - **Overrepresented sequences**: copy one sequence and search it with BLAST (https://blast.ncbi.nlm.nih.gov). Which gene does it come from?
4. Compare the report with a neighbour who has the **other tissue of the same dog**. Do the tissues look different?
5. **One dog is different.** Find the classmates who have samples of dog **1002** and compare their reports with those of another dog. Look at **Per base sequence quality** (the first base), **Sequence Length Distribution** and **Adapter Content** (which adapter type?). What differs? What could explain it, given that all 12 samples come from the same study?

## 6. Clean the reads with fastp

**fastp** removes adapter sequences, cuts off poly-G tails and throws away reads of poor quality. It also writes its own HTML report, with the quality **before and after** cleaning. Documentation: https://github.com/OpenGene/fastp

Save this as `fastp_clean.slurm` in `scripts`, with **your own** sample name:

```shell
#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --partition=pibu_el8
#SBATCH --account=class-26980-ws-2026
#SBATCH --reservation=class-26980-ws-2026
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=fastp.out

SAMPLE=1003_DRG

module load fastp/0.23.4-GCC-10.3.0

cd ~/course/dataPreProcess
fastp \
    -i raw_data/${SAMPLE}_R1.fastq.gz \
    -I raw_data/${SAMPLE}_R2.fastq.gz \
    -o clean_fastq/${SAMPLE}_R1.clean.fastq.gz \
    -O clean_fastq/${SAMPLE}_R2.clean.fastq.gz \
    --detect_adapter_for_pe \
    --trim_poly_g \
    -q 15 \
    -l 36 \
    -w 4 \
    -h clean_fastq/${SAMPLE}_fastp.html \
    -j clean_fastq/${SAMPLE}_fastp.json
```

| Option | Meaning |
|---|---|
| `-i`, `-I` | the input R1 and R2 files |
| `-o`, `-O` | the cleaned R1 and R2 files |
| `--detect_adapter_for_pe` | find the adapter sequence and cut it off the reads |
| `--trim_poly_g` | cut off runs of G at the end of reads. On the NextSeq, "no signal" is read as G |
| `-q 15` | a base with quality below Q15 counts as "low quality". A read is thrown away if more than 40% of its bases are low quality |
| `-l 36` | throw away reads shorter than 36 bases after cleaning |
| `-w 4` | use 4 CPUs |
| `-h`, `-j` | the reports, as HTML (for us) and JSON (for programs) |
| `\` at the end of a line | the command continues on the next line |

If one read of a pair is thrown away, its mate is thrown away too, so R1 and R2 stay in the same order.

```shell
sbatch fastp_clean.slurm
squeue -u $USER
```

This also takes about 10 to 20 minutes. The cleaned reads take a few GB; we each have 20 GB in our home directory. Let's check how much we use:

```shell
du -sh ~/course/dataPreProcess/clean_fastq
du -sh ~
```

Open `clean_fastq/<your sample>_fastp.html` with Live Preview.

Questions:

1. How many reads were there before and after cleaning? What percentage was removed?
2. How many reads had adapters? How many had poly-G tails?
3. Compare the quality plots before and after filtering. What changed?
4. Write the sample name, number of reads and percentage removed in the class table. Which sample has the most reads removed? Is there a difference between spinal cord and DRG, or between male and female dogs?
5. For dog 1002: open the **Adapters** section of the fastp report. Which adapter sequence did fastp detect? Is it the same one FastQC reported?

**For the brave.** Combine the FastQC and fastp reports into one page with **MultiQC**. Find the module name with `module avail MultiQC`, then run in an interactive session:

```shell
cd ~/course/dataPreProcess
module load <MultiQC module>
multiqc raw_fastqc clean_fastq -o multiqc
```

Open `multiqc/multiqc_report.html` with Live Preview. MultiQC becomes really useful with many samples: at the end of the day we will look at a MultiQC report of the whole class.

## 7. Save our work with git

The FASTQ files and reports should not go into git. Add these lines to `~/course/.gitignore` with VS Code:

```text
# sequencing data and QC reports
*.fastq
*.fastq.gz
*.fq.gz
*_fastqc.html
*_fastqc.zip
*_fastp.html
*_fastp.json
```

```shell
cd ~/course
git add .gitignore dataPreProcess/scripts
git status
git commit -m "Add FASTQ quality control and cleaning scripts"
```

Questions:

1. Which files does `git status` show as untracked? Why is `raw_data` not listed?
2. Why should large data files not go into git? Where would we store them instead?

### What did we learn?

1. A FASTQ file has 4 lines per read; the quality line gives the error probability of each base.
2. Always check the quality of the data before analysing it, and check the metadata too.
3. Some FastQC warnings are normal for RNA-seq. Understand them before "fixing" them.
4. Clean reads with a job script, so the exact settings are saved for the methods section of the paper.

Later in the course, we map these RNA-seq reads to the dog reference genome and find out which genes are active in the spinal cord and in the DRG.
