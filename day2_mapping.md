<!--
INSTRUCTOR NOTES (not rendered on the page)
Reference: /data/courses/courseB/UnixEx/dog_chr5.fa = chromosome 5 (NC_049226.1, 89,562,946 bp) of UU_Cfam_GSD_1.0 from NCBI (same file as Day 1).
Reads: /data/courses/courseB/EFNB3/fastq_chr5/<dog>_chr5_R1/R2.fastq.gz for WE0025, WE0028, WE0029
(chr5 read pairs extracted from the deposited BAMs, see prep/extract_chr5_reads.slurm; one read group per dog, lane information is lost).
Load per student: indexing 1 CPU; mapping job array 3 tasks x 4 CPUs x 12 GB. 10 students = 120 CPUs, fits the reserved node (128).
Run times are estimates, not measured on the cluster: bwa index on chr5 a few minutes; mapping roughly 15-30 min per dog with 4 threads.
The EFNB3 variant position is deliberately NOT used here; students discover it on Day 3.
Module names are from the 2025 course (SAMtools/1.13-GCC-10.3.0 confirmed); check BWA and GATK with module avail.
-->

# Mapping: Where Do the Reads Belong?

Each FASTQ file is a huge pile of short DNA pieces, without any information on where in the genome they come from. **Mapping** finds the place in the reference genome where each read fits best. Only then can we compare the puppies' DNA with the reference and find their variants on Day 3.

We map the reads of the three affected puppies (WE0025, WE0028, WE0029) to **chromosome 5** of the dog reference genome (UU_Cfam_GSD_1.0), the chromosome you explored on Day 1. Using one chromosome keeps the run times short.

## 1. Set up your working directory

```shell
cd ~/course
mkdir -p variantCalling/refIdx variantCalling/mapping variantCalling/scripts
cd variantCalling
```

The reads are in a shared course directory. As yesterday, we link them instead of copying:

```shell
ln -s /data/courses/courseB/EFNB3/fastq_chr5 fastq_data
ls -lh fastq_data/
```

The reference is the chromosome 5 file from Day 1. We copy it, because we will add index files next to it:

```shell
cp /data/courses/courseB/UnixEx/dog_chr5.fa refIdx/
ls -lh refIdx/
```

At the end of the day your directory will look like this:

```text
course
└── variantCalling
    ├── fastq_data -> /data/courses/courseB/EFNB3/fastq_chr5
    ├── refIdx       dog_chr5.fa and its index files
    ├── mapping      one sorted BAM file per dog, with index and statistics
    └── scripts
```

Questions:

1. How many FASTQ files are there? How many dogs?
2. What is the name of the chromosome in `dog_chr5.fa`? (Hint: `grep ">" refIdx/dog_chr5.fa`)

## 2. Index the reference genome

To find the right place for millions of reads quickly, the mapping program needs an **index** of the reference, much like the index at the back of a book. You build it once and can use it for every mapping to the same reference.

We make three indexes, each for different programs:

| Command | Makes | Used by |
|---|---|---|
| `bwa index` | `.amb .ann .bwt .pac .sa` | BWA, to map the reads |
| `samtools faidx` | `.fai` | samtools and IGV, to jump to any position |
| `gatk CreateSequenceDictionary` | `.dict` | GATK, for variant calling on Day 3 |

Save this as `index_ref.slurm` in `scripts`:

```shell
#!/bin/bash
#SBATCH --job-name=indexRef
#SBATCH --partition=pibu_el8
#SBATCH --account=class-26980-ws-2026
#SBATCH --reservation=class-26980-ws-2026
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=indexRef.out

module load BWA/0.7.17-GCC-10.3.0
module load SAMtools/1.13-GCC-10.3.0
module load GATK/4.2.6.1-GCCcore-10.3.0-Java-11

cd ~/course/variantCalling/refIdx
bwa index dog_chr5.fa
samtools faidx dog_chr5.fa
gatk CreateSequenceDictionary -R dog_chr5.fa -O dog_chr5.dict
```

```shell
cd ~/course/variantCalling/scripts
sbatch index_ref.slurm
squeue -u $USER
```

When the job has finished (a few minutes):

```shell
ls -lh ../refIdx/
cat ../refIdx/dog_chr5.fa.fai
```

Questions:

1. Which files were created? Which one is the biggest?
2. The `.fai` file has one line per sequence. What do the first two columns tell you?

## 3. Map the reads

We use **BWA-MEM**, the most widely used program to map short DNA reads to a genome. Mapping usually takes three steps:

1. **map** the reads (`bwa mem`), which writes a text file in **SAM** format,
2. **convert** the SAM file to the compressed binary **BAM** format (`samtools view`),
3. **sort** the reads by their position on the chromosome (`samtools sort`).

Writing the big intermediate files to disk is slow. With a **pipe** (`|`) the output of `bwa mem` goes straight into `samtools sort`, which writes a sorted BAM file. Only the final file is written to disk.

We map all three dogs at the same time with a **job array**, as on Day 1. Save this as `mapping.slurm` in `scripts`:

```shell
#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH --partition=pibu_el8
#SBATCH --account=class-26980-ws-2026
#SBATCH --reservation=class-26980-ws-2026
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --array=0-2
#SBATCH --output=mapping_%a.out

module load BWA/0.7.17-GCC-10.3.0
module load SAMtools/1.13-GCC-10.3.0

# the three dogs: task 0 is WE0025, task 1 is WE0028, task 2 is WE0029
DOGS=(WE0025 WE0028 WE0029)
DOG=${DOGS[$SLURM_ARRAY_TASK_ID]}

cd ~/course/variantCalling

bwa mem -t 4 \
    -R "@RG\tID:${DOG}\tSM:${DOG}\tLB:${DOG}_LIB\tPL:ILLUMINA" \
    refIdx/dog_chr5.fa \
    fastq_data/${DOG}_chr5_R1.fastq.gz \
    fastq_data/${DOG}_chr5_R2.fastq.gz \
  | samtools sort -@ 4 -m 1G -o mapping/${DOG}.sorted.bam -

samtools index mapping/${DOG}.sorted.bam
```

| Part | Meaning |
|---|---|
| `bwa mem -t 4` | map with 4 CPUs |
| `-R "@RG\t..."` | add a **read group**: a label on every read that says which dog it comes from (see below) |
| `refIdx/dog_chr5.fa` | the reference; BWA finds the index files next to it |
| `\| samtools sort … -` | sort the mapped reads by position; `-` means "read from the pipe" |
| `-m 1G` | memory per sorting thread |
| `samtools index` | make a `.bai` index, so programs can jump to any position in the BAM |

**Read groups** matter as soon as you analyse several samples together. On Day 3, GATK uses the `SM` (sample) name to know which reads belong to which dog:

| Tag | Meaning | Here |
|---|---|---|
| `ID` | unique name of this read group | the dog's name |
| `SM` | sample: the individual the DNA comes from | the dog's name |
| `LB` | library: the DNA preparation that was sequenced | `<dog>_LIB` |
| `PL` | sequencing platform | `ILLUMINA` |

```shell
git -C ~/course add variantCalling/scripts
git -C ~/course commit -m "Add reference indexing and mapping scripts"
sbatch mapping.slurm
squeue -u $USER
```

Mapping takes roughly 15 to 30 minutes per dog. The three tasks run at the same time.

Questions:

1. Why is the time limit (`--time`) much longer than the expected run time?
2. While you wait: open the BWA manual (http://bio-bwa.sourceforge.net/bwa.shtml) and find what `-t` does. What else could you change with options?

## 4. Check the mapping

When the jobs have finished, look at the result. `samtools view` shows the reads in a BAM file as text, `-H` shows only the header:

```shell
cd ~/course/variantCalling/mapping
module load SAMtools/1.13-GCC-10.3.0
ls -lh
samtools view -H WE0025.sorted.bam | grep -v "^@PG"
samtools view WE0025.sorted.bam | head -n 3
```

Questions:

1. Where in the header do you find your read group?
2. Each read line has 11 columns plus extra tags. Find the columns for the chromosome, the position, the mapping quality and the CIGAR string (for example `100M`). What does `M` mean? (Look up the SAM format: https://samtools.github.io/hts-specs/SAMv1.pdf, section 1.4.)

### Mapping statistics

`samtools flagstat` counts how many reads mapped, how many pairs mapped properly, and more. Run it for all three dogs in one go:

```shell
for DOG in WE0025 WE0028 WE0029
do
    samtools flagstat ${DOG}.sorted.bam > ${DOG}_flagstat.txt
done
cat WE0025_flagstat.txt
```

`samtools coverage` tells you how much of the chromosome is covered by reads, and how deeply:

```shell
samtools coverage WE0025.sorted.bam
```

Questions:

1. What percentage of reads mapped in each dog? What percentage are **properly paired**?
2. Why could some reads not map, even though they all came from chromosome 5?
3. flagstat shows 0 duplicates. Does that mean there are no duplicate reads? (Hint: we haven't marked them yet. That's on Day 3.)
4. What is the mean **depth** (column `meandepth`) for each dog? Is that enough to call a heterozygous variant reliably?
5. What percentage of chromosome 5 is covered (`coverage`)? Why not 100%?

## 5. Look at the alignments

### samtools tview

`samtools tview` shows the alignment as text, right in the terminal. Start an interactive session first:

```shell
srun --partition=pibu_el8 --account=class-26980-ws-2026 --reservation=class-26980-ws-2026 --mem=4G --time=01:00:00 --pty /bin/bash
cd ~/course/variantCalling
module load SAMtools/1.13-GCC-10.3.0
samtools tview mapping/WE0025.sorted.bam refIdx/dog_chr5.fa
```

| Key | Action |
|---|---|
| arrow keys | move left and right |
| `g` | go to a position, e.g. `NC_049226.1:10000000` |
| `.` | switch between dots (same as reference) and letters |
| `?` | help |
| `q` | quit |

Press `g`, type `NC_049226.1:10000000` and press Enter. Scroll right slowly.

Questions:

1. In the dot view, what do `.` and `,` stand for?
2. Find a position where some reads show a different base from the reference. Is it in **all** reads or in about **half** of them? What would each mean for the dog's genotype?
3. Find a place where reads are missing, or where all reads disagree with the reference. What could be the reason?

Type `q` to quit, and `exit` to leave the compute node.

### IGV

**IGV** (Integrative Genomics Viewer) shows alignments graphically on your own computer. The BAM files are large, so first cut out a smaller region of 1 Mb on the cluster:

```shell
cd ~/course/variantCalling/mapping
module load SAMtools/1.13-GCC-10.3.0
for DOG in WE0025 WE0028 WE0029
do
    samtools view -b -o ${DOG}.region.bam ${DOG}.sorted.bam NC_049226.1:10000000-11000000
    samtools index ${DOG}.region.bam
done
ls -lh *.region.bam*
```

Download to your computer, with right-click → **Download** in the VS Code Explorer:

- `refIdx/dog_chr5.fa` and `refIdx/dog_chr5.fa.fai`
- the three `.region.bam` files and their `.region.bam.bai` files

Then in IGV:

1. **Genomes → Load Genome from File…** and choose `dog_chr5.fa`. (IGV's built-in dog genomes name the chromosome `chr5`, but our reference calls it `NC_049226.1`, so we load our own.)
2. **File → Load from File…** and choose the three `.region.bam` files. Each needs its `.bai` file in the same folder.
3. Type `NC_049226.1:10,000,000-10,001,000` in the search box and press Go.

Questions:

1. What does the grey histogram track above the reads show?
2. Find a variant position. Do all three dogs have the same genotype there?
3. Right-click on a read and choose **Color alignments by → read strand**. Are the variant bases on both strands?

## 6. Save your work with git

BAM files and index files are far too large for git. Add these lines to `~/course/.gitignore`:

```text
# reference, index and alignment files
*.fa
*.fa.*
*.dict
*.bam
*.bai
*.sam
```

```shell
cd ~/course
git add .gitignore variantCalling/scripts
git status
git commit -m "Ignore reference and alignment files"
```

Question: Would you track `WE0025_flagstat.txt` in git? Why or why not?

### What did we learn?

1. Mapping finds where each read belongs; an index of the reference makes it fast.
2. Pipes avoid writing large intermediate files.
3. Read groups label every read with its sample, which GATK needs on Day 3.
4. Always check the mapping (flagstat, coverage) before you call variants.

On Day 3 you use your three BAM files to find the variant that makes the puppies hop.
