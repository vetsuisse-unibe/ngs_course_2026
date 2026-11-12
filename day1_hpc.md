<!--
INSTRUCTOR NOTES (not rendered on the page)
Students need small paired FASTQ files for the three affected dogs in /data/courses/courseB/EFNB3/fastq_subset/:
  WE0025_R1.fastq.gz WE0025_R2.fastq.gz
  WE0028_R1.fastq.gz WE0028_R2.fastq.gz
  WE0029_R1.fastq.gz WE0029_R2.fastq.gz
A few million read pairs each is enough (e.g. the reads extracted for the chr5 region, see variant calling notes).
Check before the course: clusterstate.sh still exists. Course reservation (one node): account and reservation class-26980-ws-2026, valid during and after class hours.
Facts for question 5.3/5.4 (ENA, checked 2026-10-15): WE0028 run ERR12042749 = 467,919,136 reads, ~70.6 Gb, 21.5 GB of FASTQ; ~29x coverage.
-->

# HPC Exercises: Putting the Cluster to Work

The **cluster** has many computers (nodes) that can share the work, and a scheduler called **SLURM** decides who runs what, where and when.

## 1. Connect to the cluster

Log in to the Bioinformatics server (login8.hpc.binf.unibe.ch) with the Remote SSH extension of Visual Studio Code, as in the Unix exercises, and open a terminal.

The computer you are now on is the **login node**. It is shared by everybody and only meant for small tasks: editing files, copying, submitting jobs. **Real work runs on the compute nodes.**

## 2. Explore the cluster

To see all compute nodes and how busy they are:

```shell
clusterstate.sh
```

| Column | Description |
|---------|-------------|
| `NODELIST` | Node name (machine identifier) |
| `CPUS(A/I/O/T)` | CPUs **Allocated / Idle / Other / Total** |
| `CPU_LOAD` | Average CPU load on that node |
| `FREEMEMORY` | Approximate available memory (in MB) |
| `STATE` | Current node status, e.g. `idle`, `mixed`, `down`, `completing` |

| State | Meaning |
|--------|----------|
| `idle` | Node is available for new jobs |
| `mixed` | Node is running some jobs but still has free CPUs |
| `alloc` | Node is fully allocated to running jobs |
| `completing` | Node is finishing current jobs |
| `down` | Node is offline or under maintenance |

For this course, one node has been **reserved** for us. Every job you submit this week must use the course **account** and **reservation**, during and after class hours. To see the reservation and its node:

```shell
scontrol show reservation class-26980-ws-2026
```

Questions:

1. Which nodes are **idle** and which are **down**?
2. How many CPUs does the course node have? How many students are in this course? How many CPUs is that per person if everybody runs a job at the same time?
3. Why do some nodes show as *mixed* instead of *idle* or *alloc*?

### Software modules

Bioinformatics programs are installed as **modules**. You load the ones you need, in the version you need.

```shell
module avail FastQC
module avail BWA
module avail SAMtools
module avail GATK
```

```shell
module load SAMtools/1.13-GCC-10.3.0
module list
samtools --version
module purge
module list
```

Questions:

1. Which versions of BWA are installed?
2. Why is it important to write down the exact version of every program you use? (Think about the methods section of a paper.)

## 3. Set up git to track your scripts

Every analysis script you write this week will be saved with **git**. Git keeps every version of your scripts, so you can always see what you changed and go back to a version that worked.

```shell
git -v
```

Git commands follow the structure `git verb options`: the verb tells git what to do (`add`, `commit`, …) and the options give the details.

Tell git who you are. You only do this once per computer. Use your own name, and the email address you will use for your GitHub account later.

```shell
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
git config --list --show-origin
```

### Create the course repository

```shell
cd
mkdir course
cd course
git init
mkdir -p hpc-exercises/scripts
cd hpc-exercises
```

`git init` creates a hidden `.git` directory where git stores the history of the project.

Some files should **never** go into git: big data files, and the log files that SLURM writes for every job. Open a new file in VS Code, add these lines and save it as `.gitignore` in `hpc-exercises`:

```text
# SLURM log files
*.out
*.err
slurm-*.out

# sequencing data
*.fastq.gz
*.fq.gz
```

```shell
git add .gitignore
git commit -m "Add .gitignore for SLURM logs and sequencing data"
```

At the end of the day your directory will look like this:

```text
/home/<student>
└── course
    └── hpc-exercises
        └── scripts
```

## 4. Your first job

Open a new file in VS Code, paste these lines and save it as `first_job.sh` in `course/hpc-exercises/scripts`:

```shell
#!/bin/bash
hostname
date
sleep 30
date
```

- `#!/bin/bash` (the "shebang") tells the system this is a bash script.
- `hostname` prints the name of the computer the script runs on.
- `date` prints the date and time; `sleep 30` waits 30 seconds.

```shell
cd ~/course/hpc-exercises/scripts
chmod +x first_job.sh
git add first_job.sh
git commit -m "Add first test job"
```

Run it on the login node (only because it is tiny):

```shell
./first_job.sh
```

Now let SLURM run it on a compute node:

```shell
sbatch --account=class-26980-ws-2026 --reservation=class-26980-ws-2026 -t 00:05:00 first_job.sh
squeue -u $USER
```

- `sbatch` submits a job script to SLURM.
- `--account` and `--reservation` send the job to the course node reserved for us; `-t 00:05:00` is the time limit (5 minutes).
- `squeue -u $USER` lists your jobs. `$USER` is replaced by your user name.

Questions:

1. What is the job ID?
2. Where did the output go? Compare the `hostname` in it with the one you got on the login node.
3. Run `git status`. Is the new output file tracked? Why not?

## 5. Count the puppies' reads

Now a real job. How many reads were sequenced for each puppy? We have a small subset of each dog's reads for the course:

```shell
ls -lh /data/courses/courseB/EFNB3/fastq_subset/
```

Each dog has two files, R1 and R2: the two ends of the same DNA fragments (**paired-end** sequencing).

Look at the first read of one file. A FASTQ file has **4 lines per read**: a name, the sequence, a `+` line and the quality of each base.

```shell
zcat /data/courses/courseB/EFNB3/fastq_subset/WE0025_R1.fastq.gz | head -n 4
```

So the number of reads in a file is its number of lines divided by 4. Counting the lines of big files takes a while, so we let a compute node do it.

This time we write the SLURM options **inside** the script, on lines starting with `#SBATCH`. Then you never have to remember them again. Save this as `count_reads.slurm` in `scripts`:

```shell
#!/bin/bash
#SBATCH --job-name=countReads
#SBATCH --account=class-26980-ws-2026
#SBATCH --reservation=class-26980-ws-2026
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --output=countReads.out

for FILE in /data/courses/courseB/EFNB3/fastq_subset/*.fastq.gz
do
    echo $FILE
    zcat $FILE | wc -l
done
```

| Line | Meaning |
|---|---|
| `--job-name` | a name you will recognise in `squeue` |
| `--account`, `--reservation` | run on the course node reserved for us |
| `--mem=1G` | 1 GB of memory |
| `--time=00:30:00` | stop the job after 30 minutes |
| `--output=countReads.out` | write everything the script prints into this file |
| `for FILE in … do … done` | repeat the two lines in between for every FASTQ file, one after the other |
| `zcat $FILE \| wc -l` | unpack the file and count its lines |

```shell
git add count_reads.slurm
git commit -m "Add read counting job"
sbatch count_reads.slurm
squeue -u $USER
```

When the job has finished:

```shell
cat countReads.out
```

Questions:

1. How many reads does each file have? (Divide the number of lines by 4.)
2. Do R1 and R2 of the same dog have the same number of reads? Why must they?
3. How long is one read? Count the letters of the sequence line you looked at above.
4. The dog genome is about 2.4 billion bases long. The full data set of WE0028 has 70 billion bases. How many times was each position of its genome sequenced on average? This is called the **coverage**.

## 6. Do it in parallel: job arrays

The job above counts the files **one after the other**. With 3 dogs that's fine. With the 1,489 control dogs used in the paper, it would take a long time. A **job array** runs the same script several times **at the same time**, once for each dog. SLURM gives each copy (each **task**) its own number, in the variable `$SLURM_ARRAY_TASK_ID`.

Save this as `count_reads_array.slurm` in `scripts`:

```shell
#!/bin/bash
#SBATCH --job-name=countArray
#SBATCH --account=class-26980-ws-2026
#SBATCH --reservation=class-26980-ws-2026
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --array=0-2
#SBATCH --output=countArray_%a.out

# the three dogs: task 0 is WE0025, task 1 is WE0028, task 2 is WE0029
DOGS=(WE0025 WE0028 WE0029)
DOG=${DOGS[$SLURM_ARRAY_TASK_ID]}

echo "Task $SLURM_ARRAY_TASK_ID is counting the lines of dog $DOG"
zcat /data/courses/courseB/EFNB3/fastq_subset/${DOG}_R1.fastq.gz | wc -l
```

| Line | Meaning |
|---|---|
| `--array=0-2` | run 3 tasks, numbered 0, 1 and 2 |
| `--output=countArray_%a.out` | one output file per task; `%a` is replaced by the task number |
| `DOGS=(WE0025 WE0028 WE0029)` | a list of the three dogs. Bash counts from 0: `${DOGS[0]}` is WE0025 |
| `DOG=${DOGS[$SLURM_ARRAY_TASK_ID]}` | each task picks its own dog from the list |

```shell
git add count_reads_array.slurm
git commit -m "Add job array for read counting"
sbatch count_reads_array.slurm
squeue -u $USER
```

When the tasks have finished:

```shell
cat countArray_*.out
```

Questions:

1. How many tasks do you see in `squeue`? What do their job IDs look like?
2. Does each task print the same number of lines for its dog's R1 file as `countReads.out` in section 5?
3. To analyse 100 dogs, which two lines of the script would you change?

## 7. Interactive jobs

Sometimes you want to work **on** a compute node yourself, for example to look at alignments with `samtools tview` on Day 2. `srun --pty` gives you a terminal on a compute node:

```shell
srun --account=class-26980-ws-2026 --reservation=class-26980-ws-2026 --cpus-per-task=1 --mem=2G --time=00:30:00 --pty /bin/bash
hostname
module load SAMtools/1.13-GCC-10.3.0
samtools --version | head -n 1
exit
```

Question: What changed in your prompt after `srun`, and after `exit`?

## 8. Cancel a job

Made a mistake? `scancel` stops a job. Save this as `long_job.slurm`:

```shell
#!/bin/bash
#SBATCH --job-name=longJob
#SBATCH --account=class-26980-ws-2026
#SBATCH --reservation=class-26980-ws-2026
#SBATCH --mem=1G
#SBATCH --time=00:20:00
#SBATCH --output=longJob.out

hostname
date
sleep 600
date
```

```shell
sbatch long_job.slurm
squeue -u $USER
```

Wait until the job is running (`R` in the `ST` column), then cancel it with its job ID:

```shell
scancel <job_id>
squeue -u $USER
```

`scancel -u $USER` cancels **all** your jobs. Use it with care!

Question: Look at `longJob.out`. Which `date` lines were printed?

## 9. How much did your jobs use?

On a shared cluster, every CPU and gigabyte you request is blocked for everybody else, even if your job doesn't use it. `sacct` shows what your jobs requested and what they really used:

```shell
sacct -u $USER --format=JobID,JobName,State,Elapsed,ReqMem,MaxRSS,AllocCPUS,NodeList
```

- `Elapsed`: how long the job ran
- `ReqMem`: memory you requested
- `MaxRSS`: the most memory the job actually used (shown on the `.batch` line)

Questions:

1. Look at the `Elapsed` time of `countReads` (all six files, one after the other) and of the three `countArray` tasks. Why did the array finish sooner?
2. You requested 1 GB for the read counting. How much was used? What would be a sensible request next time?
3. The course node has 128 CPUs and 512 GB of memory for about 10 students. What happens to the others if each student asks for 100 GB of memory for a job that needs 200 MB?
4. What is the `State` of the job you cancelled?

## 10. Save your work

```shell
cd ~/course
git status
git log --oneline
```

Questions:

1. How many commits have you made today?
2. Which files does `git status` list as untracked? Why are the `.out` files not listed?

### Good habits

1. Put the SLURM options in the script (`#SBATCH`), not on the command line.
2. Request what your job needs, not more. Check with `sacct`.
3. Never run heavy work on the login node.
4. Commit your scripts before you run them, with a clear message.
5. Don't put data or log files in git.

Tomorrow you will run real bioinformatics programs as SLURM jobs: quality control and mapping of the puppies' reads.
