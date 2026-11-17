<!--
INSTRUCTOR NOTES (not rendered on the page)
Data: RNA-seq of 6 healthy dogs (1001, 1003, 1005 male; 1002, 1004, 1006 female), spinal cord and DRG, ENA PRJNA906635.
Shared files needed (made with the scripts in prep/):
  /data/courses/courseB/RNA-seq/reference/hisat2/UU_Cfam_GSD_1.0_ROSY.*.ht2   rnaseq_reference.slurm
  /data/courses/courseB/RNA-seq/reference/UU_Cfam_GSD_1.0_ROSY.sorted.gtf            ROSY GTF: gene_id = symbol_N (e.g. ENPP1_3); the "gene" attribute is the plain symbol, so featureCounts uses -g gene
  data/dog_SC_DRG_counts.tsv and data/dog_SC_DRG_samples.tsv in THIS repository (copied from
  /data/courses/courseB/RNA-seq/counts after rnaseq_map_full.slurm + rnaseq_featurecounts.slurm)
Part 1-2: each student selects one of the 12 samples and maps its full read pair (R1 + R2, 40-73 million pairs) from fastq_full, then counts it.
These are the same raw files as in rnaseq_map_full.slurm, so a student's counts match that sample's column in the class count table.
Estimates, not measured: HISAT2 ~30-60 min with 8 threads, ~8 GB memory; sorted BAM ~3-5 GB per student (home quota 20 GB).
Part 3: DESeq2 in RStudio on the students' laptops, with the full count table read from GitHub.
Library strandedness is not given in the metadata; students test it in Part 2, and rnaseq_featurecounts.slurm picks it automatically.
The reference is a female German Shepherd: no Y chromosome, and XIST is not annotated under that name, so sex genes cannot be used as a check.
Module names from the 2025 course; check: HISAT2/2.2.1-gompi-2021a, SAMtools/1.13-GCC-10.3.0, Subread/2.0.3-GCC-10.3.0.
-->

# RNA-seq: What Makes the Spinal Cord Different?

On Day 2 you checked and cleaned the RNA-seq reads of your sample: spinal cord or dorsal root ganglia (DRG) from one of six healthy dogs. Today we find out **which genes are active**. First you map a part of your own reads to the dog genome and count how many reads fall on each gene. Then we compare all 12 samples in R to find the genes that differ between spinal cord and DRG.

Two quick reminders of the biology:

- The **spinal cord** belongs to the central nervous system. It contains motor neurons, interneurons and their support cells, such as oligodendrocytes, which make the myelin around nerve fibres in the brain and spinal cord.
- The **dorsal root ganglia** sit just outside the spinal cord. They contain the cell bodies of the **sensory** neurons (touch, pain, temperature, position) and Schwann cells, which make the myelin in the peripheral nervous system.

Before we start, let's predict:

1. Which kinds of genes do we expect to be more active in the spinal cord? Which in the DRG?
2. The puppies' disease gene, *EFNB3*, matters most while the spinal cord develops. These dogs are adults. Do we expect it to be strongly expressed?

## Part 1: Map your reads with HISAT2

RNA is spliced: the introns are removed, so one read can span two exons that lie thousands of bases apart in the genome. BWA, which we used for DNA, would fail on such reads. **HISAT2** is a **splice-aware** aligner: it can split a read across an intron.

```
genome:  ──[ exon 1 ]────────── intron ──────────[ exon 2 ]──
read:          [ ──────── ]                     [ ──── ]
               one read, split over two exons
```

The 12 samples are in the shared directory `/data/courses/courseB/RNA-seq/fastq_full/`. **Choose any one read pair (one sample) from this table** and maps both of its files, R1 and R2:

| Sample | Dog | Sex | Tissue | R1 file | R2 file | Read pairs |
|---|---|---|---|---|---|---|
| `1001_SpinalCord` | 1001 | male | spinal cord | `1001_SpinalCord_R1.fastq.gz` | `1001_SpinalCord_R2.fastq.gz` | 73.1 M |
| `1001_DRG` | 1001 | male | DRG | `1001_DRG_R1.fastq.gz` | `1001_DRG_R2.fastq.gz` | 63.7 M |
| `1002_SpinalCord` | 1002 | female | spinal cord | `1002_SpinalCord_R1.fastq.gz` | `1002_SpinalCord_R2.fastq.gz` | 43.9 M |
| `1002_DRG` | 1002 | female | DRG | `1002_DRG_R1.fastq.gz` | `1002_DRG_R2.fastq.gz` | 45.3 M |
| `1003_SpinalCord` | 1003 | male | spinal cord | `1003_SpinalCord_R1.fastq.gz` | `1003_SpinalCord_R2.fastq.gz` | 56.9 M |
| `1003_DRG` | 1003 | male | DRG | `1003_DRG_R1.fastq.gz` | `1003_DRG_R2.fastq.gz` | 52.8 M |
| `1004_SpinalCord` | 1004 | female | spinal cord | `1004_SpinalCord_R1.fastq.gz` | `1004_SpinalCord_R2.fastq.gz` | 52.1 M |
| `1004_DRG` | 1004 | female | DRG | `1004_DRG_R1.fastq.gz` | `1004_DRG_R2.fastq.gz` | 39.9 M |
| `1005_SpinalCord` | 1005 | male | spinal cord | `1005_SpinalCord_R1.fastq.gz` | `1005_SpinalCord_R2.fastq.gz` | 55.0 M |
| `1005_DRG` | 1005 | male | DRG | `1005_DRG_R1.fastq.gz` | `1005_DRG_R2.fastq.gz` | 60.4 M |
| `1006_SpinalCord` | 1006 | female | spinal cord | `1006_SpinalCord_R1.fastq.gz` | `1006_SpinalCord_R2.fastq.gz` | 45.3 M |
| `1006_DRG` | 1006 | female | DRG | `1006_DRG_R1.fastq.gz` | `1006_DRG_R2.fastq.gz` | 59.4 M |

Write down the chosen sample; we need its name in all scripts below. In Part 3 we analyse all 12 samples together.

```shell
cd ~/course
mkdir -p rnaseq/scripts rnaseq/mapping
cd rnaseq
```

Save this as `hisat2_map.slurm` in `scripts`, with the name of **the sample you selected**:

```shell
#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH --partition=pibu_el8
#SBATCH --account=class-26980-ws-2026
#SBATCH --reservation=class-26980-ws-2026
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=hisat2.out

SAMPLE=1003_DRG

module load HISAT2/2.2.1-gompi-2021a
module load SAMtools/1.13-GCC-10.3.0

INDEX=/data/courses/courseB/RNA-seq/reference/hisat2/UU_Cfam_GSD_1.0_ROSY
FASTQ=/data/courses/courseB/RNA-seq/fastq_full

cd ~/course/rnaseq/mapping

hisat2 -p 8 -x $INDEX \
    -1 $FASTQ/${SAMPLE}_R1.fastq.gz \
    -2 $FASTQ/${SAMPLE}_R2.fastq.gz \
    --summary-file ${SAMPLE}.hisat2.txt \
  | samtools sort -@ 4 -m 1G -o ${SAMPLE}.sorted.bam -
samtools index ${SAMPLE}.sorted.bam
```

| Part | Meaning |
|---|---|
| `hisat2 -p 8` | map with 8 CPUs |
| `-x $INDEX` | the index of the whole dog genome, made in advance |
| `-1`, `-2` | R1 and R2 of the read pairs; HISAT2 reads the gzipped files directly |
| `--summary-file` | save the mapping statistics |
| `\| samtools sort` | sort the mapped reads by position, as on Day 2 |

Anyone who selected their own Day 2 sample could also map its cleaned reads instead: `FASTQ=~/course/dataPreProcess/clean_fastq`, and add `.clean` before `.fastq.gz` in the two file names. The quality of this data is high, so the result is almost the same.

The sorted BAM file of a full sample takes a few GB, and you have 20 GB in your home directory. Check the free space first, and if needed delete the chr5 BAM files from Day 2. The dedup BAM files from Day 3 contain the same reads:

```shell
du -sh ~
rm ~/course/variantCalling/mapping/*.sorted.bam*
```

```shell
cd ~/course/rnaseq/scripts
sbatch hisat2_map.slurm
squeue -u $USER
```

When the job has finished (roughly 30 to 60 minutes; a good moment to read about featureCounts below):

```shell
cat ~/course/rnaseq/mapping/*.hisat2.txt
```

Questions:

1. What is the **overall alignment rate**? How does it compare with the DNA mapping on Day 2?
2. How many pairs aligned **concordantly exactly 1 time**, and how many **more than once**? What could make a read map to several places?
3. Let's compare the numbers between neighbours. Are they similar for all samples?

## Part 2: Count reads per gene

To measure how active each gene is, we count how many read pairs fall on its exons. **featureCounts** does this using the gene annotation (GTF file) of the dog genome.

One thing we don't know is whether the library is **stranded**. In a stranded library, every read tells us which DNA strand the RNA came from; in an unstranded library it doesn't. featureCounts needs to know (`-s 0` unstranded, `-s 1` stranded, `-s 2` reversely stranded). Let's find out by trying all three. Start an interactive session:

```shell
srun --partition=pibu_el8 --account=class-26980-ws-2026 --reservation=class-26980-ws-2026 --cpus-per-task=4 --mem=8G --time=01:00:00 --pty /bin/bash
SAMPLE=1003_DRG
cd ~/course/rnaseq/mapping
module load Subread/2.0.3-GCC-10.3.0
GTF=/data/courses/courseB/RNA-seq/reference/UU_Cfam_GSD_1.0_ROSY.sorted.gtf

for S in 0 1 2
do
    featureCounts -T 4 -p --countReadPairs -s $S -t exon -g gene \
        -a $GTF -o ${SAMPLE}_counts_s$S.txt ${SAMPLE}.sorted.bam
done

ls -lh ${SAMPLE}_counts_s*
grep -E "Status|Assigned" ${SAMPLE}_counts_s*.txt.summary
```

`ls` should list six files, for example for 1003_DRG: `1003_DRG_counts_s0.txt`, `1003_DRG_counts_s1.txt`, `1003_DRG_counts_s2.txt` and a `.summary` file for each. If they are missing, scroll up and read the featureCounts message. The usual reasons are a wrong `SAMPLE` name (the BAM file is not found) or the mapping job not having finished yet.

| Option | Meaning |
|---|---|
| `-p --countReadPairs` | the reads are paired; count each **pair** once |
| `-s` | strandedness: 0, 1 or 2 |
| `-t exon -g gene` | count reads on exons, and add up the exons of each gene; `gene` is the gene symbol, e.g. `EFNB3` |
| `-a` | the gene annotation |

Questions:

1. How many pairs were **assigned** to a gene with `-s 0`, `-s 1` and `-s 2`?
2. In a stranded library, one of `-s 1` or `-s 2` assigns about as many pairs as `-s 0`, and the other almost none. In an unstranded library, `-s 1` and `-s 2` each assign about half. Which is it here?
3. Look at the `.summary` file of your sample, for example `1003_DRG_counts_s0.txt.summary` (or the one that fits). Why are some pairs **Unassigned_NoFeatures**? Why are some **Unassigned_Ambiguity**?

Now we look at a few genes in your counts. We are still in `~/course/rnaseq/mapping` on the compute node. Column 7 is the count:

```shell
cd ~/course/rnaseq/mapping
grep -w -E "^(MOG|OLIG2|MPZ|PRPH|SCN10A|EFNB3)" ${SAMPLE}_counts_s0.txt | cut -f 1,7
exit
```

Question: Write the counts in the class table. Do spinal cord and DRG samples differ for these genes?

## Part 3: Differential expression in R

To compare the tissues properly, we use all 12 samples, mapped and counted in advance with the same steps as above. We analyse them with **DESeq2** in RStudio on your own computer.

### Install the packages (once)

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Cf.eg.db"))
install.packages(c("ggplot2", "pheatmap"))
```

`org.Cf.eg.db` is the gene database for the dog (*Canis familiaris*); we need it for the GO enrichment at the end.

### Read the data

Open a new R script in RStudio, paste the code step by step and run it.

```r
library(DESeq2)
library(ggplot2)
library(pheatmap)

url <- "https://raw.githubusercontent.com/vetsuisse-unibe/ngs_course_2026/main/data/"
counts  <- read.delim(paste0(url, "dog_SC_DRG_counts.tsv"), row.names = 1, check.names = FALSE)
samples <- read.delim(paste0(url, "dog_SC_DRG_samples.tsv"), row.names = 1)

samples$dog    <- factor(samples$dog)
samples$sex    <- factor(samples$sex)
samples$tissue <- factor(samples$tissue, levels = c("DRG", "SpinalCord"))

dim(counts)
head(counts)
samples
all(colnames(counts) == rownames(samples))
```

Questions:

1. How many genes and how many samples are there?
2. The gene names are gene **symbols** (like `MOG`), because the NCBI annotation uses them. Why is that convenient?
3. Why must `all(colnames(counts) == rownames(samples))` be `TRUE`?

### A paired design

Each dog gave **both** tissues. Dogs differ from each other (age, genetics, how the tissue was collected), and we don't want those differences to hide the tissue effect. So we tell DESeq2 that samples from the same dog belong together: the design `~ dog + tissue` means "compare the tissues **within** each dog".

```r
dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~ dog + tissue)

# keep genes with at least 10 reads in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
nrow(dds)

dds <- DESeq(dds)
res <- results(dds, contrast = c("tissue", "SpinalCord", "DRG"), alpha = 0.05)
summary(res)
```

Questions:

1. How many genes were removed by the filter? Why remove genes with almost no reads?
2. `log2FoldChange` > 0 means higher in which tissue?
3. How many genes are significantly higher in spinal cord, and how many in DRG (adjusted p-value < 0.05)? Is that more or less than we expected?

### Do the samples group by tissue? (PCA)

```r
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = c("tissue", "sex", "dog"), returnData = TRUE)
var <- round(100 * attr(pca, "percentVar"))

ggplot(pca, aes(PC1, PC2, colour = tissue, shape = sex, label = dog)) +
    geom_point(size = 4) +
    geom_text(vjust = -1, size = 3) +
    labs(x = paste0("PC1: ", var[1], "% variance"),
         y = paste0("PC2: ", var[2], "% variance")) +
    theme_bw()
```

Questions:

1. What separates the samples along PC1?
2. Do the two tissues of the same dog lie close together on PC2? Do males and females separate?

### Check the result with known markers

Before we trust a gene list, we check genes whose answer we know:

| Gene | Cell type | Expected in |
|---|---|---|
| `MOG`, `OLIG2` | oligodendrocytes (central nervous system myelin) | spinal cord |
| `MPZ` | Schwann cells (peripheral nervous system myelin) | DRG |
| `PRPH`, `SCN10A`, `POU4F1` | sensory neurons | DRG |
| `EFNB3` | the puppies' disease gene | ? |

```r
markers <- c("MOG", "OLIG2", "MPZ", "PRPH", "SCN10A", "POU4F1", "EFNB3")
markers <- markers[markers %in% rownames(dds)]   # genes removed by the filter are skipped

norm <- counts(dds, normalized = TRUE)[markers, ]
df <- data.frame(gene   = rep(markers, ncol(norm)),
                 sample = rep(colnames(norm), each = length(markers)),
                 count  = as.vector(norm))
df$tissue <- samples[df$sample, "tissue"]

ggplot(df, aes(tissue, count + 1, colour = tissue)) +
    geom_jitter(width = 0.1, size = 2) +
    scale_y_log10() +
    facet_wrap(~ gene, scales = "free_y") +
    labs(y = "normalised count + 1 (log scale)") +
    theme_bw()

res[markers, c("baseMean", "log2FoldChange", "padj")]
```

Questions:

1. Do the marker genes behave as expected? If not, what could be wrong?
2. Is *EFNB3* expressed in adult spinal cord? In DRG? Does that fit what we know about its role during development?

### Volcano plot and heatmap

```r
resdf <- as.data.frame(res)
resdf$gene <- rownames(resdf)
resdf$significant <- !is.na(resdf$padj) & resdf$padj < 0.05 & abs(resdf$log2FoldChange) > 1

ggplot(resdf, aes(log2FoldChange, -log10(padj), colour = significant)) +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_colour_manual(values = c("grey70", "firebrick")) +
    geom_text(data = subset(resdf, gene %in% markers), aes(label = gene), colour = "black", vjust = -0.8) +
    labs(x = "log2 fold change (spinal cord vs DRG)", y = "-log10 adjusted p-value") +
    theme_bw()

top <- head(rownames(res)[order(res$padj)], 40)
pheatmap(assay(vsd)[top, ], scale = "row",
         annotation_col = samples[, c("tissue", "sex")],
         show_rownames = TRUE, fontsize_row = 7)
```

Questions:

1. In the volcano plot, what does a point at the top right mean? At the top left?
2. In the heatmap, do the samples cluster by tissue? Pick three genes from the top of the heatmap and look them up in GeneCards (https://www.genecards.org). Do they make sense for the tissue where they are higher?

### Which biological processes differ? (GO enrichment)

A list of hundreds of genes is hard to read. **Gene Ontology (GO) enrichment** asks which biological processes are over-represented among them.

```r
library(clusterProfiler)
library(org.Cf.eg.db)

universe <- rownames(resdf)
up_sc  <- subset(resdf, padj < 0.05 & log2FoldChange >  2)$gene
up_drg <- subset(resdf, padj < 0.05 & log2FoldChange < -2)$gene

go_sc <- enrichGO(up_sc, OrgDb = org.Cf.eg.db, keyType = "SYMBOL", ont = "BP", universe = universe)
go_drg <- enrichGO(up_drg, OrgDb = org.Cf.eg.db, keyType = "SYMBOL", ont = "BP", universe = universe)

dotplot(go_sc, showCategory = 15) + ggtitle("Higher in spinal cord")
dotplot(go_drg, showCategory = 15) + ggtitle("Higher in DRG")
```

Questions:

1. Which processes come up for the spinal cord? Which for the DRG? Do they fit the biology at the top of this page?
2. Why do we give a `universe` (all tested genes) instead of all dog genes?

### Save your results

```r
resdf <- resdf[order(resdf$padj), ]
write.csv(resdf, "SpinalCord_vs_DRG_DESeq2.csv", row.names = FALSE)
```

**For the brave.** Compare males and females: change the design to `~ dog + tissue` → `~ sex + tissue` and get `results(dds, contrast = c("sex", "male", "female"))`. Why can't we use `~ dog + sex + tissue`? (Hint: every dog has only one sex.) How many genes differ between the sexes? The reference genome is from a female dog: why can't Y-chromosome genes show up?

## Part 4: Save your work with git

```shell
cd ~/course
git add rnaseq/scripts
git commit -m "Add RNA-seq mapping script"
```

We save the R script into the course repository too (for example `rnaseq/scripts/deseq2_SC_vs_DRG.R`), then commit and push it.

### What did we learn?

1. RNA-seq reads need a splice-aware aligner such as HISAT2.
2. featureCounts turns alignments into a count per gene; check the strandedness first.
3. A paired design (`~ dog + tissue`) compares tissues within each animal.
4. Check known marker genes before trusting a list of differentially expressed genes.
5. GO enrichment turns a gene list into biology.
