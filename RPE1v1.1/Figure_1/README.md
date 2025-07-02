# Figure 1 - Quality control of the RPE-1 genome.
Workflow followed to assess diploid assembly quality and correcteness.

## Prerequisites

- **Data**
  - PacBio HiFi reads (FASTQ)
  - Oxford Nanopore reads (FASTQ)
  - Hi-C reads (FASTQ)

- **Tools**
  - R (ggplot2, dplyr, readr and stringr packages)
  - Verkko
  - meryl
  - Merqury
  - CRAQ
  - Flagger
  - Bandage
  
## Genome assembly with *Verkko v1.4*

```
/lustre/home/evolpe/.conda/envs/lastverkko/bin/verkko -d /lustre/home/evolpe/PacBio/verkko_phased --slurm \
--hifi [RPE1_HiFi_1.fastq.gz] [RPE1_HiFi_2.fastq.gz] [RPE1_HiFi_3.fastq.gz] [RPE1_HiFi_4.fastq.gz] \
--nano [RPE1_ULK114_20221211.fastq.gz] [RPE1_ULK114_20221215.fastq.gz] [RPE1_ULK114_20221219.fastq.gz] [RPE1_ULK114_20230124.fastq.gz] \
       [RPE1_ULK114_20230303.fastq.gz] [RPE1_ULK114_20230307.fastq.gz] [RPE1_ULK114_20230309.fastq.gz] [RPE1_ULK114_20230313.fastq.gz] [RPE1_ULK114_20230728.fastq.gz] \
--hic1 [RPE1_HiC_R1.fastq.gz] \
--hic2 [RPE1_HiC_R2.fastq.gz]
```

This step produces the assembly graph file:

```
/lustre/home/evolpe/PacBio/verkko_phased/assembly.gfa
```

The graph can be visualized with Bandage:

```
Bandage GUI
# Open: /lustre/home/evolpe/PacBio/verkko_phased/assembly.gfa
```

This corresponds to **panel A** in Figure 1. 

## Genome quality evaluation with *Merqury*

#### Step 1: k-mer database generation with meryl (k=31)

```
meryl count k=31 output [HiFi_1_31.meryl] [RPE1_HiFi_1.fastq.gz]
meryl count k=31 output [HiFi_2_31.meryl] [RPE1_HiFi_2.fastq.gz]
meryl count k=31 output [HiFi_3_31.meryl] [RPE1_HiFi_3.fastq.gz]
meryl count k=31 output [HiFi_4_31.meryl] [RPE1_HiFi_4.fastq.gz]

# Merge multiple datasets
meryl union-sum output [readdb_31.meryl] [HiFi_1_31.meryl] [HiFi_2_31.meryl] [HiFi_3_31.meryl] [HiFi_4_31.meryl]
```

#### Step 2: Assembly quality assessment with Merqury

```
merqury.sh readdb_31.meryl [assembly.hap1.fa] [assembly.hap2.fa] [meryl_output_prefix]
```

This command allows to generate the assembly spectrum plots shown in **panel B** of Figure 1.

## Identification of errors with *CRAQ*

```
craq -mgs 1000 -q 10 -avgl 45 -avgs 60 -pl T -b T -D [output_directory] -g [reference_index] -sms [hifi_reads] -x map-hifi 
```

The BED files obtained using this tool were used for **panel C** of Figure 1.

## Genome quality evaluation with *Flagger*

#### Step 1: k-mer database generation with meryl (k=15)

```
meryl count k=15 output [merylDB] [assembly.fa]
meryl print greater-than distinct=0.9998 [merylDB] > [repetitive_k15.txt]
```

#### Step 2: alignment of HiFi reads to the genome assembly using winnowmap
```
winnowmap -W [repetitive_k15.txt] -ax map-pb -Y -L --eqx --cs -I8g [reference_index] [hifi_reads] | \
samtools view -hb | \
samtools sort -@ 128  > [reads.bam]
samtools index [reads.bam]
```

#### Step 3: evaluation of read depth and coverage along the genome

```
samtools depth -aa -Q 0 [reads.bam] > [reads.bam.depth]
depth2cov -d [reads.bam.depth] -f [assembly.fa.fai] -o [reads.bam.depth.cov]
cov2counts -i [reads.bam.depth.cov] -o [reads_alignments.counts]
```

#### Step 4: Coverage based evaluation to extract haploid, error, duplicated and collpased regions of the genome

```
fit_gmm.py --counts [reads_alignments.counts]  --output [reads_alignments.table]
find_blocks_from_table [reads_alignments.table]

```
This workflow outputs four BED files, each corresponding to different genomic region classifications:

- **Erroneous** – Regions with low read coverage, potentially unreliable.  
  *File:* `rpe1.v1.1.error.bed`

- **Collapsed** – Regions where multiple copies are collapsed into a single sequence.  
  *File:* `rpe1.v1.1.collapsed.bed`

- **Haploid** – Regions with expected coverage, representing correctly assembled haploid sequence.  
  *File:* `rpe1.v1.1.haploid.bed`

- **Duplicated** – No falsely duplicated regions were detected.  

Regions flagged as **Erroneous** or **Collapsed** can be considered as "low confidence".

The `flagger_barplot.R` script was used to generate **panel D** of Figure 1, based on the BED files produced in the previous steps.
