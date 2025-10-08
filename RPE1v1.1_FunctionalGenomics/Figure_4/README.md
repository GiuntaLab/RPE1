# RNA-seq analysis with matched and unmatched reference genomes

This folder contains a description of analyses performed on **bulk RNA-seq from RPE-1 cells** to investigate the **impact of the chosen reference genome** on downstream interpretation. Multiple reference genomes were used, including **matched** (RPE-1) and unmatched (hg38). This approach allowed us to isolate the effects of reference genome choice on apparent gene expression, as any differences observed would be solely due to mapping discrepancies rather than biological variation.

### Prerequisites:
- TrimGalore v0.6.10, FastQC v0.11.9 and MultiQC v1.24 for quality control and adapter trimming;
- STAR v2.7.11b for read alignments against the RPE-1 hap1, RPE-1 hap2, and hg38 separately;
- samtools v1.21 for data processing and filtering;
- featureCounts v2.0.1 for the quantification of aligned reads;
- edgeR v4.0.16 to perform the differential expression analysis on R. 

## Workflow for RNA-seq analysis

### 1) Quality control and adapter trimming
Once the publicly available datasets were downloaded locally (SRA accession numbers: SRR23924195, SRR23924196, and SRR23924197), TrimGalore, FastQC, and MultiQC were used with default parameters to assess the quality of the paired-end reads and to remove adapters.

```commandline
#TrimGalore with default parameters
trim_galore --paired --basename <fastq_name> --gzip <fastq_name>_R1.fastq.gz <fastq_name>_R2.fastq.gz -output_dir fastq_val
```

### 2) Read alignment against haploid genomes
Paired-end reads were aligned to the respective haploid reference genome using STAR. Each genome was indexed using the respective GTF file, which contains the same genes annotation (kept only common genes) for the three genomes.

```
#!/bin/sh

STAR --runThreadN 2 \
     --genomeDir genome_index_folder/ \
     --readFilesIn trimmed_reads/${sample}_val_1.fq.gz \
                   trimmed_reads/${sample}_val_2.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix aligned_reads/${sample}_${genome_name}_ \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes Standard
```

### 3) Quantification of aligned reads on annotated genes
Aligned reads were quantified at the gene level using featureCounts. The corresponding annotation file (.gtf) for each genome was provided to assign reads to annotated genes, generating a count matrix suitable for differential expression analysis.
```
#!/bin/sh

featureCounts -T 20 \
                -a ${genome_name}.gtf \
                -p -o feature_counts/${genome_name}_feature_counts.txt aligned_reads/*bam
```


### 4) Differential expression analysis
Differential expression analysis was performed in R using the custom script *dge_analysis.R* and the edgeR package (v4.0.16). \
This script imports raw count data, normalizes expression values, and identifies significantly differentially expressed genes between conditions (triplicates of the same condition mapped to 3 different genomes). \
Standard visualization and statistical outputs were generated to evaluate the impact of the chosen reference genome on RNA-seq interpretation.