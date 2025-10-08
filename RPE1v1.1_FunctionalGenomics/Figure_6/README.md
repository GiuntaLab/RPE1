# Centromeric proteins enrichment at centromeres

This folder contains a description of analyses performed on **CUT&RUN and ChIP-seq datasets** to investigate protein enrichment and spread across centromeric arrays, as well as to study the **impact of the chosen reference genome** on downstream interpretation. Multiple reference genomes were used, including **matched** (RPE-1) and unmatched (CHM13, hg38, HG002, and YAO).

### Prerequisites:
- TrimGalore v0.6.10, FastQC v0.11.9 and MultiQC v1.24 for quality control and adapter trimming;
- bowtie2 v2.4.4 for diploid (RPE-1, YAO, and HG002) and haploid (CHM13, hg38) read alignments;
- samtools v1.21 for data processing and filtering;
- MACS3 v3.0.1 for peak calling;
- karyoploteR v1.34.2 and ComplexHeatmap v2.24 for data visualization (R packages).

## Workflow for centromeric peak calling

### 1) Quality control and adapter trimming
Once the publicly available datasets were downloaded locally, TrimGalore, FastQC, and MultiQC were used with default parameters to assess the quality of the paired-end reads and to remove adapters.

```commandline
#TrimGalore with default parameters
trim_galore --paired --basename <fastq_name> --gzip <fastq_name>_R1.fastq.gz <fastq_name>_R2.fastq.gz -output_dir fastq_val
```

### 2) Read alignment against haploid and diploid reference genomes
Read alignment was performed against multiple chromosome-level genomes: the RPE1v1.1 phased haplotypes (diploid mapping), T2T-CHM13v2.0 (haploid mapping), hg38 (GRCh38.p14, haploid mapping), the HG002v1.0.1 paternal and maternal haplotypes (diploid mapping) and the YAO paternal and maternal haplotypes (diploid mapping).

```
#!/bin/sh

Set='--end-to-end --sensitive --no-mixed --threads 8'

#Generating the BAM file for each sample and genome discarding unmapped reads
bowtie2 -x <reference_genomes> -1 fastq_val/<fastq_name>_val_1.fq.gz -2 fastq_val/<fastq_name>_val_2.fq.gz $Set | \
    samtools view -b -F 12 - | \
    samtools sort -@ 8 - -O BAM -o bam/<fastq_name>.<reference_genomes>.bam

#Generating BAM files removing alignments with MAPQ < 20
samtools view -bh -q 20 bam/<fastq_name>.<reference_genomes>.bam | \ 
    samtools sort -@ 8 - -O BAM -o bam/<fastq_name>.<reference_genomes>.q20.bam

#Indexing bam files
samtools index -@ 8 bam/<fastq_name>.<reference_genomes>.bam
samtools index -@ 8 bam/<fastq_name>.<reference_genomes>.q20.bam
```

Next, samtools tool was used to merge replicates of the same sample into a single file:
```commandline
samtools merge bam_merge/<output_merged_files>.bam <sample1_file1>.bam <sample1_file2>.bam <sample1_file3>.bam
```

### 3) Peak calling 
Unfiltered and filtered alignment files were used to perform the peak calling step using MACS3. The tool was used to call broad peaks by comparing immunoprecipitated samples against the corresponding controls with default parameters. 

```commandline
#!/bin/sh
#CALLPEAK
#Diploid
macs3 callpeak -t IP_<fastq_name>.<reference_genomes>.bam -c INPUT_<fastq_name>.<reference_genomes>.bam -f BAMPE -B -g 6.00e9 -q 0.00001 -n macs3_peaks/<output_file_name> --broad --tempdir tmp_dir/ 
macs3 callpeak -t IP_<fastq_name>.<reference_genomes>.q20.bam -c INPUT_<fastq_name>.<reference_genomes>.q20.bam -f BAMPE -B -g 6.00e9 -q 0.00001 -n macs3_peaks/<output_file_name> --broad --tempdir tmp_dir/

#Haploid
macs3 callpeak -t IP_<fastq_name>.<reference_genomes>.bam -c INPUT_<fastq_name>.<reference_genomes>.bam -f BAMPE -B -g 3.00e9 -q 0.00001 -n macs3_peaks/<output_file_name> --broad --tempdir tmp_dir/ 
macs3 callpeak -t IP_<fastq_name>.<reference_genomes>.q20.bam -c INPUT_<fastq_name>.<reference_genomes>.q20.bam -f BAMPE -B -g 3.00e9 -q 0.00001 -n macs3_peaks/<output_file_name> --broad --tempdir tmp_dir/
```

This command will output a plain text containing peak data, two peak files (broadPeak and gappedPeak), and two bedraph files (treat and control). To obtain an enrichment signal for downstream visualization, we used the MACS3 *bdgcmp* module to compare the *_treat_pileup.bdg and the *_control_lambda.bdg tracks produced by MACS3 callpeak. The -m parameter was set to qpois to calculate, for each genomic bin, the log10 of the Poisson q-value.

```
#!/bin/sh
#BedGraph comparison
macs3 bdgcmp -t <output_file_name>_treat_pileup.bdg -c <output_file_name>_control_lambda.bdg -m qpois -o bdg_compare/<output_file_name>_comp_qval.bdg
```

The karyoploteR v1.34.2 R package was used to visualize protein enrichment tracks across all chromosomes and genomes with the kpBars function. The code is stored in the *plot_bdg_signals.R* file which contains an example for plotting **Fig. 6c, d**.

## Evaluation of the haplotype-specificity of CENP-A peaks
To assess the haplotype-specificity and uniqueness of the genomic region containing CENP-A broad peaks identified on the RPE1v1.1 genome using filtered alignment files, we divided the LHOR regions into non-overlapping 500 nucleotide bins using the GenomicRanges v1.60 R package.
For each haplotype, two FASTA files were created using bedtools *getfasta* and *intersect* modules, one containing bins overlapping with at least one peak and the other the remaining bins. The resulting FASTA files were aligned to the other haplotype reference sequence to evaluate sequence similarity using minimap2 with parameters *-x asm20 -c --cs=long*, generating PAF files, and bowtie2 v2.4.4 with parameters *-end-to-end --sensitive* generating BAM files.
- For the alignments produced by minimap2, **percent identity** was calculated as (Number of matching bases / Alignment block length)×100. 
- For the alignments produced by bowtie2, **percent identity** was calculated as (Alignmentblock length - NM value)/ Alignment block length)×100.

BAM files produced using bowtie2 were first postprocessed using the script *p_identity_bam.py* to obtain the percentage identity of each mapped read. The resulting TSV files were then used to plot the values shown in **Supplementary Fig. S16a**.\
PAF files were entirely postprocessed using the script *plot_paf.py* both for computing the percentage of identity and for plotting the results shown in **Supplementary Fig. 16b**.  

## Evaluation of CDR occupancy
The coordinates of the broad peaks identified by MACS3 were intersected with those of the centromere dip regions (CDRs) using the bedtools *intersect* module. For each CDR, the average of -log10(q-value) was calculated from the broadPeak file for each centromeric protein and it was plotted using the ComplexHeatmap v2.24 R package (**Supplementary Fig. 8c, d**). 
The custom code *cdr_occupancy.R* was used to generate the two heatmaps.
