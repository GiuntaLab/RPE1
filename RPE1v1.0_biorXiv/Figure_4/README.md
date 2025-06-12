# Figure 4 - The isogenomic reference genome identifies high-confidence epigenomic landscapes
This folder contains all the workflow used to perform the epigenetic analysis.  

## Prerequisites:
- SRA Toolkit, samtools, awk
- Bowtie2
- MACS3
- dorado
- modkit
- karyoploteR
- regioneR
  
## Chromatin phased profiles 
CUT&RUN reads from the publicly available dataset of hTERT RPE-1 cells [GSE132193](https://doi.org/10.15252/embj.2019102924) were used to assess the epigenetic landscape of human centromeres comparing the results within our *de novo* phased genome assembly of the RPE-1 cells with the HG002v1.0, the CHM13v2.0, and GRCH38.p14.

### Read mapping
We performed the diploid read mapping for the RPE-1 and HG002 genomes while the haploid read mapping for the CHM13 and HG38 genomes. We used the following command from [Bowtie2](https://github.com/BenLangmead/bowtie2):

```
bowtie2 --end-to-end -x [index_reference] [read1.fastq] [read2.fastq]
```
The following rules were used to distinguish the enrichment profile and the high-confidence peaks:

- Enrichment profile: we removed FLAG 2308 and we retained peaks with q-value ≤0.00001;
- High-confidence profile: we removed FLAG 2308, filtered by MAPQ >20, and we retained peaks with q-value ≤0.00001;

### Peak calling
CENP-A peaks were determined using [MACS3](https://github.com/macs3-project/MACS), calculating the ratio between immunoprecipitated samples and background with these parameters:

```
macs3 callpeak -t [read_IP] -c [read_INPUT] -f BAMPE -B -g 3.03e9 -q 0.00001 -n CHM13_CENPA #haploid mapping
macs3 callpeak -t [read_IP] -c [read_INPUT] -f BAMPE -B -g 6.06e9 -q 0.00001 -n RPE1v1.0dip_CENPA #diploid mapping
```
Among all the outputs of this command, we processed the bedGraph files filtering by q-value ≤0.001 (Supplementary Fig. 12), q-value ≤0.00001 (Fig. 4a, b (small panel n.2), c), and q-value ≤0.0000001 (Fig. 4b (small panel n.3-4)).
The bedGraph files were used in karyoploteR to plot each chromatin density profile.

## Methylation phased profiles using *Dorado* and *Modkit*
Using [Dorado](https://github.com/nanoporetech/dorado), we downloaded the simplex basecalling model (dna_r10.4.1_e8.2_400bps_sup@v4.2.0) and we called the methylation profile (5mC) on Oxford Nanopore reads. We processed the output using [Modkit](https://github.com/nanoporetech/modkit).

### Read mapping and modified-base calling
Dorado outputs directly the modified 5mC in the SAM/BAM files mapping in parallel all the ONT reads against the diploid RPE1v1.0 reference genome (Hap1 + Hap2). We used the following command:

```
/data/Packages/dorado-0.4.3-linux-x64/bin/./dorado basecaller -v
/data/Packages/dorado-0.4.3-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v4.2.0
/data/FAST5_POD5/
--reference rpe1.v1.0.fasta
--modified-bases 5mC
--modified-bases-threshold 0.08
-c 43000
--secondary yes
-N 5 > rpe1.v1.0.fasta_5mC.bam
```

### Data processing
The output is a modBAM file with the evaluation of each 5-methylcytosine (5mC) performed using raw data (POD5) and the matched reference (rpe1.v1.0.fasta). After sorting the modBAM file, we used Modkit to convert the modBAM to bedMethyl file using the following command: 

```
modkit pileup rpe1.v1.0.fasta_5mC.bam.sort.bam
RPE1.v1.0.modkit.bed
--ref /rpe1.v1.0.fasta
--cpg
--bedgraph
--prefix RPE1.v1.0
--filter-threshold 0.80
```

Then, we used *awk* to filter the methylation profile in the plus strand and the modified fraction (N<sub>mod</sub> / N<sub>valid_cov</sub>) over 60. 

### Graphical representation
R scripts and example files for Figure 4 panel C (i.e. chromosome 9).
