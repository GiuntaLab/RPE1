# Giunta Lab Human Diploid Cell Lines Project
We have sequenced experimentally-ameanable human diploid laboratory cell lines and generated complete phased assemblies to use as matched reference genomes to analyze sequencing data generated from the same cell line, an approach we refer to as *isogenomic reference genome*. The improvement in alignment quality using matched reads-reference enables high-precision mapping for profiling phased epigenome and methylome. This proof-of-concept calls for a comprehensive catalog of complete assemblies for commonly used cells for a widespread application of isogenomic reference genomes to enable high-precision multi-omics analyses.
## Latest assembly release
### RPE1v1.0
#### The complete diploid reference genome of RPE-1 identifies human phased epigenetic landscapes
Volpe et al., [preprint](https://pubmed.ncbi.nlm.nih.gov/38168337/)

#### Sequencing data
We generated the complete phased genome assembly of one of the most widely used non-cancer cell lines (RPE-1) with a stable diploid karyotype. We produced state-of-the-art sequencing data using third generation sequencing, Pacific Biosciences (PacBio) high-fidelity (HiFi) and Oxford Nanopore Technologies (ONT). ONT long and ultra-long (UL) reads were exclusively generated using R10.14 (pore chemistry V14 yielding 99.9% base accuracy). We aimed at above-average reads depth coverage with 46x sequence coverage of HiFi reads, 80x of ONT, 30x ONT-UL reads (>100 kb). During DNA extraction, we assessed the length of high molecular weight DNA from hTERT RPE-1 cells using Femto Pulse with a yield of native DNA size distribution around 116 kb (main peak) and up to 1 Mb in length (smaller peaks) with centrifugation below 1000 RPM. For NGS, we generated 100x Illumina, 60x PRC-free Illumina and 30x reads for the Hi-C data (Arima Genomics) for haplotypes phasing.

Files:
- PacBio HiFi raw reads (fastq)
- ONT reads (fastq)
- Illumina (fastq)
- Hi-C Arima reads (fastq)

#### Assemblers
We tested the capabilities of two genome assemblers, [Hifiasm](https://github.com/chhylp123/hifiasm) (v. 0.19.8-r603) and [Verkko](https://github.com/marbl/verkko) (v.1.4) under the following specification: 1 processor, 128 threads and 2 GB of memory per thread, for a total of 100 hours jobs launched on the [Sapienza TERASTAT2 server](https://www.dss.uniroma1.it/it/HPCTerastat2). These resources were managed by the Slurm system, and various combinations of threads and memory allocations were tested to achieve the optimal balance between time and memory efficiency. Following the integration of Hi-C data, we proceeded using Verkko for the final assembly (Supplementary Note 2). 

Files:
- Verkko Unphased (GFA) with tangled t(10q:Xq) visualized with [Bandage](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4595904/)
- Partially phased (fasta)

#### Phasing of haplotypes
In absence of parental information to support [Trio-binning](https://www.nature.com/articles/nbt.4277), we obtained fully phased haplotypes for the RPE-1 genome using contact maps. We used the Vertebrate Genome Project (VGP) [pipeline](https://github.com/VGP/vgp-assembly). Hi-C raw reads were aligned against the merged assembly composed of both RPE-1 haplotypes and unassigned reads generated from Verkko. The alignment step was performed using the short-read aligner [BWA](https://github.com/lh3/bwa). All reads were retained, including those classified as supplementary, with low mapping quality or having multiple alignments. The final aligned file was converted into a 3D Contact Map viewable through [PretextView](https://github.com/wtsi-hpag/PretextView). PretextView allows the modification of the contact map file by changing contig positions along the diagonal and finding the correct chromatin interaction path. The final diploid Hi-C contact map was based on Nadolina Brajuka [RapidCuration2.0](https://github.com/Nadolina/Rapid-curation-2.0). The unassigned reads over 300 kb were assigned to the chromosome merged to the assembly file, yet a remaining 667 contigs <100 kb could not be aligned manually due to short size. In the final contact map, each scaffold was assigned to a specific haplotype Meta Data Tag (Hap 1 and Hap 2). This latest contact map was then converted into an A Golden Path (AGP) file, which was used as input for the subsequent separation of both haplotypes and generating the two fully phased haploid FASTA files after running [Curation2.0_pipe.sh](https://github.com/Nadolina/Rapid-curation-2.0/blob/main/curation_2.0_pipe.sh). RPE1v1.0 base accuracy quality score (Phred) was QV 64.1 for Hap 1 and QV 61.8 for Hap 2.

Files:
- Contact map (Pretext view BAM file)
- AGP
- RPE1v1.0 assembly (diploid fasta)

#### Genome Quality
- [CRAQ](https://github.com/JiaoLaboratory/CRAQ) (Clipping Reveals Assembly Quality):
We mapped PacBio HiFi and Illumina reads to evaluate the quality of the genome using CRAQ. Files obtained:
  - out_regional.Report
  - locER_out/out_final.CRH.bed
  - locER_out/out_final.CRE.bed
  - strER_out (and subfolders were empty due to lack of structural variants (SVs) found in the RPE1v1.0 assembly)
- [Flagger](https://github.com/mobinasri/flagger):
We used PacBio HiFi reads to detect misassemblies. Files obtained:
  - rpe1.error.bed
  - rpe1.collapsed.bed
  - rpe1.haploid.bed
  - rpe1.diplicated.bed (file empty as no duplicated regions were found in the RPE1v1.0 assembly)

#### Genome-to-genome comparison
- [SyRI](https://github.com/schneebergerlab/syri) (Synteny and Rearrangement Identifier):
We used Minimap2.0 for genome-to-genome alignment to obtain the SAM files. SyRI gave the following output files:
  - chm13.hap1.out
  - chm13.hap2.out
  - hap1.hap2.out
  - syri.summary
  - VCF
where the .out includes not aligned, SVs, SNPs, *invOut* for inversions; *TLOut* for translocations; *invTL* for inverted translocations; *invDupOut* for inverted duplications; *dupOut* for duplications; *ctxOut* for cross-chromosomal translocations; *synOut* for synthenic regions.

#### Identification and curation of RPE-1 specific structural variants
Multi-step pipeline for the manual curation of the RPE-1 specific structural variant identified as 46,XX,dup(10q),t(Xq;10q),del(Xq28)
1)	Reads alignment: The *de novo* diploid genome assembly RPE1v1.0 was used as a reference to map HiFi and ONT reads with Minimap2.0.
2)	Visualizaion: Long-read alignments were visualized on IGV, revealing an increase in reads coverage on chromosome 10q (long arm). Reads interruption was found, with ~100 bp difference in mapping position between chromosome 10 of Hap 1 and 2.
3)	Read alignment quality between haplotypes: Chromosome 10 Hap 1 showed only reads with mapping quality of 0, while chromosome 10 Hap 2 showed reads with mapping quality between 10-60 and supplementary alignments in the telomeric region of chromosome X Hap 1, corresponding to the translocation breakpoint.
4)	Manual curation (translocation): The sequence of the duplicated and translocated long arm of chromosome 10 Hap 2 was added to the telomeric region of chromosome X Hap 1. This addition was done merging the two previously mentioned sequences into the existing FASTA file of the RPE-1 assembly.
5)	Read alignment verification: Minimap2.0 was used to align RPE-1 HiFi and ONT reads against the modified FASTA file in the telomeric region of chromosome X Hap 1. The IGV visualization of the fusion point position on chromosome X revealed a microdeletion of 3,603 bp in read alignment, suggesting the loss of these bases during the rearrangement between chromosome X and 10.
6)	Manual curation (deletion): The bases were deleted from the FASTA file of chromosome X and the modified genome was aligned against the RPE-1 HiFi and ONT reads.
7)	Final verification: The RPE1.v1.0 diploid genome shows reads that are completely aligned to the fusion point on chromosome X.

Files:
- BAM dup(10q) Hap 1 - ONT
- BAM dup(10q) Hap 2 - ONT
- BAM dup(10q) Hap 1 - HiFi
- BAM dup(10q) Hap 2 - HiFi
- BAM t(Xq;10q) - ONT
- BAM t(Xq;10q),del(Xq28) - ONT

#### Repeat annotation
We obtained RPE-1 monomers, monomers organization and genome-wide alpha-satellite DNA annotation using:  
- [RepeatMasker](https://www.repeatmasker.org/)
- [HumAS-HMMER_for_AnVIL](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL?tab=readme-ov-file)
RPE-1 HOR Structural Variant (StV) prediction using HOR-monomer annotation using:  
- [StV](https://github.com/fedorrik/stv)
Sequence identity heatmap were obtined using StainedGlass:
- [StainedGlass0.5](https://github.com/mrvollger/StainedGlass)

#### Alignment
Reads alignment was performed on diploid RPE1v1.0 or CHM13v2.0 using [NucFreq](https://github.com/mrvollger/NucFreq) v0.1 and visualized with Nucplot.py. Alignments were done whole-genome. NM and mapQ were extracted from the RPE1v1.0 Hap1 BAM, RPE1v1.0 Hap2 BAM, and CHM13 BAM whole-genome or from SyRI coordinates of highly-diverged regions (HDR).  

*Information related to [Figure 3](https://www.biorxiv.org/content/10.1101/2023.11.01.565049v2.full.pdf+html)*. [See linked scripts]

#### Epigenetic analysis
Centromere chromatin phased landscapes were obtained from RPE-1 CUT&RUN CENP-A [dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132193) (GSE132193). Comparison in short reads alignment mapping using the following reference genomes:
- RPE1v1.0 Hap 1
- RPE1v1.0 Hap 2
- CHM13v2.0
- GRCh38.p14
- HG002v1.0 maternal
- HG002v1.0 paternal
Data were aligned using [Bowtie](https://github.com/BenLangmead/bowtie) and CENP-A high-confidence peaks were determined with [MACS3](https://github.com/macs3-project/MACS).

Files:
- rpe1.v1.0.bedgraph
- CHM13v2.0.bedgraph
- GRCh38.p14.bedgraph
- HG002v1.0.bedgraph

#### Methylation profiles
Methylation profiles for 5-methylcytosine (5mc) were generated from the ONT RPE-1 POD5 (3.6 TB) using [Dorado v4.2.0 basecalling model](https://github.com/nanoporetech/dorado/releases) and the output processed with Modkit. Following the evaluation of reads coverage values N<sub>canonical</sub> and N<sub>mod</sub>, we selected the filter (N<sub>mod</sub> / N<sub>valid_cov</sub>) >60 applied to the bedMethyl output. 

Files:
- Dorado output: BAM cytosines with tags MM (5mC) and ML (probability score)
- Modkit output: bedMethyl with position, score and coverage for Hap 1
- Modkit output: bedMethyl with position, score and coverage for Hap 2
  
*Information related to [Figure 4](https://www.biorxiv.org/content/10.1101/2023.11.01.565049v2.full.pdf+html)*. [See linked scripts]
  





