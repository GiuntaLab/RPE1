# Giunta Lab Human Diploid Cell Lines Project
We have sequenced experimentally-ameanable human diploid laboratory cell lines and generated complete phased assemblies to use as matched reference genomes for sequencing data generated from the same cell line, an approach we refer to as *isogenomic reference genome*. This proof-of-concept calls for a comprehensive catalog of complete genome assemblies for commonly used cells for a widespread application of isogenomic reference genomes to enable high-precision multi-omics analyses.
## Latest assembly release
### RPE1v1.0
#### Sequencing data
We generated the complete phased genome assembly of one of the most widely used non-cancer cell lines (RPE-1) with a stable diploid karyotype. We produced state-of-the-art sequencing data using third generation sequencing, Pacific Biosciences (PacBio) high-fidelity (HiFi) and Oxford Nanopore Technologies (ONT). ONT long and ultra-long (UL; > 100 kb) reads were exclusively generated using the pore chemistry R10.14 (V14). We aimed at above-average reads depth coverage with 46x sequence coverage of HiFi reads, 80x of ONT, 30x ONT-UL reads over 100 kb with a QV score >20. During DNA extraction, we assessed using Femto Pulse the legnth of high molecular weight DNA from hTERT RPE-1 cells with a yield of native DNA size distribution around 116 kb (main peak) and up to 1 Mb in length (smaller peaks) when centrifugation was kept below 1000 RPM. For NGS, we generated 100x Illumina, 60x PRC-free Illumina and 30x reads for the Hi-C data (Arima Genomics) - see "Phasing of haplotypes".
Files:
- PacBio HiFi raw reads (fastq)
- ONT reads (fastq)
- Illumina (fastq)
- Hi-C Arima reads (fastq)

#### Assemblers
We tested the capabilities of two genome assemblers, Hifiasm (v. 0.19.8-r603) and Verkko (v.1.4). We utilized both tools under the following specifications: 1 processor, 128 threads and 2 GB of memory per thread, for a total of 100 hours launching the jobs on the Sapienza TERASTAT2 server. These resources were managed by the Slurm system, and various combinations of threads and memory allocations were tested to achieve the optimal balance between time and memory efficiency. Following the integration of Hi-C data, we proceeded using Verkko (See 'Supplementary Note 2). 
Files:
- Verkko Unphased (GFA) with tangled t(10q:Xq) visualized with [Bandage](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4595904/)
- Partially phased (fasta)

#### Phasing of haplotypes
In absence of parental information to support [Trio-binning](https://www.nature.com/articles/nbt.4277), we obtained fully phased haplotypes for the RPE-1 genome using contact maps. We used the Vertebrate Genome Project (VGP) [pipeline](https://github.com/VGP/vgp-assembly). Briefly, Hi-C raw reads were aligned against the merged assembly composed of both RPE-1 haplotypes and unassigned reads generated from Verkko. The alignment step was performed using the short-read aligner [BWA](https://github.com/lh3/bwa). All reads were retained, including those classified as supplementary, with low mapping quality or having multiple alignments. The final aligned file was converted into a 3D Contact Map viewable through [PretextView](https://github.com/wtsi-hpag/PretextView). PretextView allows the modification of the contact map file by changing contig positions along the diagonal and finding the correct chromatin interaction path, thus generating a final diploid Hi-C contact map based on Nadolina Brajuka [RapidCuration2.0](https://github.com/Nadolina/Rapid-curation-2.0). The unassigned reads over 300 kb were assigned to the chromosome merged to the assembly file, yet a remaining 667 contigs <100 kb could not be aligned manually because of their short size. In the final contact map, each scaffold was assigned to a specific haplotype Meta Data Tag (Hap 1 and Hap 2). This latest contact map was then converted into an A Golden Path (AGP) file, which was used as input for the subsequent separation of both haplotypes and generating the two fully phased haploid FASTA files after running [Curation2.0_pipe.sh](https://github.com/Nadolina/Rapid-curation-2.0/blob/main/curation_2.0_pipe.sh).
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
were the .out includes not aligned, SVs, SNPs, *invOut* for inversions; *TLOut* for translocations; *invTL* for inverted translocations; *invDupOut* for inverted duplications; *dupOut* for duplications; *ctxOut* for cross-chromosomal translocations; *synOut* for synthenic regions.
  - VCF

#### Structural rearrangements
Multi-step pipeline to identify the RPE-1 specific structural variant t(Xq;dup10q)
Resolution of the karyotypically stable translocation of the RPE-1 genome, which involved chromosome 10 and chromosome X, was conducted through a multi-step pipeline we developed:
1)	The *de novo* diploid genome assembly RPE1v1.0 was used as a reference to map HiFi and ONT reads with Minimap2.
2)	Long-read alignments were visualized on IGV, revealing an interruption in read coverage on chromosome 10. The interruption position was slightly different between chromosome 10 of Hap 1 and 2, consistently occurring on the long arm.
3)	Read alignment at this point showed different features. Chromosome 10 Hap 1 showed only reads with mapping quality of 0, while chromosome 10 Hap 2 showed reads with mapping quality between 10-60 and supplementary alignments in the telomeric region of chromosome X Hap 1, corresponding to the breakpoint position of chromosome 10 Hap 2.
4)	The sequence of the duplicated and translocated long arm of chromosome 10 Hap 2 was added to the telomeric region of chromosome X Hap 1. This addition was done merging the two previously mentioned sequences into the existing FASTA file of the RPE-1 assembly.
5)	Minimap2 was used to align RPE-1 HiFi and ONT reads against the modified FASTA file in the telomeric region of chromosome X Hap 1. The IGV visualization of the fusion point position on chromosome X revealed a microdeletion of 3,603 bp in read alignment, suggesting the loss of these bases during the rearrangement between chromosome X and 10.
6)	These bases were deleted from the FASTA file of chromosome X, and the modified genome was aligned against the RPE-1 HiFi and ONT reads.
7)	The RPE1.v1.0 diploid genome shows reads that are completely aligned to the fusion point on chromosome X.
Files:

#### Repeat annotation
- RepeatMasker [https://www.repeatmasker.org/]
- HumAS-HMMER_for_AnVIL [https://github.com/fedorrik/HumAS-HMMER_for_AnVIL?tab=readme-ov-file]
...    
- HOR Structural Variant (StV) prediction using HOR-monomer annotation [https://github.com/fedorrik/stv]
   
#### Epigenetic analysis
We used publicly available CUT&RUN CENP-A dataset from RPE-1 cells (GSE132193) and compared the alignments of RPE-1 reads to different reference genomes:
- RPE1v1.0 Hap 1
- RPE1v1.0 Hap 2
- CHM13v2.0
- GRCh38.p14
- HG002v1.0 maternal
- HG002v1.0 paternal
Data were aligned using Bowtie [https://github.com/BenLangmead/bowtie] and CENP-A peaks were determined with MACS3 [https://github.com/macs3-project/MACS]. We obtained:
-rpe1.v1.0.bedgraph
-CHM13v2.0.bedgraph
-GRCh38.p14.bedgraph
-HG002v1.0.bedgraph

- RPE-1 5mC methylation profile from ONT reads 
 
 
  





