# The reference genome of the human diploid cell line RPE-1
We have sequenced the diploid genome of the human diploid laboratory cell line RPE-1. This near complete, fully phased can be used as matched reference genomes to analyze sequencing data generated from the same cell line, an approach we refer to as *isogenomic reference genome*. The improvement in alignment quality using matched reads-reference enables high-precision mapping for profiling phased epigenome and methylome. This proof-of-concept calls for a comprehensive catalog of complete assemblies for commonly used cells for a widespread application of isogenomic reference genomes to enable high-precision multi-omics analyses.

## RPE1 Genome Assembly Versions

### RPE1v1.1  
This is the latest release, incorporating manual curation of gaps present in RPE1v1.0 and the recovery of previously unassigned telomeric sequences.  
**The reference genome of the human diploid cell line RPE-1**  
Volpe, Colantoni et al., *Manuscript Under Revision*  
Documentation, scripts, files and figures relative to the main figures are available in the `RPE1v1.1/` folder. Supplemental information is available in the **Zenodo Repository**.

#### Genome availabiliy
The RPE1v1.1 genome has been deposited in NCBI GenBank under the accession numbers [JBJQNK000000000] (https://www.ncbi.nlm.nih.gov/nuccore/JBJQNK000000000) (Hap1) and [JBJQNL000000000](https://www.ncbi.nlm.nih.gov/nuccore/JBJQNL000000000) (Hap2), with links to BioProject accession numbers [PRJNA1193286](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1193286/) (Hap1) and [PRJNA1193302](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1193302/) (Hap2), both under the umbrella BioProject accession number [PRJNA1195024](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1195024/).

### RPE1v1.0  
**The complete diploid reference genome of RPE-1 identifies human phased epigenetic landscapes**  
Volpe et al., [preprint](https://pubmed.ncbi.nlm.nih.gov/38168337/)
Documentation, scripts, files and figures relative to the main figures are available in the `RPE1v1.0_biorXiv/` folder. 

## Sequencing data
We generated the near complete phased genome assembly of one of the most widely used non-cancer cell lines (RPE-1) with a stable diploid karyotype. We produced state-of-the-art sequencing data using third generation sequencing: Pacific Biosciences (PacBio) high-fidelity (HiFi) and Oxford Nanopore Technologies (ONT). ONT long and ultra-long (UL, >100 kb) reads were exclusively generated using R10.14 (pore chemistry V14 yielding 99.9% base accuracy). We aimed at above-average reads depth coverage with 46x sequence coverage of HiFi reads, 95x of ONT long, and 30x ONT-UL reads. During DNA extraction, we assessed the length of high molecular weight DNA from hTERT RPE-1 cells using Femto Pulse with a yield of native DNA size distribution around 116 kb (main peak) and up to 1 Mb in length (smaller peaks) with centrifugation below 1000 RPM. For haplotype phasing, we generated 30x Hi-C reads (Arima Genomics).

#### Sequencing data availability

The HiFi, ONT, and Hi-C sequencing data used to generate the genome have been deposited in the SRA under BioProject [PRJNA1193286](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1193286/). Accession numbers are as follows:

- **PacBio HiFi raw reads**:  
  [SRR33464826](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464826), [SRR33464827](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464827),
[SRR33464828](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464828),
[SRR33464829](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464829)

- **ONT reads**:  
  [SRR33464817](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464817), [SRR33464818](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464818), [SRR33464819](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464819), [SRR33464820](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464820), [SRR33464821](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464821), [SRR33464822](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464822), [SRR33464823](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464823), [SRR33464824](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464824), [SRR33464830](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464830), [SRR33464831](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464831)

- **Hi-C Arima reads**:  
  [SRR33464825](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33464825)


## Genome Assembly procedure

### Automated genome assembly  
The first draft of the genome was generated using [Verkko](https://github.com/marbl/verkko) v1.4 with HiFi, ONT and Hi-C reads. Verkko was run on the [Sapienza TERASTAT2 cluster](https://www.dss.uniroma1.it/it/HPCTerastat2).  

### Dual manual curation  
In absence of parental information to support [Trio-binning](https://www.nature.com/articles/nbt.4277), we obtained fully phased haplotypes for the RPE-1 genome using contact maps. We used the Vertebrate Genome Project (VGP) [pipeline](https://github.com/VGP/vgp-assembly). Hi-C raw reads were aligned against the merged assembly composed of both RPE-1 haplotypes and unassigned reads generated from Verkko. The alignment step was performed using the short-read aligner [BWA](https://github.com/lh3/bwa). All reads were retained, including those classified as supplementary, with low mapping quality or having multiple alignments. The final aligned file was converted into a 3D Contact Map viewable through [PretextView](https://github.com/wtsi-hpag/PretextView). PretextView allows the modification of the contact map file by changing contig positions along the diagonal and finding the correct chromatin interaction path. The final diploid Hi-C contact map was based on Nadolina Brajuka [RapidCuration2.0](https://github.com/Nadolina/Rapid-curation-2.0). Contigs over 300 kb were assigned to the chromosome merged to the assembly file, yet a remaining 667 contigs <100 kb could not be aligned manually due to short size. In the final contact map, each scaffold was assigned to a specific haplotype Meta Data Tag (Hap 1 and Hap 2). This latest contact map was then converted into an A Golden Path (AGP) file, which was used as input for the subsequent separation of both haplotypes and generating the two fully phased haploid FASTA files after running [Curation2.0_pipe.sh](https://github.com/Nadolina/Rapid-curation-2.0/blob/main/curation_2.0_pipe.sh). 

## Further Manual Curation: v1.0 → v1.1

To address gaps remaining after dual manual curation, we applied two complementary gap-closing strategies:

- **Long-read alignment-based gap filling**  
  ONT reads >100 kb were aligned to the assembly using [minimap2](https://github.com/lh3/minimap2).  
  Gaps were identified when reads spanned at least 40 kb on both flanking regions. These gaps were temporarily filled with ONT read sequences.  
  HiFi reads were then aligned to the patched genome, and final gap-filling was performed using the consensus sequence from these alignments.

- **Assembly graph-based gap filling**  
  ONT reads were aligned to the assembly graph using [gfalign](https://github.com/vgl-hub/gfalign).  
  When alignments supported a single unambiguous path through a graph bubble, the corresponding sequence of unitigs was extracted and used to fill the gap.

Additionally, the p-arm telomeres of chromosomes 16 and X, and the q-arm telomere of chromosome 4, were identified among the unassigned sequences generated by Verkko and subsequently recovered through manual curation.

## Curation of RPE-1-Specific Structural Variant

Multi-step pipeline for the manual curation of the RPE-1-specific structural variant identified as *46,XX,dup(10q),t(Xq;10q),del(Xq28)*:

1. **Reads alignment**  
   The de novo diploid genome assembly RPE1v1.0 was used as a reference to map HiFi and ONT reads using [minimap2](https://github.com/lh3/minimap2).

2. **Visualization**  
   Long-read alignments were visualized using IGV, revealing an increased read coverage on the long arm of chromosome 10 (10q). A read interruption was observed, with ~100 bp difference in mapping position between chromosome 10 of Hap1 and Hap2.

3. **Assessment of haplotypic differences in read mapping quality**  
   Chromosome 10 Hap1 displayed only reads with mapping quality 0. In contrast, Hap2 exhibited reads with mapping qualities between 10–60 and supplementary alignments in the telomeric region of chromosome X Hap1, corresponding to the translocation breakpoint.

4. **Manual curation (translocation)**  
   The sequence of the duplicated and translocated 10q arm from chromosome 10 Hap2 was appended to the telomeric region of chromosome X Hap1. This was achieved by merging the sequences into the existing FASTA file of the RPE-1 assembly.

5. **Read alignment verification**  
   RPE-1 HiFi and ONT reads were re-aligned to the modified FASTA using minimap2. IGV inspection of the fusion point on chromosome X revealed a 3,603 bp microdeletion in the read alignments, suggesting sequence loss during the X–10 rearrangement.

6. **Manual curation (deletion)**  
   The deleted bases were manually removed from the FASTA sequence of chromosome X. The modified genome was then aligned again to the HiFi and ONT reads for confirmation.

7. **Final verification**  
   In the final RPE1v1.0 diploid genome, long reads align fully across the fusion point on chromosome X, confirming successful correction of the structural variant.

## Genome Quality

The quality of the RPE1v1.1 genome was assessed using a variety of tools:

- **[NucFreq](https://github.com/mrvollger/NucFreq)**  
  Primary alignments of HiFi reads to the RPE1 diploid genome were used to generate a NucFreq plot, which displays the frequencies of the most and second most common bases at each genomic position. The `HetDetection.R` script was used to identify possible errors in the assembly.

- **[Merqury](https://github.com/marbl/merqury)**  
  By comparing the k-mers in the HiFi reads to the k-mers found in the assembly, we obtained a quality value (QV) of 64.1 for Hap1 and 61.8 for Hap2. Furthermore, k-mer spectra revealed a multiplicity profile consistent with a near-complete, haplotype-resolved assembly.

- **[CRAQ](https://github.com/JiaoLaboratory/CRAQ)**  
  Using HiFi read clipping information, CRAQ detected only ~2 Mb of potential assembly errors, primarily corresponding to regions with unresolved gaps and specific telomeric sequences. The AQI score exceeded 98, underscoring the quality of the assembly.

- **[Flagger](https://github.com/mobinasri/flagger)**  
  Based on the alignment of HiFi reads to the RPE1v1.1 assembly, Flagger identified 2.4 Mb as assembly errors and 55 Mb as collapsed regions. The genomic coordinates of these low-confidence regions are provided as BED files and can be used as reference tracks in downstream analyses.

- **[Compleasm](https://github.com/huangnengCSU/compleasm)**  
  Compleasm identified 99.71% complete genes (with 3.28% duplicated) in Hap1, and 99.73% (with 0.7% duplicated) in Hap2. The missing genes were 0.21% and 0.09%, while the fragmented genes were 0.08% and 0.09% for Hap1 and Hap2, respectively.

- **[SecPhase](https://github.com/mobinasri/secphase)**  
  SecPhase was used to identify HiFi reads aligning to the incorrect haplotype. No read needed to be relocalized, indicating that Hi-C data were sufficient to phase the assembly into two haplotype blocks.

- **[breakpointR](https://github.com/daewoooo/breakpointR) R package**  
  Strand-seq reads from [Sanders et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7612647/), generated from 80 RPE-1 cells, were aligned using [BWA-MEM](https://github.com/lh3/bwa) to the RPE1v1.1 diploid genome. `breakpointR` was used to detect strand state patterns in each cell that supported the correct haplotype phasing of the assembly.

## Genome-to-genome comparison  

*RPE1 Hap1 vs. RPE1 Hap2*, *CHM13 vs. RPE1 Hap1*, and *CHM13 vs. RPE1 Hap2* comparisons were performed using two approaches:

- **[SyRI](https://github.com/schneebergerlab/syri)**  
Based on the  genome-to-genome alignments obtained using [minimap2](https://github.com/lh3/minimap2), SyRI was used to identify synthenic regions (*SYN*), inversions (*INV*), insertions (*INS*), translocations (*TRANS*), inverted translocations (*INVTR*), duplications (*DUP*), inverted duplications (*INVDP*), unaligned regions (*NOTAL*), highly diverged regions (*HDR*), copy gains (*CPG*) and SNPs (*SNP*). 

- **[NUCmer and dnadiff](https://github.com/mummer4/mummer)**  
On average, RPE1 haplotypes aligned over segments of 471 kb with a sequence identity of 99.83%. 393 Mb of Hap1 did not align to Hap2, and 323 Mb of Hap2 did not align to Hap1, reflecting the different length of the two haplotypes.

## Repeat annotation  

- **[RepeatMasker](https://www.repeatmasker.org/)** and **[HumAS-HMMER_for_AnVIL](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL?tab=readme-ov-file)**  
Used for annotation of repeats, including alpha-satellite monomers.

- **[StV](https://github.com/fedorrik/stv)**  
Used for HOR Structural Variant (StV) prediction based on HOR-monomer annotation.

- **[StainedGlass0.5](https://github.com/mrvollger/StainedGlass)**  
Generated sequence identity heatmaps.

## Gene annotation  
[Liftoff](https://github.com/agshumate/Liftoff) was used to map genes from the human transcriptome annotation of ENSEMBL 112 (GRCh38.p14 genome) to the RPE1v1.1 genome.

## Epigenetic analysis  
Centromere chromatin phased landscapes were obtained from RPE-1 CUT&RUN CENP-A [dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132193) (GSE132193). Short reads were aligned using [Bowtie](https://github.com/BenLangmead/bowtie) to the following reference genomes:
- RPE1.0 Hap1
- RPE1.0 Hap2
- CHM13v2.0
- GRCh38.p14
- HG002v1.0 maternal
- HG002v1.0 paternal
CENP-A high-confidence peaks were determined with [MACS3](https://github.com/macs3-project/MACS).

## Methylation profiles  
Methylation profiles for 5-methylcytosine (5mc) were generated from the ONT RPE-1 POD5 (3.6 TB) using [Dorado v4.2.0 basecalling model](https://github.com/nanoporetech/dorado/releases) and the output processed with Modkit. Following the evaluation of reads coverage values N<sub>canonical</sub> and N<sub>mod</sub>, we selected the filter (N<sub>mod</sub> / N<sub>valid_cov</sub>) >60 applied to the bedMethyl output. 
