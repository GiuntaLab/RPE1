# Giunta Lab Human Diploid Cell Lines Project
We have sequenced experimentally-ameanable human diploid laboratory cell lines and generated complete phased assemblies to use as matched reference genomes for sequencing data generated from the same cell line, an approach we refer to as *isogenomic reference genome*. This proof-of-concept calls for a comprehensive catalog of complete genome assemblies for commonly used cells for a widespread application of isogenomic reference genomes to enable high-precision multi-omics analyses.
## Latest assembly release
### RPE1v1.0
#### Sequencing data
We generated the complete phased genome assembly of one of the most widely used non-cancer cell lines (RPE-1) with a stable diploid karyotype. We produced state-of-the-art sequencing data using third generation sequencing, Pacific Biosciences (PacBio) high-fidelity (HiFi) and Oxford Nanopore Technologies (ONT). ONT long and ultra-long (UL; > 100 kb) reads were exclusively generated using the pore chemistry R10.14 (V14). We aimed at above-average reads depth coverage with 46x sequence coverage of HiFi reads, 80x of ONT, 30x ONT-UL reads over 100 kb with a QV score >20. During DNA extraction, we assessed using Femto Pulse the legnth of high molecular weight DNA from hTERT RPE-1 cells with a yield of native DNA size distribution around 116 kb (main peak) and up to 1 Mb in length (smaller peaks) when centrifugation was kept below 1000 RPM. For NGS, we generated 100x Illumina, 60x PRC-free Illumina and 30x reads for the Hi-C data (Arima Genomics) - see "Phasing of haplotypes".

#### Phasing of haplotypes
In absence of parental information to support [Trio-binning](https://www.nature.com/articles/nbt.4277), we obtained fully phased haplotypes for the RPE-1 genome using contact maps. We used the Vertebrate Genome Project (VGP) [pipeline](https://github.com/VGP/vgp-assembly). Briefly, Hi-C raw reads were aligned against the merged assembly composed of both RPE-1 haplotypes and unassigned reads generated from Verkko. The alignment step was performed using the short-read aligner [BWA](https://github.com/lh3/bwa). All reads were retained, including those classified as supplementary, with low mapping quality or having multiple alignments. The final aligned file was converted into a 3D Contact Map viewable through [PretextView](https://github.com/wtsi-hpag/PretextView). PretextView allows the modification of the contact map file by changing contig positions along the diagonal and finding the correct chromatin interaction path, thus generating a final diploid Hi-C contact map based on Nadolina Brajuka [RapidCuration2.0](https://github.com/Nadolina/Rapid-curation-2.0). The unassigned reads over 300 kb were assigned to the chromosome merged to the assembly file, yet a remaining 667 contigs <100 kb could not be aligned manually because of their short size. In the final contact map, each scaffold was assigned to a specific haplotype Meta Data Tag (Hap 1 and Hap 2). This latest contact map was then converted into an A Golden Path (AGP) file, which was used as input for the subsequent separation of both haplotypes and generating the two fully phased haploid FASTA files after running [Curation2.0_pipe.sh](https://github.com/Nadolina/Rapid-curation-2.0/blob/main/curation_2.0_pipe.sh).

#### Genome Quality
- CRAQ (Clipping Reveals Assembly Quality): We use PacBio HiFi and Illumina reads to evaluate the quality of the genome using [CRAQ}(https://github.com/JiaoLaboratory/CRAQ). We obtained:
- 






