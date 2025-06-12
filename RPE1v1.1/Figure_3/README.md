# Figure 3 - Variation between RPE-1 haplotypes
Workflow followed to calculate and plot SNP density and HDRs starting from RPE1v1.1 Hap1 vs. Hap2 comparison.

## Prerequisites
- python3
- R
- karyoploteR R package
- minimap2
- SyRI

#### Alignment of RPE1v1.1 Hap2 (query) to Hap1 (reference) using *minimap2*

```
minimap2 -ax asm5 --eqx [refgenome] [qrygenome] > [out.sam]
```
#### The output SAM file is provided as input to SyRI, along with reference and query sequence.

```
syri -c [out.sam] -r [refgenome] -q [qrygenome] -k -F S
```

#### SyRI output files are available in the **Zenodo Repository**. The `hap1.hap2.syri.out` file is used to create highly diverged regions (HDRs) and to calculate SNP density across 10 kb genomic bins using the snp_density_hdr.py Python script.

```
python snp_density_hdr.py [hap1.hap2.syri.out] --bin 10000 --output_density hap1.hap2.snp.density --output_hdr hap1.hap2.HDR.txt
```
#### The output files, along with the custom `RPE1_hap1_cenblack_stain.txt` file that specifies the positions of the live HORs detected by *HumAS-HMMER_for_AnVIL* on RPE1v1.1 Hap1, are used to generate the final karyoplot with the `SNP_density_karyoploteR.R` R script.


