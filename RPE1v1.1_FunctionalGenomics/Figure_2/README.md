# Figure 2 - Intra-haplotype variation between RPE-1 centromere pairs

(a) Example of a chromosome pair showing the two haplotypes, structural rearrangements between them identified by **SyRI**, and **StainedGlass** plots depicting sequence similarity within the centromeric live HOR array of each haplotype.  
(b) The same visualization as in (a) is shown for chromosomes 2, 4, 5, 6, 9, 11, 14, 16 and X of the RPE1v1.1 genome, along with haplotype-specific HOR structural variants annotated beside the StainedGlass plots.  
(c) Live α-satellite array length for each RPE1v1.1 chromosome; Hap1 (red) and Hap2 (orange) are shown.  
(d) Ratio of live α-satellite array lengths between the two haplotypes for each chromosome.  

**Chr:** chromosome **Hap:** haplotype **HOR:** higher-order repeats  

---

## Prerequisites

### Data
- RPE1v1.1 Hap1 assembly (FASTA)  
- RPE1v1.1 Hap2 assembly (FASTA)  
- CHM13 reference genome (FASTA)  
- *HumAS-HMMER_for_AnVIL* output on RPE1v1.1 (both haplotypes)  
- `stv_processed_with_mer_colors.tsv` (for exclusive HOR track plotting)

### Tools
- [**StainedGlass v6.7.0**](https://github.com/mrvollger/StainedGlass)  
- **SyRI** and **plotsr**  
- **minimap2**  
- [**HumAS-HMMER_for_AnVIL**](https://github.com/fedorrik/HumAS-HMMER_for_AnVIL)  
- [**StV**](https://github.com/fedorrik/stv)  
- **R** (ggplot2, dplyr, readr, stringr)  
- **GraphPad Prism v10.3.1**  
- **Adobe Illustrator**

---

## Genome similarity analysis with *StainedGlass* v6.7.0

StainedGlass was run with the following parameters:

```bash
window=1000
mm_f=10000
mm_s=5000
n=batch 4
```

Executed using Snakemake:
```bash
snakemake --cores 24 make_figures
```

This step produced the sequence-similarity plots shown in panels (a) and (b) of Figure 2.

---

## Structural rearrangement analysis with *SyRI*

Alignment of RPE1v1.1 Hap1 or Hap2 (query) to CHM13 (reference) using minimap2:

```bash
minimap2 -ax asm5 --eqx [refgenome] [qrygenome] > [out.sam]
```

SyRI run:
```bash
syri -c [out.sam] -r [refgenome] -q [qrygenome] -k -F S
```

Visualization with plotsr:
```bash
plotsr --sr [syri.out] --genomes [genomes.txt] --tracks [centromeres_track.txt]        --chrord [chord.txt] -H 6 -W 4 -b pdf -v > [plotsr.pdf]
```

SyRI output files are archived in the **Zenodo Repository**.  
The resulting plots were refined with **Adobe Illustrator** to compose panels (a) and (b).

---

## HOR structural-variant annotation with *StV*

Positions of live HORs detected by *HumAS-HMMER_for_AnVIL* on RPE1v1.1 (both Hap1 and Hap2) were used as input for *StV*.  
From the StV output, only **haplotype-specific structural variants** were retained and plotted beside the StainedGlass plots (panels a–b).

---

## Exclusive HOR tracks plotting

The script `plot_stv_exclusive_tracks.py` was used to generate per-chromosome exclusive HOR tracks from `stv_processed_with_mer_colors.tsv`.  
This script outputs PDFs in the `exclusive_tracks_only` directory, one for each chromosome–haplotype pair.

Run:
```bash
python plot_stv_exclusive_tracks.py
```

---

## Computation of n-mer size (*n* equal to number of monomers)

Repeat-unit sizes (monomer lengths) were derived from the *HumAS-HMMER* output and used as input for *StV*.  
Only structural variants (STVs) specific to each haplotype were plotted.

---

## Live HOR array length and ratio calculation

The final length of the live HOR (LHOR) regions and the ratio of lengths between Hap1 and Hap2 were plotted using **GraphPad Prism v10.3.1**.  
Panels (c) and (d) represent these data for each chromosome.

---

## Notes

- All plots use consistent coordinate spaces for Hap1 and Hap2 assemblies.  
- Final layout and labeling were composed in **Adobe Illustrator**.
