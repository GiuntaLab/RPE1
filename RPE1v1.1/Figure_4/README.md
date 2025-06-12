# Figure 4 - Variation between RPE-1 and CHM13 genomes.
Workflow followed to compare RPE1v1.1 haplotypes to CHM13 genome using SyRI.

## Prerequisites
- minimap2
- SyRI
- plotsr

#### Alignment of RPE1v1.1 Hap1 or Hap2 (query) to CHM13 (reference) using *minimap2*

```
minimap2 -ax asm5 --eqx [refgenome] [qrygenome] > [out.sam]
```

#### The output SAM file is provided as input to SyRI, along with reference and query sequence. plotsr is used to visualize the structural rearrangements. SyRI output files are available in the **Zenodo Repository**. 

```
syri -c [out.sam] -r [refgenome] -q [qrygenome] -k -F S
plotsr --sr [syri.out] --genomes [genomes.txt] --tracks [centromeres_track.txt] --chrord [chord.txt] -H 6 -W 4 -b pdf -v > [plotsr.pdf]
```
The resulting plot was edited with Adobe Illustrator to generate Figure 4.
