# Figure 3 - Whole-genome sequencing read alignment analysis

This figure illustrates the analysis of long- and short-read alignments of RPE-1 sequencing data against multiple reference genomes (RPE1v1.1, CHM13 v2.0, T2T-YAO, and HG002 v1.0.1).  
HiFi and Illumina reads were aligned, filtered for primary alignments, and analyzed for heterozygosity and alignment quality metrics.

---

## Prerequisites

### Data
- PacBio HiFi reads: `rpe1.cell1.fastq`, `rpe1.cell2.fastq`, `rpe1.cell3.fastq`, `rpe1.cell4.fastq`
- Illumina short reads (FASTQ)
- Reference genomes: `RPE1v1.1` (diploid), `CHM13 v2.0`, `T2T-YAO`, `HG002 v1.0.1`

### Tools
- [**minimap2 v2.17-r941**](https://github.com/lh3/minimap2)
- [**samtools v1.17**](https://github.com/samtools/samtools)
- [**NucFreq v0.1**](https://github.com/mrvollger/NucFreq)
- **HetDetection.R** (custom section by G. Logsdon)
- [**karyoploteR v1.28.0**](https://github.com/bernatgel/karyoploteR)
- [**cutadapt v4.0**](https://github.com/marcelm/cutadapt)
- [**BWA-MEM v0.7.17-r1188**](https://github.com/lh3/bwa)
- [**bedtools**](https://github.com/arq5x/bedtools2)
- **GraphPad Prism v10.3.1**

---

## HiFi read alignment

PacBio HiFi reads were aligned against diploid RPE1v1.1, CHM13, T2T-YAO, and HG002 references using minimap2:

```bash
minimap2 -ax map-hifi index_reference *hifi_fastq.gz --MD --secondary=no | samtools view -b - | samtools sort -@ 128 - -O BAM -o hifi_reads_primary.bam
samtools index hifi_reads_primary.bam
```

---

## NucFreq plotting

Primary alignments were extracted using samtools (`-F 2308`) and split by chromosome.  
Nucleotide frequency plots were generated with the **NucPlot.py** script from the *NucFreq* toolkit using the following bash pipeline:

```bash
#!/bin/bash

# Directory containing BAM files
BAM_DIR="/data/YAO_ALIGNEMENT/split"

# Path to NucPlot script
NUCPLOT="/data/NucFreq-0.1/NucPlot.py"

# Parameters
YMAX=300
MINOBED=2
DPI=500

# Loop through all BAM files in the directory
for BAM in "$BAM_DIR"/*.bam; do
    BASE=$(basename "$BAM" .bam)
    OBED="${BASE}.obed.bed"
    OUTPUT="${BASE}.png"
    echo "Processing $BAM..."
    python3 "$NUCPLOT" --obed "$OBED" --minobed "$MINOBED" --dpi "$DPI" -y "$YMAX" "$BAM" "$OUTPUT"
done
```

The plotting script above is stored in the `scripts/` directory and can be executed with:
```bash
bash nucplot.sh
```

For **CHM13 alignments**, heatmaps of heterozygosity density were plotted **below the NucFreq tracks** using **karyoploteR**.  
These heatmaps use as input the **output tables from `hetDetection.R`**, which identifies regions where the second most common base is present in at least 10% of reads in ≥5 positions within a 500 bp window.

---

## Heterozygosity detection and heatmap visualization

The *HetDetection.R* script was applied to the NucFreq output to identify heterozygous regions.  
Example snippet from the R workflow:

```r
library(karyoploteR)

pdf('/data/HetDetection_nucfreq/chr17_hap2_het.pdf', onefile = TRUE, paper = 'a4r', width = 30, height = 30)
pp <- getDefaultPlotParams(plot.type = 1)
pp$ideogramheight <- 20
pp$data1height <- 200
custom.genome <- toGRanges(data.frame(chr = c('chr17'), start = c(1), end = c(90000000)))

kp <- plotKaryotype(plot.type = 1, genome = custom.genome, chromosomes = 'chr17', cex = 0.5, main = 'Het chr17 hap2', plot.params = pp)

chr17_het <- read.table(file = "/data/HetDetection_nucfreq/chr17.hap2.HiFi.tbl", header = TRUE, sep = "	", stringsAsFactors = FALSE)
chr17_het <- setNames(chr17_het, c("chr", "start", "end", "het_ratio"))
chr17_het <- toGRanges(chr17_het[, c("chr", "start", "end", "het_ratio")])

kp <- kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 4, tick.col = 'black', units = 'Mb', add.units = TRUE)
kp <- kpHeatmap(kp, data.panel = 1, chr17_het, y = chr17_het$het_ratio, ymin = 0, ymax = max(chr17_het$het_ratio),
                r0 = 0.8, r1 = 1, color = c("darkorange", "white"))

legend(x = "right",
       legend = rep(NA, 101),
       fill = colorRampPalette(colors = c("darkorange", "white"))(101),
       border = NA, box.lwd = NA, y.intersp = 0.07, cex = 0.4)

dev.off()
```

This workflow was used to produce the heterozygosity heatmaps shown below the NucFreq plots in Figure 3.

---

## Short-read alignment and quality metric extraction

Illumina short reads were aligned to the same reference genomes using BWA-MEM, while HiFi reads were aligned using minimap2.

```bash
# Short-read alignment
bwa index ref.fa
bwa mem -k25 -p ref.fa *illumina.fastq | samtools view -b -F 2308 - | samtools sort -@ 128 - -O BAM -o hifi_reads_primary.bam
samtools index hifi_reads_primary.bam

# HiFi read alignment
minimap2 -ax map-hifi index_reference *hifi_fastq.gz --MD --secondary=no | samtools view -b - | samtools sort -@ 128 - -O BAM -o hifi_reads_primary.bam
samtools index hifi_reads_primary.bam
```

---

## Extraction of HDR and alignment metrics

### Extract HDR regions from SyRI output
```bash
grep 'HDR.*SYN\|SYN.*HDR' Syri.out > final_bed_file.bed
```

### Subset BAM by HDR coordinates
```bash
bedtools intersect -abam *.primary.bam -b *syri.hdr.syn.bed > hap1.chm13.subsample.primary.final.hdr.bam
```

### Extract NM (edit distance) and MAPQ (mapping quality)
```bash
# Extract NM value per read
ls */*bam | while read f; do
    NAME=$(basename $f .bam)
    samtools view $f | awk '{for(i=12;i<=NF;i++){if($i ~ /^NM:/){print $1 "\t" $i}}}' |     sed 's/NM:i://' | gzip -9 > $NAME.NM.tsv.gz
done

# Extract mapping quality
for bam_file in *bam; do
    output_file="${bam_file%.bam}.mapq.gz"
    samtools view "$bam_file" | awk '{print $1,$5}' | gzip > "$output_file"
done
```

The resulting files were used to create box plots and dot plots representing NM value and mapping quality distributions across different reference genomes.

---

## Statistical analysis and visualization

All alignment statistics (NM values, MAPQ scores, and per-chromosome HDR sums) were visualized in **GraphPad Prism v10.3.1**.  
Differences between reference genomes were assessed using **Student’s t-test**.  
Box plots, dot plots, and heterozygosity heatmaps correspond to panels (c–h) of Figure 3.

---

## Notes

- Primary alignments were consistently used for all analyses (`-F 2308`).  
- NucFreq plots were used to identify regions of heterozygosity across chromosomes.  
- Heterozygosity heatmaps (generated with *karyoploteR*) visualize the local density of heterozygous SNPs below the NucFreq tracks for CHM13.  
- Statistical visualization was performed in **GraphPad Prism v10.3.1**.
