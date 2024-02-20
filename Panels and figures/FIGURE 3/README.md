# Figure 3- The isogenomic reference genome improves reads alignment
Workflow used to test and evaluate reference matched genome alignment improvements

## Prerequisites
- Minimap2
- python3
- Nucfreq
- bedtools

## HiFi reads alignment using diploid and haploid genome as a reference
### Long reads alignment

```
minimap2  -ax map-hifi [index_reference] [*/*hifi_fastq.gz] | samtools view -b  - | samtools sort -@ 128 - -O BAM -o  [hifi_reads.bam]
samtools index [hifi_reads.bam]
```
### Extract each chromosome from BAM file and visualize read-depth profile using *Nucfreq*

```
samtools view -b [hifi_reads.bam] [list.txt] > [chr*.bam]
NucFreq-0.1/NucPlot.py  [chr*.bam] --ylim 300 dpi 500 [chr*.png]
```
## RPE-1 short and long reads sequencing from two different batches aligned against from RPE-1 Hap 1, Hap 2, and CHM13 genomes
### Short reads alignment keeping only primary alignment

```
bwa index [ref.fa]
bwa mem -k25 -p [ref.fa] [*/*R1.illumina.fastq] [*/*R2.illumina.fastq]| samtools view -b -F 2308 - | samtools sort -@ 128 - -O BAM -o  [illumina_reads_primary.bam]
samtools index [illumina_reads_primary.bam]
```

### Long reads alignment keeping only primary alignment

```
minimap2  -ax map-hifi [index_reference] [*/*hifi_fastq.gz]  --MD --secondary=no| samtools view -b  - | samtools sort -@ 128 - -O BAM -o  [hifi_reads_primary.bam]
samtools index [hifi_reads_primary.bam]
```

## Extract HDR in syntenic regions from Syri output

```
grep 'HDR.*SYN\|SYN.*HDR' [syri.out] > [hdr.syntenic.bed]
```
### Extract edit distance and mapping quality from HDR BAM files
Extract bed file coordinates from BAM file

```
bedtools intersect -abam [*.primary.bam] -b [*/*syri.hdr.syn.bed] > [*/*primaryhdr.bam]
```
### Extract edit distance and mapping quality values from each BAM file

```
ls */*bam | while read f; do echo $f; NAME=$(basename $f .bam); samtools view $f |  awk '{for(i=12;i<=NF;i++){if($i ~ /^NM:/){print $1 "\t" $i}}}' | sed 's/NM:i://' | gzip -9 > $NAME.NM.tsv.gz; done

for bam_file in *bam; do
output_file="${bam_file%.bam}.output.gz"  
samtools view "$bam_file" | awk '{print $1,$5}' | gzip > "$output_file"
done
```
The output files were used to create boxplot graphs showed in Figure 3 panel C-H
