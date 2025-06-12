# Figure 1 - Quality evaluation of the RPE-1 diploid genome assembly
Workflow followed to assess diploid assembly quality and correcteness.

## Prerequisites:
- Meryl
- Winnowmap
- samtools
- python3
- minimap2

## Mis-assemblies detection with *Flagger*
### Read mapping

PacBio HiFi reads were aligned with winnowmap against the final diploid genome following *Flagger* workflow:

```
winnowmap -W repetitive_k15.txt -ax map-pb -Y -L --eqx --cs -I8g [index_reference] [reads_hifi] | \
samtools view -hb | \
samtools sort -@ 128  > [reads.bam]
samtools index [reads.bam]
```

### Reads depth and coverage evaluation along the genome:

```
samtools depth -aa -Q 0 [reads.bam] > [reads.bam.depth]
depth2cov -d [reads.bam.depth] -f [assembly.fasta.fai] -o [reads.bam.depth.cov]
cov2counts -i [reads.bam.depth.cov] -o [reads_alignments.counts]
```
### Coverage based evaluation to extract haploid, error, duplicated and collpased regions of the genome:

```
fit_gmm.py --counts [reads_alignments.counts]  --output [reads_alignments.table]
find_blocks_from_table [reads_alignments.table]
```
This workflow outputs four bed files with haploid, erroneous, (falsely) duplicated and collapsed regions:
- RPE-1 haploid (6 Gb)
- RPE-1 errors (2 Mb)
- RPE-1 duplicated (0 B)
- RPE-1 collapsed (53 Mb)

Bed files were used for Figure 1 panel C. 

## Errors identification with *CRAQ*
Clipped information derived from raw reads mapped back to the assembly identify regional and structural assembly errors (https://github.com/JiaoLaboratory/CRAQ).
```
craq -mgs 1000 -q 10 -avgl 45 -avgs 60 -pl T -b T -D [output_directory] -g [index_reference] -sms [reads_hifi] -ngs [reads_illumina_pcrfree] -x map-hifi 
```
Bed files obtained using this tool were used for Figure 1 panel B.

# Figure 1 - Syntenic and structural rearrangments in RPE-1 genome:
Workflow used to identies specific structural rearrangements comparing RPE-1 haplotypes and haploid CHM13 genome.

## Prerequisites:
- python3
- plotsr
- minimap2

## Genome-to-genome alignment using *minimap2*

```
minimap2 -ax asm5 --eqx [refgenome] [qrygenome] > [out.sam]
```
## The SAM was used to identifies and plot synteny and rearrangements between CHM13/Hap 1-CHM13/Hap 2:

```
syri -c [out.sam] -r [refgenome] -q [qrygenome] -k -F S
plotsr --sr [syri.out] --genomes [genomes.txt] --tracks [centromeres_track.txt] --chrord [chord.txt] -H 6 -W 4 -b pdf -v > [plotsr.pdf]

```
Syntenic and structural rearrangements are shown in Figure 1 panel D.




