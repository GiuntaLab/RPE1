# Figure 5 - Comparison between HPRC assemblies and hTERT RPE-1 assembly.
Workflow followed to compare hTERT RPE-1 assembly with HPRC assemblies using pangenome graph construction and similarity analysis.

## Prerequisites
- wget
- samtools
- bgzip
- pggb
- odgi

#### Download and prepare hTERT RPE-1 assembly

```shell
cd /scratch
wget -c https://www.dropbox.com/scl/fo/i5v31h6zkc0l0x67efwnd/h/GIUNTAlab_RPE1v1.1.zip?rlkey=zmv4kz1unkqg28nljerwjnkuh&dl=0
unzip GIUNTAlab_RPE1v1.1.zip\?rlkey=zmv4kz1unkqg28nljerwjnkuh\&dl=0
sed '/>/ s/>chr\([^_]*\)_hap1/>rpe1#1#chr\1_hap1/; />/ s/>chr\([^_]*\)_hap2/>rpe1#2#chr\1_hap2/' GIUNTAlab_RPE1v1.1/rpe1v1.1.fasta | bgzip -@ 48 -l 9 > rpe1v1.1.fa.gz && samtools faidx rpe1v1.1.fa.gz
```

#### Download and process HPRC assemblies by chromosome

```shell
(seq 1 22; echo X; echo Y) | while read c; do
    chr=chr$c
    gfa_gz="$chr.hprc-v1.0-pggb.gfa.gz"
    wget -c "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/$gfa_gz"
    gunzip "$gfa_gz"
    cat \
        <(samtools faidx $dir_base/assemblies/rpe1v1.1.fa.gz $(grep "${chr}_" $dir_base/assemblies/rpe1v1.1.fa.gz.fai | cut -f 1)) \
        <(odgi paths -i $chr.hprc-v1.0-pggb.gfa -f -t 24 | sed -e 's/^>chm13#chr/>chm13#1#chr/' -e 's/^>grch38#chr/>grch38#1#chr/') | \
        bgzip -@ 16 -l 9 > rpe1+hprc.$chr.fa.gz && samtools faidx rpe1+hprc.$chr.fa.gz
done
```

#### Construct pangenome graphs with PGGB for each chromosome

```shell
(seq 1 22; echo X; echo Y) | while read c; do
    chr=chr$c
    pggb -i $dir_base/assemblies/rpe1+hprc.$chr.fa.gz -o pggb.rpe1+hprc.$chr -p 98 -s 10k -G 1300 -D /scratch/$chr -t 96
done
```

#### Compute similarity metrics for individual chromosomes using ODGI

```shell
(seq 1 22; echo X; echo Y) | while read c; do
    chr=chr$c
    odgi similarity -i $dir_base/graphs/pggb.rpe1+hprc.$chr/rpe1+hprc.$chr.fa.gz.*.smooth.final.og -d -D '#' -p 2 -t 96 > pggb.rpe1+hprc.$chr.similarity.by-haplotype.tsv
done
```

#### Combine all graphs and compute genome-wide similarity

```shell
ls $dir_base/graphs/*/*.og | sort -V > graphs-to-squeeze.txt
odgi squeeze -f graphs-to-squeeze.txt -o rpe1+hprc.all.og -t 96 -P
odgi similarity -i rpe1+hprc.all.og -d -D '#' -p 2 -t 96 > pggb.rpe1+hprc.all.similarity.by-haplotype.tsv
```

The resulting similarity matrices and pangenome graphs were visualized and edited with Adobe Illustrator to generate Figure 5.
