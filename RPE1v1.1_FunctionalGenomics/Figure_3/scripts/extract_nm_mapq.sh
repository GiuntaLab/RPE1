#!/bin/bash
# Purpose: Extract alignment quality metrics (NM and MAPQ) from BAM files.
# Notes:
#   - Replace the glob patterns as needed.
#   - Produces compressed TSVs per BAM with NM and MAPQ values.
#   - Assumes samtools is available in PATH.

# NM extraction
# For every BAM found under subfolders: produce NAME.NM.tsv.gz with read_id and NM tag
ls */*bam | while read f; do
    echo "$f"
    NAME=$(basename "$f" .bam)
    samtools view "$f" |       awk '{for(i=12;i<=NF;i++){if($i ~ /^NM:/){print $1 "\t" $i}}}' |       sed 's/NM:i://' | gzip -9 > "$NAME.NM.tsv.gz"
done

#MAPQ extraction
# For every BAM in current directory: produce <bam>.output.gz with read_id and MAPQ
for bam_file in *bam; do
    output_file="${bam_file%.bam}.output.gz"
    samtools view "$bam_file" | awk '{print $1,$5}' | gzip > "$output_file"
done
