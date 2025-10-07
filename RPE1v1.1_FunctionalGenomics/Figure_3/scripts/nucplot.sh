#!/bin/bash
# Purpose: Generate NucFreq (NucPlot.py) nucleotide frequency plots per chromosome BAM.
# Inputs:
#   - BAM_DIR: directory containing per-chromosome BAM files (primary alignments only; e.g., filtered with -F 2308)
#   - NUCPLOT: path to NucFreq's NucPlot.py
# Parameters:
#   - YMAX: y-axis max for the plot
#   - MINOBED: minimum number of observations for OBED output
#   - DPI: image resolution
# Process (pseudocode):
#   for each BAM in BAM_DIR:
#     derive base name
#     define OBED (NucFreq bed-like output)
#     define PNG output path
#     run NucPlot.py with the chosen parameters to create the plot and the OBED bed
#
# Edit these paths:
BAM_DIR="/path/to/split"
NUCPLOT="/path/to/NucFreq-0.1/NucPlot.py"
# Parameters (same as your original)
YMAX=300
MINOBED=2
DPI=500

# Loop over BAMs 
for BAM in "$BAM_DIR"/*.bam; do
    BASE=$(basename "$BAM" .bam)
    OBED="${BASE}.obed.bed"
    OUTPUT="${BASE}.png"

    echo "Processing $BAM..."
    python3 "$NUCPLOT" --obed "$OBED" --minobed "$MINOBED" --dpi "$DPI" -y "$YMAX" "$BAM" "$OUTPUT"
done
