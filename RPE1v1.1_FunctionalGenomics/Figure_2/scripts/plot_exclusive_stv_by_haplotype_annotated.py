import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import re

# -----------------------------------------------------------------------------
# Purpose
#   Plot haplotype-exclusive HOR (Higher-Order Repeat) tracks per chromosome.
#   For each chromosome, identify HORs present in one haplotype and absent
#   in the other, then draw a horizontal track with segments colored by
#   repeat-unit size (e.g., k-mer length) and export a PDF per haplotype.
#
# Inputs
#   - stv_processed_with_mer_colors.tsv (TSV)
#       Required columns (at minimum):
#         chrom               : chromosome + haplotype label (e.g., chr17_hap1)
#         start, end          : genomic coordinates for the HOR block
#         name                : string containing HOR base + unit (parsed below)
#         repeat_unit_size    : integer, used to map consistent colors
#         color_rgb           : tuple-like string '(R,G,B)' (overwritten here)
#
# Outputs
#   - one PDF per (chromosome, haplotype) in ./exclusive_tracks_only/
#     file name: {chromosome}_{hap}_exclusive_HORs_only.pdf
#


#PARAMETERS 
figsize = (12, 1.5)
output_dir = "exclusive_tracks_only"
os.makedirs(output_dir, exist_ok=True)

#LOAD DATA
df = pd.read_csv("stv_processed_with_mer_colors.tsv", sep="\t")
df['color_rgb'] = df['color_rgb'].apply(eval)

#ARSE HOR_base and HOR_unit FROM 'name' 
def extract_base_and_unit(name):
    # Expect pattern like ... :Sxxx.<unit>; capture S... as HOR_base and the rest as HOR_unit
    match = re.search(r":(S[^.]+)\.(.+)$", name)
    if match:
        return match.group(1), match.group(2)
    return None, None

df[['HOR_base', 'HOR_unit']] = df['name'].apply(lambda x: pd.Series(extract_base_and_unit(x)))


# Palette is cycled if there are more unique sizes than provided colors.
color_list = [
    (31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40),
    (148, 103, 189), (140, 86, 75), (227, 119, 194), (127, 127, 127),
    (188, 189, 34), (23, 190, 207)
]
unique_sizes = sorted(df['repeat_unit_size'].dropna().unique())
strong_color_map = {size: color_list[i % len(color_list)] for i, size in enumerate(unique_sizes)}
df['color_rgb'] = df['repeat_unit_size'].map(strong_color_map)

# PLOT ONLY EXCLUSIVE TRACKS PER HAPLOTYPE 
def plot_exclusive_track(df_excl, chrom_label, output_pdf):
    # Skip if no exclusive intervals
    if df_excl.empty:
        return
    fig, ax = plt.subplots(figsize=figsize)

    # Draw each interval as a thick horizontal segment at y=0
    for _, row in df_excl.iterrows():
        color = tuple(c / 255 for c in row['color_rgb'])
        ax.plot([row['start'], row['end']], [0, 0], lw=6, color=color)

    ax.set_title(f"{chrom_label} HORs exclusive only")
    ax.set_xlabel("Genomic Position (bp)")
    ax.set_yticks([])
    ax.set_ylim(-0.5, 0.5)

    legend_elements = [
        mpatches.Patch(color=tuple(c / 255 for c in strong_color_map[size]), label=f"{int(size)}-mer")
        for size in sorted(df_excl['repeat_unit_size'].dropna().unique())
    ]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.01, 1), loc='upper left', title="Repeat unit size")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, output_pdf))
    plt.close()

# Derive unique chromosome base names 
chrom_bases = sorted(set(c.split("_hap")[0] for c in df['chrom'].unique()))
for chrom_base in chrom_bases:
    # Partition by haplotype
    df_hap1 = df[df['chrom'] == f"{chrom_base}_hap1"].copy()
    df_hap2 = df[df['chrom'] == f"{chrom_base}_hap2"].copy()

    for hap_df, other_df, hap_label in [(df_hap1, df_hap2, 'hap1'), (df_hap2, df_hap1, 'hap2')]:
        if hap_df.empty:
            continue

        # Exclusive keys = present in current haplotype but not in the other
        excl_keys = set(zip(hap_df['HOR_base'], hap_df['HOR_unit'])) - set(zip(other_df['HOR_base'], other_df['HOR_unit']))

        # Filter rows to those exclusive HORs only
        df_excl = hap_df[hap_df.apply(lambda x: (x['HOR_base'], x['HOR_unit']) in excl_keys, axis=1)].copy()

        # Compose output file name and plot
        chrom_full = f"{chrom_base}_{hap_label}"
        plot_exclusive_track(df_excl, chrom_full, f"{chrom_full}_exclusive_HORs_only.pdf")
