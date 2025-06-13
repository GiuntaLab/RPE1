# This script parses a SyRI output TSV file to identify SNP and HDR positions,
# computes SNP density in fixed-size bins along chromosomes,
# and writes the resulting densities to a BED-like output file and
# HDR positions to a text file.

import argparse
from collections import defaultdict

def get_chrom_base(chrom):
    return chrom.split('_')[0]

def parse_syri_file(filepath, hdroutpath):
    chrom_lengths = defaultdict(int)
    snp_positions = defaultdict(set)
    out_file = open(hdroutpath, "w")
    with open(filepath) as f:
        for line in f:
            if line.startswith("-") or line.strip() == "":
                continue
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            chrA, startA, endA = fields[0], int(fields[1]), int(fields[2])
            chrB = fields[5]
            event_type = fields[-2]

            # Update chromosome lengths
            baseA, baseB = get_chrom_base(chrA), get_chrom_base(chrB)            
            chrom_lengths[baseA] = max(chrom_lengths[baseA], endA)

            if baseA != baseB:
                continue


            # Store SNP positions and write HDRs to file
            if event_type == "SNP":
                snp_positions[baseA].add(startA)
            elif event_type == "HDR":
                out_file.write("\t".join([baseA,fields[1],fields[2]])+"\n")

    out_file.close()
    return chrom_lengths, snp_positions

def bin_snps(chrom_lengths, snp_positions, bin_size, output_file):
    with open(output_file, 'w') as out:
        for chrom in sorted(chrom_lengths):
            print(chrom)
            chrom_len = chrom_lengths[chrom]
            bins = range(0, chrom_len + bin_size, bin_size)
            snps = snp_positions.get(chrom, set())

            for start in bins:
                if start < chrom_len:
                       end = min(start + bin_size, chrom_len)
                       count = sum(1 for pos in snps if start <= pos < end)
                       density = count / (end - start) if (end - start) > 0 else 0
                       out.write(f"{chrom}\t{start}\t{end}\t{density:.6f}\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("syri_file", help="Input syri.out file")
    parser.add_argument("--bin", type=int, default=10000, help="Bin size (default: 10000)")
    parser.add_argument("--output_density", default="snp_density.bed", help="Output SNP density file name")
    parser.add_argument("--output_hdr", default="HDR.txt", help="Output HDR file name")    
    args = parser.parse_args()

    chrom_lengths, snp_positions = parse_syri_file(args.syri_file, args.output_hdr)
    bin_snps(chrom_lengths, snp_positions, args.bin, args.output_density)

if __name__ == "__main__":
    main()

