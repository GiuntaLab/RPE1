###COMPUTING PERCENTAGE IDENTITY FOR BAM FILES FROM BOWTIE###
#bowtie2 -x "$hap1_genome" -U $read_hap2 --end-to-end --sensitive --threads 10 | samtools view -bS | samtools sort -@ 8 - -O BAM -o "hap1_genome_$read_hap2_bowtie.bam"

import pysam
import sys

bam_file = sys.argv[1]
output_file = sys.argv[2]

bam = pysam.AlignmentFile(bam_file, "rb")
with open(output_file, "w") as out:
    out.write("percent_identity\tsample\n")
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        # NM = number of mismatches + indels
        nm = read.get_tag("NM") if read.has_tag("NM") else None
        if nm is None:
            continue
        ref_len = read.reference_length  # length on the reference (matches + mismatches + deletions)
        if ref_len == 0:
            continue
        pid = 100 * (ref_len - nm) / ref_len
        out.write(f"{pid:.2f}\t{bam_file}\n")
bam.close()
