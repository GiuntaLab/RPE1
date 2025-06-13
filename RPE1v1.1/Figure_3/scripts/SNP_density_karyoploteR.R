library(karyoploteR)

pp <- getDefaultPlotParams(plot.type=1)
pp$data1height <- 350
pp$data1inmargin <- 20
pp$data1outmargin <- 100
pp$leftmargin <- 0.1


custom.genome <- toGRanges(data.frame(chr=c("chr1",'chr2', 'chr3', 'chr4', 'chr5', 'chr6','chr7','chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX'), start=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), end=c(245130294, 242676988, 201036031, 192521913, 182259134, 172045469, 163647842, 145680301, 137521496, 134988382, 134265560, 133236650, 101016140, 97561019, 95272790, 89871686, 82991831, 76653175, 60904190, 65674320, 38797553, 45250002, 227216509)))

custom.genome

# Load SNP density data, produced using the snp_density_calc.py script
snp_bed <- read.table("hap1.hap2.snp.density", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(snp_bed) <- c("chr", "start", "end", "snp_density")

snp_bed <- snp_bed[order(snp_bed$chr, snp_bed$start), ]

snp_bed <- subset(snp_bed, end > start)

# Convert to GRanges
snp_gr <- GRanges(seqnames = snp_bed$chr,
                  ranges = IRanges(start = snp_bed$start+1, end = snp_bed$end),
                  score = snp_bed$snp_density)

# Import HDR
hdr_table <- read.table("hap1.hap2.HDR.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(hdr_table) <- c("chr", "start", "end")

hdr_gr <- GRanges(seqnames = hdr_table$chr,
                  ranges = IRanges(start = hdr_table$start, end = hdr_table$end))

scaling.value <- 1.5


kp <- plotKaryotype(genome=custom.genome, 
                    cytobands = toGRanges("RPE1_hap1_cenblack_stain.txt"),
                    chromosomes='all', 
                    cex = 1.2, main='RPE-1 Hap1 vs Hap2 SNP density', 
                    plot.params = pp,
                    plot.type = 1)

kpAddBaseNumbers(kp, tick.dist=10000000, tick.len=6, tick.col='black',  units='Mb', add.units=TRUE, cex=0.7)

# Plot HDR
kpPlotRegions(kp, data=hdr_gr, avoid.overlapping = FALSE, col=transparent("lightgrey"))

# Plot SNP density as a line plot
kpLines(kp, data = snp_gr, y = snp_gr$score, col = "#5459AC", lwd = 0.5, ymin=0, ymax=0.08, cex =0.7)



# Optional: add axis
kpAxis(kp, ymin = 0, ymax = 0.08, cex = .8, tick.len=300000)





