#!/usr/bin/env Rscript
# Purpose: Plot heterozygosity heatmaps using karyoploteR from hetDetection output tables.
#
# Input table format (TSV):
#   chr    start    end    het_ratio
#
# Example paths (placeholders):
#   Input:  /path/to/HetDetection_nucfreq/chr${chromosome}.${assembly}.HiFi.tbl
#   Output: /path/to/HetDetection_nucfreq/chr${chromosome}_${assembly}_het.pdf
#
# Notes:
#   - Replace /path/to/ with your working directory.
#   - Works for any chromosome (e.g. 1–22, X) and any assembly (hap1 or hap2).
#   - Example: chr2.hap2.HiFi.tbl → chr2_hap2_het.pdf
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(karyoploteR)
  library(GenomicRanges)
})

# Define input/output paths (generic placeholders)
input_tbl  <- "/path/to/HetDetection_nucfreq/chr${chromosome}.${assembly}.HiFi.tbl"
output_pdf <- "/path/to/HetDetection_nucfreq/chr${chromosome}_${assembly}_het.pdf"

# Open PDF output
pdf(output_pdf, onefile = TRUE, paper = 'a4r', width = 30, height = 30)

# Plot parameters
pp <- getDefaultPlotParams(plot.type = 1)
pp$ideogramheight <- 20
pp$data1height <- 200

# Define a simple genome range for the selected chromosome (adjust end as needed)
custom.genome <- toGRanges(data.frame(chr = c("chr${chromosome}"), start = c(1), end = c(90000000)))

# Initialize karyotype plot
kp <- plotKaryotype(plot.type = 1, genome = custom.genome, chromosomes = "chr${chromosome}", cex = 0.5,
                    main = "Het chr${chromosome} ${assembly}", plot.params = pp)

# Load hetDetection table and convert to GRanges
het_tbl <- read.table(file = input_tbl, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(het_tbl) <- c("chr", "start", "end", "het_ratio")
gr <- toGRanges(het_tbl[, c("chr","start","end","het_ratio")])

# Add ticks and heatmap
kp <- kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 4, tick.col = 'black', units = 'Mb', add.units = TRUE)
kp <- kpHeatmap(kp, data.panel = 1, gr, y = gr$het_ratio, ymin = 0, ymax = max(gr$het_ratio, na.rm = TRUE),
                r0 = 0.8, r1 = 1, color = c("darkorange", "white"))

# Legend
lgd_ <- rep(NA, 101)
lgd_[c(101, 50, 1)] <- c(0, 20, 50)
legend(x = "right",
       legend = lgd_,
       fill = colorRampPalette(colors = c("darkorange", "white"))(101),
       border = NA, box.lwd = NA,
       y.intersp = 0.07,
       cex = 0.4)

dev.off()
