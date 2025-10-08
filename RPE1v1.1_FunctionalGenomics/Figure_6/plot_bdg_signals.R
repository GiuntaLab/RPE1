library(regioneR)
library(karyoploteR)
library(GenomicRanges)
library(rtracklayer)

#####LOAD INPUT FILES#####
#Custom genome
final_custom_genome #Custom genome of interest (bed file), with chr names (e.g. chr1_chm13, chr1_hap1, chr1_hap2..), start and end coordinates of each region to plot

#Centromere annotation
centr_annot_hap1 #Centromere annotation provided by HumAS https://github.com/fedorrik/HumAS-HMMER_for_AnVIL
centr_annot_hap2 #Centromere annotation provided by HumAS https://github.com/fedorrik/HumAS-HMMER_for_AnVIL
centr_annot_chm13 #Centromere annotation provided by HumAS https://github.com/fedorrik/HumAS-HMMER_for_AnVIL

#Methylation
meth_hap1 #5mC bedgraph file of Hap1
meth_hap1 #5mC bedgraph file of Hap1
meth_chm13 #5mC bedgraph file of Hap1

#Bedgraph files of CENP-A peaks produced by bdgcmp (MACS3 command)
cenpa_chm13_hap2 #Signals of CENP-A from CHM13 cell line mapped to RPE1 haplotypes
cenpa_chm13_q20_hap2 #Signals of CENP-A from CHM13 cell line mapped to RPE1 haplotypes and without alignments with MAPQ<20

cenpa_rpe1_hap2 #Signals of CENP-A from RPE-1 cell line mapped to RPE1 haplotypes
cenpa_rpe1_q20_hap2 #Signals of CENP-A from RPE-1 cell line mapped to RPE1 haplotypes and without alignments with MAPQ<20

cenpa_rpe1_chm13 #Signals of CENP-A from RPE-1 cell line mapped to CHM13
cenpa_rpe1_q20_chm13 #Signals of CENP-A from RPE-1 cell line mapped to CHM13 and without alignments with MAPQ<20

cenpa_chm13_chm13 #Signals of CENP-A from CHM13 cell line mapped to CHM13
cenpa_chm13_q20_chm13 #Signals of CENP-A from CHM13 cell line mapped to CHM13 and without alignments with MAPQ<20

cenpa_chm13_hap1 #Signals of CENP-A from CHM13 cell line mapped to RPE1 haplotypes
cenpa_chm13_q20_hap1 #Signals of CENP-A from CHM13 cell line mapped to RPE1 haplotypes and without alignments with MAPQ<20

cenpa_rpe1_hap1 #Signals of CENP-A from RPE-1 cell line mapped to RPE1 haplotypes
cenpa_rpe1_q20_hap1 #Signals of CENP-A from RPE-1 cell line mapped to RPE1 haplotypes and without alignments with MAPQ<20

#Colors
color_chm13_rpe1 <- "darkorchid4"
color_chm13 <- "#009A67"
color_rpe1hap1 <- "#E86e0F"
color_rpe1hap2 <- "#F39200"
background <- "#d3d5d3"
meth_col <- "black"

####PLOTS####
#loop_plot_for_cenpa_rpe1_chm13_q20####
chr_names <- c(paste0("chr", 1:22), "chrX")
#pdf("/home/lcorda/chr19_plot.pdf", paper = "a4")
for (chr in chr_names) {

  #Define the path and file to save
  path <- "/data/isogenomic_plots"
  jpeg(file.path(path, paste0(chr, "_plot.jpeg")), res = 1200, width = 8000, height = 8000, quality = 100)

  #Define graphical parameters
  pp <- getDefaultPlotParams(plot.type=1)
  pp$ideogramheight <- 0
  pp$data2max <- 0
  pp$data1max <- 100
  pp$data1height <- 400
  pp$bottommargin <- 100
  pp$topmargin <- 100
  pp$data1inmargin <- 0
  pp$leftmargin <- 0.35
  pp$rightmargin <- 0.15

  #Define the chromosome to consider
  chr_zoom <- c(paste0(chr, "_hap1"), paste0(chr, "_chm13"), paste0(chr, "_hap2"))
  final_custom_genome_f <- final_custom_genome[seqnames(final_custom_genome) %in% chr_zoom]

  #Plot ideograms
  kp <- plotKaryotype(plot.type=1, genome = final_custom_genome_f, plot.params=pp, cex=0.8, chromosomes = chr_zoom, cytobands = NULL)
  kpAddMainTitle(kp, main = "CENP-A profile within multiple genomes (enrichment = grey, Q20 = colors, meth = black)", cex = 0.8)
  kpAddBaseNumbers(kp, tick.dist = 1000000, tick.len = 2, tick.col="black", cex=0.5, minor.tick.dist = 500000, minor.tick.col = "black", minor.tick.len = 1.5, units = "Mb")

  #Plot centromeric arrays annotation
  #Hap1
  centr_annot_hap1_f <- subsetByOverlaps(centr_annot_hap1, final_custom_genome_f, type = "within")
  kp <- kpPlotRegions(kp, data.panel=1,r0=3, r1=5, data = centr_annot_hap1_f, col=centr_annot_hap1_f$gieStain, avoid.overlapping = F)
  #Hap2
  centr_annot_chm13_f <- subsetByOverlaps(centr_annot_chm13, final_custom_genome_f, type = "within")
  kp <- kpPlotRegions(kp, data.panel=1,r0=3, r1=5, data = centr_annot_chm13_f, col=centr_annot_chm13_f$gieStain, avoid.overlapping = F)
  #CHM13
  centr_annot_hap2_f <- subsetByOverlaps(centr_annot_hap2, final_custom_genome_f, type = "within")
  kp <- kpPlotRegions(kp, data.panel=1,r0=3, r1=5, data = centr_annot_hap2_f, col=centr_annot_hap2_f$gieStain, avoid.overlapping = F)

  #Plot bdg signals across datasets and genomes
  #chm13-rpe1-hap2-cenpa (chm13 dataset mapped against rpe1 hap2 genome)
  cenpa_chm13_hap2_f <- subsetByOverlaps(cenpa_chm13_hap2, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_hap2"), data.panel=1,r0=8, r1=20, data = cenpa_chm13_hap2_f, y1 = cenpa_chm13_hap2_f$V4, ymax = max(cenpa_chm13_hap2_f$V4), col=background, border = background)
  kpAxis(kp, chr = paste0(chr, "_hap2"), data = cenpa_chm13_hap2_f, data.panel=1,r0=8, r1=20, ymax=max(cenpa_chm13_hap2_f$V4), cex=0.4, side = 2)
  cenpa_chm13_q20_hap2_f <- subsetByOverlaps(cenpa_chm13_q20_hap2, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_hap2"), data.panel=1,r0=8, r1=20, data = cenpa_chm13_q20_hap2_f, y1 = cenpa_chm13_q20_hap2_f$V4, ymax = max(cenpa_chm13_q20_hap2_f$V4), col=color_chm13_rpe1, border = color_chm13_rpe1)
  kpAxis(kp, chr = paste0(chr, "_hap2"), data = cenpa_rpe1_q20_hap1_f, data.panel=1,r0=8, r1=20, ymax=max(cenpa_rpe1_q20_hap1_f$V4), cex=0.4, side = 1)

  #rpe1-rpe1-hap2-cenpa (rpe1 dataset mapped against rpe1 genome)
  cenpa_rpe1_hap2_f <- subsetByOverlaps(cenpa_rpe1_hap2, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_hap2"), data.panel=1,r0=24, r1=36, data = cenpa_rpe1_hap2_f, y1 = cenpa_rpe1_hap2_f$V4, ymax = max(cenpa_rpe1_hap2_f$V4), col=background, border = background)
  kpAxis(kp, chr = paste0(chr, "_hap2"), data = cenpa_rpe1_hap2_f, data.panel=1,r0=24, r1=36, ymax=max(cenpa_rpe1_hap2_f$V4), cex=0.4, side = 2)
  cenpa_rpe1_q20_hap2_f <- subsetByOverlaps( cenpa_rpe1_q20_hap2, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_hap2"), data.panel=1,r0=24, r1=36, data = cenpa_rpe1_q20_hap2_f, y1 = cenpa_rpe1_q20_hap2_f$V4, ymax = max(cenpa_rpe1_q20_hap2_f$V4), col=color_rpe1hap2, border = color_rpe1hap2)
  kpAxis(kp, chr = paste0(chr, "_hap2"), data = cenpa_rpe1_q20_hap2_f, data.panel=1,r0=24, r1=36, ymax=max(cenpa_rpe1_q20_hap2_f$V4), cex=0.4, side = 1)

  #hap2-meth (5mC signals of hap2)
  meth_hap2_f <- subsetByOverlaps(meth_hap2, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_hap2"), data.panel=1,r0=40, r1=52, data = meth_hap2_f, y1 = meth_hap2_f$V4, ymax = 1, col=meth_col)
  kpAxis(kp, chr = paste0(chr, "_hap2"), data = meth_hap2_f, data.panel=1,r0=40, r1=52, ymax=1, cex=0.4, side = 1)
  #kpAddLabels(kp, labels="5mC RPE1 hap2",r0=40, r1=52, data.panel = 1, cex = 0.5, side = "left", label.margin = 0.07, font = 1)

  #rpe1-chm13-cenpa (rpe1 dataset mapped against chm13 genome)
  cenpa_rpe1_chm13_f <- subsetByOverlaps(cenpa_rpe1_chm13, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_chm13"), data.panel=1,r0=8, r1=20, data = cenpa_rpe1_chm13_f, y1 = cenpa_rpe1_chm13_f$V4, ymax = max(cenpa_rpe1_chm13_f$V4), col=background, border = background)
  kpAxis(kp, chr = paste0(chr, "_chm13"), data = cenpa_rpe1_chm13_f, data.panel=1,r0=8, r1=20, ymax=max(cenpa_rpe1_chm13_f$V4), cex=0.4, side = 2)
  cenpa_rpe1_q20_chm13_f <- subsetByOverlaps(cenpa_rpe1_q20_chm13, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_chm13"), data.panel=1,r0=8, r1=20, data = cenpa_rpe1_q20_chm13_f, y1 = cenpa_rpe1_q20_chm13_f$V4, ymax = max(cenpa_rpe1_q20_chm13_f$V4), col=color_chm13, border = color_chm13)
  kpAxis(kp, chr = paste0(chr, "_chm13"), data = cenpa_rpe1_q20_chm13_f, data.panel=1,r0=8, r1=20, ymax=max(cenpa_rpe1_q20_chm13_f$V4), cex=0.4, side = 1)

  #chm13-chm13-cenpa (chm13 dataset mapped against chm13 genome)
  cenpa_chm13_chm13_f <- subsetByOverlaps(cenpa_chm13_chm13, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_chm13"), data.panel=1,r0=24, r1=36, data = cenpa_chm13_chm13_f, y1 = cenpa_chm13_chm13_f$V4, ymax = max(cenpa_chm13_chm13_f$V4), col=background, border = background)
  kpAxis(kp, chr = paste0(chr, "_chm13"), data = cenpa_chm13_chm13_f, data.panel=1,r0=24, r1=36, ymax=max(cenpa_chm13_chm13_f$V4), cex=0.4, side = 2)
  cenpa_chm13_q20_chm13_f <- subsetByOverlaps(cenpa_chm13_q20_chm13, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_chm13"), data.panel=1, r0=24, r1=36, data = cenpa_chm13_q20_chm13_f, y1 = cenpa_chm13_q20_chm13_f$V4, ymax = max(cenpa_chm13_q20_chm13_f$V4), col=color_chm13_rpe1, border = color_chm13_rpe1)
  kpAxis(kp, chr = paste0(chr, "_chm13"), data = cenpa_chm13_q20_chm13_f, data.panel=1,r0=24, r1=36, ymax=max(cenpa_chm13_q20_chm13_f$V4), cex=0.4, side = 1)

  #chm13-meth (5mC signals of chm13)
  meth_chm13_f <- subsetByOverlaps(meth_chm13, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_chm13"), data.panel=1,r0=40, r1=52, data = meth_chm13_f, y1 = meth_chm13_f$V4, ymax = 1, col=meth_col)
  kpAxis(kp, chr = paste0(chr, "_chm13"), data = meth_chm13_f, data.panel=1,r0=40, r1=52, ymax=1, cex=0.4, side = 1)
  #kpAddLabels(kp, labels="5mC CHM13",r0=40, r1=52, data.panel = 1, cex = 0.5, side = "left", label.margin = 0.07, font = 1)

  #chm13-rpe1-hap1-cenpa (chm13 dataset mapped against rpe1 hap1 genome)
  cenpa_chm13_hap1_f <- subsetByOverlaps(cenpa_chm13_hap1, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_hap1"), data.panel=1,r0=8, r1=20, data = cenpa_chm13_hap1_f, y1 = cenpa_chm13_hap1_f$V4, ymax = max(cenpa_chm13_hap1_f$V4), col=background, border = background)
  kpAxis(kp, chr = paste0(chr, "_hap1"), data = cenpa_chm13_hap1_f, data.panel=1,r0=8, r1=20, ymax=max(cenpa_chm13_hap1_f$V4), cex=0.4, side = 2)
  cenpa_chm13_q20_hap1_f <- subsetByOverlaps(cenpa_chm13_q20_hap1, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_hap1"), data.panel=1,r0=8, r1=20, data = cenpa_chm13_q20_hap1_f, y1 = cenpa_chm13_q20_hap1_f$V4, ymax = max(cenpa_chm13_q20_hap1_f$V4), col=color_chm13_rpe1, border = color_chm13_rpe1)
  kpAxis(kp, chr = paste0(chr, "_hap1"), data = cenpa_rpe1_q20_hap1_f, data.panel=1,r0=8, r1=20, ymax=max(cenpa_rpe1_q20_hap1_f$V4), cex=0.4, side = 1)

  #rpe1-rpe1-hap1-cenpa (rpe1 dataset mapped against rpe1 hap1 genome)
  cenpa_rpe1_hap1_f <- subsetByOverlaps(cenpa_rpe1_hap1, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_hap1"), data.panel=1,r0=24, r1=36, data = cenpa_rpe1_hap1_f, y1 = cenpa_rpe1_hap1_f$V4, ymax = max(cenpa_rpe1_hap1_f$V4), col=background, border = background)
  kpAxis(kp, chr = paste0(chr, "_hap1"), data = cenpa_rpe1_hap1_f, data.panel=1,r0=24, r1=36, ymax=max(cenpa_rpe1_hap1_f$V4), cex=0.4, side = 2)
  cenpa_rpe1_q20_hap1_f <- subsetByOverlaps( cenpa_rpe1_q20_hap1, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_hap1"), data.panel=1,r0=24, r1=36, data = cenpa_rpe1_q20_hap1_f, y1 = cenpa_rpe1_q20_hap1_f$V4, ymax = max(cenpa_rpe1_q20_hap1_f$V4), col=color_rpe1hap1, border = color_rpe1hap1)
  kpAxis(kp, chr = paste0(chr, "_hap1"), data = cenpa_rpe1_q20_hap1_f, data.panel=1,r0=24, r1=36, ymax=max(cenpa_rpe1_q20_hap1_f$V4), cex=0.4, side = 1)

  #hap1-meth (5mC signals of hap1)
  meth_hap1_f <- subsetByOverlaps(meth_hap1, final_custom_genome_f, type = "within")
  kp <- kpBars(kp, chr = paste0(chr, "_hap1"), data.panel=1,r0=40, r1=52, data = meth_hap1_f, y1 = meth_hap1_f$V4, ymax = 1, col=meth_col)
  kpAxis(kp, chr = paste0(chr, "_hap1"), data = meth_hap1_f, data.panel=1,r0=40, r1=52, ymax=1, cex=0.4, side = 1)
  #kpAddLabels(kp, labels="5mC RPE1 hap1",r0=40, r1=52, data.panel = 1, cex = 0.5, side = "left", label.margin = 0.07, font = 1)

  dev.off()
}