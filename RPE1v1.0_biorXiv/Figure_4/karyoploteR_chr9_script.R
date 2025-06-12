#clear the current workspace or environment
rm(list = ls())

setwd("")

library(karyoploteR)

pdf("FIG4_chr9_cenpa_meth_over5_diploid.pdf", onefile = T, paper = "a4")

#Loading the methylation profile of each haplotypes (+ strand and score => 60) 
meth1 <- toGRanges("RPE1.v1.0.modkit.bedMethyl.+.over60.hap1_chr9.bdg")
meth2 <- toGRanges("RPE1.v1.0.modkit.bedMethyl.+.over60.hap2_chr9.bdg")

#Loading the CENP-A profile of each haplotypes
#Enrichment peaks: qvalue <= 0.00001
RPE1v1.0dip_CENPA_hap1_over5_chr9.bdg <- toGRanges("RPE1v1.0dip_CENPA_hap1_over5_chr9.bdg")
RPE1v1.0dip_CENPA_hap2_over5_chr9.bdg <- toGRanges("RPE1v1.0dip_CENPA_hap2_over5_chr9.bdg")
HG002v1.0dip_CENPA_mat_over5_chr9.bdg <- toGRanges("HG002v1.0dip_CENPA_maternal_over5_chr9.bdg")
HG002v1.0dip_CENPA_pat_over5_chr9.bdg <- toGRanges("HG002v1.0dip_CENPA_paternal_over5_chr9.bdg")
CHM13_CENPA_over5_chr9.bdg <- toGRanges("CHM13_CENPA_over5_chr9.bdg")
HG38_CENPA_over5_chr9.bdg <- toGRanges("HG38_CENPA_over5_chr9.bdg")

#High-confidence peaks: qvalue <= 0.00001 and MAPQ => 20
RPE1v1.0dip_CENPAmapQ20_hap1_over5_chr9.bdg <- toGRanges("RPE1v1.0dip_CENPAmapQ20_hap1_over5_chr9.bdg")
RPE1v1.0dip_CENPAmapQ20_hap2_over5_chr9.bdg <- toGRanges("RPE1v1.0dip_CENPAmapQ20_hap2_over5_chr9.bdg")
HG002v1.0dip_CENPAmapQ20_mat_over5_chr9.bdg <- toGRanges("HG002v1.0dip_CENPAmapQ20_maternal_over5_chr9.bdg")
HG002v1.0dip_CENPAmapQ20_pat_over5_chr9.bdg <- toGRanges("HG002v1.0dip_CENPAmapQ20_paternal_over5_chr9.bdg")
CHM13_CENPA_mapQ20_over5_chr9.bdg <- toGRanges("CHM13_CENPA_mapQ20_over5_chr9.bdg")
HG38_CENPA_mapQ20_over5_chr9.bdg <- toGRanges("HG38_CENPA_mapQ20_over5_chr9.bdg")

#set plot parameters 6Mb
pp <- getDefaultPlotParams(plot.type=2)
pp$ideogramheight <- 15
pp$ideogramlateralmargin <- 0.008
pp$data2height <- 1200
pp$data1inmargin <- 10
pp$data2inmargin <- 30
pp$data1outmargin <-10
pp$data1height <- 1200
pp$leftmargin <- 0.3
pp$rightmargin <- 0.3

#chr9 custom ideogram
customGenome <- toGRanges(data.frame(chr = "9" , start = 1, end = 150617247)) #chr9 length in CHM13v2.0

#creating karyoplote ideogram
kp <- plotKaryotype(plot.type=2, genome = customGenome, plot.params=pp, cex=0.5, zoom = GRanges("9:44000000-50000000"))
kpAddBaseNumbers(kp, tick.dist = 2000000, tick.len = 10, tick.col="black", cex=0.5, minor.tick.dist = 25000000, minor.tick.col = "black")

#CENPA RPE1 HAP2 diploid alignment
kp <- kpPlotDensity(kp, data.panel=1, data = RPE1v1.0dip_CENPA_hap2_over5_chr9.bdg, r0=0, r1=0.1, window.size = 250, col= "#DCDCDC", border = "#DCDCDC")
kp <- kpPlotDensity(kp, data.panel=1, data = RPE1v1.0dip_CENPAmapQ20_hap2_over5_chr9.bdg, r0=0, r1=0.1, window.size = 50, col= "#f7941d", border = "#f7941d")
kpAxis(kp, data.panel=1, ymax=kp$latest.plot$computed.values$max.density, tick.len = 0.02e6, r0=0, r1=0.1, cex=0.25, numticks = 2)

#5mC RPE1 HAP2 diploid alignment
kp <- kpPlotDensity(kp, data.panel=1, data = meth2, r0=0.11, r1=0.21, window.size = 5000, border = "#000000", col= "#000000")
kpAxis(kp, data.panel=1, ymax=kp$latest.plot$computed.values$max.density, tick.len = 0.02e6, r0=0.11, r1=0.21, cex=0.25, numticks = 2)

#CENPA CHM13 haploid alignment
kp <- kpPlotDensity(kp, data.panel=1, data = CHM13_CENPA_over5_chr9.bdg, r0=0.22, r1=0.32, window.size = 250, col= "#DCDCDC", border = "#DCDCDC")
kp <- kpPlotDensity(kp, data.panel=1, data = CHM13_CENPA_mapQ20_over5_chr9.bdg, r0=0.22, r1=0.32, window.size = 50, col= "#00a875", border = "#00a875")
kpAxis(kp, data.panel=1, ymax=kp$latest.plot$computed.values$max.density, tick.len = 0.02e6, r0=0.22, r1=0.32, cex=0.25, numticks = 2)

#CENPA HG002v1.0 maternal diploid alignment
kp <- kpPlotDensity(kp, data.panel=1, data = HG002v1.0dip_CENPA_mat_over5_chr9.bdg, r0=0.33, r1=0.43, window.size = 250,  col= "#DCDCDC", border = "#DCDCDC")
kp <- kpPlotDensity(kp, data.panel=1, data = HG002v1.0dip_CENPAmapQ20_mat_over5_chr9.bdg, r0=0.33, r1=0.43, window.size = 50,  col= "#0072bc", border = "#0072bc")
kpAxis(kp, data.panel=1, ymax=kp$latest.plot$computed.values$max.density, tick.len = 0.02e6, r0=0.33, r1=0.43, cex=0.25, numticks = 2)

#CENPA HG002v1.0 paternal diploid alignment
kp <- kpPlotDensity(kp, data.panel=1, data = HG002v1.0dip_CENPA_pat_over5_chr9.bdg,, r0=0.44, r1=0.54, window.size = 250,  border = "#DCDCDC", col="#DCDCDC")
kp <- kpPlotDensity(kp, data.panel=1, data = HG002v1.0dip_CENPAmapQ20_pat_over5_chr9.bdg, r0=0.44, r1=0.54, window.size = 50,  border = "#00b9f2", col="#00b9f2")
kpAxis(kp, data.panel=1, ymax=kp$latest.plot$computed.values$max.density, tick.len = 0.02e6,  r0=0.44, r1=0.54, cex=0.25, numticks = 2)

#CENPA HG38 haploid alignment
kp <- kpPlotDensity(kp, data.panel=1, data = HG38_CENPA_over5_chr9.bdg, r0=0.55, r1=0.65, window.size = 250, col= "#DCDCDC", border = "#DCDCDC")
kp <- kpPlotDensity(kp, data.panel=1, data = HG38_CENPA_mapQ20_over5_chr9.bdg, r0=0.55, r1=0.65, window.size = 50, col= "#ecde38", border = "#ecde38")
kpAxis(kp, data.panel=1, ymax=kp$latest.plot$computed.values$max.density, tick.len = 0.02e6, r0=0.55, r1=0.65, cex=0.25, numticks = 2)

#CENPA RPE1 HAP1 diploid alignment
kp <- kpPlotDensity(kp, data.panel=1, data = RPE1v1.0dip_CENPA_hap1_over5_chr9.bdg, r0=0.66, r1=0.76, window.size = 250, col= "#DCDCDC", border = "#DCDCDC")
kp <- kpPlotDensity(kp, data.panel=1, data = RPE1v1.0dip_CENPAmapQ20_hap1_over5_chr9.bdg, r0=0.66, r1=0.76, window.size = 50, col= "#f15a22", border = "#f15a22")
kpAxis(kp, data.panel=1, ymax=kp$latest.plot$computed.values$max.density, tick.len = 0.02e6, r0=0.66, r1=0.76, cex=0.25, numticks = 2)

#5mC RPE1 HAP1 diploid alignment
kp <- kpPlotDensity(kp, data.panel=1, data = meth1 ,r0=0.77, r1=0.87, window.size = 5000, border = "#000000", col= "#000000")
kpAxis(kp, data.panel=1, ymax=kp$latest.plot$computed.values$max.density, tick.len = 0.02e6, r0=0.77, r1=0.87, cex=0.25, numticks = 2)

dev.off()
