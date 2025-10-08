####COMPUTING PERCENTAGE IDENTITY FROM PAF FILE (MINIMAP ALIGNMENT) AND PLOTTING THE RESULT####
#minimap2 -x asm20 -c --cs=long -t 10 "$hap2_genome" $read_hap1 > "hap2_genome_$n.paf"

#Loading libraries
library(ggplot2)
library(dplyr)

read_paf_with_pid <- function(paf_file) {
  #Read PAF file considering only the first 12 columns
  paf <- read.table(paf_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, colClasses = c(rep("character", 12)))[, 1:12]
  
  #Check
  if (ncol(paf) < 11) stop("The file ", paf_file, " doesn't have enough columns")
  
  #Computing percentage identity (Number of matching bases / Alignment block length)×100.
  paf$V10 <- as.numeric(paf$V10)
  paf$V11 <- as.numeric(paf$V11)
  paf$percent_identity <- (paf$V10 / paf$V11) * 100
  
  #Fixing names
  full_name <- basename(paf_file)
  if (grepl("broadP", full_name)) {
    sample_name <- sub("(.*broadP).*", "\\1", full_name)
  } else if (grepl("nobroadP", full_name)) {
    sample_name <- sub("(.*nobroadP).*", "\\1", full_name)
  } else {
    sample_name <- tools::file_path_sans_ext(full_name)
  }
  paf$sample <- sample_name
  
  return(paf[, c("percent_identity", "sample")])
}

#My PAF files
paf_files <- c(
  "hap1_genome_broadP.500bp.f025.AS-HOR-vs-RPE1v1.1.onlylive.hap2.startend.bed.fasta.paf",
  "hap1_genome_nobroadP.500bp.f025.AS-HOR-vs-RPE1v1.1.onlylive.hap2.startend.bed.fasta.paf",
  "hap2_genome_broadP.500bp.f025.AS-HOR-vs-RPE1v1.1.onlylive.hap1.startend.bed.fasta.paf",
  "hap2_genome_nobroadP.500bp.f025.AS-HOR-vs-RPE1v1.1.onlylive.hap1.startend.bed.fasta.paf"
)

#Read all files
all_data <- do.call(rbind, lapply(paf_files, read_paf_with_pid))

#Plot with outlier
p <- ggplot(all_data, aes(x = sample, y = percent_identity, fill = sample)) +
  geom_boxplot(outlier.alpha = 0.2, width = 0.6) +
  theme_minimal(base_size = 14) +
  labs(title = "Percent Identity per Mapping", x = "Sample (PAF file)", y = "Percent Identity (%)") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2") +
  theme_linedraw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

jpeg("PAF_percent_identity_boxplot.jpeg", width = 2000, height = 1500, res = 300)
print(p)
dev.off()

#Plot without outlier
p_no_outlier <- ggplot(all_data, aes(x = sample, y = percent_identity, fill = sample)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) + # outlier.shape=NA rimuove gli outlier dal plot
  theme_minimal(base_size = 14) + ylim(90,100) +
  labs(title = "Percent Identity per Mapping (No Outliers)", x = "Sample (PAF file)", y = "Percent Identity (%)") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2") +
  theme_linedraw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

jpeg("PAF_percent_identity_boxplot_nooutlier.jpeg", width = 2000, height = 1500, res = 300)
print(p_no_outlier)
dev.off()

#Comparison 1: hap1 broadP vs noBroadP
t.test(percent_identity ~ sample, 
       data = subset(all_data, sample %in% c("hap1_genome_broadP", "hap1_genome_nobroadP")))

#Comparison 2: hap2 broadP vs noBroadP
t.test(percent_identity ~ sample, 
       data = subset(all_data, sample %in% c("hap2_genome_broadP", "hap2_genome_nobroadP")))
