#PLOTTING PERCENTAGE IDENTITY FROM BAM FILES (BOWTIE2) COMPUTED USING 'p_identity_bam.py' script ####

library(ggplot2)
library(dplyr)

#List TSV files
tsv_files <- c(
  "hap1_genome_broadP.500bp.f025.AS-HOR-vs-RPE1v1.1.onlylive.hap2.startend.bed.fq_bowtie.bam.tsv",
  "hap1_genome_nobroadP.500bp.f025.AS-HOR-vs-RPE1v1.1.onlylive.hap2.startend.bed.fq_bowtie.bam.tsv",
  "hap2_genome_broadP.500bp.f025.AS-HOR-vs-RPE1v1.1.onlylive.hap1.startend.bed.fq_bowtie.bam.tsv",
  "hap2_genome_nobroadP.500bp.f025.AS-HOR-vs-RPE1v1.1.onlylive.hap1.startend.bed.fq_bowtie.bam.tsv"
)

#Read and combine all TSV files into one dataframe
all_data <- do.call(rbind, lapply(tsv_files, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  #Fixing names
  full_name <- basename(f)
  
  #Simplify sample names
  if (grepl("broadP", full_name)) {
    sample_name <- sub("(.*broadP).*", "\\1", full_name)
  } else if (grepl("nobroadP", full_name)) {
    sample_name <- sub("(.*nobroadP).*", "\\1", full_name)
  } else {
    sample_name <- tools::file_path_sans_ext(full_name)
  }
  
  df$sample <- sample_name
  return(df)
}))

#Plot with outliers
p <- ggplot(all_data, aes(x = sample, y = percent_identity, fill = sample)) +
  geom_boxplot(outlier.alpha = 0.2, width = 0.6) +
  theme_minimal(base_size = 14) +
  labs(title = "Percent Identity per Mapping", x = "Sample (BAM)", y = "Percent Identity (%)") +
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

jpeg("BAM_percent_identity_boxplot.jpeg", width = 2000, height = 1500, res = 300)
print(p)
dev.off()

#Plot without outliers
p_no_outlier <- ggplot(all_data, aes(x = sample, y = percent_identity, fill = sample)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  theme_minimal(base_size = 14) +
  ylim(90, 100) +
  labs(title = "Percent Identity per Mapping (No Outliers)", x = "Sample (BAM)", y = "Percent Identity (%)") +
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

jpeg("BAM_percent_identity_boxplot_nooutlier.jpeg", width = 2000, height = 1500, res = 300)
print(p_no_outlier)
dev.off()

#Comparison 1: hap1 broadP vs noBroadP
t.test(percent_identity ~ sample, 
       data = subset(all_data, sample %in% c("hap1_genome_broadP", "hap1_genome_nobroadP")))

#Comparison 2: hap2 broadP vs noBroadP
t.test(percent_identity ~ sample, 
       data = subset(all_data, sample %in% c("hap2_genome_broadP", "hap2_genome_nobroadP")))