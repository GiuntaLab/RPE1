####SCRIPT FOR COMPUTING THE AVERAGE log10(qvalue) FOR EACH CDR AND PROTEIN (both from q20 filtered and unfiltered bdg files)

library(dplyr)
library(ggplot2)
library(gridExtra)
library(data.table)
library(ComplexHeatmap)

#log10(qvalue) average for all proteins from unfiltered alignment files####

#CDR annotation for both haplotypes of RPE-1 genome
data_cdr <- fread("CDR_RPE1.bed")
data_cdr <- data_cdr[,c(1,4)]
data_cdr$values <- NA
colnames(data_cdr) <- c("chr", "cdr", "values")

#CENP-A bedgraph from unfiltered BAM files intersected with CDRs
#This file processing was applied to all tables of all proteins, both from filtered and unfiltered alignments.
#After the averaging and normalization (min-max scaling) process, the tables of different proteins are merged and plotted as a heatmap.
data_cenpa <- fread("CENPA_rep1_rpe1_rpe1_comp_qval.bdg.CDRs")
data_cenpa$protein <- "cenpa"
data_f <- data_cenpa %>%
  group_by(V1, V8, protein) %>%
  mutate(qvalue_mean = mean(V4)) %>% #computing average of q-values for each chromosome, protein, and cdr
  ungroup() %>%
  distinct(V1, V8, protein, .keep_all = TRUE) %>%
  dplyr::select(V1, V8, protein, qvalue_mean)

colnames(data_f) <- c("chr", "cdr", "protein", "qvalue_mean")

#Reorder columns for merging
data_f_subset <- data_f %>%
  dplyr::select(chr, cdr, qvalue_mean, protein)

#Merge CENP-A q-values with CDRs table
merged <- data_cdr %>%
  left_join(data_f_subset, by = c("chr", "cdr"))

#Replace missing q-values with 0 and fill missing protein labels
merged_cenpa <- merged %>%
  mutate(qvalue_mean = replace_na(qvalue_mean, 0)) %>%
  mutate(protein = replace_na(protein, "cenpa"))

#Scale q-values per protein to 0-1 range for visualization
#Min-max scaling formula
merged_cenpa <- merged_cenpa %>%
  mutate(qvalue_minmax = (qvalue_mean - min(qvalue_mean, na.rm = TRUE)) /
           (max(qvalue_mean, na.rm = TRUE) - min(qvalue_mean, na.rm = TRUE))) %>%
  dplyr::select(chr, protein, cdr, qvalue_minmax)

#Convert table to wide format: CDRs as columns, q-value minmax as values
merged_cenpa <- merged_cenpa %>%
  pivot_wider(names_from = cdr, values_from = qvalue_minmax)

#CENP-B bedgraph from unfiltered BAM files intersected with CDRs
data_cenpb <- fread("CENPB_merge_rpe1_rpe1_comp_qval.bdg.CDRs")
data_cenpb$protein <- "cenpb"
data_f <- data_cenpb %>%
  group_by(V1, V8, protein) %>%
  mutate(qvalue_mean = mean(V4)) %>%
  ungroup() %>%
  distinct(V1, V8, protein, .keep_all = TRUE) %>%
  dplyr::select(V1, V8, protein, qvalue_mean)

colnames(data_f) <- c("chr", "cdr", "protein", "qvalue_mean")

data_f_subset <- data_f %>%
  dplyr::select(chr, cdr, qvalue_mean, protein)

merged <- data_cdr %>%
  left_join(data_f_subset, by = c("chr", "cdr"))

merged_cenpb <- merged %>%
  mutate(qvalue_mean = replace_na(qvalue_mean, 0)) %>%
  mutate(protein = replace_na(protein, "cenpb"))

merged_cenpb <- merged_cenpb %>%
  mutate(qvalue_minmax = (qvalue_mean - min(qvalue_mean, na.rm = TRUE)) /
           (max(qvalue_mean, na.rm = TRUE) - min(qvalue_mean, na.rm = TRUE))) %>%
  dplyr::select(chr, protein, cdr, qvalue_minmax)

merged_cenpb <- merged_cenpb %>%
  pivot_wider(names_from = cdr, values_from = qvalue_minmax)

#CENP-C bedgraph from unfiltered BAM files intersected with CDRs
data_cenpc <- fread("CENPC_merge_rpe1_rpe1_comp_qval.bdg.CDRs")
data_cenpc$protein <- "cenpc"
data_f <- data_cenpc %>%
  group_by(V1, V8, protein) %>%
  mutate(qvalue_mean = mean(V4)) %>%
  ungroup() %>%
  distinct(V1, V8, protein, .keep_all = TRUE) %>%
  dplyr::select(V1, V8, protein, qvalue_mean)

colnames(data_f) <- c("chr", "cdr", "protein", "qvalue_mean")

data_f_subset <- data_f %>%
  dplyr::select(chr, cdr, qvalue_mean, protein)

merged <- data_cdr %>%
  left_join(data_f_subset, by = c("chr", "cdr"))

merged_cenpc <- merged %>%
  mutate(qvalue_mean = replace_na(qvalue_mean, 0)) %>%
  mutate(protein = replace_na(protein, "cenpc"))

merged_cenpc <- merged_cenpc %>%
  mutate(qvalue_minmax = (qvalue_mean - min(qvalue_mean, na.rm = TRUE)) /
           (max(qvalue_mean, na.rm = TRUE) - min(qvalue_mean, na.rm = TRUE))) %>%
  dplyr::select(chr, protein, cdr, qvalue_minmax)

merged_cenpc <- merged_cenpc %>%
  pivot_wider(names_from = cdr, values_from = qvalue_minmax)

#CENP-T bedgraph from unfiltered BAM files intersected with CDRs
data_cenpt <- fread("CENPT_untr_rpe1_rpe1_comp_qval.bdg.CDRs")
data_cenpt$protein <- "cenpt"
data_f <- data_cenpt %>%
  group_by(V1, V8, protein) %>%
  mutate(qvalue_mean = mean(V4)) %>%
  ungroup() %>%
  distinct(V1, V8, protein, .keep_all = TRUE) %>%
  dplyr::select(V1, V8, protein, qvalue_mean)

colnames(data_f) <- c("chr", "cdr", "protein", "qvalue_mean")

data_f_subset <- data_f %>%
  dplyr::select(chr, cdr, qvalue_mean, protein)

merged <- data_cdr %>%
  left_join(data_f_subset, by = c("chr", "cdr"))

merged_cenpt <- merged %>%
  mutate(qvalue_mean = replace_na(qvalue_mean, 0)) %>%
  mutate(protein = replace_na(protein, "cenpt"))

merged_cenpt <- merged_cenpt %>%
  mutate(qvalue_minmax = (qvalue_mean - min(qvalue_mean, na.rm = TRUE)) /
           (max(qvalue_mean, na.rm = TRUE) - min(qvalue_mean, na.rm = TRUE))) %>%
  dplyr::select(chr, protein, cdr, qvalue_minmax)

merged_cenpt <- merged_cenpt %>%
  pivot_wider(names_from = cdr, values_from = qvalue_minmax)

#Combine all proteins' CDR occupancy data into a single data frame
merged_final <- bind_rows(merged_cenpa, merged_cenpb, merged_cenpc, merged_cenpt)

# 1. Pivot to long format: one row per chr/protein/CDR
long_df <- merged_final %>%
  pivot_longer(cols = starts_with("cdr_"), names_to = "cdr", values_to = "occupancy")  #all columns representing CDRs, name of new column storing CDR identifiers, name of new column storing values

# 2. Desired chromosome order (interleaved hap1/hap2)
chromosomes <- c(paste0("chr", 1:22, "_hap1"), paste0("chr", 1:22, "_hap2"))
chromosomes <- as.vector(rbind(chromosomes[1:22], chromosomes[23:44]))
chromosomes <- c(chromosomes, "chrX_hap1", "chrX_hap2")

# 3. Get unique proteins (columns)
proteins <- sort(unique(long_df$protein))

# 4. Split long_df into list of data.frames by CDR, where each list element corresponds to one CDR
cdr_list <- split(long_df, long_df$cdr)

 5. Define color scale for heatmaps
col_fun <- colorRamp2(c(0, 0.5, 1), c("navy","white","firebrick"))

# 6. Build heatmap list (1 per CDR)
heatmaps3 <- lapply(names(cdr_list), function(cdr_name) {
  df <- cdr_list[[cdr_name]]

  # Create matrix: rows = chr, columns = protein
  mat <- df %>%
    dplyr::select(chr, protein, occupancy) %>%
    pivot_wider(names_from = protein, values_from = occupancy) %>%
    complete(chr = chromosomes) %>%  #ensure all chromosomes are included
    arrange(factor(chr, levels = chromosomes)) %>% #enforce desired row order
    column_to_rownames("chr") %>% #row names = chromosomes
    as.matrix()

  # Make sure all proteins are present as columns (add missing ones as NA)
  for (prot in proteins) {
    if (!prot %in% colnames(mat)) {
      mat <- cbind(mat, setNames(rep(NA, nrow(mat)), prot))
    }
  }
  mat <- mat[, proteins, drop = FALSE]  # order columns

  Heatmap(
    mat,
    name = paste0("Occupancy_", cdr_name),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = chromosomes,
    column_title = toupper(cdr_name), # CDR name as column title
    col = col_fun,
    na_col = "grey",
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    border = TRUE,
    heatmap_legend_param = list(title = "Scaled -log10(qvalue)"),
  )
})

# 7. Combine individual CDR heatmaps into a single heatmap object
p3 <- Reduce(`+`, heatmaps3)

#log10(qvalue) average for all proteins from filtered (MAPQ 20) alignment files####

data_cdr <- fread("CDR_RPE1.bed")
data_cdr <- data_cdr[,c(1,4)]
data_cdr$values <- NA
colnames(data_cdr) <- c("chr", "cdr", "values")

#CENP-A bedgraph from filtered BAM (MAPQ20) files intersected with CDRs
data_cenpa <- fread("CENPA_rep1_rpe1_rpe1_q20_comp_qval.bdg.CDRs")
data_cenpa$protein <- "cenpa"
data_f <- data_cenpa %>%
  group_by(V1, V8, protein) %>%
  mutate(qvalue_mean = mean(V4)) %>%
  ungroup() %>%
  distinct(V1, V8, protein, .keep_all = TRUE) %>%
  dplyr::select(V1, V8, protein, qvalue_mean)

colnames(data_f) <- c("chr", "cdr", "protein", "qvalue_mean")

data_f_subset <- data_f %>%
  dplyr::select(chr, cdr, qvalue_mean, protein)

merged <- data_cdr %>%
  left_join(data_f_subset, by = c("chr", "cdr"))

merged_cenpa <- merged %>%
  mutate(qvalue_mean = replace_na(qvalue_mean, 0)) %>%
  mutate(protein = replace_na(protein, "cenpa"))

merged_cenpa <- merged_cenpa %>%
  mutate(qvalue_minmax = (qvalue_mean - min(qvalue_mean, na.rm = TRUE)) /
           (max(qvalue_mean, na.rm = TRUE) - min(qvalue_mean, na.rm = TRUE))) %>%
  dplyr::select(chr, protein, cdr, qvalue_minmax)

merged_cenpa <- merged_cenpa %>%
  pivot_wider(names_from = cdr, values_from = qvalue_minmax)

#CENP-B bedgraph from filtered BAM (MAPQ20) files intersected with CDRs
data_cenpb <- fread("CENPB_merge_rpe1_rpe1_q20_comp_qval.bdg.CDRs")
data_cenpb$protein <- "cenpb"
data_f <- data_cenpb %>%
  group_by(V1, V8, protein) %>%
  mutate(qvalue_mean = mean(V4)) %>%
  ungroup() %>%
  distinct(V1, V8, protein, .keep_all = TRUE) %>%
  dplyr::select(V1, V8, protein, qvalue_mean)

colnames(data_f) <- c("chr", "cdr", "protein", "qvalue_mean")

data_f_subset <- data_f %>%
  dplyr::select(chr, cdr, qvalue_mean, protein)

merged <- data_cdr %>%
  left_join(data_f_subset, by = c("chr", "cdr"))

merged_cenpb <- merged %>%
  mutate(qvalue_mean = replace_na(qvalue_mean, 0)) %>%
  mutate(protein = replace_na(protein, "cenpb"))

merged_cenpb <- merged_cenpb %>%
  mutate(qvalue_minmax = (qvalue_mean - min(qvalue_mean, na.rm = TRUE)) /
           (max(qvalue_mean, na.rm = TRUE) - min(qvalue_mean, na.rm = TRUE))) %>%
  dplyr::select(chr, protein, cdr, qvalue_minmax)

merged_cenpb <- merged_cenpb %>%
  pivot_wider(names_from = cdr, values_from = qvalue_minmax)

#CENP-C bedgraph from filtered BAM (MAPQ20) files intersected with CDRs
data_cenpc <- fread("CENPC_merge_rpe1_rpe1_q20_comp_qval.bdg.CDRs")
data_cenpc$protein <- "cenpc"
data_f <- data_cenpc %>%
  group_by(V1, V8, protein) %>%
  mutate(qvalue_mean = mean(V4)) %>%
  ungroup() %>%
  distinct(V1, V8, protein, .keep_all = TRUE) %>%
  dplyr::select(V1, V8, protein, qvalue_mean)

colnames(data_f) <- c("chr", "cdr", "protein", "qvalue_mean")

data_f_subset <- data_f %>%
  dplyr::select(chr, cdr, qvalue_mean, protein)

merged <- data_cdr %>%
  left_join(data_f_subset, by = c("chr", "cdr"))

merged_cenpc <- merged %>%
  mutate(qvalue_mean = replace_na(qvalue_mean, 0)) %>%
  mutate(protein = replace_na(protein, "cenpc"))

merged_cenpc <- merged_cenpc %>%
  mutate(qvalue_minmax = (qvalue_mean - min(qvalue_mean, na.rm = TRUE)) /
           (max(qvalue_mean, na.rm = TRUE) - min(qvalue_mean, na.rm = TRUE))) %>%
  dplyr::select(chr, protein, cdr, qvalue_minmax)

merged_cenpc <- merged_cenpc %>%
  pivot_wider(names_from = cdr, values_from = qvalue_minmax)

#CENP-T bedgraph from filtered BAM (MAPQ20) files intersected with CDRs
data_cenpt <- fread("CENPT_untr_rpe1_rpe1_q20_comp_qval.bdg.CDRs")
data_cenpt$protein <- "cenpt"
data_f <- data_cenpt %>%
  group_by(V1, V8, protein) %>%
  mutate(qvalue_mean = mean(V4)) %>%
  ungroup() %>%
  distinct(V1, V8, protein, .keep_all = TRUE) %>%
  dplyr::select(V1, V8, protein, qvalue_mean)

colnames(data_f) <- c("chr", "cdr", "protein", "qvalue_mean")

data_f_subset <- data_f %>%
  dplyr::select(chr, cdr, qvalue_mean, protein)

merged <- data_cdr %>%
  left_join(data_f_subset, by = c("chr", "cdr"))

merged_cenpt <- merged %>%
  mutate(qvalue_mean = replace_na(qvalue_mean, 0)) %>%
  mutate(protein = replace_na(protein, "cenpt"))

merged_cenpt <- merged_cenpt %>%
  mutate(qvalue_minmax = (qvalue_mean - min(qvalue_mean, na.rm = TRUE)) /
           (max(qvalue_mean, na.rm = TRUE) - min(qvalue_mean, na.rm = TRUE))) %>%
  dplyr::select(chr, protein, cdr, qvalue_minmax)

merged_cenpt <- merged_cenpt %>%
  pivot_wider(names_from = cdr, values_from = qvalue_minmax)


merged_final <- bind_rows(merged_cenpa, merged_cenpb, merged_cenpc, merged_cenpt)

# 1. Pivot to long format: one row per chr/protein/CDR
long_df <- merged_final %>%
  pivot_longer(cols = starts_with("cdr_"), names_to = "cdr", values_to = "occupancy")

# 2. Desired chromosome order (interleaved hap1/hap2)
chromosomes <- c(paste0("chr", 1:22, "_hap1"), paste0("chr", 1:22, "_hap2"))
chromosomes <- as.vector(rbind(chromosomes[1:22], chromosomes[23:44]))
chromosomes <- c(chromosomes, "chrX_hap1", "chrX_hap2")

# 3. Get unique proteins (columns)
proteins <- sort(unique(long_df$protein))

# 4. Split long_df into list of data.frames by CDR
cdr_list <- split(long_df, long_df$cdr)

col_fun <- colorRamp2(c(0, 0.5, 1), c("navy","white","firebrick"))
# 6. Build heatmap list (1 per CDR)
heatmaps4 <- lapply(names(cdr_list), function(cdr_name) {
  df <- cdr_list[[cdr_name]]

  # Create matrix: rows = chr, columns = protein
  mat <- df %>%
    dplyr::select(chr, protein, occupancy) %>%
    pivot_wider(names_from = protein, values_from = occupancy) %>%
    complete(chr = chromosomes) %>%  # Ensure all chromosomes
    arrange(factor(chr, levels = chromosomes)) %>%
    column_to_rownames("chr") %>%
    as.matrix()

  # Make sure all proteins are present
  for (prot in proteins) {
    if (!prot %in% colnames(mat)) {
      mat <- cbind(mat, setNames(rep(NA, nrow(mat)), prot))
    }
  }
  mat <- mat[, proteins, drop = FALSE]  # order proteins

  Heatmap(
    mat,
    name = paste0("Occupancy_", cdr_name),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = chromosomes,
    column_title = toupper(cdr_name),
    col = col_fun,
    na_col = "grey",
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    border = TRUE,
    heatmap_legend_param = list(title = "Scaled -log10(qvalue)")
  )
})

p4 <- Reduce(`+`, heatmaps4)

# Define a helper function to convert a ComplexHeatmap object into a grob
# This allows further manipulation with grid graphics or inclusion in ggplot2 layouts
heatmap_to_grob <- function(ht_list, title) {
  grid.grabExpr(draw(
    ht_list,
    heatmap_legend_side = "right", # position legend to the right
    column_title = title,  # add a title above the heatmap
    column_title_gp = gpar(fontsize = 14, fontface = "bold") # customize title font
  ))
}

#Convert heatmaps into a grab for plotting
g3 <- heatmap_to_grob(p3, "q-value distribution across CDRs")
g4 <- heatmap_to_grob(p4, "q-value distribution from MAPQ20 across CDRs")

#Plotting and exporting the two heatmaps in the same jpeg file
jpeg("qvalue_combined_heatmaps.jpeg",
     width = 16,          # inches
     height = 12,        # inches
     units = "in",
     quality = 100,
     res = 300)          # 300 dpi
grid.arrange(g3, g4, ncol = 2)
dev.off()