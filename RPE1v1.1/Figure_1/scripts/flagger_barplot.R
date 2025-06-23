library(ggplot2)
library(dplyr)
library(readr)
library(stringr)


read_bed_lengths <- function(path, category) {
  df <- read_tsv(path, col_names = FALSE, comment = "#", col_types = cols(
    X1 = col_character(), X2 = col_integer(), X3 = col_integer()
  )) %>%
    mutate(length = X3 - X2,
           Group = case_when(
             str_ends(X1, "_hap1") ~ "hap1",
             str_ends(X1, "_hap2") ~ "hap2",
             TRUE ~ NA_character_
           ),
           Chromosome = str_replace(X1, "_hap[12]$", "")) %>%
    filter(!is.na(Group)) %>%
    group_by(Group, Chromosome) %>%
    summarise(value = sum(length) / 1e6, .groups = "drop") %>%
    mutate(Category = category)
  return(df)
}


error_df <- read_bed_lengths("rpe1.v1.1.error.bed", "Error")
collapsed_df <- read_bed_lengths("rpe1.v1.1.collapsed.bed", "Collapsed")
haploid_df <- read_bed_lengths("rpe1.v1.1.haploid.bed", "Haploid")


all_df <- bind_rows(error_df, collapsed_df, haploid_df)


color_map <- c(
  "Error" = rgb(255, 95, 31, maxColorValue = 255),
  "Collapsed" = rgb(69, 99, 170, maxColorValue = 255),
  "Haploid" = rgb(9, 157, 54, maxColorValue = 255)
)


save_group_plot <- function(df, group_name, output_file) {
  p <- df %>%
    filter(Group == group_name) %>%
    ggplot(aes(x = Chromosome, y = value, fill = Category)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = color_map) +
    labs(title = paste("Regions per Chromosome (", group_name, ")", sep = ""),
         y = "Length (Mb)", x = "Chromosome") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    )

  ggsave(output_file, plot = p, width = 12, height = 5, dpi = 300)
}



save_group_plot(all_df, "hap1", "error_collapsed_haploid_hap1.png")
save_group_plot(all_df, "hap2", "error_collapsed_haploid_hap2.png")

