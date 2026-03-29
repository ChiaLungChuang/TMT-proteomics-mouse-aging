# ============================================================
# Project: TMT Proteomics — Multi-organ Mouse Aging
# Script:  insolubility_protein_half-life.R
# Purpose: Correlates protein half-life with age-related
#          insolubility fold change in 4 tissues (Heart,
#          Cerebellum, Neocortex, Liver); generates faceted
#          and non-faceted scatter plots with optional labels
#          and R² annotations
# Input:   protein half-life.xlsx,
#          combined_plot_data_mouseTMT_<tissue>.csv
# Output:  mouseTMT_<tissue>_Halflife_*.pdf,
#          combined_half_<tissue>.csv
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# NOTE:    Half-life data available for 4 tissues only:
#          Heart, Cerebellum, Neocortex, Liver
# ============================================================

library(tidyverse)
library(ggrepel)
library(readxl)
library(ggpmisc)

# NOTE: Set your working directory to the folder containing input files
# setwd("path/to/your/data")

# ---- 1. Load Data ----

halflife_data <- read_excel("protein half-life.xlsx")
write.csv(halflife_data, "protein half-life.csv", row.names = FALSE)
halflife_data <- read.csv("protein half-life.csv")

# Load per-tissue combined plot data
combined_plot_data_He <- read.csv("combined_plot_data_mouseTMT_He.csv")
combined_plot_data_CB <- read.csv("combined_plot_data_mouseTMT_CB.csv")
combined_plot_data_NC <- read.csv("combined_plot_data_mouseTMT_NC.csv")
combined_plot_data_Li <- read.csv("combined_plot_data_mouseTMT_Li.csv")

# ---- 2. Extract UniProt ID from GN_ID for joining ----
# GN_ID format contains a UniProt accession between pipe characters (e.g. |P12345|)
# Half-life data uses UniProt accession as the key (column: Proteins)

extract_uniprot_id <- function(df) {
  df %>% mutate(ID = str_extract(GN_ID, "(?<=\\|)[A-Z0-9]+(\\-\\d+)?(?=\\|)"))
}

combined_plot_data_He <- extract_uniprot_id(combined_plot_data_He)
combined_plot_data_CB <- extract_uniprot_id(combined_plot_data_CB)
combined_plot_data_NC <- extract_uniprot_id(combined_plot_data_NC)
combined_plot_data_Li <- extract_uniprot_id(combined_plot_data_Li)

# ---- 3. Join Half-Life Data (tissue-specific columns) ----

combined_half_He <- combined_plot_data_He %>%
  left_join(halflife_data %>%
              select(Proteins, protein.half.life..heart.) %>%
              rename(ID = Proteins, Protein_Half_Life = protein.half.life..heart.),
            by = "ID")

combined_half_CB <- combined_plot_data_CB %>%
  left_join(halflife_data %>%
              select(Proteins, protein.half.life..cerebellum.) %>%
              rename(ID = Proteins, Protein_Half_Life = protein.half.life..cerebellum.),
            by = "ID", relationship = "many-to-many")

combined_half_NC <- combined_plot_data_NC %>%
  left_join(halflife_data %>%
              select(Proteins, protein.half.life..neoortex.) %>%
              rename(ID = Proteins, Protein_Half_Life = protein.half.life..neoortex.),
            by = "ID", relationship = "many-to-many")

combined_half_Li <- combined_plot_data_Li %>%
  left_join(halflife_data %>%
              select(Proteins, protein.half.life..liver.) %>%
              rename(ID = Proteins, Protein_Half_Life = protein.half.life..liver.),
            by = "ID", relationship = "many-to-many")

# Save joined data
write.csv(combined_half_He, "combined_half_He.csv", row.names = FALSE)
write.csv(combined_half_CB, "combined_half_CB.csv", row.names = FALSE)
write.csv(combined_half_NC, "combined_half_NC.csv", row.names = FALSE)
write.csv(combined_half_Li, "combined_half_Li.csv", row.names = FALSE)

# ---- 4. Tissue-Specific Label Thresholds ----
# Filter NA genes and label proteins with extreme half-life values

combined_half_CB <- combined_half_CB %>%
  filter(!is.na(Gene_Name)) %>%
  mutate(Label_half = case_when(
    Comparison == "Old vs. Middle" & (Protein_Half_Life > 5)   ~ Gene_Name,
    Comparison != "Old vs. Middle" & (Protein_Half_Life > 6)   ~ Gene_Name,
    Protein_Half_Life < 0                                       ~ Gene_Name,
    TRUE ~ NA_character_))

combined_half_He <- combined_half_He %>%
  filter(!is.na(Gene_Name)) %>%
  mutate(Label_half = case_when(
    Comparison == "Old vs. Young" & (Protein_Half_Life > 4.5) ~ Gene_Name,
    Comparison != "Old vs. Young" & (Protein_Half_Life > 4)   ~ Gene_Name,
    Protein_Half_Life < 0.2                                    ~ Gene_Name,
    TRUE ~ NA_character_))

combined_half_NC <- combined_half_NC %>%
  filter(!is.na(Gene_Name)) %>%
  mutate(Label_half = case_when(
    Comparison == "Old vs. Middle" & (Protein_Half_Life > 4 | Protein_Half_Life < 1) ~ Gene_Name,
    Comparison != "Old vs. Middle" & Protein_Half_Life > 5                           ~ Gene_Name,
    Protein_Half_Life < 0.2                                                           ~ Gene_Name,
    TRUE ~ NA_character_))

combined_half_Li <- combined_half_Li %>%
  filter(!is.na(Gene_Name)) %>%
  mutate(Label_half = case_when(
    Comparison == "Old vs. Middle" & (Protein_Half_Life > 5 | Protein_Half_Life < 0) ~ Gene_Name,
    Comparison != "Old vs. Middle" & Protein_Half_Life > 5                           ~ Gene_Name,
    Protein_Half_Life < -0.2                                                          ~ Gene_Name,
    TRUE ~ NA_character_))

# ---- 5. Shared Plot Elements ----

half_theme <- theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12, face = "italic"),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title   = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text   = element_text(size = 14, face = "bold")
  )

comparison_colors <- scale_color_manual(
  values = c("Old vs. Young" = "purple4", "Middle vs. Young" = "maroon2", "Old vs. Middle" = "darkgreen")
)

size_scale <- scale_size_continuous(
  range  = c(1, 9), name = "P-Value Size (-log10)",
  breaks = c(1.3, 2, 3, 4, 5), limit = c(1.3, 6),
  labels = c("1.30(p<0.05)", "2", "3", "4", "5")
)

# ---- 6. Plot Function ----
# Generates 6 plot variants per tissue:
#   facet_nolabel, facet_label, facet_nolabel_R,
#   nolabel, label, nolabel_R

plot_halflife_tissue <- function(data, tissue_name, file_prefix) {
  labs_base <- labs(
    title   = paste0("MouseTMT Protein Half-Life vs. Insolubility Fold Change in ", tissue_name),
    x       = "Protein Half-Life",
    y       = "Log2 fold change in Insoluble/Soluble",
    caption = "Point size based on statistical significance of Age difference in insolubility"
  )
  
  # Faceted — no labels
  p_f_nolabel <- ggplot(data, aes(x = Protein_Half_Life, y = Ratio_FC, size = Size, color = Comparison)) +
    geom_point(alpha = 0.4) + size_scale + comparison_colors +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    labs_base + half_theme + facet_wrap(~Comparison)
  
  print(p_f_nolabel)
  ggsave(paste0(file_prefix, "_facet_nolabel.pdf"), width = 7200, height = 6400, units = "px", dpi = 600)
  
  # Faceted — with labels
  p_f_label <- ggplot(data, aes(x = Protein_Half_Life, y = Ratio_FC, size = Size,
                                label = Label_half, color = Comparison)) +
    geom_point(alpha = 0.4) +
    geom_text_repel(size = 3, segment.color = "grey50", max.overlaps = Inf,
                    box.padding = 0.5, point.padding = 0.2, min.segment.length = 0.1) +
    size_scale + comparison_colors +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    labs_base + half_theme + facet_wrap(~Comparison) +
    guides(color = guide_legend(override.aes = list(label = "")))
  
  print(p_f_label)
  ggsave(paste0(file_prefix, "_facet_label.pdf"), width = 7200, height = 6400, units = "px", dpi = 600)
  
  # Faceted — no labels + R²
  p_f_nolabel_R <- p_f_nolabel +
    stat_poly_eq(aes(label = paste(after_stat(rr.label))), rr.digits = 3,
                 label.x = 1, label.y = 1.4, parse = TRUE, color = "black", show.legend = FALSE)
  
  print(p_f_nolabel_R)
  ggsave(paste0(file_prefix, "_facet_nolabel_R.pdf"), width = 7200, height = 6400, units = "px", dpi = 600)
  
  # No facet — no labels
  p_nolabel <- ggplot(data, aes(x = Protein_Half_Life, y = Ratio_FC, size = Size, color = Comparison)) +
    geom_point(alpha = 0.4) + size_scale + comparison_colors +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    labs_base + half_theme
  
  print(p_nolabel)
  ggsave(paste0(file_prefix, "_nolabel.pdf"), width = 7200, height = 6400, units = "px", dpi = 600)
  
  # No facet — with labels
  p_label <- ggplot(data, aes(x = Protein_Half_Life, y = Ratio_FC, size = Size,
                              label = Label_half, color = Comparison)) +
    geom_point(alpha = 0.4) +
    geom_text_repel(size = 3, segment.color = "grey50", max.overlaps = Inf,
                    box.padding = 0.5, point.padding = 0.2, min.segment.length = 0.1) +
    size_scale + comparison_colors +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    labs_base + half_theme +
    guides(color = guide_legend(override.aes = list(label = "")))
  
  print(p_label)
  ggsave(paste0(file_prefix, "_label.pdf"), width = 7200, height = 6400, units = "px", dpi = 600)
  
  # No facet — no labels + R²
  p_nolabel_R <- p_nolabel +
    stat_poly_eq(aes(label = paste(after_stat(rr.label))), rr.digits = 3,
                 label.x = 1, label.y = 1.4, parse = TRUE, color = "black", show.legend = FALSE)
  
  print(p_nolabel_R)
  ggsave(paste0(file_prefix, "_nolabel_R.pdf"), width = 7200, height = 6400, units = "px", dpi = 600)
}

# ---- 7. Generate Per-Tissue Plots ----

plot_halflife_tissue(combined_half_CB, "Cerebellum", "mouseTMT_CB_Halflife")
plot_halflife_tissue(combined_half_He, "Heart",      "mouseTMT_He_Halflife")
plot_halflife_tissue(combined_half_NC, "Neocortex",  "mouseTMT_NC_Halflife")
plot_halflife_tissue(combined_half_Li, "Liver",      "mouseTMT_Li_Halflife")