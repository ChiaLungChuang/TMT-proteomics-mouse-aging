# ============================================================
# Project: TMT Proteomics — Multi-organ Mouse Aging
# Script:  Insolubility_meltTemp.R
# Purpose: Correlates protein melting temperature (thermal
#          stability) with age-related insolubility fold change
#          across 8 tissues; generates faceted and non-faceted
#          scatter plots with optional gene labels and R²
# Input:   meltTemperature.xlsx,
#          combined_plot_data_mouseTMT_<tissue>.csv
# Output:  mouseTMT_<tissue>_melTemp_*.pdf,
#          mouseTMT_meltTemp_insolubility_All.csv
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# ============================================================

library(tidyverse)
library(ggrepel)
library(readxl)
library(ggpmisc)

# NOTE: Set your working directory to the folder containing input files
# setwd("path/to/your/data")

# ---- 1. Load Data ----

melt_data <- read_excel("meltTemperature.xlsx")
write.csv(melt_data, "meltTemperature.csv", row.names = FALSE)
melt_data <- read.csv("meltTemperature.csv")

# Load per-tissue combined plot data (from 01_per-tissue-analysis scripts)
combined_plot_data_TA  <- read.csv("combined_plot_data_mouseTMT_TA.csv")
combined_plot_data_He  <- read.csv("combined_plot_data_mouseTMT_He.csv")
combined_plot_data_CB  <- read.csv("combined_plot_data_mouseTMT_CB.csv")
combined_plot_data_HC  <- read.csv("combined_plot_data_mouseTMT_HC.csv")
combined_plot_data_NC  <- read.csv("combined_plot_data_mouseTMT_NC.csv")
combined_plot_data_So  <- read.csv("combined_plot_data_mouseTMT_So.csv")
combined_plot_data_Li  <- read.csv("combined_plot_data_mouseTMT_Li.csv")
combined_plot_data_WAT <- read.csv("combined_plot_data_mouseTMT_WAT.csv")

# ---- 2. Join Melting Temperature Data to Each Tissue ----

# meltTemperature.xlsx uses column "gene_name"; rename for join
join_meltTemp <- function(plot_data) {
  plot_data %>%
    left_join(melt_data %>% rename(Gene_Name = gene_name), by = "Gene_Name")
}

combined_mT_TA  <- join_meltTemp(combined_plot_data_TA)
combined_mT_He  <- join_meltTemp(combined_plot_data_He)
combined_mT_CB  <- join_meltTemp(combined_plot_data_CB)
combined_mT_HC  <- join_meltTemp(combined_plot_data_HC)
combined_mT_NC  <- join_meltTemp(combined_plot_data_NC)
combined_mT_So  <- join_meltTemp(combined_plot_data_So)
combined_mT_Li  <- join_meltTemp(combined_plot_data_Li)
combined_mT_WAT <- join_meltTemp(combined_plot_data_WAT)

# ---- 3. Tissue-Specific Label Thresholds ----
# Label genes with melting points outside normal range (tissue-specific)

combined_mT_TA <- combined_mT_TA %>%
  filter(!is.na(Gene_Name)) %>%
  mutate(Label_mT = case_when(
    meltPoint > 59 | meltPoint < 42 ~ Gene_Name, TRUE ~ NA_character_))

combined_mT_He <- combined_mT_He %>%
  filter(!is.na(Gene_Name)) %>%
  mutate(Label_mT = case_when(
    meltPoint > 57 | meltPoint < 45 ~ Gene_Name, TRUE ~ NA_character_))

combined_mT_CB <- combined_mT_CB %>%
  filter(!is.na(Gene_Name)) %>%
  mutate(Label_mT = case_when(
    Comparison == "Old vs. Middle" & (meltPoint > 58 | meltPoint < 45) ~ Gene_Name,
    Comparison != "Old vs. Middle" & (meltPoint > 60 | meltPoint < 40) ~ Gene_Name,
    TRUE ~ NA_character_))

combined_mT_HC <- combined_mT_HC %>%
  filter(!is.na(Gene_Name)) %>%
  mutate(Label_mT = case_when(
    Comparison == "Old vs. Young" & (meltPoint > 60 | meltPoint < 40) ~ Gene_Name,
    Comparison != "Old vs. Young" & (meltPoint > 60 | meltPoint < 43) ~ Gene_Name,
    TRUE ~ NA_character_))

combined_mT_NC <- combined_mT_NC %>%
  filter(!is.na(Gene_Name)) %>%
  mutate(Label_mT = case_when(
    Comparison == "Old vs. Young" & (meltPoint > 58 | meltPoint < 42) ~ Gene_Name,
    Comparison != "Old vs. Young" & (meltPoint > 57 | meltPoint < 42) ~ Gene_Name,
    TRUE ~ NA_character_))

combined_mT_Li <- combined_mT_Li %>%
  filter(!is.na(Gene_Name)) %>%
  mutate(Label_mT = case_when(
    Comparison == "Old vs. Middle" & (meltPoint > 60 | meltPoint < 45) ~ Gene_Name,
    Comparison != "Old vs. Middle" & (meltPoint > 60 | meltPoint < 42) ~ Gene_Name,
    TRUE ~ NA_character_))

combined_mT_So <- combined_mT_So %>%
  filter(!is.na(Gene_Name)) %>%
  mutate(Label_mT = case_when(
    meltPoint > 60 | meltPoint < 42 ~ Gene_Name, TRUE ~ NA_character_))

combined_mT_WAT <- combined_mT_WAT %>%
  filter(!is.na(Gene_Name)) %>%
  mutate(Label_mT = case_when(
    Comparison == "Old vs. Middle" & (meltPoint > 60 | meltPoint < 45) ~ Gene_Name,
    Comparison != "Old vs. Middle" & (meltPoint > 60 | meltPoint < 43) ~ Gene_Name,
    TRUE ~ NA_character_))

# ---- 4. Plot Function ----

# Shared theme and scale elements
mT_theme <- theme_minimal() +
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

# Generate three plot variants per tissue:
#   _facet_label: faceted by comparison, with gene labels
#   _facet_nolabel: faceted, no labels
#   _facet_nolabel_R: faceted, no labels, with R² annotation
#   _nolabel: combined (no facet), no labels
#   _nolabel_R: combined, no labels, with R²
plot_meltTemp_tissue <- function(data, tissue_name, file_prefix) {
  labs_base <- labs(
    title   = paste0("MouseTMT Protein Melting Temperature vs. Insolubility Fold Change in ", tissue_name),
    x       = "Melting Temperature",
    y       = "Log2 fold change in Insoluble/Soluble",
    caption = "Point size based on statistical significance of Age difference in insolubility"
  )
  
  # Faceted with labels
  p_label <- ggplot(data, aes(x = meltPoint, y = Ratio_FC, size = Size,
                              label = Label_mT, color = Comparison)) +
    geom_point(alpha = 0.4) +
    geom_text_repel(size = 3, segment.color = "grey50", max.overlaps = Inf,
                    box.padding = 0.5, point.padding = 0.2, min.segment.length = 0.1) +
    size_scale + comparison_colors +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    labs_base + mT_theme +
    facet_wrap(~Comparison) +
    guides(color = guide_legend(override.aes = list(label = "")))
  
  print(p_label)
  ggsave(paste0(file_prefix, "_facet_label.pdf"), width = 7200, height = 6400, units = "px", dpi = 600)
  
  # Faceted, no labels
  p_nolabel_f <- ggplot(data, aes(x = meltPoint, y = Ratio_FC, size = Size, color = Comparison)) +
    geom_point(alpha = 0.4) +
    size_scale + comparison_colors +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    labs_base + mT_theme +
    facet_wrap(~Comparison)
  
  print(p_nolabel_f)
  ggsave(paste0(file_prefix, "_facet_nolabel.pdf"), width = 7200, height = 6400, units = "px", dpi = 600)
  
  # Faceted, no labels + R²
  p_nolabel_f_R <- p_nolabel_f +
    stat_poly_eq(aes(label = paste(after_stat(rr.label))), rr.digits = 3,
                 label.x = 1, label.y = 1.4, parse = TRUE, color = "black", show.legend = FALSE)
  
  print(p_nolabel_f_R)
  ggsave(paste0(file_prefix, "_facet_nolabel_R.pdf"), width = 7200, height = 6400, units = "px", dpi = 600)
  
  # No facet, no labels
  p_nolabel <- ggplot(data, aes(x = meltPoint, y = Ratio_FC, size = Size, color = Comparison)) +
    geom_point(alpha = 0.4) +
    size_scale + comparison_colors +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    labs_base + mT_theme
  
  print(p_nolabel)
  ggsave(paste0(file_prefix, "_nolabel.pdf"), width = 7200, height = 6400, units = "px", dpi = 600)
  
  # No facet, no labels + R²
  p_nolabel_R <- p_nolabel +
    stat_poly_eq(aes(label = paste(after_stat(rr.label))), rr.digits = 3,
                 label.x = 1, label.y = 1.4, parse = TRUE, color = "black", show.legend = FALSE)
  
  print(p_nolabel_R)
  ggsave(paste0(file_prefix, "_nolabel_R.pdf"), width = 7200, height = 6400, units = "px", dpi = 600)
}

# ---- 5. Generate Per-Tissue Plots ----

plot_meltTemp_tissue(combined_mT_TA,  "Tibialis Anterior",    "mouseTMT_TA_melTemp")
plot_meltTemp_tissue(combined_mT_He,  "Heart",                "mouseTMT_He_melTemp")
plot_meltTemp_tissue(combined_mT_CB,  "Cerebellum",           "mouseTMT_CB_melTemp")
plot_meltTemp_tissue(combined_mT_HC,  "Hippocampus",          "mouseTMT_HC_melTemp")
plot_meltTemp_tissue(combined_mT_NC,  "Neocortex",            "mouseTMT_NC_melTemp")
plot_meltTemp_tissue(combined_mT_Li,  "Liver",                "mouseTMT_Li_melTemp")
plot_meltTemp_tissue(combined_mT_So,  "Soleus",               "mouseTMT_So_melTemp")
plot_meltTemp_tissue(combined_mT_WAT, "White Adipose Tissue", "mouseTMT_WAT_melTemp")

# ---- 6. Whole-Proteome Combined Plot ----

# Shared organ color palette
organ_colors <- c(
  "Cerebellum"           = "#c7a6db",
  "Heart"                = "#c93b92",
  "Hippocampus"          = "#55397f",
  "Liver"                = "#4d9fdb",
  "Neocortex"            = "#8a68ca",
  "Soleus"               = "#a8192a",
  "Tibialis Anterior"    = "#e991b7",
  "White Adipose Tissue" = "#dfc41b"
)

# Add organ labels and combine
combined_mT_TA$Organ  <- "Tibialis Anterior"
combined_mT_He$Organ  <- "Heart"
combined_mT_CB$Organ  <- "Cerebellum"
combined_mT_HC$Organ  <- "Hippocampus"
combined_mT_NC$Organ  <- "Neocortex"
combined_mT_So$Organ  <- "Soleus"
combined_mT_Li$Organ  <- "Liver"
combined_mT_WAT$Organ <- "White Adipose Tissue"

combined_mT_All <- bind_rows(
  combined_mT_TA, combined_mT_He, combined_mT_CB, combined_mT_HC,
  combined_mT_NC, combined_mT_So, combined_mT_Li, combined_mT_WAT
)

write.csv(combined_mT_All, "mouseTMT_meltTemp_insolubility_All.csv", row.names = FALSE)

# Whole-proteome scatter: meltPoint vs Ratio_FC, colored by organ, faceted by comparison
mouseTMT_meltTemp_insolubility_All <- ggplot(combined_mT_All,
                                             aes(x = meltPoint, y = Ratio_FC, size = Size, color = Organ)) +
  geom_point(alpha = 0.3) +
  size_scale +
  scale_color_manual(values = organ_colors) +
  labs(title   = "MouseTMT Whole Proteome: Melting Temperature vs. Ratio Fold Change",
       x       = "Melting Temperature",
       y       = "Log2 fold change in Insoluble/Soluble",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  facet_wrap(~Comparison) +
  xlim(35, 70) +
  mT_theme

# With R² annotation
mouseTMT_meltTemp_insolubility_All_R <- mouseTMT_meltTemp_insolubility_All +
  stat_poly_eq(aes(label = paste(after_stat(rr.label))), rr.digits = 3,
               label.x = 1, label.y = 1.4, parse = TRUE, color = "black", show.legend = FALSE)

print(mouseTMT_meltTemp_insolubility_All_R)
ggsave("mouseTMT_meltTemp_insolubility_All_R.pdf", width = 7200, height = 6400, units = "px", dpi = 600)