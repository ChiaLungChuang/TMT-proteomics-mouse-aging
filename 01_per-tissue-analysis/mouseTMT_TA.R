# ============================================================
# Project: TMT Proteomics — Multi-organ Mouse Aging
# Script:  mouseTMT_TA.R
# Purpose: Per-tissue analysis for tibialis anterior (ta) — data reshaping,
#          ANOVA + Tukey HSD across age groups, Log2 fold
#          change calculation, and scatter plot visualization
#          (Abundance FC vs. Insolubility Ratio FC)
# Input:   FD_TA.xlsx (TMT proteomics raw data)
# Output:  mouseTMT_TA.csv, scatter plots as PDF
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# ============================================================

library(tidyverse)
library(ggrepel)
library(readxl)
library(broom)
library(extrafont)
library(ggpmisc)

# NOTE: Set your working directory to the folder containing FD_Liver.xlsx
# setwd("path/to/your/data")

# ---- 1. Data Import and Preparation ----

data <- read_excel("FD_TA.xlsx")
write.csv(data, "FD_TibialisAnterior(TA).csv", row.names = FALSE)
data_TA <- read.csv("FD_TibialisAnterior(TA).csv")

# Clean protein accession by removing "_DROME" suffix
data_TA <- data_TA %>%
  mutate(Protein.ID_cleaned = sub("_DROME$", "", Protein.Accession..))

# Create a unique GN_ID combining gene name and protein ID to distinguish isoforms
data_TA <- data_TA %>%
  mutate(GN_ID = paste(GN, Protein.ID_cleaned, sep = "_"))

write.csv(data_TA, "data_mouse_TibialisAnterior(TA).csv", row.names = FALSE)

# ---- 2. Data Reshaping ----

# Helper function to normalize replicate numbering within each Age/Solubility group
normalize_replicates <- function(df) {
  df %>%
    group_by(Age, Solubility) %>%
    mutate(
      Replicate = as.integer((dense_rank(Replicate) - 1) %% 3 + 1)
    ) %>%
    ungroup()
}

# Reshape from wide to long format; parse Age, Replicate, Solubility from column names
data_TA_long <- data_TA %>%
  pivot_longer(
    cols = starts_with("sig"),
    names_to = c("Age", "Replicate", "Solubility"),
    names_pattern = "sig\\d+[NC]?..TA_(\\d+M)(\\d+)_([SI])",
    values_to = "Measurement"
  ) %>%
  mutate(
    Age = case_when(
      Age == "6M"  ~ "Young",
      Age == "24M" ~ "Middle",
      Age == "30M" ~ "Old"
    ),
    Solubility = case_when(
      Solubility == "S" ~ "Soluble",
      Solubility == "I" ~ "Insoluble"
    ),
    Replicate   = as.integer(Replicate),
    Measurement = as.numeric(Measurement)
  ) %>%
  normalize_replicates

# Pivot to wide format to get Soluble and Insoluble as separate columns,
# then compute Ratio (Insoluble/Soluble) and total Abundance
data_TA_wider <- data_TA_long %>%
  pivot_wider(
    names_from  = Solubility,
    values_from = Measurement,
    id_cols     = c(Age, GN_ID, Replicate)
  ) %>%
  mutate(
    Ratio     = Insoluble / Soluble,
    Abundance = Soluble + Insoluble
  )

# ---- 3. Statistical Analysis: One-Way ANOVA + Tukey HSD ----

# Fit ANOVA models per protein for Abundance and Ratio
anova_results_abundance <- data_TA_wider %>%
  group_by(GN_ID) %>%
  do(aov_result = aov(Abundance ~ Age, data = .)) %>%
  ungroup()

anova_results_ratio <- data_TA_wider %>%
  group_by(GN_ID) %>%
  do(aov_result = aov(Ratio ~ Age, data = .)) %>%
  ungroup()

# Extract overall ANOVA p-values
anova_p_values_abundance <- anova_results_abundance %>%
  mutate(p.value = map_dbl(aov_result, ~ summary(.)[[1]][["Pr(>F)"]][1])) %>%
  select(GN_ID, p.value)

anova_p_values_ratio <- anova_results_ratio %>%
  mutate(p.value = map_dbl(aov_result, ~ summary(.)[[1]][["Pr(>F)"]][1])) %>%
  select(GN_ID, p.value)

# Tukey HSD post-hoc test for pairwise age comparisons
tukey_results_abundance <- anova_results_abundance %>%
  mutate(tukey_result = map(aov_result, ~ tidy(TukeyHSD(.))))

tukey_results_ratio <- anova_results_ratio %>%
  mutate(tukey_result = map(aov_result, ~ tidy(TukeyHSD(.))))

# Tidy Tukey results into wide format with named p-value columns
tidy_tukey_results_abundance <- tukey_results_abundance %>%
  unnest(tukey_result) %>%
  select(GN_ID, comparison = contrast, p.value = adj.p.value) %>%
  pivot_wider(names_from = comparison, values_from = p.value, names_prefix = "p.value_") %>%
  rename(
    ab.p.value_middle_young = `p.value_Young-Middle`,
    ab.p.value_old_middle   = `p.value_Old-Middle`,
    ab.p.value_old_young    = `p.value_Young-Old`
  )

tidy_tukey_results_ratio <- tukey_results_ratio %>%
  unnest(tukey_result) %>%
  select(GN_ID, comparison = contrast, p.value = adj.p.value) %>%
  pivot_wider(names_from = comparison, values_from = p.value, names_prefix = "p.value_") %>%
  rename(
    ra.p.value_middle_young = `p.value_Young-Middle`,
    ra.p.value_old_middle   = `p.value_Old-Middle`,
    ra.p.value_old_young    = `p.value_Young-Old`
  )

# ---- 4. Fold Change Calculation ----

# Log2 fold change for each pairwise age comparison (Abundance and Ratio)
ratio_data_TA_OY <- data_TA_wider %>%
  group_by(GN_ID) %>%
  summarize(
    Abundance_mean_old   = mean(Abundance[Age == "Old"],   na.rm = TRUE),
    Abundance_mean_young = mean(Abundance[Age == "Young"], na.rm = TRUE),
    Ratio_mean_old       = mean(Ratio[Age == "Old"],       na.rm = TRUE),
    Ratio_mean_young     = mean(Ratio[Age == "Young"],     na.rm = TRUE)
  ) %>%
  mutate(
    Abundance_FC_old_young = log2(Abundance_mean_old / Abundance_mean_young),
    Ratio_FC_old_young     = log2(Ratio_mean_old / Ratio_mean_young)
  )

ratio_data_TA_MY <- data_TA_wider %>%
  group_by(GN_ID) %>%
  summarize(
    Abundance_mean_middle = mean(Abundance[Age == "Middle"], na.rm = TRUE),
    Abundance_mean_young  = mean(Abundance[Age == "Young"],  na.rm = TRUE),
    Ratio_mean_middle     = mean(Ratio[Age == "Middle"],     na.rm = TRUE),
    Ratio_mean_young      = mean(Ratio[Age == "Young"],      na.rm = TRUE)
  ) %>%
  mutate(
    Abundance_FC_middle_young = log2(Abundance_mean_middle / Abundance_mean_young),
    Ratio_FC_middle_young     = log2(Ratio_mean_middle / Ratio_mean_young)
  )

ratio_data_TA_OM <- data_TA_wider %>%
  group_by(GN_ID) %>%
  summarize(
    Abundance_mean_middle = mean(Abundance[Age == "Middle"], na.rm = TRUE),
    Abundance_mean_old    = mean(Abundance[Age == "Old"],    na.rm = TRUE),
    Ratio_mean_middle     = mean(Ratio[Age == "Middle"],     na.rm = TRUE),
    Ratio_mean_old        = mean(Ratio[Age == "Old"],        na.rm = TRUE)
  ) %>%
  mutate(
    Abundance_FC_old_middle = log2(Abundance_mean_old / Abundance_mean_middle),
    Ratio_FC_old_middle     = log2(Ratio_mean_old / Ratio_mean_middle)
  )

# Merge all fold change and p-value data into one plotting dataframe
ratio_data_TA <- ratio_data_TA_OY %>%
  left_join(ratio_data_TA_OM, by = "GN_ID") %>%
  left_join(ratio_data_TA_MY, by = "GN_ID")

plot_data_TA <- ratio_data_TA %>%
  left_join(anova_p_values_ratio,       by = "GN_ID") %>%
  left_join(anova_p_values_abundance,   by = "GN_ID") %>%
  left_join(tidy_tukey_results_ratio,   by = "GN_ID") %>%
  left_join(tidy_tukey_results_abundance, by = "GN_ID") %>%
  mutate(
    # Point size scaled to -log10(Tukey p-value) for each comparison
    Size_ra_old_young    = ifelse(is.na(ra.p.value_old_young),    1, -log10(ra.p.value_old_young)),
    Size_ra_middle_young = ifelse(is.na(ra.p.value_middle_young), 1, -log10(ra.p.value_middle_young)),
    Size_ra_old_middle   = ifelse(is.na(ra.p.value_old_middle),   1, -log10(ra.p.value_old_middle)),
    Size_ab_old_young    = ifelse(is.na(ab.p.value_old_young),    1, -log10(ab.p.value_old_young)),
    Size_ab_middle_young = ifelse(is.na(ab.p.value_middle_young), 1, -log10(ab.p.value_middle_young)),
    Size_ab_old_middle   = ifelse(is.na(ab.p.value_old_middle),   1, -log10(ab.p.value_old_middle))
  ) %>%
  mutate(Gene_Name = sub("_.*", "", GN_ID)) %>%          # Extract gene name from GN_ID
  filter(!grepl("co|CON_", GN_ID))                        # Remove contaminant proteins

write.csv(plot_data_TA, "mouseTMT_TA.csv", row.names = FALSE)
plot_data_TA <- read.csv("mouseTMT_TA.csv")

# Load fonts for PDF export
extrafont::font_import(prompt = FALSE)
loadfonts(device = "pdf")

# ---- 5. Visualization ----
# Each comparison (OY, MY, OM) generates three plot variants:
#   _top10  : labels top 10 most significant proteins
#   _label  : labels proteins exceeding FC thresholds
#   _nolabel: clean plot without labels
# A combined facet plot and overlapping plot are also generated.

##### Old vs. Young #####

plot_data_TA_filtered_OY <- plot_data_TA %>%
  filter(ra.p.value_old_young <= 0.05)

top_genes_OY <- plot_data_TA_filtered_OY %>%
  arrange(ra.p.value_old_young) %>%
  slice_head(n = 10)

plot_data_TA_filtered_OY <- plot_data_TA_filtered_OY %>%
  mutate(Label_top = case_when(
    Gene_Name == "NA"              ~ NA_character_,
    GN_ID %in% top_genes_OY$GN_ID ~ Gene_Name,
    TRUE                           ~ NA_character_
  ))

plot_data_TA_com_filtered_OY <- plot_data_TA_filtered_OY %>%
  mutate(Comparison = "Old vs. Young") %>%
  select(GN_ID, Gene_Name, Abundance_FC_old_young, Ratio_FC_old_young, Size_ra_old_young, Comparison, Label_top) %>%
  rename(Abundance_FC = Abundance_FC_old_young, Ratio_FC = Ratio_FC_old_young, Size = Size_ra_old_young) %>%
  mutate(Label = case_when(
    Gene_Name == "NA" ~ NA_character_,
    Comparison == "Old vs. Young" & (Ratio_FC > 0.9 | Abundance_FC > 1.5 | Abundance_FC < -1.3 | Ratio_FC < -1.2) ~ Gene_Name,
    TRUE ~ NA_character_
  ))

# Shared theme for all plots
tmT_theme <- theme_minimal(base_family = "Arial") +
  theme(
    axis.text.x    = element_text(size = 14, face = "bold"),
    axis.text.y    = element_text(size = 14, face = "bold"),
    legend.text    = element_text(size = 12, face = "italic"),
    legend.title   = element_text(size = 14, face = "bold"),
    plot.title     = element_text(size = 16, face = "bold"),
    axis.title.x   = element_text(size = 14, face = "bold"),
    axis.title.y   = element_text(size = 14, face = "bold"),
    strip.text     = element_text(size = 14, face = "bold", family = "Arial")
  )

# Shared size scale
size_scale <- scale_size_continuous(
  range  = c(1, 9),
  name   = "P-Value Size (-log10)",
  breaks = c(1.3, 2, 3, 4, 5),
  limit  = c(1.3, 8),
  labels = c("1.30(p<0.05)", "2", "3", "4", "5")
)

# Top 10 labels — Old vs. Young
mouseTMT_TA_OY_top10 <- ggplot(plot_data_TA_filtered_OY,
                               aes(x = Abundance_FC_old_young, y = Ratio_FC_old_young,
                                   size = Size_ra_old_young, label = Label_top)) +
  geom_point(alpha = 0.4, color = "purple4") +
  geom_text_repel(vjust = 0.5, hjust = 0.5, size = 5, fontface = "bold",
                  family = "Arial", segment.color = "grey50", max.overlaps = Inf,
                  nudge_x = 0.2, nudge_y = 0.2,
                  box.padding = unit(0.1, "lines"), point.padding = unit(0.1, "lines"),
                  min.segment.length = 0.1, na.rm = TRUE) +
  size_scale +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(title   = "MouseTMT Tibialis Anterior (TA): Abundance FC vs. Ratio FC (Old/Young)",
       x       = "Log2 Abundance fold change (Old/Young)",
       y       = "Log2 fold change in Insoluble/Soluble (Old/Young)",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  tmT_theme

print(mouseTMT_TA_OY_top10)
ggsave("mouseTMT_TA_OY_top10.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

# FC-threshold labels — Old vs. Young
mouseTMT_TA_OY_label <- ggplot(plot_data_TA_com_filtered_OY,
                               aes(x = Abundance_FC, y = Ratio_FC, size = Size, label = Label)) +
  geom_point(alpha = 0.4, color = "purple4") +
  geom_text_repel(vjust = 0.5, hjust = 0.5, size = 5, fontface = "bold",
                  family = "Arial", segment.color = "grey50", max.overlaps = Inf,
                  nudge_x = 0.15, nudge_y = 0.15, box.padding = 0.5, point.padding = 0.5,
                  min.segment.length = 0.1, na.rm = TRUE) +
  size_scale +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(title   = "MouseTMT Tibialis Anterior (TA): Abundance FC vs. Ratio FC (Old/Young)",
       x       = "Log2 Abundance fold change (Old/Young)",
       y       = "Log2 fold change in Insoluble/Soluble (Old/Young)",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  tmT_theme

print(mouseTMT_TA_OY_label)
ggsave("mouseTMT_TA_OY_label.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

# Add R-squared annotation
mouseTMT_TA_OY_label_R <- mouseTMT_TA_OY_label +
  stat_poly_eq(aes(label = paste(after_stat(rr.label))),
               rr.digits = 3, parse = TRUE, color = "black", show.legend = FALSE)

print(mouseTMT_TA_OY_label_R)
ggsave("mouseTMT_TA_OY_label_R.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

# No labels — Old vs. Young
mouseTMT_TA_OY_nolabel <- ggplot(plot_data_TA_filtered_OY,
                                 aes(x = Abundance_FC_old_young, y = Ratio_FC_old_young, size = Size_ra_old_young)) +
  geom_point(alpha = 0.4, color = "purple4") +
  size_scale +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(title   = "MouseTMT Tibialis Anterior (TA): Abundance FC vs. Ratio FC (Old/Young)",
       x       = "Log2 Abundance fold change (Old/Young)",
       y       = "Log2 fold change in Insoluble/Soluble (Old/Young)",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  tmT_theme

print(mouseTMT_TA_OY_nolabel)
ggsave("mouseTMT_TA_OY_nolabel.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

##### Middle vs. Young #####

plot_data_TA_filtered_MY <- plot_data_TA %>%
  filter(ra.p.value_middle_young <= 0.05)

top_genes_MY <- plot_data_TA_filtered_MY %>%
  arrange(ra.p.value_middle_young) %>%
  slice_head(n = 10)

plot_data_TA_filtered_MY <- plot_data_TA_filtered_MY %>%
  mutate(Label_top = case_when(
    Gene_Name == "NA"              ~ NA_character_,
    GN_ID %in% top_genes_MY$GN_ID ~ Gene_Name,
    TRUE                           ~ NA_character_
  ))

plot_data_TA_com_filtered_MY <- plot_data_TA_filtered_MY %>%
  mutate(Comparison = "Middle vs. Young") %>%
  select(GN_ID, Gene_Name, Abundance_FC_middle_young, Ratio_FC_middle_young, Size_ra_middle_young, Comparison, Label_top) %>%
  rename(Abundance_FC = Abundance_FC_middle_young, Ratio_FC = Ratio_FC_middle_young, Size = Size_ra_middle_young) %>%
  mutate(Label = case_when(
    Gene_Name == "NA"   ~ NA_character_,
    Gene_Name == "Hrnr" ~ "Hrnr",
    Comparison == "Middle vs. Young" & (Ratio_FC > 0.8 | Abundance_FC > 0.8 | Abundance_FC < -1 | Ratio_FC < -1) ~ Gene_Name,
    TRUE ~ NA_character_
  ))

mouseTMT_TA_MY_top10 <- ggplot(plot_data_TA_filtered_MY,
                               aes(x = Abundance_FC_middle_young, y = Ratio_FC_middle_young,
                                   size = Size_ra_middle_young, label = Label_top)) +
  geom_point(alpha = 0.4, color = "maroon2") +
  geom_text_repel(vjust = 0.5, hjust = 0.5, size = 5, fontface = "bold",
                  family = "Arial", segment.color = "grey50", max.overlaps = Inf,
                  nudge_x = 0.15, nudge_y = 0.15, box.padding = 0.5, point.padding = 0.5,
                  min.segment.length = 0.1) +
  size_scale +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(title   = "MouseTMT Tibialis Anterior (TA): Abundance FC vs. Ratio FC (Middle/Young)",
       x       = "Log2 Abundance fold change (Middle/Young)",
       y       = "Log2 fold change in Insoluble/Soluble (Middle/Young)",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  tmT_theme

print(mouseTMT_TA_MY_top10)
ggsave("mouseTMT_TA_MY_top10.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

mouseTMT_TA_MY_label <- ggplot(plot_data_TA_com_filtered_MY,
                               aes(x = Abundance_FC, y = Ratio_FC, size = Size, label = Label)) +
  geom_point(alpha = 0.4, color = "maroon2") +
  geom_text_repel(vjust = 0.5, hjust = 0.5, size = 5, fontface = "bold",
                  family = "Arial", segment.color = "grey50", max.overlaps = Inf,
                  nudge_x = 0.15, nudge_y = 0.15, box.padding = 0.5, point.padding = 0.5,
                  min.segment.length = 0.1) +
  scale_size_continuous(range = c(1,9), name = "P-Value Size (-log10)",
                        breaks = c(1.3,2,3,4,5), limit = c(1.3,6),
                        labels = c("1.30(p<0.05)","2","3","4","5")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(title   = "MouseTMT Tibialis Anterior (TA): Abundance FC vs. Ratio FC (Middle/Young)",
       x       = "Log2 Abundance fold change (Middle/Young)",
       y       = "Log2 fold change in Insoluble/Soluble (Middle/Young)",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  tmT_theme

print(mouseTMT_TA_MY_label)
ggsave("mouseTMT_TA_MY_label.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

mouseTMT_TA_MY_label_R <- mouseTMT_TA_MY_label +
  stat_poly_eq(aes(label = paste(after_stat(rr.label))),
               rr.digits = 3, parse = TRUE, label.x = 1, label.y = 1,
               color = "black", show.legend = FALSE)

print(mouseTMT_TA_MY_label_R)
ggsave("mouseTMT_TA_MY_label_R.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

mouseTMT_TA_MY_nolabel <- ggplot(plot_data_TA_filtered_MY,
                                 aes(x = Abundance_FC_middle_young, y = Ratio_FC_middle_young, size = Size_ra_middle_young)) +
  geom_point(alpha = 0.4, color = "maroon2") +
  scale_size_continuous(range = c(1,9), name = "P-Value Size (-log10)",
                        breaks = c(1.3,2,3,4,5), limit = c(1.3,8),
                        labels = c("1.30(p<0.05)","2","3","4","5")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(title   = "MouseTMT Tibialis Anterior (TA): Abundance FC vs. Ratio FC (Middle/Young)",
       x       = "Log2 Abundance fold change (Middle/Young)",
       y       = "Log2 fold change in Insoluble/Soluble (Middle/Young)",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  tmT_theme

print(mouseTMT_TA_MY_nolabel)
ggsave("mouseTMT_TA_MY_nolabel.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

##### Old vs. Middle #####

plot_data_TA_filtered_OM <- plot_data_TA %>%
  filter(ra.p.value_old_middle <= 0.05)

top_genes_OM <- plot_data_TA_filtered_OM %>%
  arrange(ra.p.value_old_middle) %>%
  slice_head(n = 10)

plot_data_TA_filtered_OM <- plot_data_TA_filtered_OM %>%
  mutate(Label_top = case_when(
    Gene_Name == "NA"              ~ NA_character_,
    GN_ID %in% top_genes_OM$GN_ID ~ Gene_Name,
    TRUE                           ~ NA_character_
  ))

plot_data_TA_com_filtered_OM <- plot_data_TA_filtered_OM %>%
  mutate(Comparison = "Old vs. Middle") %>%
  select(GN_ID, Gene_Name, Abundance_FC_old_middle, Ratio_FC_old_middle, Size_ra_old_middle, Comparison, Label_top) %>%
  rename(Abundance_FC = Abundance_FC_old_middle, Ratio_FC = Ratio_FC_old_middle, Size = Size_ra_old_middle) %>%
  mutate(Label = case_when(
    Gene_Name == "NA" ~ NA_character_,
    Comparison == "Old vs. Middle" & (Ratio_FC > 0.7 | Abundance_FC > 0.9 | Abundance_FC < -0.5 | Ratio_FC < -0.5) ~ Gene_Name,
    TRUE ~ NA_character_
  ))

mouseTMT_TA_OM_top10 <- ggplot(plot_data_TA_filtered_OM,
                               aes(x = Abundance_FC_old_middle, y = Ratio_FC_old_middle,
                                   size = Size_ra_old_middle, label = Label_top)) +
  geom_point(alpha = 0.4, color = "darkgreen") +
  geom_text_repel(vjust = 0.5, hjust = 0.5, size = 5, fontface = "bold",
                  family = "Arial", segment.color = "grey50", max.overlaps = Inf,
                  nudge_x = 0.15, nudge_y = 0.15, box.padding = 0.5, point.padding = 0.5,
                  min.segment.length = 0.1) +
  scale_size_continuous(range = c(1,9), name = "P-Value Size (-log10)",
                        breaks = c(1.3,2,3,4,5), limit = c(1.3,6),
                        labels = c("1.30(p<0.05)","2","3","4","5")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(title   = "MouseTMT Tibialis Anterior (TA): Abundance FC vs. Ratio FC (Old/Middle)",
       x       = "Log2 Abundance fold change (Old/Middle)",
       y       = "Log2 fold change in Insoluble/Soluble (Old/Middle)",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  tmT_theme

print(mouseTMT_TA_OM_top10)
ggsave("mouseTMT_TA_OM_top10.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

mouseTMT_TA_OM_label <- ggplot(plot_data_TA_com_filtered_OM,
                               aes(x = Abundance_FC, y = Ratio_FC, size = Size, label = Label)) +
  geom_point(alpha = 0.4, color = "darkgreen") +
  geom_text_repel(vjust = 0.5, hjust = 0.5, size = 5, fontface = "bold",
                  family = "Arial", segment.color = "grey50", max.overlaps = Inf,
                  nudge_x = 0.15, nudge_y = 0.15, box.padding = 0.5, point.padding = 0.5,
                  min.segment.length = 0.1) +
  scale_size_continuous(range = c(1,9), name = "P-Value Size (-log10)",
                        breaks = c(1.3,2,3,4,5), limit = c(1.3,6),
                        labels = c("1.30(p<0.05)","2","3","4","5")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(title   = "MouseTMT Tibialis Anterior (TA): Abundance FC vs. Ratio FC (Old/Middle)",
       x       = "Log2 Abundance fold change (Old/Middle)",
       y       = "Log2 fold change in Insoluble/Soluble (Old/Middle)",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  tmT_theme

print(mouseTMT_TA_OM_label)
ggsave("mouseTMT_TA_OM_label.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

mouseTMT_TA_OM_label_R <- mouseTMT_TA_OM_label +
  stat_poly_eq(aes(label = paste(after_stat(rr.label))),
               rr.digits = 3, parse = TRUE, color = "black", show.legend = FALSE)

print(mouseTMT_TA_OM_label_R)
ggsave("mouseTMT_TA_OM_label_R.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

mouseTMT_TA_OM_nolabel <- ggplot(plot_data_TA_filtered_OM,
                                 aes(x = Abundance_FC_old_middle, y = Ratio_FC_old_middle, size = Size_ra_old_middle)) +
  geom_point(alpha = 0.4, color = "darkgreen") +
  scale_size_continuous(range = c(1,9), name = "P-Value Size (-log10)",
                        breaks = c(1.3,2,3,4,5), limit = c(1.3,6),
                        labels = c("1.30(p<0.05)","2","3","4","5")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(title   = "MouseTMT Tibialis Anterior (TA): Abundance FC vs. Ratio FC (Old/Middle)",
       x       = "Log2 Abundance fold change (Old/Middle)",
       y       = "Log2 fold change in Insoluble/Soluble (Old/Middle)",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  tmT_theme

print(mouseTMT_TA_OM_nolabel)
ggsave("mouseTMT_TA_OM_nolabel.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

# ---- 6. Combined Three-Comparison Plots ----

combined_plot_data_mouseTMT_TA <- bind_rows(
  plot_data_TA_com_filtered_OM,
  plot_data_TA_com_filtered_MY,
  plot_data_TA_com_filtered_OY
)

write.csv(combined_plot_data_mouseTMT_TA, "combined_plot_data_mouseTMT_TA.csv", row.names = FALSE)

comparison_colors <- scale_color_manual(
  values = c("Old vs. Young" = "purple4", "Middle vs. Young" = "maroon2", "Old vs. Middle" = "darkgreen")
)

# Faceted — with labels
mouseTMT_TA_combined_facet_label <- ggplot(combined_plot_data_mouseTMT_TA,
                                           aes(x = Abundance_FC, y = Ratio_FC, size = Size, label = Label, color = Comparison)) +
  geom_point(alpha = 0.4) +
  geom_text_repel(vjust = 0.5, hjust = 0.5, size = 5, fontface = "bold",
                  family = "Arial", segment.color = "grey50", max.overlaps = Inf,
                  nudge_x = 0.15, nudge_y = 0.15, box.padding = 0.5, point.padding = 0.5,
                  min.segment.length = 0.1) +
  scale_size_continuous(range = c(1,9), name = "P-Value Size (-log10)",
                        breaks = c(1.3,2,3,4,5), limit = c(1.3,6),
                        labels = c("1.30(p<0.05)","2","3","4","5")) +
  comparison_colors +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(title   = "MouseTMT Tibialis Anterior (TA): Abundance FC vs. Ratio FC (All Comparisons)",
       x       = "Log2 Abundance fold change",
       y       = "Log2 fold change in Insoluble/Soluble",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  tmT_theme +
  facet_wrap(~Comparison) +
  guides(color = guide_legend(override.aes = list(label = "")))

print(mouseTMT_TA_combined_facet_label)
ggsave("mouseTMT_TA_combined_facet_label.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

# Faceted — no labels
mouseTMT_TA_combined_facet_nolabel <- ggplot(combined_plot_data_mouseTMT_TA,
                                             aes(x = Abundance_FC, y = Ratio_FC, size = Size, color = Comparison)) +
  geom_point(alpha = 0.4) +
  scale_size_continuous(range = c(1,9), name = "P-Value Size (-log10)",
                        breaks = c(1.3,2,3,4,5), limit = c(1.3,6),
                        labels = c("1.30(p<0.05)","2","3","4","5")) +
  comparison_colors +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(title   = "MouseTMT Tibialis Anterior (TA): Abundance FC vs. Ratio FC (All Comparisons)",
       x       = "Log2 Abundance fold change",
       y       = "Log2 fold change in Insoluble/Soluble",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  tmT_theme +
  facet_wrap(~Comparison)

print(mouseTMT_TA_combined_facet_nolabel)
ggsave("mouseTMT_TA_combined_facet_nolabel.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

# Overlapping (no facet) — no labels
mouseTMT_TA_combined_nofacet_nolabel <- ggplot(combined_plot_data_mouseTMT_TA,
                                               aes(x = Abundance_FC, y = Ratio_FC, size = Size, color = Comparison)) +
  geom_point(alpha = 0.3) +
  scale_size_continuous(range = c(1,9), name = "P-Value Size (-log10)",
                        breaks = c(1.3,2,3,4,5), limit = c(1.3,6),
                        labels = c("1.30(p<0.05)","2","3","4","5")) +
  comparison_colors +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(title   = "MouseTMT Tibialis Anterior (TA): Abundance FC vs. Ratio FC (All Comparisons)",
       x       = "Log2 Abundance fold change",
       y       = "Log2 fold change in Insoluble/Soluble",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  tmT_theme

print(mouseTMT_TA_combined_nofacet_nolabel)
ggsave("mouseTMT_TA_combined_nofacet_nolabel.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

# Overlapping (no facet) — with labels
mouseTMT_TA_combined_nofacet_label <- ggplot(combined_plot_data_mouseTMT_TA,
                                             aes(x = Abundance_FC, y = Ratio_FC, size = Size, label = Label, color = Comparison)) +
  geom_point(alpha = 0.3) +
  geom_text_repel(vjust = 0.5, hjust = 0.5, size = 5, fontface = "bold",
                  family = "Arial", segment.color = "grey50", max.overlaps = Inf,
                  nudge_x = 0.15, nudge_y = 0.15, box.padding = 0.5, point.padding = 0.5,
                  min.segment.length = 0.1) +
  scale_size_continuous(range = c(1,9), name = "P-Value Size (-log10)",
                        breaks = c(1.3,2,3,4,5), limit = c(1.3,6),
                        labels = c("1.30(p<0.05)","2","3","4","5")) +
  comparison_colors +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(title   = "MouseTMT Tibialis Anterior (TA): Abundance FC vs. Ratio FC (All Comparisons)",
       x       = "Log2 Abundance fold change",
       y       = "Log2 fold change in Insoluble/Soluble",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  tmT_theme +
  guides(color = guide_legend(override.aes = list(label = "")))

print(mouseTMT_TA_combined_nofacet_label)
ggsave("mouseTMT_TA_combined_nofacet_label.pdf", width = 7200, height = 6400, units = "px", dpi = 600)