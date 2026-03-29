# ============================================================
# Project: TMT Proteomics — Multi-organ Mouse Aging
# Script:  Organs_mRNA.R
# Purpose: Cross-organ network graphs with nodes colored by
#          mRNA expression levels (average TPM per tissue),
#          for all three age comparisons. Proteins undetectable
#          at mRNA level are colored black.
# Input:   combined_plot_data_mouseTMT_All.csv,
#          combined_abundance_All.csv (from Organs_abundance.R),
#          mRNA levels.xlsx (multi-sheet: one sheet per organ)
# Output:  mRNA_colored_network_<comparison>_kk.pdf (labeled and no-label)
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# ============================================================

library(readxl)
library(ggnewscale)
library(dplyr)
library(tidyr)
library(ggraph)
library(igraph)
library(tidygraph)
library(scales)
library(ggplot2)
library(stringr)
library(ggrepel)

# NOTE: Set your working directory to the folder containing input CSVs
# and the mRNA levels.xlsx file
# setwd("path/to/your/data")

# ---- 1. Load Data ----

combined_plot_data_mouseTMT_All <- read.csv("combined_plot_data_mouseTMT_All.csv")
combined_abundance_All          <- read.csv("combined_abundance_All.csv")

# Split by age comparison
combined_All_old_vs_young         <- combined_plot_data_mouseTMT_All %>% filter(Comparison == "Old vs. Young")
combined_All_old_vs_middle        <- combined_plot_data_mouseTMT_All %>% filter(Comparison == "Old vs. Middle")
combined_All_middle_vs_young      <- combined_plot_data_mouseTMT_All %>% filter(Comparison == "Middle vs. Young")
combined_abundance_All_old_vs_young    <- combined_abundance_All %>% filter(Comparison == "Old vs. Young")
combined_abundance_All_old_vs_middle   <- combined_abundance_All %>% filter(Comparison == "Old vs. Middle")
combined_abundance_All_middle_vs_young <- combined_abundance_All %>% filter(Comparison == "Middle vs. Young")

# ---- 2. Load mRNA Data (multi-sheet Excel) ----

excel_file <- "mRNA levels.xlsx"

# Map sheet names (lowercase) to standardized organ names
organ_mapping <- c(
  "white adipose"  = "White Adipose Tissue",
  "heart"          = "Heart",
  "liver"          = "Liver",
  "neocortex"      = "Neocortex",
  "hippocampus"    = "Hippocampus",
  "cerebellum"     = "Cerebellum",
  "soleus"         = "Soleus",
  "tibialis anterior" = "Tibialis Anterior"
)

# Read and combine all sheets; add Organ column from sheet name
mRNA_All <- lapply(excel_sheets(excel_file), function(sheet) {
  df <- read_excel(excel_file, sheet = sheet)
  df %>% mutate(Organ = organ_mapping[tolower(sheet)])
}) %>%
  bind_rows() %>%
  rename(Gene_Name = Symbol, average_TPMs = `average TPMs`)

# ---- 3. Graph Builder Function ----

# Builds a tbl_graph for one comparison's mRNA-annotated network
build_mRNA_graph <- function(combined_abundance_comp, mRNA_data) {
  # Join mRNA levels to the abundance data
  combined_mRNA <- combined_abundance_comp %>%
    left_join(mRNA_data, by = c("Gene_Name", "Organ")) %>%
    select(GN_ID, Gene_Name, average_TPMs, Organ, Comparison) %>%
    filter(!is.na(Gene_Name))
  
  # Per-gene summary: take max TPM across organs; classify inter/intra-organ
  mRNA_summary <- combined_mRNA %>%
    mutate(average_TPMs = ifelse(is.na(average_TPMs), 0, average_TPMs)) %>%
    group_by(GN_ID) %>%
    summarise(
      average_TPMs = max(average_TPMs, na.rm = TRUE),
      Organ        = if (n_distinct(Organ) > 1) "Inter-organ" else first(Organ),
      gene_type    = if (n_distinct(Organ) > 1) "inter_organ"  else "intra_organ",
      .groups = "drop"
    )
  
  # Edges: gene → organ
  edges <- combined_mRNA %>%
    select(from = GN_ID, to = Organ) %>%
    distinct()
  
  # Vertices: genes with mRNA annotation + label for TPM == 0 or < 1
  vertices <- combined_mRNA %>%
    distinct(name = GN_ID, Gene_Name) %>%
    arrange(name) %>%
    left_join(mRNA_summary, by = c("name" = "GN_ID")) %>%
    mutate(type = "gene") %>%
    bind_rows(
      combined_mRNA %>%
        distinct(name = Organ) %>%
        mutate(Gene_Name = NA, average_TPMs = NA, Organ = name, gene_type = NA, type = "organ")
    ) %>%
    group_by(name) %>%
    mutate(label = case_when(
      any(average_TPMs == 0, na.rm = TRUE) & average_TPMs == 0 ~ Gene_Name,
      !any(average_TPMs == 0, na.rm = TRUE) & average_TPMs < 1 ~ Gene_Name,
      TRUE ~ NA_character_
    )) %>%
    mutate(label = if_else(row_number() == 1, label, NA_character_)) %>%
    ungroup()
  
  # Build graph and annotate node categories
  g <- tbl_graph(nodes = vertices, edges = edges, directed = FALSE) %>%
    activate(nodes) %>%
    mutate(
      degree        = centrality_degree(),
      node_category = case_when(
        type == "organ"            ~ "Organ",
        Organ == "Inter-organ"     ~ "Inter-organ Gene",
        gene_type == "intra_organ" ~ "Intra-organ Gene",
        TRUE                       ~ "Other"
      )
    )
  return(list(graph = g, data = combined_mRNA))
}

# ---- 4. mRNA Network Plot Function ----

# Generates kk layout plots: with organ labels + dot labels, with organ labels only,
# and no labels; saves all three as PDFs
plot_mRNA_network <- function(g, comparison_label, file_prefix) {
  set.seed(42)
  
  # Shared node size scale
  node_size_scale <- scale_size_manual(
    values = c("Organ" = 10, "Inter-organ Gene" = 3, "Intra-organ Gene" = 3),
    name   = "Node Type"
  )
  
  # Shared mRNA color gradient
  tpm_color_scale <- scale_color_gradientn(
    colors   = c("blue", "grey", "red", "green"),
    name     = "log10(levels of mRNA)",
    na.value = "grey50",
    guide    = guide_colorbar(title.position = "top", barwidth = 1, barheight = 8)
  )
  
  shared_theme <- theme_void() +
    labs(title    = paste("Inter-Organ Gene Network (", comparison_label, ")", sep = ""),
         subtitle = "Large nodes are organs;\nGenes are colored by levels of mRNA.") +
    theme(legend.position = "right", legend.box = "vertical",
          plot.title    = element_text(hjust = 0.5, face = "bold", size = 24),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          plot.margin   = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          legend.title  = element_text(size = 18),
          legend.text   = element_text(size = 14))
  
  # Base plot: organ nodes (grey) + TPM-colored gene nodes
  base_plot <- ggraph(g, layout = 'kk') +
    geom_edge_link(alpha = 0.1) +
    # Organ nodes — fixed grey color
    geom_node_point(data = function(x) filter(x, type == "organ"),
                    aes(size = node_category), color = "grey50") +
    # Genes with TPM > 0: colored by log10(TPM)
    ggnewscale::new_scale_color() +
    geom_node_point(data = function(x) filter(x, type == "gene" & average_TPMs > 0),
                    aes(color = log10(average_TPMs), size = node_category)) +
    tpm_color_scale +
    # Genes with TPM == 0: black
    ggnewscale::new_scale_color() +
    geom_node_point(data = function(x) filter(x, type == "gene" & average_TPMs == 0),
                    aes(color = "Zero value", size = node_category)) +
    scale_color_manual(values = c("Zero value" = "black"), name = "mRNA undetectable and NA") +
    node_size_scale +
    shared_theme
  
  # Version 1: organ name labels + dot labels (TPM == 0 and TPM < 1)
  p_labeldots <- base_plot +
    geom_node_text(data = function(x) filter(x, type == "organ"),
                   aes(label = name), repel = TRUE, size = 5, fontface = "bold") +
    geom_text_repel(data = function(x) filter(x, !is.na(label)),
                    aes(x = x, y = y, label = label),
                    size = 3, color = "black", max.overlaps = Inf)
  
  print(p_labeldots)
  ggsave(paste0(file_prefix, "_labeldots_kk.pdf"),
         plot = p_labeldots, width = 14400, height = 12800, units = "px", dpi = 600)
  
  # Version 2: organ name labels only (no dot labels)
  p_label <- base_plot +
    geom_node_text(data = function(x) filter(x, type == "organ"),
                   aes(label = name), repel = TRUE, size = 5, fontface = "bold")
  
  print(p_label)
  ggsave(paste0(file_prefix, "_kk.pdf"),
         plot = p_label, width = 14400, height = 12800, units = "px", dpi = 600)
  
  # Version 3: no labels at all
  print(base_plot)
  ggsave(paste0(file_prefix, "_nolabel_kk.pdf"),
         plot = base_plot, width = 14400, height = 12800, units = "px", dpi = 600)
}

# ---- 5. Generate All Three Comparison Networks ----

result_my <- build_mRNA_graph(combined_abundance_All_middle_vs_young, mRNA_All)
write.csv(result_my$graph %>% activate(nodes) %>% as_tibble(), "mRNA_middle_young.csv", row.names = FALSE)
plot_mRNA_network(result_my$graph, "Middle vs. Young", "mRNA_colored_network_middle_young")

result_oy <- build_mRNA_graph(combined_abundance_All_old_vs_young, mRNA_All)
write.csv(result_oy$graph %>% activate(nodes) %>% as_tibble(), "mRNA_old_young.csv", row.names = FALSE)
plot_mRNA_network(result_oy$graph, "Old vs. Young", "mRNA_colored_network_old_young")

result_om <- build_mRNA_graph(combined_abundance_All_old_vs_middle, mRNA_All)
write.csv(result_om$graph %>% activate(nodes) %>% as_tibble(), "mRNA_old_middle.csv", row.names = FALSE)
plot_mRNA_network(result_om$graph, "Old vs. Middle", "mRNA_colored_network_old_middle")