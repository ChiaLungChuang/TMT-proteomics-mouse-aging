# ============================================================
# Project: TMT Proteomics — Multi-organ Mouse Aging
# Script:  Organs_abundance.R
# Purpose: Cross-organ network graphs with nodes colored by
#          baseline protein abundance (log-transformed
#          Abundance_mean_young) for all three age comparisons
# Input:   combined_plot_data_mouseTMT_<tissue>.csv,
#          ratio_data_<tissue>_OY.csv (from per-tissue scripts)
# Output:  combined_abundance_All.csv,
#          inter_organ_abundance_colored_network_<comparison>_kk.pdf
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# ============================================================

library(dplyr)
library(tidyr)
library(ggraph)
library(igraph)
library(tidygraph)
library(scales)
library(ggplot2)
library(ggnewscale)

# NOTE: This script expects combined_plot_data_mouseTMT_<tissue>.csv and
# ratio_data_<tissue>_OY.csv files in the working directory,
# produced by the per-tissue scripts in 01_per-tissue-analysis/
# setwd("path/to/your/data")

# ---- 1. Load Per-Tissue Data and Attach Abundance ----

# Extract Abundance_mean_young (baseline abundance at 6M) from each tissue
ratio_selected_TA  <- read.csv("ratio_data_TA_OY.csv")  %>% select(GN_ID, Abundance_mean_young)
ratio_selected_He  <- read.csv("ratio_data_He_OY.csv")  %>% select(GN_ID, Abundance_mean_young)
ratio_selected_CB  <- read.csv("ratio_data_CB_OY.csv")  %>% select(GN_ID, Abundance_mean_young)
ratio_selected_HC  <- read.csv("ratio_data_HC_OY.csv")  %>% select(GN_ID, Abundance_mean_young)
ratio_selected_NC  <- read.csv("ratio_data_NC_OY.csv")  %>% select(GN_ID, Abundance_mean_young)
ratio_selected_So  <- read.csv("ratio_data_So_OY.csv")  %>% select(GN_ID, Abundance_mean_young)
ratio_selected_Li  <- read.csv("ratio_data_Li_OY.csv")  %>% select(GN_ID, Abundance_mean_young)
ratio_selected_WAT <- read.csv("ratio_data_WAT_OY.csv") %>% select(GN_ID, Abundance_mean_young)

combined_plot_data_TA  <- read.csv("combined_plot_data_mouseTMT_TA.csv")
combined_plot_data_He  <- read.csv("combined_plot_data_mouseTMT_He.csv")
combined_plot_data_CB  <- read.csv("combined_plot_data_mouseTMT_CB.csv")
combined_plot_data_HC  <- read.csv("combined_plot_data_mouseTMT_HC.csv")
combined_plot_data_NC  <- read.csv("combined_plot_data_mouseTMT_NC.csv")
combined_plot_data_So  <- read.csv("combined_plot_data_mouseTMT_So.csv")
combined_plot_data_Li  <- read.csv("combined_plot_data_mouseTMT_Li.csv")
combined_plot_data_WAT <- read.csv("combined_plot_data_mouseTMT_WAT.csv")

# Join abundance to combined plot data and label each tissue
combined_abundance_TA  <- combined_plot_data_TA  %>% left_join(ratio_selected_TA,  by = "GN_ID") %>% mutate(Organ = "Tibialis Anterior")
combined_abundance_He  <- combined_plot_data_He  %>% left_join(ratio_selected_He,  by = "GN_ID") %>% mutate(Organ = "Heart")
combined_abundance_CB  <- combined_plot_data_CB  %>% left_join(ratio_selected_CB,  by = "GN_ID") %>% mutate(Organ = "Cerebellum")
combined_abundance_HC  <- combined_plot_data_HC  %>% left_join(ratio_selected_HC,  by = "GN_ID") %>% mutate(Organ = "Hippocampus")
combined_abundance_NC  <- combined_plot_data_NC  %>% left_join(ratio_selected_NC,  by = "GN_ID") %>% mutate(Organ = "Neocortex")
combined_abundance_So  <- combined_plot_data_So  %>% left_join(ratio_selected_So,  by = "GN_ID") %>% mutate(Organ = "Soleus")
combined_abundance_Li  <- combined_plot_data_Li  %>% left_join(ratio_selected_Li,  by = "GN_ID") %>% mutate(Organ = "Liver")
combined_abundance_WAT <- combined_plot_data_WAT %>% left_join(ratio_selected_WAT, by = "GN_ID") %>% mutate(Organ = "White Adipose Tissue")

combined_abundance_All <- bind_rows(
  combined_abundance_TA, combined_abundance_He, combined_abundance_CB,
  combined_abundance_HC, combined_abundance_NC, combined_abundance_So,
  combined_abundance_Li, combined_abundance_WAT
)

write.csv(combined_abundance_All, "combined_abundance_All.csv", row.names = FALSE)
combined_abundance_All <- read.csv("combined_abundance_All.csv")

# Split by comparison
combined_abundance_All_old_vs_young    <- combined_abundance_All %>% filter(Comparison == "Old vs. Young")
combined_abundance_All_old_vs_middle   <- combined_abundance_All %>% filter(Comparison == "Old vs. Middle")
combined_abundance_All_middle_vs_young <- combined_abundance_All %>% filter(Comparison == "Middle vs. Young")

# ---- 2. Graph Builder Function ----

# Builds a tbl_graph for one comparison with Abundance_mean_young as node attribute
build_abundance_graph <- function(data) {
  data <- data %>% filter(!is.na(Gene_Name))
  
  # Per-gene summary: take max abundance; classify inter/intra-organ
  gene_abundance <- data %>%
    group_by(GN_ID) %>%
    summarise(
      Abundance_mean_young = max(Abundance_mean_young, na.rm = TRUE),
      Organ     = if (n_distinct(Organ) > 1) "Inter-organ" else first(Organ),
      gene_type = if (n_distinct(Organ) > 1) "inter_organ"  else "intra_organ",
      .groups = "drop"
    )
  
  edges <- data %>% select(from = GN_ID, to = Organ) %>% distinct()
  
  vertices <- data %>%
    distinct(name = GN_ID, Gene_Name) %>%
    arrange(name) %>%
    left_join(gene_abundance, by = c("name" = "GN_ID")) %>%
    mutate(type = "gene") %>%
    bind_rows(
      data %>%
        distinct(name = Organ) %>%
        mutate(Gene_Name = NA, Abundance_mean_young = NA, Organ = name, gene_type = NA, type = "organ")
    )
  
  g <- tbl_graph(nodes = vertices, edges = edges, directed = FALSE) %>%
    activate(nodes) %>%
    mutate(
      degree        = centrality_degree(),
      node_category = case_when(
        type == "organ"            ~ "Organ",
        Organ == "Inter-organ"     ~ "Inter-organ Gene",
        gene_type == "intra_organ" ~ "Intra-organ Gene"
      )
    )
  return(g)
}

# ---- 3. Abundance Network Plot Function ----

# Generates kk layout plots with and without organ labels; saves as PDFs
plot_abundance_network <- function(g, comparison_label, file_prefix) {
  set.seed(42)
  
  # Shared layers
  base_layers <- list(
    geom_edge_link(alpha = 0.1),
    geom_node_point(
      data = function(x) filter(x, type == "organ"),
      aes(size = node_category), color = "grey50"
    ),
    ggnewscale::new_scale_color(),
    geom_node_point(
      data = function(x) filter(x, type == "gene"),
      aes(color = log(Abundance_mean_young), size = node_category)
    ),
    scale_color_gradientn(
      colors   = c("blue", "grey", "red", "green"),
      name     = "Log(Abundance)",
      na.value = "grey50",
      guide    = guide_colorbar(title.position = "top", barwidth = 1, barheight = 8)
    ),
    scale_size_manual(
      values = c("Organ" = 10, "Inter-organ Gene" = 3, "Intra-organ Gene" = 3),
      name   = "Node Type"
    ),
    theme_void(),
    labs(
      title    = paste("Inter-Organ Gene Network (", comparison_label, ")", sep = ""),
      subtitle = "Large nodes are organs;\nGenes are colored by total protein levels at 6 months (Young)."
    ),
    theme(legend.position = "right", legend.box = "vertical",
          plot.title    = element_text(hjust = 0.5, face = "bold", size = 24),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          plot.margin   = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          legend.title  = element_text(size = 18),
          legend.text   = element_text(size = 14))
  )
  
  # With organ labels
  p_label <- ggraph(g, layout = 'kk')
  for (l in base_layers) p_label <- p_label + l
  p_label <- p_label +
    geom_node_text(data = function(x) filter(x, type == "organ"),
                   aes(label = name), repel = TRUE, size = 5, fontface = "bold")
  
  print(p_label)
  ggsave(paste0(file_prefix, "_kk.pdf"),
         plot = p_label, width = 14400, height = 12800, units = "px", dpi = 600)
  
  # Without organ labels
  p_nolabel <- ggraph(g, layout = 'kk')
  for (l in base_layers) p_nolabel <- p_nolabel + l
  
  print(p_nolabel)
  ggsave(paste0(file_prefix, "_nolabel_kk.pdf"),
         plot = p_nolabel, width = 14400, height = 12800, units = "px", dpi = 600)
}

# ---- 4. Generate All Three Comparison Networks ----

graph_old_vs_middle   <- build_abundance_graph(combined_abundance_All_old_vs_middle   %>% filter(!is.na(Gene_Name)))
graph_old_vs_young    <- build_abundance_graph(combined_abundance_All_old_vs_young    %>% filter(!is.na(Gene_Name)))
graph_middle_vs_young <- build_abundance_graph(combined_abundance_All_middle_vs_young %>% filter(!is.na(Gene_Name), !is.na(Abundance_mean_young)))

plot_abundance_network(graph_old_vs_middle,   "Old vs. Middle",   "inter_organ_abundance_colored_network_old_middle")
plot_abundance_network(graph_old_vs_young,    "Old vs. Young",    "inter_organ_abundance_colored_network_old_young")
plot_abundance_network(graph_middle_vs_young, "Middle vs. Young", "inter_organ_abundance_colored_network_middle_young")