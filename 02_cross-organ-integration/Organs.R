# ============================================================
# Project: TMT Proteomics — Multi-organ Mouse Aging
# Script:  Organs.R
# Purpose: Cross-organ integration — combines per-tissue outputs,
#          classifies intra- vs. inter-organ proteins, generates
#          whole-proteome scatter plots, and builds organ-gene
#          network graphs (Fruchterman-Reingold and Kamada-Kawai
#          layouts) for all three age comparisons
# Input:   combined_plot_data_mouseTMT_<tissue>.csv (from per-tissue scripts)
# Output:  combined_plot_data_mouseTMT_All.csv, network graph PDFs,
#          protein distribution CSVs
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# ============================================================

library(ggplot2)
library(dplyr)
library(ggrepel)
library(igraph)
library(ggraph)
library(tidygraph)
library(extrafont)
library(ggpmisc)

# NOTE: Set your working directory to the folder containing all
# combined_plot_data_mouseTMT_*.csv files produced by per-tissue scripts
# setwd("path/to/your/data")

# ---- 1. Load Per-Tissue Combined Plot Data ----

# NOTE: combined_plot_data_mouseTMT_He.csv and _Se.csv (Heart and Serum)
# are loaded here; ensure these files exist in your working directory.
# Heart was generated with variable suffix "HE" in the per-tissue script
# but saved as "He" here — verify filename consistency in your environment.

combined_plot_data_TA  <- read.csv("combined_plot_data_mouseTMT_TA.csv")
combined_plot_data_He  <- read.csv("combined_plot_data_mouseTMT_He.csv")
combined_plot_data_CB  <- read.csv("combined_plot_data_mouseTMT_CB.csv")
combined_plot_data_HC  <- read.csv("combined_plot_data_mouseTMT_HC.csv")
combined_plot_data_NC  <- read.csv("combined_plot_data_mouseTMT_NC.csv")
combined_plot_data_So  <- read.csv("combined_plot_data_mouseTMT_So.csv")
combined_plot_data_Li  <- read.csv("combined_plot_data_mouseTMT_Li.csv")
combined_plot_data_Se  <- read.csv("combined_plot_data_mouseTMT_Se.csv")
combined_plot_data_WAT <- read.csv("combined_plot_data_mouseTMT_WAT.csv")

ratio_data_TA_OY  <- read.csv("ratio_data_TA_OY.csv")
ratio_data_He_OY  <- read.csv("ratio_data_He_OY.csv")
ratio_data_CB_OY  <- read.csv("ratio_data_CB_OY.csv")
ratio_data_HC_OY  <- read.csv("ratio_data_HC_OY.csv")
ratio_data_NC_OY  <- read.csv("ratio_data_NC_OY.csv")
ratio_data_So_OY  <- read.csv("ratio_data_So_OY.csv")
ratio_data_Li_OY  <- read.csv("ratio_data_Li_OY.csv")
ratio_data_WAT_OY <- read.csv("ratio_data_WAT_OY.csv")

# ---- 2. Add Organ Labels and Combine ----

combined_plot_data_TA  <- combined_plot_data_TA  %>% mutate(Organ = "Tibialis Anterior")
combined_plot_data_He  <- combined_plot_data_He  %>% mutate(Organ = "Heart")
combined_plot_data_CB  <- combined_plot_data_CB  %>% mutate(Organ = "Cerebellum")
combined_plot_data_HC  <- combined_plot_data_HC  %>% mutate(Organ = "Hippocampus")
combined_plot_data_NC  <- combined_plot_data_NC  %>% mutate(Organ = "Neocortex")
combined_plot_data_So  <- combined_plot_data_So  %>% mutate(Organ = "Soleus")
combined_plot_data_Li  <- combined_plot_data_Li  %>% mutate(Organ = "Liver")
combined_plot_data_Se  <- combined_plot_data_Se  %>% mutate(Organ = "Serum")
combined_plot_data_WAT <- combined_plot_data_WAT %>% mutate(Organ = "White Adipose Tissue")

combined_plot_data_mouseTMT_All <- bind_rows(
  combined_plot_data_TA, combined_plot_data_He, combined_plot_data_CB,
  combined_plot_data_HC, combined_plot_data_NC, combined_plot_data_So,
  combined_plot_data_Li, combined_plot_data_Se, combined_plot_data_WAT
)

write.csv(combined_plot_data_mouseTMT_All, "combined_plot_data_mouseTMT_All.csv", row.names = FALSE)
combined_plot_data_mouseTMT_All <- read.csv("combined_plot_data_mouseTMT_All.csv")

# Strip FC and label columns for a leaner all-organs table
plot_data_mouseTMT_All <- combined_plot_data_mouseTMT_All[,
                                                          !names(combined_plot_data_mouseTMT_All) %in% c("Abundance_FC", "Ratio_FC", "Label_top", "Label")
]
write.csv(plot_data_mouseTMT_All, "plot_data_mouseTMT_All.csv", row.names = FALSE)

# Remove serum for tissue-level analyses
combined_plot_data_mouseTMT_no_serum <- combined_plot_data_mouseTMT_All %>%
  filter(Organ != "Serum")

# ---- 3. Protein Distribution Across Organs ----

# Count how many distinct organs each protein appears in per comparison group
summarized_data <- combined_plot_data_mouseTMT_All %>%
  group_by(GN_ID, Gene_Name, Comparison) %>%
  summarise(
    Organ_Count = n_distinct(Organ),
    Organs      = paste(unique(Organ), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(Gene_Type = if_else(Organ_Count == 1, "intra-organ", "inter-organ"))

# Filter to inter-organ proteins (>= 2 organs) per comparison
data_old_vs_young    <- summarized_data %>% filter(Comparison == "Old vs. Young"    & Organ_Count >= 2)
data_old_vs_middle   <- summarized_data %>% filter(Comparison == "Old vs. Middle"   & Organ_Count >= 2)
data_middle_vs_young <- summarized_data %>% filter(Comparison == "Middle vs. Young" & Organ_Count >= 2)

# All proteins per comparison (including intra-organ)
data_old_vs_young_v2    <- summarized_data %>% filter(Comparison == "Old vs. Young")
data_old_vs_middle_v2   <- summarized_data %>% filter(Comparison == "Old vs. Middle")
data_middle_vs_young_v2 <- summarized_data %>% filter(Comparison == "Middle vs. Young")

write.csv(data_old_vs_young,       "data_old_vs_young.csv",       row.names = FALSE)
write.csv(data_old_vs_middle,      "data_old_vs_middle.csv",      row.names = FALSE)
write.csv(data_middle_vs_young,    "data_middle_vs_young.csv",    row.names = FALSE)
write.csv(summarized_data,         "data_all_comparison.csv",     row.names = FALSE)
write.csv(data_old_vs_young_v2,    "data_old_vs_young_v2.csv",    row.names = FALSE)
write.csv(data_old_vs_middle_v2,   "data_old_vs_middle_v2.csv",   row.names = FALSE)
write.csv(data_middle_vs_young_v2, "data_middle_vs_young_v2.csv", row.names = FALSE)

# Dot plots: proteins by number of organs they appear in
create_dot_plot <- function(data, title) {
  ggplot(data, aes(x = Gene_Name, y = Organ_Count)) +
    geom_point(aes(size = Organ_Count), color = "blue", alpha = 0.7) +
    geom_text(aes(label = Organs), vjust = -0.5, size = 3, angle = 90, hjust = 1) +
    labs(x = "Gene Name", y = "Number of Organs", title = title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(size = "none")
}

dot_old_vs_young    <- create_dot_plot(data_old_vs_young,    "Old vs. Young: Organ Count >= 2")
dot_old_vs_middle   <- create_dot_plot(data_old_vs_middle,   "Old vs. Middle: Organ Count >= 2")
dot_middle_vs_young <- create_dot_plot(data_middle_vs_young, "Middle vs. Young: Organ Count >= 2")

print(dot_old_vs_young)
ggsave("plot_old_vs_young.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

# Bar plots: proteins by organ count, colored by comparison
create_bar_plot <- function(data, title) {
  ggplot(data, aes(x = reorder(Gene_Name, -Organ_Count), y = Organ_Count, fill = Comparison)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Organs), vjust = -0.5, size = 3, angle = 90, hjust = 1) +
    labs(x = "Gene Name", y = "Number of Organs", title = title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = c(
      "Old vs. Young"    = "purple",
      "Middle vs. Young" = "maroon2",
      "Old vs. Middle"   = "green"
    ))
}

bar_old_vs_young    <- create_bar_plot(data_old_vs_young,    "Old vs. Young: Organ Count")
bar_old_vs_middle   <- create_bar_plot(data_old_vs_middle,   "Old vs. Middle: Organ Count")
bar_middle_vs_young <- create_bar_plot(data_middle_vs_young, "Middle vs. Young: Organ Count")

print(bar_old_vs_young)
ggsave("bar_old_vs_young.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

# ---- 4. Whole-Proteome Scatter Plots (All Organs) ----

# Shared organ color palette used throughout network graphs
organ_colors <- c(
  "Cerebellum"          = "#c7a6db",
  "Heart"               = "#c93b92",
  "Hippocampus"         = "#55397f",
  "Liver"               = "#4d9fdb",
  "Neocortex"           = "#8a68ca",
  "Soleus"              = "#a8192a",
  "Tibialis Anterior"   = "#e991b7",
  "White Adipose Tissue"= "#dfc41b"
)

# With serum (all organs)
mouseTMT_whole_proteome <- ggplot(combined_plot_data_mouseTMT_All,
                                  aes(x = Abundance_FC, y = Ratio_FC, size = Size, color = Organ)) +
  geom_point(alpha = 0.3) +
  scale_size_continuous(range = c(1,9), name = "P-Value Size (-log10)",
                        breaks = c(1.3,2,3,4,5), limit = c(1.3,6),
                        labels = c("1.30(p<0.05)","2","3","4","5")) +
  labs(title   = "MouseTMT Whole Proteome: Abundance FC vs. Ratio FC",
       x       = "Log2 Abundance fold change",
       y       = "Log2 fold change in Insoluble/Soluble",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  facet_wrap(~Comparison)

print(mouseTMT_whole_proteome)
ggsave("mouseTMT_whole_proteome.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

# Without serum, organ-specific colors
mouseTMT_whole_proteome_no_serum <- ggplot(combined_plot_data_mouseTMT_no_serum,
                                           aes(x = Abundance_FC, y = Ratio_FC, size = Size, color = Organ)) +
  geom_point(alpha = 0.3) +
  scale_size_continuous(range = c(1,9), name = "P-Value Size (-log10)",
                        breaks = c(1.3,2,3,4,5), limit = c(1.3,6),
                        labels = c("1.30(p<0.05)","2","3","4","5")) +
  scale_color_manual(values = organ_colors) +
  labs(title   = "MouseTMT Whole Proteome: Abundance FC vs. Ratio FC (no serum)",
       x       = "Log2 Abundance fold change",
       y       = "Log2 fold change in Insoluble/Soluble",
       caption = "Point size based on statistical significance of Age difference in insolubility") +
  facet_wrap(~Comparison)

# Add R-squared annotation per facet
mouseTMT_whole_proteome_no_serum_R <- mouseTMT_whole_proteome_no_serum +
  stat_poly_eq(aes(label = paste(after_stat(rr.label))),
               rr.digits = 3, label.x = 1, label.y = 1.4,
               parse = TRUE, color = "black", show.legend = FALSE)

print(mouseTMT_whole_proteome_no_serum_R)
ggsave("mouseTMT_whole_proteome_no_serum_R.pdf", width = 7200, height = 6400, units = "px", dpi = 600)

# ---- 5. Split Data by Comparison for Network Graphs ----

data_no_serum_old_vs_middle   <- combined_plot_data_mouseTMT_no_serum %>% filter(Comparison == "Old vs. Middle")
data_no_serum_old_vs_young    <- combined_plot_data_mouseTMT_no_serum %>% filter(Comparison == "Old vs. Young")
data_no_serum_middle_vs_young <- combined_plot_data_mouseTMT_no_serum %>% filter(Comparison == "Middle vs. Young")

# ---- 6. Network Graph Helper Functions ----

# Build a tbl_graph from a single-comparison dataframe
build_organ_gene_graph <- function(data) {
  data <- data %>% filter(!is.na(GN_ID), !is.na(Organ), !is.na(Gene_Name))
  
  gene_colors <- data %>%
    group_by(GN_ID) %>%
    summarise(
      Organ     = if (n_distinct(Organ) > 1) "Inter-organ" else first(Organ),
      gene_type = if (n_distinct(Organ) > 1) "inter_organ"  else "intra_organ",
      .groups = "drop"
    )
  
  edges <- data %>% select(from = GN_ID, to = Organ)
  
  vertices <- data %>%
    distinct(name = GN_ID, Gene_Name) %>%
    arrange(name) %>%
    left_join(gene_colors, by = c("name" = "GN_ID")) %>%
    mutate(type = "gene") %>%
    bind_rows(
      data %>%
        distinct(name = Organ) %>%
        mutate(Organ = name, type = "organ", gene_type = NA, Gene_Name = NA)
    )
  
  g <- tbl_graph(nodes = vertices, edges = edges, directed = FALSE) %>%
    activate(nodes) %>%
    mutate(
      degree    = centrality_degree(),
      node_type = case_when(
        type == "organ"            ~ "Organ",
        Organ == "Inter-organ"     ~ "Inter-organ",
        TRUE                       ~ "Intra-organ"
      )
    )
  return(g)
}

# Shared network plot layers
network_layers <- function(g, organ_col_palette, title_str, layout_str, layout_args = list()) {
  layout_call <- c(list(graph = g, layout = layout_str), layout_args)
  p <- do.call(ggraph, layout_call) +
    geom_edge_link(alpha = 0.1) +
    geom_node_point(aes(color = Organ, size = node_type), show.legend = TRUE) +
    scale_color_manual(values = c(organ_col_palette, "Inter-organ" = "grey15")) +
    scale_size_manual(values = c("Organ" = 10, "Inter-organ" = 6, "Intra-organ" = 3), name = "Node Type") +
    theme_void() +
    labs(title    = title_str,
         subtitle = "Large nodes are organs, medium grey nodes are inter-organ genes,\nsmall colored nodes are intra-organ genes.") +
    theme(legend.position  = "right",
          legend.box        = "vertical",
          plot.title        = element_text(hjust = 0.5, face = "bold", size = 24),
          plot.subtitle     = element_text(hjust = 0.5, size = 16),
          plot.margin       = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          legend.title      = element_text(size = 18),
          legend.text       = element_text(size = 14))
  return(p)
}

# ---- 7. All-Comparisons Network (combined_plot_data_mouseTMT_no_serum) ----

# Build graph across all comparisons combined
unique_genes_all <- combined_plot_data_mouseTMT_no_serum %>%
  filter(!is.na(GN_ID), !is.na(Organ)) %>%
  group_by(GN_ID, Gene_Name, Organ) %>%
  summarise(Comparisons = paste(unique(Comparison), collapse = ", "), .groups = "drop")

gene_colors_all <- unique_genes_all %>%
  group_by(GN_ID) %>%
  summarise(
    Organ     = if (n_distinct(Organ) > 1) "Inter-organ" else first(Organ),
    gene_type = if (n_distinct(Organ) > 1) "inter_organ"  else "intra_organ",
    .groups = "drop"
  )

edges_all <- unique_genes_all %>% select(from = GN_ID, to = Organ)

vertices_all <- unique_genes_all %>%
  distinct(name = GN_ID, Gene_Name) %>%
  arrange(name) %>%
  left_join(gene_colors_all, by = c("name" = "GN_ID")) %>%
  mutate(type = "gene") %>%
  bind_rows(
    unique_genes_all %>%
      distinct(name = Organ) %>%
      mutate(Organ = name, type = "organ", gene_type = NA, Gene_Name = NA, Gene_Name = NA)
  )

graph_all <- tbl_graph(nodes = vertices_all, edges = edges_all, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(
    degree    = centrality_degree(),
    node_type = case_when(
      type == "organ"        ~ "Organ",
      Organ == "Inter-organ" ~ "Inter-organ",
      TRUE                   ~ "Intra-organ"
    )
  )

set.seed(42)

# fr — with inter-organ gene labels
inter_organ_genes_colored_network_fr2 <- ggraph(graph_all, layout = 'fr', niter = 1000) +
  geom_edge_link(alpha = 0.1) +
  geom_node_point(aes(color = Organ, size = node_type), show.legend = TRUE) +
  scale_color_manual(values = c(organ_colors, "Inter-organ" = "grey15")) +
  scale_size_manual(values = c("Organ" = 10, "Inter-organ" = 6, "Intra-organ" = 3), name = "Node Type") +
  theme_void() +
  labs(title    = "Inter-Organ Gene Network (All comparisons)",
       subtitle = "Large nodes are organs, medium grey nodes are inter-organ genes,\nsmall colored nodes are intra-organ genes.") +
  theme(legend.position = "right", legend.box = "vertical",
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 24),
        plot.subtitle = element_text(hjust = 0.5, size = 16),
        plot.margin   = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title  = element_text(size = 18), legend.text = element_text(size = 14)) +
  geom_text_repel(aes(x = x, y = y, label = ifelse(node_type == "Inter-organ", Gene_Name, "")),
                  size = 3, color = "black", max.overlaps = Inf,
                  box.padding = 0.5, point.padding = 0.2, segment.color = "grey50")

print(inter_organ_genes_colored_network_fr2)
ggsave("inter_organ_genes_colored_network_fr2.pdf", width = 14400, height = 12800, units = "px", dpi = 600)

# fr — no labels
inter_organ_genes_colored_network_fr_nolabel2 <- ggraph(graph_all, layout = 'fr', niter = 10000) +
  geom_edge_link(alpha = 0.1) +
  geom_node_point(aes(color = Organ, size = node_type), show.legend = TRUE, alpha = 0.5) +
  scale_color_manual(values = c(organ_colors, "Inter-organ" = "grey15")) +
  scale_size_manual(values = c("Organ" = 6, "Inter-organ" = 3, "Intra-organ" = 1), name = "Node Type") +
  theme_void() +
  labs(title    = "Inter-Organ Gene Network (All comparisons)",
       subtitle = "Large nodes are organs, medium grey nodes are inter-organ genes,\nsmall colored nodes are intra-organ genes.") +
  theme(legend.position = "right", legend.box = "vertical",
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 24),
        plot.subtitle = element_text(hjust = 0.5, size = 16),
        plot.margin   = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title  = element_text(size = 18), legend.text = element_text(size = 14))

print(inter_organ_genes_colored_network_fr_nolabel2)
ggsave("inter_organ_genes_colored_network_fr_nolabel2.pdf", width = 14400, height = 12800, units = "px", dpi = 600)

# kk — with labels
inter_organ_genes_colored_network_kk <- ggraph(graph_all, layout = 'kk', epsilon = 1e-10) +
  geom_edge_link(alpha = 0.1) +
  geom_node_point(aes(color = Organ, size = node_type), show.legend = TRUE) +
  scale_color_manual(values = c(organ_colors, "Inter-organ" = "grey15")) +
  scale_size_manual(values = c("Organ" = 10, "Inter-organ" = 6, "Intra-organ" = 3), name = "Node Type") +
  theme_void() +
  labs(title    = "Inter-Organ Gene Network (All comparisons)",
       subtitle = "Large nodes are organs, medium grey nodes are inter-organ genes,\nsmall colored nodes are intra-organ genes.") +
  theme(legend.position = "right", legend.box = "vertical",
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 24),
        plot.subtitle = element_text(hjust = 0.5, size = 16),
        plot.margin   = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title  = element_text(size = 18), legend.text = element_text(size = 14)) +
  geom_text_repel(aes(x = x, y = y, label = ifelse(node_type == "Inter-organ", Gene_Name, "")),
                  size = 3, color = "black", max.overlaps = Inf,
                  box.padding = 0.5, point.padding = 0.2, segment.color = "grey50")

print(inter_organ_genes_colored_network_kk)
ggsave("inter_organ_genes_colored_network_kk.pdf", width = 14400, height = 12800, units = "px", dpi = 600)

# kk — no labels
inter_organ_genes_colored_network_kk_nolabel2 <- ggraph(graph_all, layout = 'kk', epsilon = 1e-10) +
  geom_edge_link(alpha = 0.1) +
  geom_node_point(aes(color = Organ, size = node_type), show.legend = TRUE) +
  scale_color_manual(values = c(organ_colors, "Inter-organ" = "grey15")) +
  scale_size_manual(values = c("Organ" = 6, "Inter-organ" = 3, "Intra-organ" = 1), name = "Node Type") +
  theme_void() +
  labs(title    = "Inter-Organ Gene Network (All comparisons)",
       subtitle = "Large nodes are organs, medium grey nodes are inter-organ genes,\nsmall colored nodes are intra-organ genes.") +
  theme(legend.position = "right", legend.box = "vertical",
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 24),
        plot.subtitle = element_text(hjust = 0.5, size = 16),
        plot.margin   = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title  = element_text(size = 18), legend.text = element_text(size = 14))

print(inter_organ_genes_colored_network_kk_nolabel2)
ggsave("inter_organ_genes_colored_network_kk_nolabel2.pdf", width = 14400, height = 12800, units = "px", dpi = 600)

# ---- 8. Per-Comparison Network Graphs ----
# Build graphs for Old vs. Middle, Old vs. Young, Middle vs. Young

graph_old_vs_middle   <- build_organ_gene_graph(data_no_serum_old_vs_middle)
graph_old_vs_young    <- build_organ_gene_graph(data_no_serum_old_vs_young)
graph_middle_vs_young <- build_organ_gene_graph(data_no_serum_middle_vs_young)

write.csv(
  graph_old_vs_middle   %>% activate(nodes) %>% as_tibble(),
  "data_old_vs_middle_V2.csv", row.names = FALSE
)
write.csv(
  graph_old_vs_young    %>% activate(nodes) %>% as_tibble(),
  "data_old_vs_young_V2.csv", row.names = FALSE
)
write.csv(
  graph_middle_vs_young %>% activate(nodes) %>% as_tibble(),
  "data_middle_vs_young_V2.csv", row.names = FALSE
)

# Load fonts for PDF export
extrafont::font_import(prompt = FALSE)
loadfonts(device = "pdf")

# Helper: generate fr and kk plots (with and without labels) for one comparison
generate_network_plots <- function(g, organ_palette, comparison_label, file_prefix) {
  set.seed(42)
  
  for (layout_str in c("fr", "kk")) {
    layout_args <- if (layout_str == "kk") list(epsilon = 1e-8) else list(niter = 1000)
    
    # With inter-organ gene labels
    p_label <- do.call(ggraph, c(list(graph = g, layout = layout_str), layout_args)) +
      geom_edge_link(alpha = 0.1) +
      geom_node_point(aes(color = Organ, size = node_type), show.legend = TRUE) +
      scale_color_manual(values = c(organ_palette, "Inter-organ" = "grey15")) +
      scale_size_manual(values = c("Organ" = 8, "Inter-organ" = 5, "Intra-organ" = 3), name = "Node Type") +
      theme_void() +
      labs(title    = paste("Inter-Organ Gene Network (", comparison_label, ")", sep = ""),
           subtitle = "Large nodes are organs, medium grey nodes are inter-organ genes,\nsmall colored nodes are intra-organ genes.") +
      theme(legend.position = "right", legend.box = "vertical",
            plot.title    = element_text(hjust = 0.5, face = "bold", size = 24),
            plot.subtitle = element_text(hjust = 0.5, size = 16),
            plot.margin   = margin(0.5, 0.5, 0.5, 0.5, "cm"),
            legend.title  = element_text(size = 18), legend.text = element_text(size = 14)) +
      geom_text_repel(aes(x = x, y = y, label = ifelse(node_type == "Inter-organ", Gene_Name, "")),
                      size = 3, color = "black", family = "Arial", max.overlaps = Inf,
                      box.padding = 0.5, point.padding = 0.2, segment.color = "grey50")
    
    print(p_label)
    ggsave(paste0(file_prefix, "_", layout_str, ".pdf"),
           width = 14400, height = 12800, units = "px", dpi = 600)
    
    # Without labels
    p_nolabel <- p_label + guides(text = "none")  # Remove text_repel
    p_nolabel$layers <- p_nolabel$layers[!sapply(p_nolabel$layers, function(l)
      inherits(l$geom, "GeomTextRepel"))]
    
    print(p_nolabel)
    ggsave(paste0(file_prefix, "_", layout_str, "_nolabel.pdf"),
           width = 14400, height = 12800, units = "px", dpi = 600)
  }
}

generate_network_plots(graph_old_vs_middle,   organ_colors, "Old vs. Middle",   "inter_organ_genes_colored_network_old_middle")
generate_network_plots(graph_old_vs_young,    organ_colors, "Old vs. Young",    "inter_organ_genes_colored_network_old_young")
generate_network_plots(graph_middle_vs_young, organ_colors, "Middle vs. Young", "inter_organ_genes_colored_network_middle_young")

# ---- 9. Summary Statistics ----

summary_stats <- graph_all %>%
  activate(nodes) %>%
  as_tibble() %>%
  group_by(type, gene_type) %>%
  summarise(count = n(), avg_degree = mean(degree), .groups = "drop")

print("Summary Statistics (all comparisons combined):")
print(summary_stats)