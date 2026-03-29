# ============================================================
# Project: TMT Proteomics — Multi-organ Mouse Aging
# Script:  Organs_abundance_upset_plot.R
# Purpose: UpSet plots showing overlap of significant proteins
#          across organs for each age comparison; also exports
#          inter-organ intersection tables as CSV
# Input:   combined_abundance_All.csv (from Organs_abundance.R)
# Output:  upset_plot_<comparison>_colored.pdf,
#          <comparison>_intersections.csv
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# ============================================================

library(UpSetR)
library(dplyr)

# NOTE: Requires combined_abundance_All.csv in the working directory
# setwd("path/to/your/data")

# ---- 1. Load Data ----

combined_abundance_All <- read.csv("combined_abundance_All.csv")

combined_abundance_All_old_vs_young    <- combined_abundance_All %>%
  filter(Comparison == "Old vs. Young",    !is.na(Gene_Name))
combined_abundance_All_old_vs_middle   <- combined_abundance_All %>%
  filter(Comparison == "Old vs. Middle",   !is.na(Gene_Name))
combined_abundance_All_middle_vs_young <- combined_abundance_All %>%
  filter(Comparison == "Middle vs. Young", !is.na(Gene_Name), !is.na(Abundance_mean_young))

# ---- 2. Shared Organ Color Palette ----

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

# ---- 3. Helper Functions ----

# Build a named list of GN_IDs per organ for one comparison's data
create_organ_lists <- function(data) {
  organs <- unique(data$Organ)
  organ_lists <- setNames(
    lapply(organs, function(organ) {
      data %>% filter(Organ == organ) %>% pull(GN_ID) %>% unique()
    }),
    organs
  )
  return(organ_lists)
}

# Generate an UpSet plot with per-organ color queries
custom_upset_plot <- function(data_lists, title) {
  ordered_colors <- organ_colors[names(data_lists)]
  
  upset(
    fromList(data_lists),
    nsets           = length(data_lists),
    order.by        = "freq",
    main.bar.color  = "#4a4a4a",
    matrix.color    = "black",
    text.scale      = 1.5,
    point.size      = 3,
    line.size       = 1,
    mainbar.y.label = "Intersection Size",
    sets.x.label    = "Proteins per Organ",
    set_size.show   = FALSE,
    set_size.angles = 0,
    show.numbers    = TRUE,
    shade.color     = "#f0f0f0",
    # Color bars by organ
    queries = lapply(names(data_lists), function(organ) {
      list(
        query       = intersects,
        params      = list(organ),
        color       = ordered_colors[organ],
        active      = TRUE,
        show.legend = TRUE
      )
    })
  )
}

# Save intersection details (proteins in > 1 organ) to CSV
save_intersection_details <- function(data_lists, filename) {
  if (length(data_lists) == 0 || all(sapply(data_lists, length) == 0)) {
    warning(paste("Skipping", filename, "- No data available."))
    return(NULL)
  }
  
  # Build binary membership matrix
  all_proteins <- unique(unlist(data_lists))
  upset_matrix <- do.call(cbind, lapply(data_lists, function(organ_genes) {
    as.integer(all_proteins %in% organ_genes)
  }))
  rownames(upset_matrix) <- all_proteins
  
  intersections <- apply(upset_matrix, 1, function(row) {
    paste(colnames(upset_matrix)[as.logical(row)], collapse = "|")
  })
  
  result <- data.frame(
    Proteins       = rownames(upset_matrix),
    Organs         = intersections,
    Number_of_Organs = rowSums(upset_matrix, na.rm = TRUE)
  )
  
  # Keep only proteins present in more than one organ
  result_filtered <- result[result$Number_of_Organs > 1, ]
  write.csv(result_filtered, filename, row.names = FALSE)
}

# ---- 4. Generate UpSet Plots and Intersection Tables ----

old_young_lists    <- create_organ_lists(combined_abundance_All_old_vs_young)
old_middle_lists   <- create_organ_lists(combined_abundance_All_old_vs_middle)
middle_young_lists <- create_organ_lists(combined_abundance_All_middle_vs_young)

# Old vs. Young
pdf("upset_plot_old_vs_young_colored.pdf", width = 14, height = 10)
custom_upset_plot(old_young_lists, "Old vs. Young")
dev.off()

# Old vs. Middle
pdf("upset_plot_old_vs_middle_colored.pdf", width = 14, height = 10)
custom_upset_plot(old_middle_lists, "Old vs. Middle")
dev.off()

# Middle vs. Young
pdf("upset_plot_middle_vs_young_colored.pdf", width = 14, height = 10)
custom_upset_plot(middle_young_lists, "Middle vs. Young")
dev.off()

# Export intersection CSVs
save_intersection_details(old_young_lists,    "old_vs_young_intersections.csv")
save_intersection_details(old_middle_lists,   "old_vs_middle_intersections.csv")
save_intersection_details(middle_young_lists, "middle_vs_young_intersections.csv")