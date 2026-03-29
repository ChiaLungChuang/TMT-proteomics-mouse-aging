# ============================================================
# Project: TMT Proteomics — Multi-organ Mouse Aging
# Script:  Insolubility_GOterm_upset.R
# Purpose: GO enrichment analysis (BP) comparing proteins with
#          increased vs. decreased insolubility across organs
#          and age comparisons; generates UpSet plots, Venn
#          diagrams, heatmaps, and per-organ GO term summaries
# Input:   combined_abundance_All.csv (from Organs_abundance.R)
# Output:  GO_*.csv, GO_*.pdf, GO_UpSet_Shared_*.pdf,
#          Organ_Specific_GO_Analysis_Summary.txt
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# Dependencies: BiocManager packages — install once before use:
#   BiocManager::install(c("clusterProfiler", "org.Mm.eg.db",
#                          "GO.db", "ComplexHeatmap"))
# ============================================================

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(VennDiagram)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library(ggVennDiagram)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(pdftools)
library(grid)

# NOTE: Set your working directory to the folder containing combined_abundance_All.csv
# setwd("path/to/your/data")

# ---- 1. Shared Organ Color Palette ----

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

organs <- names(organ_colors)

# ---- 2. Load Data ----

combined_abundance_All <- read.csv("combined_abundance_All.csv")

# Reshape to wide format for cross-comparison analysis
df_subset_All <- combined_abundance_All[, c("Gene_Name", "GN_ID", "Organ", "Comparison", "Ratio_FC")]
df_wide_All <- reshape(df_subset_All,
                       idvar    = c("Gene_Name", "GN_ID", "Organ"),
                       timevar  = "Comparison",
                       direction = "wide")
names(df_wide_All) <- gsub("Ratio_FC.", "", names(df_wide_All))

# ---- 3. GO Analysis Helper Functions ----

# Run GO enrichment for a gene list (BP ontology, BH correction, p < 0.05)
run_go_analysis <- function(gene_list, title) {
  if (length(gene_list) < 2) {
    cat("\nNot enough genes for GO analysis in", title, "\n")
    return(NULL)
  }
  
  go_results <- enrichGO(
    gene          = unique(gene_list),
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  if (is.null(go_results) || nrow(go_results@result) == 0) {
    cat("\nNo significant GO terms found for", title, "\n")
    return(NULL)
  }
  
  write.csv(go_results@result, paste0("GO_", gsub(" ", "_", title), ".csv"), row.names = FALSE)
  
  pdf(paste0("GO_", gsub(" ", "_", title), "_plots.pdf"), width = 12, height = 8)
  print(barplot(go_results, showCategory = 20) +
          labs(title = paste("GO Enrichment -", title)) + theme_minimal())
  print(dotplot(go_results, showCategory = 20) +
          labs(title = paste("GO Enrichment -", title)) + theme_minimal())
  dev.off()
  
  return(go_results)
}

# Extract gene symbols for a given GO term from an enrichResult object
get_proteins_in_term <- function(go_result, term) {
  idx <- which(go_result@result$Description == term)
  if (length(idx) > 0) {
    return(unlist(strsplit(go_result@result$geneID[idx], "/")))
  }
  return(character(0))
}

# ---- 4. Classify Proteins by Insolubility Direction ----

# Increased: positive Ratio_FC in any comparison
proteins_increase <- df_wide_All[
  (!is.na(df_wide_All$"Old vs. Young")    & df_wide_All$"Old vs. Young"    > 0) |
    (!is.na(df_wide_All$"Middle vs. Young") & df_wide_All$"Middle vs. Young" > 0) |
    (!is.na(df_wide_All$"Old vs. Middle")   & df_wide_All$"Old vs. Middle"   > 0),
]

# Decreased: negative Ratio_FC in any comparison
proteins_decrease <- df_wide_All[
  (!is.na(df_wide_All$"Old vs. Young")    & df_wide_All$"Old vs. Young"    < 0) |
    (!is.na(df_wide_All$"Middle vs. Young") & df_wide_All$"Middle vs. Young" < 0) |
    (!is.na(df_wide_All$"Old vs. Middle")   & df_wide_All$"Old vs. Middle"   < 0),
]

# Both: proteins showing increase in some comparisons and decrease in others
proteins_both <- df_wide_All[
  paste(df_wide_All$Gene_Name, df_wide_All$Organ) %in%
    intersect(
      paste(proteins_increase$Gene_Name, proteins_increase$Organ),
      paste(proteins_decrease$Gene_Name, proteins_decrease$Organ)
    ),
]

# Exclusive groups (no mixed direction)
proteins_increase_only <- proteins_increase[
  !(paste(proteins_increase$Gene_Name, proteins_increase$Organ) %in%
      paste(proteins_both$Gene_Name, proteins_both$Organ)),
]

proteins_decrease_only <- proteins_decrease[
  !(paste(proteins_decrease$Gene_Name, proteins_decrease$Organ) %in%
      paste(proteins_both$Gene_Name, proteins_both$Organ)),
]

# Save classified protein lists
write.csv(proteins_increase_only, "Proteins_Increase_Only.csv", row.names = FALSE)
write.csv(proteins_decrease_only, "Proteins_Decrease_Only.csv", row.names = FALSE)
write.csv(proteins_both,          "Proteins_Both_Patterns.csv",  row.names = FALSE)

# ---- 5. GO Analysis: Global Categories ----

cat("\nRunning global GO analysis by insolubility category...\n")
go_increase <- run_go_analysis(unique(proteins_increase_only$Gene_Name), "Increase_Only")
go_decrease <- run_go_analysis(unique(proteins_decrease_only$Gene_Name), "Decrease_Only")
go_both     <- run_go_analysis(unique(proteins_both$Gene_Name),          "Both_Patterns")

# ---- 6. Venn Diagram: Increased vs. Decreased Proteins ----

pdf("Protein_Categories_Venn.pdf", width = 10, height = 10)
venn_plot <- draw.pairwise.venn(
  area1      = length(unique(proteins_increase_only$Gene_Name)),
  area2      = length(unique(proteins_decrease_only$Gene_Name)),
  cross.area = length(unique(proteins_both$Gene_Name)),
  category   = c("Increased Insolubility", "Decreased Insolubility"),
  fill       = c("red", "blue"),
  alpha      = 0.5
)
grid.draw(venn_plot)
dev.off()

# ---- 7. Mixed-Pattern Plots per Organ ----

for (organ in unique(df_wide_All$Organ)) {
  organ_data <- proteins_both[proteins_both$Organ == organ, ]
  if (nrow(organ_data) > 0) {
    plot_data <- data.frame(
      Gene_Name  = rep(organ_data$Gene_Name, 3),
      Comparison = rep(c("Old vs. Young", "Middle vs. Young", "Old vs. Middle"), each = nrow(organ_data)),
      Ratio_FC   = c(organ_data$"Old vs. Young", organ_data$"Middle vs. Young", organ_data$"Old vs. Middle")
    ) %>% filter(!is.na(Ratio_FC))
    
    if (nrow(plot_data) > 0) {
      p <- ggplot(plot_data, aes(x = Comparison, y = Ratio_FC, group = Gene_Name)) +
        geom_line() + geom_point() +
        geom_text(data = plot_data[plot_data$Comparison == "Old vs. Young", ],
                  aes(label = Gene_Name), hjust = -0.1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste("Proteins with Mixed Patterns in", organ))
      
      ggsave(paste0("Mixed_Patterns_", gsub(" ", "_", organ), ".pdf"), p, width = 10, height = 6)
    }
  }
}

# ---- 8. Organ-Specific GO Analysis ----

# Run enrichGO for increased and decreased proteins per organ
organ_go_results <- list()

for (organ in organs) {
  cat("\nAnalyzing GO terms for", organ, "...\n")
  
  organ_increase <- proteins_increase_only[proteins_increase_only$Organ == organ, ]
  organ_decrease <- proteins_decrease_only[proteins_decrease_only$Organ == organ, ]
  
  cat("  Increased:", nrow(organ_increase), " | Decreased:", nrow(organ_decrease), "\n")
  
  go_inc <- if (nrow(organ_increase) >= 2) {
    tryCatch(enrichGO(gene = unique(organ_increase$Gene_Name), OrgDb = org.Mm.eg.db,
                      keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                      pvalueCutoff = 0.05, readable = TRUE),
             error = function(e) { cat("  Error (increase):", e$message, "\n"); NULL })
  } else NULL
  
  go_dec <- if (nrow(organ_decrease) >= 2) {
    tryCatch(enrichGO(gene = unique(organ_decrease$Gene_Name), OrgDb = org.Mm.eg.db,
                      keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
                      pvalueCutoff = 0.05, readable = TRUE),
             error = function(e) { cat("  Error (decrease):", e$message, "\n"); NULL })
  } else NULL
  
  organ_go_results[[organ]] <- list(increase = go_inc, decrease = go_dec)
  
  if (!is.null(go_inc)) cat("  Found", nrow(go_inc@result), "terms for increased proteins\n")
  if (!is.null(go_dec)) cat("  Found", nrow(go_dec@result), "terms for decreased proteins\n")
}

# ---- 9. GO Term Heatmap Across Organs ----

all_go_terms <- unique(unlist(lapply(organ_go_results, function(x) {
  c(if (!is.null(x$increase)) x$increase@result$Description else NULL,
    if (!is.null(x$decrease)) x$decrease@result$Description else NULL)
})))

heatmap_data <- matrix(0, nrow = length(all_go_terms), ncol = length(organs),
                       dimnames = list(all_go_terms, organs))

for (organ in organs) {
  if (!is.null(organ_go_results[[organ]]$increase)) {
    inc_res <- organ_go_results[[organ]]$increase@result
    idx <- match(inc_res$Description, all_go_terms)
    heatmap_data[idx, organ] <- -log10(inc_res$p.adjust)
  }
  if (!is.null(organ_go_results[[organ]]$decrease)) {
    dec_res <- organ_go_results[[organ]]$decrease@result
    idx <- match(dec_res$Description, all_go_terms)
    # Negative sign indicates decreased insolubility direction
    heatmap_data[idx, organ] <- -(-log10(dec_res$p.adjust))
  }
}

nonzero_rows <- which(rowSums(abs(heatmap_data)) > 0)
if (length(nonzero_rows) > 0) {
  heatmap_data <- heatmap_data[nonzero_rows, ]
  pdf("GO_Terms_Across_Organs_Heatmap.pdf", width = 12, height = 20)
  print(pheatmap(heatmap_data, cluster_rows = TRUE, cluster_cols = TRUE,
                 fontsize_row = 8, main = "GO Term Enrichment Across Organs",
                 color = colorRampPalette(c("blue", "white", "red"))(100),
                 breaks = seq(-max(abs(heatmap_data)), max(abs(heatmap_data)), length.out = 101)))
  dev.off()
}

write.csv(heatmap_data, "GO_Terms_by_Organ_with_Direction.csv")

# ---- 10. Per-Organ GO Overlap: Venn Diagrams and Bar Plots ----

for (organ in organs) {
  inc_res <- organ_go_results[[organ]]$increase
  dec_res <- organ_go_results[[organ]]$decrease
  
  if (is.null(inc_res) || is.null(dec_res)) next
  
  inc_terms <- inc_res@result$Description
  dec_terms <- dec_res@result$Description
  overlap   <- intersect(inc_terms, dec_terms)
  
  # Venn diagram
  pdf(paste0("GO_Term_Venn_", gsub(" ", "_", organ), ".pdf"), width = 10, height = 10)
  venn_p <- draw.pairwise.venn(
    area1      = length(inc_terms), area2 = length(dec_terms),
    cross.area = length(overlap),
    category   = c("Increased", "Decreased"),
    fill = c("red", "blue"), alpha = 0.5,
    cat.col = c("red", "blue"), cat.cex = 1.5, cat.dist = 0.07,
    euler.d = TRUE, scaled = TRUE
  )
  grid.draw(venn_p)
  dev.off()
  
  # Bar plot of overlapping terms
  if (length(overlap) > 0) {
    overlap_df <- data.frame(
      GO_Term       = overlap,
      Increase_padj = inc_res@result$p.adjust[match(overlap, inc_terms)],
      Decrease_padj = dec_res@result$p.adjust[match(overlap, dec_terms)],
      Increase_Count = inc_res@result$Count[match(overlap, inc_terms)],
      Decrease_Count = dec_res@result$Count[match(overlap, dec_terms)]
    ) %>% arrange(Increase_padj + Decrease_padj)
    
    write.csv(overlap_df, paste0("GO_Term_Overlap_Details_", gsub(" ", "_", organ), ".csv"), row.names = FALSE)
    
    top_terms <- head(overlap_df, 10)
    plot_data <- data.frame(
      GO_Term   = rep(top_terms$GO_Term, 2),
      Direction = rep(c("Increased", "Decreased"), each = nrow(top_terms)),
      Significance = c(-log10(top_terms$Increase_padj), -log10(top_terms$Decrease_padj))
    )
    
    p <- ggplot(plot_data, aes(x = reorder(GO_Term, Significance), y = Significance, fill = Direction)) +
      geom_bar(stat = "identity", position = "dodge") + coord_flip() +
      scale_fill_manual(values = c("Increased" = "red", "Decreased" = "blue")) +
      theme_minimal() +
      labs(title = paste("Top Overlapping GO Terms in", organ),
           x = "GO Term", y = "-log10(adjusted p-value)") +
      theme(axis.text.y = element_text(size = 8), legend.position = "bottom")
    
    ggsave(paste0("GO_Term_Overlap_Barplot_", gsub(" ", "_", organ), ".pdf"), p, width = 12, height = 8)
  }
}

# ---- 11. Per-Comparison GO UpSet Plots ----

# Runs GO enrichment for each organ × age comparison, then creates UpSet plots
# showing which proteins are shared across organs within the same GO term

perform_GO_analysis <- function(data, age_comp) {
  go_results <- data.frame()
  
  for (organ in unique(data$Organ)) {
    organ_data <- data[data$Organ == organ & !is.na(data[[age_comp]]), ]
    increased  <- unique(organ_data$Gene_Name[organ_data[[age_comp]] > 0])
    decreased  <- unique(organ_data$Gene_Name[organ_data[[age_comp]] < 0])
    
    for (direction in c("Increased", "Decreased")) {
      genes <- if (direction == "Increased") increased else decreased
      
      if (length(genes) >= 2) {
        go_res <- tryCatch(
          clusterProfiler::enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                                    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.2, readable = TRUE),
          error = function(e) NULL
        )
        
        if (!is.null(go_res) && nrow(go_res@result) > 0) {
          go_results <- rbind(go_results, data.frame(
            GO_Term   = go_res@result$Description,
            GO_ID     = go_res@result$ID,
            Organ     = organ,
            Direction = direction,
            Proteins  = go_res@result$geneID,
            Pvalue    = go_res@result$pvalue,
            Padj      = go_res@result$p.adjust,
            GeneRatio = go_res@result$GeneRatio
          ))
        }
      }
    }
  }
  
  label <- gsub(" vs. ", "_vs_", age_comp)
  saveRDS(go_results,  paste0("GO_analysis_results_", label, ".rds"))
  write.csv(go_results, paste0("GO_analysis_results_", label, ".csv"), row.names = FALSE)
  return(go_results)
}

create_upset_plots <- function(go_results, age_comp) {
  filtered_terms  <- c()
  shared_proteins <- list()
  
  for (go_term in unique(go_results$GO_Term)) {
    term_data  <- go_results[go_results$GO_Term == go_term, ]
    organ_sets <- list()
    
    for (org in unique(term_data$Organ)) {
      org_data <- term_data[term_data$Organ == org, ]
      if (any(org_data$Direction == "Increased"))
        organ_sets[[paste0(org, "_increased")]] <- unlist(strsplit(
          org_data$Proteins[org_data$Direction == "Increased"], "/"))
      if (any(org_data$Direction == "Decreased"))
        organ_sets[[paste0(org, "_decreased")]] <- unlist(strsplit(
          org_data$Proteins[org_data$Direction == "Decreased"], "/"))
    }
    
    all_proteins <- unlist(organ_sets)
    shared <- names(table(all_proteins)[table(all_proteins) > 1])
    
    if (length(shared) > 0) {
      filtered_terms <- c(filtered_terms, go_term)
      shared_proteins[[go_term]] <- data.frame(
        GO_Term = go_term,
        Protein = shared,
        Organs  = sapply(shared, function(p) {
          paste(names(organ_sets)[sapply(organ_sets, function(x) p %in% x)], collapse = "|")
        })
      )
    }
  }
  
  shared_df <- do.call(rbind, shared_proteins)
  write.csv(shared_df, paste0("Shared_Proteins_", gsub(" vs. ", "_vs_", age_comp), ".csv"), row.names = FALSE)
  
  temp_dir   <- tempdir()
  temp_files <- c()
  
  for (i in seq_along(filtered_terms)) {
    pdf_name   <- file.path(temp_dir, paste0("plot_", i, ".pdf"))
    temp_files <- c(temp_files, pdf_name)
    pdf(pdf_name, width = 8, height = 6)
    
    go_term   <- filtered_terms[i]
    term_data <- go_results[go_results$GO_Term == go_term, ]
    organ_sets <- list()
    
    for (org in unique(term_data$Organ)) {
      org_data <- term_data[term_data$Organ == org, ]
      if (any(org_data$Direction == "Increased"))
        organ_sets[[paste0(org, "_increased")]] <- unlist(strsplit(
          org_data$Proteins[org_data$Direction == "Increased"], "/"))
      if (any(org_data$Direction == "Decreased"))
        organ_sets[[paste0(org, "_decreased")]] <- unlist(strsplit(
          org_data$Proteins[org_data$Direction == "Decreased"], "/"))
    }
    
    m <- make_comb_mat(organ_sets)
    set_colors <- sapply(names(organ_sets), function(set_name) {
      organ <- strsplit(set_name, "_")[[1]][1]
      organ_colors[organ]
    })
    
    ht <- UpSet(m,
                column_title      = substr(go_term, 1, 50),
                set_order         = names(organ_sets),
                comb_order        = order(comb_size(m), decreasing = TRUE),
                top_annotation    = upset_top_annotation(m, annotation_name_side = "left", annotation_name_rot = 90),
                right_annotation  = upset_right_annotation(m),
                row_names_gp      = gpar(col = set_colors),
                height            = unit(6, "cm"),
                width             = unit(9, "cm"))
    draw(ht)
    dev.off()
  }
  
  pdf_combine(temp_files,
              output = paste0("GO_UpSet_Shared_", gsub(" vs. ", "_vs_", age_comp), ".pdf"))
  unlink(temp_files)
}

# Run for all three age comparisons
cat("\nRunning per-comparison GO UpSet analysis...\n")

go_results_OY <- perform_GO_analysis(df_wide_All, "Old vs. Young")
create_upset_plots(go_results_OY, "Old vs. Young")

go_results_MY <- perform_GO_analysis(df_wide_All, "Middle vs. Young")
create_upset_plots(go_results_MY, "Middle vs. Young")

go_results_OM <- perform_GO_analysis(df_wide_All, "Old vs. Middle")
create_upset_plots(go_results_OM, "Old vs. Middle")

# ---- 12. Summary Text File ----

sink("Organ_Specific_GO_Analysis_Summary.txt")
cat("ORGAN-SPECIFIC GO ANALYSIS SUMMARY\n================================\n\n")
for (organ in organs) {
  cat("\n", organ, "\n")
  cat(paste(rep("-", nchar(organ)), collapse = ""), "\n")
  if (!is.null(organ_go_results[[organ]]$increase)) {
    cat("\nTop GO terms for increased insolubility:\n")
    print(head(organ_go_results[[organ]]$increase@result[, c("Description", "Count", "p.adjust")], 5))
  }
  if (!is.null(organ_go_results[[organ]]$decrease)) {
    cat("\nTop GO terms for decreased insolubility:\n")
    print(head(organ_go_results[[organ]]$decrease@result[, c("Description", "Count", "p.adjust")], 5))
  }
}
cat("\n\nFinal Statistics:\n----------------")
cat("\nTotal GO terms analyzed:", nrow(heatmap_data))
cat("\n\nGO terms per organ:\n")
print(colSums(abs(heatmap_data) > 0))
sink()