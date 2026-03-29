# ============================================================
# Project: TMT Proteomics — Multi-organ Mouse Aging
# Script:  mouseTMT_insolubility_subcellular.R
# Purpose: Annotates significantly regulated proteins from the
#          published TMT dataset with GO Cellular Component
#          (subcellular localization) data. Annotations sourced
#          from org.Mm.eg.db (primary), then UniProt and Ensembl
#          as fallbacks for unannotated genes. Classifies
#          proteins by direction of insolubility change
#          (increased vs. decreased) and tests for compartment
#          enrichment using chi-square tests.
# Input:   Data S2_TMT significanly-regulated proteins.xlsx
#          (multi-sheet: one sheet per organ/condition)
# Output:  Data S2_TMT significantly-regulated proteins_
#            with_subcellular_locations.xlsx
#          Data S2_TMT proteins_with_ALL_annotations.xlsx
#          Data S2_TMT_proteins_BY_DIRECTION_with_locations.xlsx
#          subcellular_localization_summary.csv
#          subcellular_localization_BY_DIRECTION.csv
#          statistical_tests_increased_vs_decreased.csv
#          genes_with_alternative_annotations.csv
# Author:  Chia-Lung Chuang
# Updated: 2026-03
# Dependencies: BiocManager packages — install once before use:
#   BiocManager::install(c("org.Mm.eg.db", "GO.db", "biomaRt"))
# ============================================================

library(org.Mm.eg.db)
library(GO.db)
library(readxl)
library(writexl)
library(tidyverse)
library(dplyr)
library(biomaRt)

# NOTE: Set your working directory to the folder containing the input Excel file
# setwd("path/to/your/data")

# ============================================================
# 1. Load Input Data
# ============================================================

file_path  <- "Data S2_TMT significanly-regulated proteins.xlsx"
sheet_names <- excel_sheets(file_path)
cat("Found sheets:", paste(sheet_names, collapse = ", "), "\n")

all_sheets <- setNames(
  lapply(sheet_names, function(s) { cat("Reading sheet:", s, "\n"); read_excel(file_path, sheet = s) }),
  sheet_names
)

# ============================================================
# 2. GO Cellular Component Annotation
# ============================================================

# Retrieves GO CC terms for a vector of mouse gene symbols
# using org.Mm.eg.db → GO.db
get_go_cellular_components <- function(gene_symbols) {
  cat("Getting GO CC annotations for", length(gene_symbols), "genes...\n")
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
  
  # Gene symbol → Entrez ID
  entrez_ids <- mapIds(org.Mm.eg.db, keys = gene_symbols,
                       column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  cat("Mapped", sum(!is.na(entrez_ids)), "genes to Entrez IDs\n")
  
  # Entrez ID → GO CC terms
  go_terms <- AnnotationDbi::select(org.Mm.eg.db,
                                    keys    = entrez_ids[!is.na(entrez_ids)],
                                    columns = c("SYMBOL", "GO", "ONTOLOGY"),
                                    keytype = "ENTREZID") %>%
    filter(ONTOLOGY == "CC")
  
  cat("Found", nrow(go_terms), "GO CC annotations\n")
  
  # GO ID → human-readable term name
  unique_go_ids  <- unique(go_terms$GO)
  go_id_to_term  <- AnnotationDbi::select(GO.db, keys = unique_go_ids,
                                          columns = c("GOID", "TERM"), keytype = "GOID")
  
  # Summarize: one row per gene with all CC terms concatenated
  go_terms %>%
    left_join(go_id_to_term, by = c("GO" = "GOID")) %>%
    group_by(SYMBOL) %>%
    summarize(All_Cellular_Components = paste(unique(TERM), collapse = "; "),
              GO_CC_IDs               = paste(unique(GO),   collapse = "; "),
              .groups = "drop")
}

# ============================================================
# 3. Binary Compartment Columns
# ============================================================

# Adds one binary column per major compartment based on
# pattern matching against the All_Cellular_Components string
create_compartment_columns <- function(df, location_column = "All_Cellular_Components") {
  compartment_patterns <- list(
    Nucleus               = c("nucleus", "nuclear", "chromosome", "nucleolus", "nucleoplasm"),
    Cytoplasm             = c("cytoplasm", "cytosol", "cytoplasmic"),
    Mitochondrion         = c("mitochondri"),
    Endoplasmic_Reticulum = c("endoplasmic reticulum", "\\ber\\b"),
    Golgi                 = c("golgi"),
    Plasma_Membrane       = c("plasma membrane", "cell membrane"),
    Membrane              = c("membrane"),
    Ribosome              = c("ribosom", "ribosomal"),
    Lysosome              = c("lysosome", "lysosomal"),
    Peroxisome            = c("peroxisome"),
    Cytoskeleton          = c("cytoskeleton", "actin", "tubulin", "microtubule"),
    Extracellular         = c("extracellular", "secreted", "cell surface")
  )
  
  for (comp_name in names(compartment_patterns)) {
    patterns    <- compartment_patterns[[comp_name]]
    loc_strings <- tolower(df[[location_column]])
    df[[comp_name]] <- as.integer(
      !is.na(loc_strings) & sapply(loc_strings, function(loc)
        any(sapply(patterns, function(p) grepl(p, loc))))
    )
  }
  return(df)
}

compartment_cols <- c("Nucleus", "Cytoplasm", "Mitochondrion", "Endoplasmic_Reticulum",
                      "Golgi", "Plasma_Membrane", "Ribosome", "Lysosome",
                      "Peroxisome", "Cytoskeleton", "Extracellular")

# ============================================================
# 4. Annotate All Genes (run once across all sheets)
# ============================================================

all_genes <- unique(unlist(lapply(all_sheets, function(s) s$Gene_Name)))
all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
cat("\nTotal unique genes across all sheets:", length(all_genes), "\n")

go_annotations <- get_go_cellular_components(all_genes)
go_annotations  <- create_compartment_columns(go_annotations)

cat("Genes with GO annotations:", nrow(go_annotations), "\n")
cat("Genes without annotations:", length(all_genes) - nrow(go_annotations), "\n")

# ============================================================
# 5. Merge Annotations into Each Sheet
# ============================================================

annotated_sheets <- lapply(sheet_names, function(sheet) {
  cat("\nProcessing sheet:", sheet, "\n")
  sheet_data <- all_sheets[[sheet]]
  
  sheet_annotated <- sheet_data %>%
    left_join(go_annotations, by = c("Gene_Name" = "SYMBOL"))
  
  annotation_cols <- setdiff(colnames(go_annotations), "SYMBOL")
  sheet_annotated <- sheet_annotated[, c(colnames(sheet_data), annotation_cols)]
  
  cat("  Annotated:", sum(!is.na(sheet_annotated$All_Cellular_Components)),
      "/", nrow(sheet_data), "\n")
  sheet_annotated
})
names(annotated_sheets) <- sheet_names

# Save primary annotated file
output_file <- "Data S2_TMT significantly-regulated proteins_with_subcellular_locations.xlsx"
write_xlsx(annotated_sheets, output_file)
cat("\nSaved:", output_file, "\n")

# ============================================================
# 6. Fallback Annotations: UniProt and Ensembl
# ============================================================

# Collect genes missing GO CC annotations
missing_info <- local({
  missing_list <- lapply(sheet_names, function(s) {
    annotated_sheets[[s]] %>%
      filter(is.na(All_Cellular_Components)) %>%
      select(Gene_Name, GN_ID) %>%
      distinct() %>%
      mutate(Organ = s)
  })
  all_missing    <- bind_rows(missing_list)
  unique_missing <- unique(all_missing$Gene_Name)
  cat("\nGenes without GO annotations:", length(unique_missing), "\n")
  list(all_missing = all_missing, unique_missing = unique_missing)
})

# UniProt fallback
get_uniprot_locations <- function(gene_symbols) {
  cat("\n=== Querying UniProt for", length(gene_symbols), "genes ===\n")
  tryCatch({
    up      <- UniProt.ws::UniProt.ws(taxId = 10090)
    results <- AnnotationDbi::select(up, keys = gene_symbols,
                                     columns = c("UNIPROTKB", "GENES", "SUBCELLULAR-LOCATIONS"),
                                     keytype = "Gene_Name")
    results %>%
      filter(!is.na(`SUBCELLULAR-LOCATIONS`)) %>%
      select(GENES, `SUBCELLULAR-LOCATIONS`) %>%
      distinct() %>%
      rename(Gene_Name = GENES, UniProt_Location = `SUBCELLULAR-LOCATIONS`)
  }, error = function(e) {
    cat("UniProt query failed:", e$message, "\n")
    data.frame(Gene_Name = character(), UniProt_Location = character())
  })
}

# Ensembl fallback
get_ensembl_locations <- function(gene_symbols) {
  cat("\n=== Querying Ensembl for", length(gene_symbols), "genes ===\n")
  tryCatch({
    ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
    results  <- getBM(attributes = c("external_gene_name", "go_id", "name_1006", "namespace_1003"),
                      filters    = "external_gene_name",
                      values     = gene_symbols, mart = ensembl) %>%
      filter(namespace_1003 == "cellular_component") %>%
      group_by(external_gene_name) %>%
      summarize(Ensembl_Location = paste(unique(name_1006), collapse = "; "),
                Ensembl_GO       = paste(unique(go_id),     collapse = "; "),
                .groups = "drop") %>%
      rename(Gene_Name = external_gene_name)
    cat("Ensembl found annotations for", nrow(results), "genes\n")
    results
  }, error = function(e) {
    cat("Ensembl query failed:", e$message, "\n")
    data.frame(Gene_Name = character(), Ensembl_Location = character(), Ensembl_GO = character())
  })
}

# Gene type check (identifies non-coding RNAs that won't have CC terms)
check_gene_type <- function(gene_symbols) {
  tryCatch({
    AnnotationDbi::select(org.Mm.eg.db, keys = gene_symbols,
                          columns = c("SYMBOL", "GENETYPE"), keytype = "SYMBOL")
  }, error = function(e) data.frame(SYMBOL = gene_symbols))
}

uniprot_results <- get_uniprot_locations(missing_info$unique_missing)
ensembl_results <- get_ensembl_locations(missing_info$unique_missing)
gene_types      <- check_gene_type(missing_info$unique_missing)

# Combine fallback sources into one unified table
alternative_annotations <- missing_info$all_missing %>%
  left_join(uniprot_results, by = "Gene_Name") %>%
  left_join(ensembl_results, by = "Gene_Name") %>%
  left_join(gene_types,      by = c("Gene_Name" = "SYMBOL")) %>%
  mutate(
    Alternative_Location = case_when(
      !is.na(UniProt_Location) ~ UniProt_Location,
      !is.na(Ensembl_Location) ~ Ensembl_Location,
      TRUE ~ NA_character_
    ),
    Annotation_Source = case_when(
      !is.na(UniProt_Location) ~ "UniProt",
      !is.na(Ensembl_Location) ~ "Ensembl",
      TRUE ~ "Not Found"
    )
  )

write.csv(alternative_annotations, "genes_with_alternative_annotations.csv", row.names = FALSE)
cat("Saved: genes_with_alternative_annotations.csv\n")

# Update annotated_sheets with fallback annotations where GO was missing
final_annotated_sheets <- lapply(sheet_names, function(sheet) {
  organ_alt <- alternative_annotations %>%
    filter(Organ == sheet) %>%
    select(Gene_Name, Alternative_Location, Annotation_Source)
  
  updated <- annotated_sheets[[sheet]] %>%
    left_join(organ_alt, by = "Gene_Name") %>%
    mutate(
      All_Cellular_Components = ifelse(
        is.na(All_Cellular_Components) & !is.na(Alternative_Location),
        Alternative_Location, All_Cellular_Components),
      Annotation_Source = ifelse(is.na(Annotation_Source), "Gene Ontology", Annotation_Source)
    ) %>%
    select(-Alternative_Location)
  
  improved <- sum(is.na(annotated_sheets[[sheet]]$All_Cellular_Components)) -
    sum(is.na(updated$All_Cellular_Components))
  if (improved > 0) cat(sheet, ": added", improved, "alternative annotations\n")
  updated
})
names(final_annotated_sheets) <- sheet_names

write_xlsx(final_annotated_sheets, "Data S2_TMT proteins_with_ALL_annotations.xlsx")
cat("Saved: Data S2_TMT proteins_with_ALL_annotations.xlsx\n")

# ============================================================
# 7. Subcellular Localization Summary per Organ
# ============================================================

summary_table <- bind_rows(lapply(sheet_names, function(sheet) {
  data <- final_annotated_sheets[[sheet]]
  row  <- data.frame(Organ = sheet, Total_Proteins = nrow(data))
  for (comp in compartment_cols) {
    if (comp %in% colnames(data)) {
      count <- sum(data[[comp]] == 1, na.rm = TRUE)
      row[[paste0(comp, "_Count")]]   <- count
      row[[paste0(comp, "_Percent")]] <- round(count / nrow(data) * 100, 1)
    }
  }
  row
}))

write.csv(summary_table, "subcellular_localization_summary.csv", row.names = FALSE)
cat("Saved: subcellular_localization_summary.csv\n")

# ============================================================
# 8. Direction Classification: Increased vs. Decreased Insolubility
# ============================================================

# Classify proteins by sign of log2FC(insoluble/soluble ratio)
classify_by_direction <- function(sheets) {
  lapply(sheets, function(sheet_data) {
    sheet_data %>%
      mutate(Direction = case_when(
        `log2FC(insoluble/soluble protein ratio)` > 0 ~ "Increased_Insolubility",
        `log2FC(insoluble/soluble protein ratio)` < 0 ~ "Decreased_Insolubility",
        TRUE ~ "No_Change"
      ))
  })
}

classified_sheets <- classify_by_direction(final_annotated_sheets)

# Per-direction compartment summaries
direction_summaries <- local({
  summarize_direction <- function(data, direction_label, direction_value) {
    sub <- data %>% filter(Direction == direction_value)
    row <- data.frame(Direction = direction_label, Total_Proteins = nrow(sub))
    for (comp in compartment_cols) {
      if (comp %in% colnames(sub)) {
        count <- sum(sub[[comp]] == 1, na.rm = TRUE)
        row[[paste0(comp, "_Count")]]   <- count
        row[[paste0(comp, "_Percent")]] <- round(count / nrow(sub) * 100, 1)
      }
    }
    row
  }
  
  bind_rows(lapply(sheet_names, function(sheet) {
    data <- classified_sheets[[sheet]]
    bind_rows(
      summarize_direction(data, "Increased", "Increased_Insolubility") %>% mutate(Organ = sheet),
      summarize_direction(data, "Decreased", "Decreased_Insolubility") %>% mutate(Organ = sheet)
    )
  }))
})

write.csv(direction_summaries %>% filter(Direction == "Increased"),
          "subcellular_localization_INCREASED_insolubility.csv", row.names = FALSE)
write.csv(direction_summaries %>% filter(Direction == "Decreased"),
          "subcellular_localization_DECREASED_insolubility.csv", row.names = FALSE)
write.csv(direction_summaries, "subcellular_localization_BY_DIRECTION.csv", row.names = FALSE)

# Print detailed comparison table per organ
cat("\n=== INCREASED vs. DECREASED INSOLUBILITY: COMPARTMENT DISTRIBUTION ===\n")
for (sheet in sheet_names) {
  data           <- classified_sheets[[sheet]]
  increased_data <- data %>% filter(Direction == "Increased_Insolubility")
  decreased_data <- data %>% filter(Direction == "Decreased_Insolubility")
  
  cat("\n", strrep("=", 60), "\n", toupper(sheet), "\n", strrep("=", 60), "\n", sep = "")
  cat(sprintf("%-25s %15s %15s %12s\n", "Compartment", "Increased", "Decreased", "Difference"))
  cat(strrep("-", 70), "\n")
  
  for (comp in compartment_cols) {
    if (comp %in% colnames(data)) {
      inc_n   <- sum(increased_data[[comp]] == 1, na.rm = TRUE)
      inc_pct <- round(inc_n / max(nrow(increased_data), 1) * 100, 1)
      dec_n   <- sum(decreased_data[[comp]] == 1, na.rm = TRUE)
      dec_pct <- round(dec_n / max(nrow(decreased_data), 1) * 100, 1)
      cat(sprintf("%-25s %6d (%5.1f%%) %6d (%5.1f%%) %+8.1f%%\n",
                  comp, inc_n, inc_pct, dec_n, dec_pct, inc_pct - dec_pct))
    }
  }
}

# ============================================================
# 9. Direction-Split Excel Output
# ============================================================

direction_sheets <- unlist(lapply(sheet_names, function(sheet) {
  data <- classified_sheets[[sheet]]
  list(
    setNames(list(data %>% filter(Direction == "Increased_Insolubility") %>%
                    arrange(desc(`log2FC(insoluble/soluble protein ratio)`))),
             paste0(sheet, "_Increased")),
    setNames(list(data %>% filter(Direction == "Decreased_Insolubility") %>%
                    arrange(`log2FC(insoluble/soluble protein ratio)`)),
             paste0(sheet, "_Decreased"))
  )
}), recursive = FALSE)

write_xlsx(direction_sheets, "Data S2_TMT_proteins_BY_DIRECTION_with_locations.xlsx")
cat("Saved: Data S2_TMT_proteins_BY_DIRECTION_with_locations.xlsx\n")

# ============================================================
# 10. Chi-Square Tests: Compartment Enrichment by Direction
# ============================================================

chi_square_results <- bind_rows(lapply(sheet_names, function(sheet) {
  data           <- classified_sheets[[sheet]]
  increased_data <- data %>% filter(Direction == "Increased_Insolubility")
  decreased_data <- data %>% filter(Direction == "Decreased_Insolubility")
  
  bind_rows(lapply(compartment_cols, function(comp) {
    if (!comp %in% colnames(data)) return(NULL)
    
    inc_yes <- sum(increased_data[[comp]] == 1, na.rm = TRUE)
    dec_yes <- sum(decreased_data[[comp]] == 1, na.rm = TRUE)
    contingency <- matrix(
      c(inc_yes, nrow(increased_data) - inc_yes,
        dec_yes, nrow(decreased_data) - dec_yes),
      nrow = 2, byrow = TRUE)
    
    if (!all(contingency >= 5)) return(NULL)
    
    test <- chisq.test(contingency)
    data.frame(Organ = sheet, Compartment = comp,
               Chi_Square = test$statistic, P_Value = test$p.value,
               Significant = test$p.value < 0.05)
  }))
}))

write.csv(chi_square_results, "statistical_tests_increased_vs_decreased.csv", row.names = FALSE)
cat("Saved: statistical_tests_increased_vs_decreased.csv\n")

# Print significant results
sig_results <- chi_square_results %>% filter(Significant)
if (nrow(sig_results) > 0) {
  cat("\n=== SIGNIFICANT COMPARTMENT ENRICHMENTS (p < 0.05) ===\n")
  print(sig_results %>% select(Organ, Compartment, P_Value) %>% arrange(P_Value))
} else {
  cat("\nNo significant compartment enrichments found.\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files created:\n")
cat("  1.", output_file, "\n")
cat("  2. Data S2_TMT proteins_with_ALL_annotations.xlsx\n")
cat("  3. Data S2_TMT_proteins_BY_DIRECTION_with_locations.xlsx\n")
cat("  4. subcellular_localization_summary.csv\n")
cat("  5. subcellular_localization_BY_DIRECTION.csv\n")
cat("  6. genes_with_alternative_annotations.csv\n")
cat("  7. statistical_tests_increased_vs_decreased.csv\n")