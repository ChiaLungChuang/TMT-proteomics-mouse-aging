# TMT Proteomics — Multi-organ Mouse Aging

Reproducible R-based analytical workflow for quantitative TMT proteomics across 8 mouse tissues, examining age-related changes in protein abundance and solubility from Young (6M), Middle (24M), and Old (30M) mice.

Published in: *Developmental Cell*, February 19, 2026.
https://doi.org/10.1016/j.devcel.2026.01.015

---

## Project Overview

This repository contains the complete analysis pipeline for a systems-level proteomics study of mouse aging. The workflow characterizes how protein **abundance** and **insolubility** (measured as the Insoluble/Soluble ratio) change with age across 8 tissues simultaneously, enabling both tissue-specific and cross-organ perspectives on protein aggregation during aging.

**Key analytical questions addressed:**
- Which proteins become more or less insoluble with age, in which tissues?
- Are certain proteins dysregulated across multiple organs (inter-organ proteins)?
- Do protein biophysical properties (thermal stability, half-life) predict age-related insolubility changes?
- Which subcellular compartments are enriched for age-related protein aggregation?

---

## Repository Structure

```
TMT-proteomics-mouse-aging/
│
├── README.md
│
├── 01_per-tissue-analysis/          # Independent analysis for each of 8 tissues
│   ├── mouseTMT_liver.R
│   ├── mouseTMT_heart.R
│   ├── mouseTMT_cerebellum.R
│   ├── mouseTMT_hippocampus.R
│   ├── mouseTMT_neocortex.R
│   ├── mouseTMT_soleus.R
│   ├── mouseTMT_TA.R
│   └── mouseTMT_WhiteAdiposeTissue.R
│
├── 02_cross-organ-integration/      # Systems-level analysis combining all tissues
│   ├── Organs.R
│   ├── Organs_mRNA.R
│   ├── Organs_abundance.R
│   └── Organs_abundance_upset_plot.R
│
└── 03_insolubility-sub-analysis/    # Deep-dive analyses of the insolubility dataset
    ├── Insolubility_meltTemp.R
    ├── insolubility_protein_half-life.R
    ├── Insolubility_GOterm_upset.R
    ├── mouseTMT_insolubility_subcellular.R
    └── mouseTMT_insolubility_with_function.R
```

---

## Analytical Workflow

### Stage 1 — Per-Tissue Analysis (`01_per-tissue-analysis/`)

Each tissue script is structurally identical and runs independently. The 8 tissues analyzed are: **Liver, Heart, Cerebellum, Hippocampus, Neocortex, Soleus, Tibialis Anterior, White Adipose Tissue**.

**Steps per tissue:**
1. Import raw TMT intensities from Excel; clean protein accession numbers
2. Reshape from wide to long format; parse Age (Young/Middle/Old), Replicate, and Solubility (Soluble/Insoluble) from TMT channel column names
3. Calculate per-protein **Insolubility Ratio** (Insoluble / Soluble) and total **Abundance** (Soluble + Insoluble)
4. One-way ANOVA across age groups followed by Tukey's HSD post-hoc test for all pairwise comparisons (Old vs. Young, Middle vs. Young, Old vs. Middle)
5. Log2 fold change calculation for each pairwise comparison
6. Scatter plots: Log2 Abundance FC (x) vs. Log2 Insolubility Ratio FC (y), with point size scaled to −log10(Tukey adjusted p-value); multiple labeling and faceting variants exported as PDF

**Statistical approach:** One-way ANOVA with Tukey's Honest Significant Difference correction. Point size in scatter plots encodes −log10(Tukey adjusted p-value), with a minimum threshold of 1.3 (p < 0.05).

**Key outputs per tissue:** `mouseTMT_<tissue>.csv`, `ratio_data_<tissue>_OY/MY/OM.csv`, `combined_plot_data_mouseTMT_<tissue>.csv`, scatter plot PDFs

---

### Stage 2 — Cross-Organ Integration (`02_cross-organ-integration/`)

**`Organs.R`** — Combines all 8 per-tissue outputs. Classifies proteins as intra-organ (significant in only one tissue) or inter-organ (significant in ≥2 tissues). Generates whole-proteome scatter plots and organ-gene network graphs using Fruchterman-Reingold and Kamada-Kawai layouts (ggraph/tidygraph).

**`Organs_mRNA.R`** — Builds the same network graphs with gene nodes colored by mRNA expression level (average TPM per tissue from a companion transcriptomics dataset), highlighting proteins whose age-related insolubility changes cannot be explained by transcriptional regulation alone.

**`Organs_abundance.R`** — Network graphs with gene nodes colored by baseline protein abundance (log-transformed Abundance_mean_young), providing context on which highly expressed proteins are most prone to age-related aggregation.

**`Organs_abundance_upset_plot.R`** — UpSet plots showing the overlap of significant proteins across the 8 organs for each age comparison, with organs colored by the study's unified color palette.

---

### Stage 3 — Insolubility Sub-Analyses (`03_insolubility-sub-analysis/`)

**`Insolubility_meltTemp.R`** — Correlates protein thermal stability (melting temperature from a separate thermal proteome profiling dataset) with age-related insolubility fold change across all 8 tissues. Tests whether thermally unstable proteins are more prone to age-related aggregation.

**`insolubility_protein_half-life.R`** — Correlates tissue-specific protein half-life with insolubility fold change in 4 tissues (Heart, Cerebellum, Neocortex, Liver) where matched half-life data are available.

**`Insolubility_GOterm_upset.R`** — GO biological process enrichment analysis (clusterProfiler, org.Mm.eg.db) comparing proteins with **increased** vs. **decreased** insolubility across organs and age comparisons. Generates UpSet plots highlighting shared proteins across GO terms and organs, per-organ Venn diagrams, and a heatmap of −log10(p.adjust) across all organs.

**`mouseTMT_insolubility_subcellular.R`** — Annotates proteins with GO Cellular Component (subcellular localization) terms using org.Mm.eg.db as primary source, with UniProt and Ensembl as fallbacks. Classifies proteins by direction of insolubility change and applies chi-square tests to identify compartments enriched for age-related protein aggregation.

**`mouseTMT_insolubility_with_function.R`** — A function-based refactor of the per-tissue pipeline that additionally computes **Soluble FC vs. Insoluble FC** scatter plots (i.e., plots the two TMT fractions separately rather than as a ratio). Includes Pearson correlation analyses between soluble and insoluble fold changes per tissue and age comparison.

---

## Requirements

**R version:** ≥ 4.5.1

**CRAN packages:**
```r
install.packages(c(
  "tidyverse", "ggplot2", "ggrepel", "readxl", "writexl",
  "broom", "scales", "ggpmisc", "ggVennDiagram",
  "VennDiagram", "gridExtra", "RColorBrewer", "pheatmap",
  "UpSetR", "pdftools", "circlize", "igraph", "ggraph",
  "tidygraph", "ggnewscale", "extrafont"
))
```

**Bioconductor packages** (install once):
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "org.Mm.eg.db",
  "GO.db",
  "ComplexHeatmap",
  "biomaRt"
))
```

---

## Input Data

Raw TMT intensity data are not included in this repository as the underlying datasets are associated with the published manuscript. Scripts expect the following input files in the working directory:

| File | Used by | Description |
|------|---------|-------------|
| `FD_<Tissue>.xlsx` | `01_per-tissue-analysis/` | Raw TMT protein intensities per tissue |
| `meltTemperature.xlsx` | `Insolubility_meltTemp.R` | Protein melting temperatures |
| `protein half-life.xlsx` | `insolubility_protein_half-life.R` | Tissue-specific protein half-lives |
| `mRNA levels.xlsx` | `Organs_mRNA.R` | Average TPM per gene per tissue |
| `Data S2_TMT significanly-regulated proteins.xlsx` | `mouseTMT_insolubility_subcellular.R` | Published significant protein list |

---

## Organ Color Palette

A consistent color palette is used across all cross-organ visualizations:

| Tissue | Hex |
|--------|-----|
| Cerebellum | `#c7a6db` |
| Heart | `#c93b92` |
| Hippocampus | `#55397f` |
| Liver | `#4d9fdb` |
| Neocortex | `#8a68ca` |
| Soleus | `#a8192a` |
| Tibialis Anterior | `#e991b7` |
| White Adipose Tissue | `#dfc41b` |

Age comparison colors: Old vs. Young = `purple4` · Middle vs. Young = `maroon2` · Old vs. Middle = `darkgreen`

---

## Usage

Set the working directory to the folder containing your input data files before running any script:

```r
setwd("path/to/your/data")
source("01_per-tissue-analysis/mouseTMT_liver.R")
```

Scripts within each stage are independent of each other within that stage, but Stage 2 and Stage 3 scripts expect the CSV outputs from Stage 1 to be present in the working directory.

**Recommended run order:**
1. Run all 8 scripts in `01_per-tissue-analysis/` (any order)
2. Run `Organs.R` → `Organs_abundance.R` → `Organs_mRNA.R` → `Organs_abundance_upset_plot.R`
3. Run scripts in `03_insolubility-sub-analysis/` (any order, after Stage 1 and 2)

<p align="center">
  <img src="image1.png" width="800">
</p>


---

## Author

**Chia-Lung Chuang**  
Postdoctoral Research Associate  
Department of Developmental Neurobiology  
St. Jude Children's Research Hospital  

---

## Citation

If you use this code, please cite the associated publication:

Stephan A*, Graca FA*, Pagala VR*, Hunt LC, Chuang C-L, Grime J,
Alford D, Wang X, High AA, Peng J, Demontis F. Insolubilome profiling
defines molecular features that influence protein insolubility with
aging. Developmental Cell. 2026 Feb 19.
https://doi.org/10.1016/j.devcel.2026.01.015

* These authors contributed equally.
