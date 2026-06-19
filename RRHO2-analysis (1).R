# RRHO analysis
#
# Study: Insulin modulates mPFC gene expression and emotional behavior in a sex-specific manner following fetal growth restriction.
#

library(tidyverse)
library(RRHO2)

# -----------------------------
# Parameters
# -----------------------------

boundary_value <- 0.03
heatmap_min <- 0
heatmap_max <- 300

# -----------------------------
# Function to read and prepare DE data
# -----------------------------

read_de_data <- function(file) {
  data <- read.table(
    file,
    sep = "\t",
    header = TRUE
  )
  
  if ("Pvalue" %in% colnames(data)) {
    data <- data %>%
      rename(PValue = Pvalue)
  }
  
  data %>%
    select(ensembl_gene_id, logFC, PValue) %>%
    mutate(
      dde = -log10(PValue) * sign(logFC)
    ) %>%
    filter(
      !is.na(ensembl_gene_id),
      !is.na(dde),
      dde != 0
    ) %>%
    select(ensembl_gene_id, dde)
}

# -----------------------------
# Read differential expression data
#
# The DE data files to be input are at the following GEO page:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295966
# 
# The example input filenames below correspond to the Male and P90 conditions.
# Adjust for other conditions with the files in the GEO page as appropriate. 
#
# -----------------------------

insulin_data <- read_de_data(
  "GSE295966_P90_Males_INS_comparison_mPFC.txt"
)

fr_data <- read_de_data(
  "GSE295966_P90_Males_FR_comparison_mPFC.txt"
)

fr_insulin_data <- read_de_data(
  "GSE295966_P90_Males_FRINS_comparison_mPFC.txt"
)

# -----------------------------
# Function to run RRHO
# -----------------------------

run_rrho <- function(data1, data2, label1, label2) {
  data1_rrho <- semi_join(
    data1,
    data2,
    by = "ensembl_gene_id"
  )
  
  data2_rrho <- semi_join(
    data2,
    data1,
    by = "ensembl_gene_id"
  )
  
  RRHO2_initialize(
    data1_rrho,
    data2_rrho,
    labels = c(label1, label2),
    boundary = boundary_value,
    log10.ind = TRUE
  )
}

# -----------------------------
# RRHO analyses
# -----------------------------

rrho_insulin_vs_fr <- run_rrho(
  insulin_data,
  fr_data,
  label1 = "INS",
  label2 = "FR"
)

rrho_fr_vs_frins <- run_rrho(
  fr_data,
  fr_insulin_data,
  label1 = "FR",
  label2 = "FRINS"
)

rrho_insulin_vs_frins <- run_rrho(
  insulin_data,
  fr_insulin_data,
  label1 = "INS",
  label2 = "FRINS"
)

# -----------------------------
# Plot RRHO maps
# -----------------------------

RRHO2_heatmap(
  rrho_insulin_vs_fr,
  maximum = heatmap_max,
  minimum = heatmap_min
)

RRHO2_heatmap(
  rrho_fr_vs_frins,
  maximum = heatmap_max,
  minimum = heatmap_min
)

RRHO2_heatmap(
  rrho_insulin_vs_frins,
  maximum = heatmap_max,
  minimum = heatmap_min
)
