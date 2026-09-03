# Insulin modulates mPFC gene expression and emotional behavior in a sex-specific manner following fetal growth restriction

This repository contains the R scripts used to reproduce the transcriptomic analyses presented in the manuscript:

**Miguel PM et al.**
*Insulin modulates mPFC gene expression and emotional behavior in a sex-specific manner following fetal growth restriction.*

Published in **Brain, Behavior, and Immunity** (2026).

DOI: https://doi.org/10.1016/j.bbi.2026.106887

---

## Repository contents

### `RRHO2-analysis.R`

This script performs **Rank-Rank Hypergeometric Overlap (RRHO2)** analyses to compare genome-wide differential expression signatures between experimental conditions.

The workflow:

* imports differential expression results;
* ranks genes according to signed differential expression (`−log10(P-value) × sign(logFC)`);
* performs pairwise RRHO analyses using the **RRHO2** package;
* generates RRHO heatmaps for each comparison.

---

### `Cell-type-enrichment-analyses.R`

This script evaluates the enrichment of differentially expressed genes within published prefrontal cortex cell-type marker gene sets.

The workflow:

* filters differentially expressed genes according to statistical significance and fold-change thresholds;
* retains genes with valid Ensembl/BioMart annotations;
* restricts the analysis to protein-coding genes;
* compares differentially expressed genes with marker genes for astrocytes, microglia, oligodendrocytes, and pyramidal neurons;
* performs enrichment analyses using the **GeneOverlap** package;
* generates heatmaps displaying enrichment statistics.

---

## Input data

The scripts require differential expression results generated from the RNA-seq analyses.

These files are publicly available through the Gene Expression Omnibus:

**GEO accession:** GSE295966

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295966

Example input files:

```text
GSE295966_P90_Males_FR_comparison_mPFC.txt
GSE295966_P90_Males_INS_comparison_mPFC.txt
GSE295966_P90_Males_FRINS_comparison_mPFC.txt
```

Equivalent files are available for the other experimental conditions included in the study.

---

## Reference files for the cell-type enrichment analyses

In addition to the differential expression files available through GEO accession **GSE295966**, the `Cell-type-enrichment-analyses.R` script uses the following reference files.

### `Rat_external_gene_name_available_in_biomart.csv`

This file was used to retain genes with valid Ensembl/BioMart annotations.

It was generated on **19 February 2024** using the **biomaRt** package in R and the Ensembl database version available on that date.

---

### `protein_coding_rat_genes.xlsx`

This file contains the Ensembl identifiers of rat protein-coding genes and was used to restrict the enrichment analyses to protein-coding genes.

---

### `Rat_genes_for_PFC.csv`

This file contains the rat cell-type marker genes used as the reference panel for the enrichment analyses.

The file was generated using the following steps:

1. Mouse cell-type marker genes were obtained from the study indexed under **PMID: 29204516**.

2. The original marker gene files were obtained from the NeuroExpresso resource:

   https://github.com/PavlidisLab/neuroExpressoAnalysis/tree/master/analysis/01.SelectGenes/Markers_1.2

3. Mouse marker genes were converted to rat orthologs separately for each cell type using the **biomaRt** package in R and the Ensembl database version available on **19 February 2024**.

4. The resulting rat orthologs were compiled into the final reference panel used for the cell-type enrichment analyses.

---

## Software requirements

The scripts were developed in **R** and require the following packages:

* tidyverse
* RRHO2
* data.table
* GeneOverlap
* qpcR
* gplots
* RColorBrewer
* readxl
* biomaRt

---

## Data availability

RNA-seq data and differential expression results are publicly available through the Gene Expression Omnibus under accession:

**GSE295966**

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295966

---

## Citation

If you use these scripts, reference files, or associated data, please cite:

Miguel PM et al.
*Insulin modulates mPFC gene expression and emotional behavior in a sex-specific manner following fetal growth restriction.*
**Brain, Behavior, and Immunity** (2026).

DOI: https://doi.org/10.1016/j.bbi.2026.106887
