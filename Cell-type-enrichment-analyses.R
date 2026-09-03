# Gene overlap analysis between differentially expressed genes
# and PFC cell-type marker genes.
#
# Study: Insulin modulates mPFC gene expression and emotional behavior
# in a sex-specific manner following fetal growth restriction.

library(data.table)
library(tidyverse)
library(GeneOverlap)
library(qpcR)
library(gplots)
library(RColorBrewer)

# -----------------------------
# Parameters
# -----------------------------

genome_size <- 21294
p_cutoff <- 0.05
logfc_cutoff <- 0.3875

# -----------------------------
# Reference files
# -----------------------------

rat_ensembl <- fread("Rat_external_gene_name_available_in_biomart.csv")

protein_coding_genes <- readxl::read_excel(
  "protein_coding_rat_genes.xlsx",
  col_names = FALSE
)

names(protein_coding_genes) <- c("ensembl_gene_id", "protein_coding")

# -----------------------------
# Cell-type marker genes
# -----------------------------

marker_genes <- read.csv("Rat_genes_for_PFC.csv")[, c(1:3, 11)] %>%
  as.data.table()

names(marker_genes) <- c(
  "Astrocyte",
  "Microglia",
  "Oligo",
  "Pyramidal"
)

# -----------------------------
# Function to read and filter DE genes
# -----------------------------

read_de_genes <- function(file, comparison_name) {
  read.table(
    file,
    sep = "\t",
    header = TRUE
  ) %>%
    dplyr::filter(PValue <= p_cutoff, abs(logFC) > logfc_cutoff) %>%
    dplyr::filter(external_gene_name %in% rat_ensembl$external_gene_name) %>%
    dplyr::filter(ensembl_gene_id %in% protein_coding_genes$ensembl_gene_id) %>%
    dplyr::select(external_gene_name) %>%
    dplyr::distinct() %>%
    dplyr::filter(
      !is.na(external_gene_name),
      external_gene_name != "NA",
      external_gene_name != ""
    ) %>%
    dplyr::rename(!!comparison_name := external_gene_name) %>%
    as.data.table()
}

# -----------------------------
# Custom heatmap function
# -----------------------------

drawHeatmap_custom <- function(object,
                               what = c("odds.ratio", "Jaccard"),
                               log.scale = FALSE,
                               adj.p = FALSE,
                               cutoff = 0.05,
                               ncolused = 5,
                               grid.col = c("Greens", "Blues", "Greys",
                                            "Oranges", "Purples", "Reds")) {
  
  what <- match.arg(what)
  grid.col <- match.arg(grid.col)
  
  pv.mat <- getMatrix(object, "pval")
  
  plot.mat <- switch(
    what,
    odds.ratio = getMatrix(object, "odds.ratio"),
    Jaccard = getMatrix(object, "Jaccard")
  )
  
  if (what == "odds.ratio" && log.scale) {
    plot.mat <- log2(plot.mat)
  }
  
  if (adj.p) {
    pv.mat <- matrix(
      p.adjust(pv.mat, method = "BH"),
      nrow = nrow(pv.mat)
    )
  }
  
  insig.val <- 1
  
  if ((what == "odds.ratio" && log.scale) || what == "Jaccard") {
    insig.val <- 0
  }
  
  plot.mat[pv.mat >= cutoff] <- insig.val
  
  note.mat <- format(pv.mat, digits = 1)
  
  note.mat[pv.mat < 0.01] <- format(
    pv.mat,
    digits = 1,
    scientific = TRUE
  )[pv.mat < 0.01]
  
  # non-significant cells as dots
  note.mat[plot.mat == insig.val] <- "•"
  
  row_sep <- 1:(nrow(plot.mat) - 1)
  col_sep <- 1:(ncol(plot.mat) - 1)
  
  longedge <- max(nrow(plot.mat), ncol(plot.mat))
  
  row_cexrc <- 0.4 + 1 / log10(longedge + 2)
  col_cexrc <- row_cexrc
  key_size <- 0.2 + 1 / log10(longedge + 4)
  
  margins_use <- c(
    max(nchar(colnames(plot.mat))) * 0.8 + 5,
    max(nchar(rownames(plot.mat))) * 0.8 + 5
  )
  
  main.txt <- switch(
    what,
    odds.ratio = ifelse(log.scale, "log2(Odds Ratio)", "Odds Ratio"),
    Jaccard = "Jaccard Index"
  )
  
  gplots::heatmap.2(
    plot.mat,
    cellnote = note.mat,
    main = main.txt,
    xlab = "",
    col = RColorBrewer::brewer.pal(ncolused, grid.col),
    notecol = "white",
    margins = margins_use,
    colsep = col_sep,
    rowsep = row_sep,
    key = TRUE,
    keysize = key_size,
    cexRow = row_cexrc,
    cexCol = col_cexrc,
    scale = "none",
    Colv = NA,
    Rowv = NA,
    trace = "none",
    dendrogram = "none",
    density.info = "none",
    sepcolor = "white",
    sepwidth = c(0.002, 0.002),
    notecex = 1.6
  )
}

  nr <- nrow(note.mat)
  nc <- ncol(note.mat)
  
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      label <- note.mat[i, j]
      
      text_col <- ifelse(label == "•", "black", "white")
      text_cex <- ifelse(label == "•", 1.8, 1.2)
      
      graphics::text(
        x = j,
        y = nr - i + 1,
        labels = label,
        col = text_col,
        cex = text_cex,
        font = 2
      )
    }
  }

# -----------------------------
# Read differential expression data
# -----------------------------

FR <- read_de_genes(
  "GSE295966_P90_Males_FR_comparison_mPFC.txt",
  "FR"
)

INS <- read_de_genes(
  "GSE295966_P90_Males_INS_comparison_mPFC.txt",
  "INS"
)

FRINS <- read_de_genes(
  "GSE295966_P90_Males_FRINS_comparison_mPFC.txt",
  "FRINS"
)

# -----------------------------
# Combine gene lists
# -----------------------------

de_gene_lists <- qpcR:::cbind.na(FR, INS, FRINS)

# -----------------------------
# Check gene list sizes
# -----------------------------

sapply(marker_genes, length)
sapply(de_gene_lists, length)

# -----------------------------
# Gene overlap analysis
# -----------------------------

gom_obj <- newGOM(
  de_gene_lists,
  marker_genes,
  genome_size
)

# -----------------------------
# Plot heatmap
# -----------------------------

drawHeatmap_custom(
  gom_obj,
  ncolused = 5,
  grid.col = "Blues"
)