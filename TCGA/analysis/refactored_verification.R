# Verification script for refactored correlation analysis
# This script tests the refactored version and compares with original results

library(dplyr)
library(tidyr)
library(here)

# Load the pre-computed original results (from memory)
# These should already be loaded in the session:
# - coeffs: correlation coefficients matrix
# - pvals: p-values matrix
# - r2: r-squared matrix

# Load necessary data
mir_gene_pairs <- readRDS(here("TCGA", "output", "rds", "TCGA_selected_miRNA_host_genes.rds"))
mir_gene_pairs_fltr <- unique(mir_gene_pairs[, c("gene_id", "pri_mir_id", "pri_mir_id_2")])

mir_gene_pairs_fltr$miR <- paste(mir_gene_pairs_fltr$pri_mir_id, "-3p", sep="")
mir_gene_pairs_fltr_5p <- mir_gene_pairs_fltr
mir_gene_pairs_fltr_5p$miR <- gsub("-3p", "-5p", mir_gene_pairs_fltr_5p$miR)
mir_gene_pairs_fltr <- rbind(mir_gene_pairs_fltr, mir_gene_pairs_fltr_5p)
mir_gene_pairs_fltr <- mir_gene_pairs_fltr[order(mir_gene_pairs_fltr$gene), ]

dfi_data <- readRDS(here("TCGA", "output", "rds", "TCGA_dfi_long.rds"))

# ===== REFACTORED VERSION =====

# Pre-load all miRNA expression files by cancer type
miR_expr_by_cancer <- sapply(
  unique(dfi_data$cancer),
  function(cancer) {
    read.table(
      here("TCGA", "data", paste0(cancer, ".miRseq_mature_RPM_log2.txt")),
      header = TRUE
    )
  },
  simplify = FALSE
)

# Helper function to process a single gene-miRNA pair
compute_single_pair <- function(i, mir_gene_pairs_fltr, dfi_data, miR_expr_by_cancer) {
  gene_id <- mir_gene_pairs_fltr$gene_id[i]
  miRNA <- mir_gene_pairs_fltr$miR[i]
  
  # Filter dfi table for the gene of interest
  dfi_gene <- dfi_data[dfi_data$gene_id == gene_id, ]
  
  if (nrow(dfi_gene) == 0) {
    return(NULL)
  }
  
  # Apply correlation computation across cancers where this gene appears
  cancer_results <- lapply(unique(dfi_gene$cancer), function(cancer) {
    miR_expr <- miR_expr_by_cancer[[cancer]]
    
    # Check if miRNA exists in the abundance data
    if (!(miRNA %in% row.names(miR_expr))) {
      return(NULL)
    }
    
    # Extract and format miRNA abundance
    miR_cancer <- as.data.frame(cbind(colnames(miR_expr), unlist(miR_expr[miRNA, ])))
    colnames(miR_cancer) <- c("sample", "miR_expr")
    miR_cancer$sample <- gsub("\\.", "-", miR_cancer$sample)
    miR_cancer <- miR_cancer[grep("-01$", miR_cancer$sample), ]
    
    # Skip if too many missing values (>20%)
    if (sum(is.na(miR_cancer$miR_expr)) / nrow(miR_cancer) > 0.2) {
      return(NULL)
    }
    
    # Combine with host gene transcription data
    cancer_data <- dfi_gene[dfi_gene$cancer == cancer, ]
    cancer_data <- cancer_data[grep("\\-01$", cancer_data$sample), ]
    cancer_data$miRNA_expr <- as.numeric(
      miR_cancer[match(row.names(cancer_data), miR_cancer$sample), ]$miR_expr
    )
    cancer_data <- cancer_data[rowSums(is.na(cancer_data)) == 0, ]
    
    if (nrow(cancer_data) == 0) {
      return(NULL)
    }
    
    # Compute correlation
    corr <- cor.test(
      as.numeric(cancer_data$dfi),
      cancer_data$miRNA_expr
    )
    
    # Return result as a data frame row
    data.frame(
      miRNA = miRNA,
      cancer = cancer,
      coefficient = corr$estimate,
      p_value = corr$p.value,
      r_squared = (corr$estimate)^2,
      stringsAsFactors = FALSE
    )
  })
  
  # Combine results for this pair
  do.call(rbind, Filter(Negate(is.null), cancer_results))
}

# Apply across all gene-miRNA pairs and combine into a single long dataframe
results_list <- lapply(
  1:nrow(mir_gene_pairs_fltr),
  compute_single_pair,
  mir_gene_pairs_fltr = mir_gene_pairs_fltr,
  dfi_data = dfi_data,
  miR_expr_by_cancer = miR_expr_by_cancer
)

corr_results_long <- do.call(rbind, Filter(Negate(is.null), results_list))
row.names(corr_results_long) <- NULL

# Convert to wide format for compatibility with downstream code
coeffs_new <- pivot_wider(
  corr_results_long,
  names_from = cancer,
  values_from = coefficient,
  values_fill = NA
) |> column_to_rownames("miRNA") |> as.matrix() |> as.data.frame()

pvals_new <- pivot_wider(
  corr_results_long,
  names_from = cancer,
  values_from = p_value,
  values_fill = NA
) |> column_to_rownames("miRNA") |> as.matrix() |> as.data.frame()

r2_new <- pivot_wider(
  corr_results_long,
  names_from = cancer,
  values_from = r_squared,
  values_fill = NA
) |> column_to_rownames("miRNA") |> as.matrix() |> as.data.frame()

# ===== VERIFICATION =====
cat("=== VERIFICATION OF REFACTORED CODE ===\n\n")

cat("Original coeffs dimensions:", nrow(coeffs), "x", ncol(coeffs), "\n")
cat("New coeffs dimensions:", nrow(coeffs_new), "x", ncol(coeffs_new), "\n\n")

cat("Original pvals dimensions:", nrow(pvals), "x", ncol(pvals), "\n")
cat("New pvals dimensions:", nrow(pvals_new), "x", ncol(pvals_new), "\n\n")

cat("Original r2 dimensions:", nrow(r2), "x", ncol(r2), "\n")
cat("New r2 dimensions:", nrow(r2_new), "x", ncol(r2_new), "\n\n")

# Check if dimensions match
if (all.equal(dim(coeffs), dim(coeffs_new)) == TRUE) {
  cat("✓ Dimensions match for coeffs\n")
} else {
  cat("✗ Dimensions DO NOT match for coeffs\n")
}

# Compare values (accounting for potential row/column ordering differences)
# Align both matrices by row and column names
coeffs_aligned <- coeffs[rownames(coeffs) %in% rownames(coeffs_new), colnames(coeffs) %in% colnames(coeffs_new)]
coeffs_new_aligned <- coeffs_new[rownames(coeffs_new) %in% rownames(coeffs), colnames(coeffs_new) %in% colnames(coeffs)]

# Reorder to match
coeffs_new_aligned <- coeffs_new_aligned[match(rownames(coeffs_aligned), rownames(coeffs_new_aligned)), 
                                         match(colnames(coeffs_aligned), colnames(coeffs_new_aligned))]

cat("\nComparing coefficient values (first 5x5 block):\n")
print(head(coeffs_aligned, 5))
cat("\nNew version (first 5x5 block):\n")
print(head(coeffs_new_aligned, 5))

# Calculate difference
max_diff_coeffs <- max(abs(coeffs_aligned - coeffs_new_aligned), na.rm = TRUE)
cat("\nMaximum difference in coefficients:", max_diff_coeffs, "\n")

max_diff_pvals <- max(abs(pvals[rownames(pvals) %in% rownames(pvals_new), colnames(pvals) %in% colnames(pvals_new)] - 
                           pvals_new[rownames(pvals_new) %in% rownames(pvals), colnames(pvals_new) %in% colnames(pvals)], 
                           na.rm = TRUE)
cat("Maximum difference in p-values:", max_diff_pvals, "\n")

max_diff_r2 <- max(abs(r2[rownames(r2) %in% rownames(r2_new), colnames(r2) %in% colnames(r2_new)] - 
                        r2_new[rownames(r2_new) %in% rownames(r2), colnames(r2_new) %in% colnames(r2)], 
                        na.rm = TRUE)
cat("Maximum difference in r2:", max_diff_r2, "\n")

cat("\nRefactored version produces a long-format dataframe:\n")
cat("Dimensions:", nrow(corr_results_long), "rows x", ncol(corr_results_long), "columns\n")
cat("Columns:", paste(colnames(corr_results_long), collapse = ", "), "\n")
cat("First few rows:\n")
print(head(corr_results_long, 10))
