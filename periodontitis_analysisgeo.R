# ANALYSIS OF GEO DATASET FOR PERIODONTITIS
# Let us create the project structure
# Set the working directory
setwd("D:/MHPE7/periodontitis")

# Create analysis folders
dir.create("HLA_Analysis", showWarnings = FALSE)
dir.create("HLA_Analysis/Results", showWarnings = FALSE)
dir.create("HLA_Analysis/Plots", showWarnings = FALSE)
dir.create("HLA_Analysis/Checkpoints", showWarnings = FALSE)

# Set analysis path
analysis_path <- file.path(getwd(), "HLA_Analysis")

###################################################################################################################

# ============================================================================
# LOAD REQUIRED PACKAGES
# ============================================================================

# Install if needed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages_needed <- c(
  # Bioconductor
  "oligo", "limma", "GEOquery",
  "hgu133plus2.db", "org.Hs.eg.db", "AnnotationDbi",
  # CRAN
  "tidyverse", "pheatmap", "RColorBrewer",
  "ggplot2", "reshape2", "corrplot"
)

for(pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE)) {
    if(pkg %in% c("oligo", "limma", "GEOquery", "hgu133plus2.db", "org.Hs.eg.db", "AnnotationDbi")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

cat("✓ All packages loaded\n")

##########################################################################################################################

# ============================================================================
# LOAD AND PROCESS GPL570 ANNOTATION
# ============================================================================

cat("\n========================================\n")
cat("Loading GPL570 annotation\n")
cat("========================================\n")

# Check if annotation file exists
if(file.exists("GPL570.soft.gz")) {
  cat("Found GPL570.soft.gz\n")
  
  # Read the SOFT file
  gpl <- getGEO(filename = "GPL570.soft.gz")
  
  # Extract annotation table
  gpl_table <- Table(gpl)
  cat("Annotation table dimensions:", dim(gpl_table), "\n")
  cat("Column names:", paste(colnames(gpl_table)[1:10], collapse=", "), "...\n")
  
  # Save for later use
  saveRDS(gpl_table, "HLA_Analysis/Checkpoints/gpl570_annotation.rds")
  cat("✓ Saved GPL570 annotation\n")
  
} else {
  cat("GPL570.soft.gz not found. Downloading from GEO...\n")
  gpl <- getGEO("GPL570", destdir = ".")
  gpl_table <- Table(gpl)
  saveRDS(gpl_table, "HLA_Analysis/Checkpoints/gpl570_annotation.rds")
}

##########################################################################################################################

# ============================================================================
# COMPILE COMPREHENSIVE HLA GENE LIST
# ============================================================================

cat("\n========================================\n")
cat("Creating HLA gene list\n")
cat("========================================\n")

# Comprehensive HLA and related genes
hla_genes <- c(
  # HLA Class I classical
  "HLA-A", "HLA-B", "HLA-C",
  # HLA Class I non-classical
  "HLA-E", "HLA-F", "HLA-G",
  # HLA Class II classical
  "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
  "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
  # HLA Class II non-classical (antigen loading)
  "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
  # Antigen processing machinery
  "TAP1", "TAP2", "TAPBP", "B2M", "CD74",
  # Immunoproteasome subunits
  "PSMB8", "PSMB9", "PSMB10",
  # MHC class II transcriptional regulators
  "CIITA", "RFX5", "RFXAP", "RFXANK", "NFYA", "NFYB", "NFYC"
)

cat("Total HLA-related genes:", length(hla_genes), "\n")
print(hla_genes)

# Save gene list
saveRDS(hla_genes, "HLA_Analysis/Checkpoints/hla_gene_list.rds")
cat("✓ Saved HLA gene list\n")

#########################################################################################################################

# ============================================================================
# MAP HLA GENES TO AFFYMETRIX PROBE IDs
# ============================================================================

cat("\n========================================\n")
cat("Mapping HLA genes to probe IDs\n")
cat("========================================\n")

# Method 1: Using annotation package (preferred) #######################################################################
library(hgu133plus2.db)

# Get all probe IDs for HLA genes
hla_probe_map <- mapIds(hgu133plus2.db,
                        keys = hla_genes,
                        column = "PROBEID",
                        keytype = "SYMBOL",
                        multiVals = "list")

# Convert to data frame for easy viewing
hla_probes_df <- data.frame(
  Gene = rep(names(hla_probe_map), lengths(hla_probe_map)),
  ProbeID = unlist(hla_probe_map),
  row.names = NULL
)

cat("\nProbe mapping summary:\n")
print(table(hla_probes_df$Gene))

# Get unique probe IDs
hla_probe_ids <- unique(unlist(hla_probe_map))
cat("\nTotal unique probe IDs:", length(hla_probe_ids), "\n")

# Method 2: Using GPL570 annotation (backup) #############################################################################
# Get gene symbols from annotation
gpl_table <- readRDS("HLA_Analysis/Checkpoints/gpl570_annotation.rds")

# Find which columns contain gene symbols
gene_cols <- grep("Gene.Symbol|GeneSymbol|gene_symbol|Gene Symbol", 
                  colnames(gpl_table), ignore.case = TRUE, value = TRUE)
cat("\nGene symbol columns in GPL570:", paste(gene_cols, collapse=", "), "\n")

# Use the first gene symbol column
gene_col <- gene_cols[1]

# Extract probe-gene mappings
probe_gene_map <- data.frame(
  ProbeID = gpl_table$ID,
  GeneSymbol = gpl_table[[gene_col]],
  stringsAsFactors = FALSE
)

# Clean up multiple gene symbols (some probes map to multiple genes)
probe_gene_map$GeneSymbol <- gsub(" /// .*", "", probe_gene_map$GeneSymbol)

# Filter for HLA genes
hla_from_gpl <- probe_gene_map[probe_gene_map$GeneSymbol %in% hla_genes, ]

cat("\nHLA probes found in GPL570:", nrow(hla_from_gpl), "\n")

# Save both mappings
saveRDS(hla_probes_df, "HLA_Analysis/Checkpoints/hla_probe_mapping.rds")
saveRDS(hla_probe_ids, "HLA_Analysis/Checkpoints/hla_probe_ids.rds")
saveRDS(hla_from_gpl, "HLA_Analysis/Checkpoints/hla_from_gpl.rds")

cat("\n✓ Saved HLA probe mappings\n")

###########################################################################################################################

# ============================================================================
# READ AND NORMALIZE CEL FILES (LONGEST STEP - CHECKPOINTED)
# ============================================================================

cat("\n========================================\n")
cat("Reading and normalizing CEL files\n")
cat("========================================\n")

# Check if normalized data already exists
if(file.exists("HLA_Analysis/Checkpoints/exprs_matrix_all.rds")) {
  cat("Loading previously normalized data...\n")
  exprs_matrix <- readRDS("HLA_Analysis/Checkpoints/exprs_matrix_all.rds")
  pheno_data <- readRDS("HLA_Analysis/Checkpoints/pheno_data.rds")
  
} else {
  # Get all CEL files
  cel_files <- list.files(pattern = "\\.CEL$", recursive = TRUE, full.names = TRUE)
  cat("Found", length(cel_files), "CEL files\n")
  
  # Read CEL files (this takes time)   #########################################################
  cat("Reading CEL files... (this will take 20-30 minutes)\n")
  start_time <- Sys.time()
  
  raw_data <- read.celfiles(cel_files)
  
  cat("Reading completed in", 
      round(difftime(Sys.time(), start_time, units = "mins"), 1), "minutes\n")
  
  # RMA Normalization
  cat("\nPerforming RMA normalization...\n")
  start_time <- Sys.time()
  
  norm_data <- rma(raw_data)
  exprs_matrix <- exprs(norm_data)
  
  cat("RMA completed in", 
      round(difftime(Sys.time(), start_time, units = "mins"), 1), "minutes\n")
  
  # Save normalized data
  saveRDS(exprs_matrix, "HLA_Analysis/Checkpoints/exprs_matrix_all.rds")
  cat("✓ Saved normalized expression matrix\n")
  
  # Clean up to save memory
  rm(raw_data, norm_data)
  gc()
}

cat("Expression matrix dimensions:", dim(exprs_matrix), "\n")

############################################################################################################################

# ============================================================================
# PARSING GEO SERIES MATRIX
# ============================================================================

cat("\n========================================\n")
cat("PARSING GEO SERIES MATRIX\n")
cat("========================================\n")

# Read the compressed file
con <- gzfile("GSE10334_series_matrix.txt.gz", "rt")
lines <- readLines(con)
close(con)

cat("Read", length(lines), "lines from compressed file\n")

# ----------------------------------------------------------------------------
# PART 1: Get sample IDs from data table header
# ----------------------------------------------------------------------------
data_start <- which(lines == "!series_matrix_table_begin") + 1
header_line <- lines[data_start]
sample_ids <- strsplit(header_line, "\t")[[1]][-1]  # Remove "ID_REF"
cat("\n✓ Found", length(sample_ids), "sample IDs in data table\n")

# ----------------------------------------------------------------------------
# PART 2: Get sample titles (all in one line, tab-separated)
# ----------------------------------------------------------------------------
title_line <- lines[grep("^!Sample_title", lines)]
cat("Found title line:", substr(title_line, 1, 50), "...\n")

# Split the title line and remove the first element ("!Sample_title")
all_titles <- strsplit(title_line, "\t")[[1]][-1]
cat("✓ Extracted", length(all_titles), "sample titles\n")

# ----------------------------------------------------------------------------
# PART 3: Verify they match
# ----------------------------------------------------------------------------
if(length(all_titles) == length(sample_ids)) {
  cat("✓ Number of titles matches number of samples\n")
} else {
  stop("CRITICAL ERROR: Title count (", length(all_titles), 
       ") != Sample count (", length(sample_ids), ")")
}

# ----------------------------------------------------------------------------
# PART 4: Create phenotype data frame
# ----------------------------------------------------------------------------
pheno_data <- data.frame(
  SampleID = sample_ids,
  Title = all_titles,
  stringsAsFactors = FALSE
)

# Extract disease status
pheno_data$Status <- ifelse(
  grepl("Affected site|Affected", pheno_data$Title, ignore.case = TRUE), 
  "Diseased",
  ifelse(grepl("Unaffected site|Unaffected|Healthy", pheno_data$Title, ignore.case = TRUE), 
         "Healthy", NA)
)

# Extract patient ID
pheno_data$PatientID <- sapply(strsplit(pheno_data$Title, " "), function(words) {
  # Find "patient" and get next word
  patient_pos <- which(words == "patient")
  if(length(patient_pos) > 0 && patient_pos < length(words)) {
    return(gsub(",", "", words[patient_pos + 1]))
  }
  return(NA)
})

# Extract sample number
pheno_data$SampleNum <- sapply(strsplit(pheno_data$Title, " "), function(words) {
  # Find "sample" and get next word
  sample_pos <- which(words == "sample")
  if(length(sample_pos) > 0 && sample_pos < length(words)) {
    return(gsub(",", "", words[sample_pos + 1]))
  }
  return(NA)
})

##########################################################################################################################
#########################################################################################################################

# ============================================================================
# Let us continue with our setup
# ============================================================================

setwd("D:/MHPE7/periodontitis")

cat("Loading existing data from checkpoints...\n")

# Load expression matrix
exprs_matrix <- readRDS("HLA_Analysis/Checkpoints/exprs_matrix_all.rds")

# Load phenotype data
if(file.exists("HLA_Analysis/Checkpoints/pheno_data_corrected.rds")) {
  pheno_data <- readRDS("HLA_Analysis/Checkpoints/pheno_data_corrected.rds")
} else {
  pheno_data <- readRDS("HLA_Analysis/Checkpoints/pheno_data.rds")
}

# Load HLA probe data
hla_probe_ids <- readRDS("HLA_Analysis/Checkpoints/hla_probe_ids.rds")
hla_probes_df <- readRDS("HLA_Analysis/Checkpoints/hla_probe_mapping.rds")

# Load best_per_gene results
if(file.exists("HLA_Analysis/Results/HLA_best_per_gene_final.csv")) {
  best_per_gene <- read.csv("HLA_Analysis/Results/HLA_best_per_gene_final.csv")
} else {
  best_per_gene <- read.csv("HLA_Analysis/Results/HLA_best_per_gene_corrected.csv")
}

# ----------------------------------------------------------------------------
# Get HLA expression matrix with proper row names
# ----------------------------------------------------------------------------

# Find probes present in expression matrix
probes_present <- intersect(hla_probe_ids, rownames(exprs_matrix))
cat("\nProbes present:", length(probes_present), "of", length(hla_probe_ids), "\n")

# Subset expression matrix to HLA probes
hla_exprs <- exprs_matrix[probes_present, ]

# Subset probe mapping to match
hla_probes_df_subset <- hla_probes_df[hla_probes_df$ProbeID %in% probes_present, ]

# IMPORTANT: Ensure they are in the SAME ORDER
hla_probes_df_subset <- hla_probes_df_subset[match(probes_present, hla_probes_df_subset$ProbeID), ]

# Check dimensions
cat("hla_exprs rows:", nrow(hla_exprs), "\n")
cat("hla_probes_df_subset rows:", nrow(hla_probes_df_subset), "\n")

# Only assign row names if dimensions match
if(nrow(hla_exprs) == nrow(hla_probes_df_subset)) {
  rownames(hla_exprs) <- paste0(hla_probes_df_subset$Gene, "_", hla_probes_df_subset$ProbeID)
  cat("✓ Row names assigned successfully\n")
} else {
  cat("⚠️ Dimension mismatch. Creating generic row names.\n")
  rownames(hla_exprs) <- paste0("HLA_", 1:nrow(hla_exprs))
  
  # Create a mapping file for reference
  probe_mapping <- data.frame(
    RowName = rownames(hla_exprs),
    ProbeID = probes_present,
    Gene = hla_probes_df_subset$Gene[match(probes_present, hla_probes_df_subset$ProbeID)]
  )
  write.csv(probe_mapping, "HLA_Analysis/Checkpoints/probe_row_mapping.csv", row.names = FALSE)
  cat("Created probe mapping file for reference\n")
}

# Filter to valid samples (with known status)
valid_idx <- !is.na(pheno_data$Status)
hla_exprs_valid <- hla_exprs[, valid_idx]
pheno_valid <- pheno_data[valid_idx, ]

cat("\nData loaded successfully!\n")
cat("Samples:", ncol(hla_exprs_valid), "(", 
    sum(pheno_valid$Status == "Diseased"), "Diseased,",
    sum(pheno_valid$Status == "Healthy"), "Healthy)\n")
cat("HLA probes:", nrow(hla_exprs_valid), "\n")

# Check if Rtsne and umap are installed, install if needed
if(!require("Rtsne")) {
  install.packages("Rtsne")
  library(Rtsne)
}

if(!require("umap")) {
  install.packages("umap")
  library(umap)
}

##########################################################################################################################
# ============================================================================
# T-SNE WITH SIGNIFICANT GENES
# ============================================================================

# Get significant genes (FDR < 0.05)
sig_genes <- best_per_gene$ProbeID[best_per_gene$adj.P.Val < 0.05]

# Find which significant genes are in our expression matrix
sig_genes_present <- intersect(sig_genes, rownames(hla_exprs_valid))
cat("Using", length(sig_genes_present), "significant probes for t-SNE\n")

# Subset to significant genes
hla_exprs_sig <- hla_exprs_valid[rownames(hla_exprs_valid) %in% sig_genes_present, ]

# Run t-SNE
set.seed(42)
tsne_result <- Rtsne(t(hla_exprs_sig), 
                     perplexity = 30,
                     theta = 0.5,
                     dims = 2,
                     pca = TRUE,
                     max_iter = 1000,
                     verbose = TRUE)

# Create dataframe
tsne_df <- data.frame(
  tSNE1 = tsne_result$Y[, 1],
  tSNE2 = tsne_result$Y[, 2],
  Status = pheno_valid$Status
)

# Plot
p_tsne <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Status)) +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(level = 0.95, linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Diseased" = "firebrick", "Healthy" = "steelblue")) +
  labs(title = "t-SNE of Significant HLA Genes",
       x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        legend.title = element_blank())

print(p_tsne)
ggsave("HLA_Analysis/Plots/TSNE_HLA.png", p_tsne, width = 8, height = 6, dpi = 300)

##########################################################################################################################
# ============================================================================
# QUANTIFY THE T-SNE SEPARATION
# ============================================================================

# Add group labels
tsne_df$Group <- pheno_valid$Status

# Calculate centroids (mean position) for each group
centroids <- aggregate(cbind(tSNE1, tSNE2) ~ Group, data = tsne_df, FUN = mean)
print("Centroids (mean positions):")
print(centroids)

# Calculate distance between centroids
dist_between <- sqrt((centroids[1,2] - centroids[2,2])^2 + 
                       (centroids[1,3] - centroids[2,3])^2)
cat("\nDistance between group centroids:", round(dist_between, 3), "\n")

# Calculate within-group spread (average distance from centroid)
library(cluster)

# For Healthy group
healthy_dist <- dist(tsne_df[tsne_df$Group == "Healthy", c("tSNE1", "tSNE2")])
healthy_spread <- mean(healthy_dist)
cat("Healthy group spread (avg pairwise distance):", round(healthy_spread, 3), "\n")

# For Diseased group
diseased_dist <- dist(tsne_df[tsne_df$Group == "Diseased", c("tSNE1", "tSNE2")])
diseased_spread <- mean(diseased_dist)
cat("Diseased group spread (avg pairwise distance):", round(diseased_spread, 3), "\n")

# Ratio of spreads
cat("Diseased/Healthy spread ratio:", round(diseased_spread/healthy_spread, 3), "\n")

# Silhouette score (measures cluster quality)
dist_matrix <- dist(tsne_df[, c("tSNE1", "tSNE2")])
sil <- silhouette(as.numeric(factor(tsne_df$Group)), dist_matrix)
mean_sil <- mean(sil[, 3])
cat("\nMean silhouette score:", round(mean_sil, 3))
cat("\n(0.25-0.5 = weak structure, 0.5-0.7 = reasonable structure, >0.7 = strong structure)\n")

###########################################################################################################################
# ============================================================================
# ENHANCED T-SNE PLOT
# ============================================================================

library(ggrepel)  # For labels if needed

# Create enhanced plot
p_tsne_enhanced <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Status)) +
  # Add density contours
  stat_density2d(aes(fill = after_stat(level)), 
                 geom = "polygon", 
                 alpha = 0.2,
                 bins = 10) +
  # Add points
  geom_point(size = 2.5, alpha = 0.8) +
  # Add ellipses
  stat_ellipse(level = 0.95, linewidth = 1, linetype = "dashed") +
  # Colors
  scale_color_manual(values = c("Diseased" = "#D55E00", "Healthy" = "#0072B2")) +
  scale_fill_gradient(low = "gray90", high = "gray50", guide = "none") +
  # Labels
  labs(title = "t-SNE Visualization of HLA Gene Expression",
       subtitle = paste0("27 significant HLA genes · 247 samples · Perplexity = 30"),
       x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2",
       caption = "Ellipses show 95% confidence regions") +
  # Theme
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
        panel.grid.minor = element_blank())

print(p_tsne_enhanced)
ggsave("HLA_Analysis/Plots/TSNE_HLA_enhanced.png", p_tsne_enhanced, 
       width = 9, height = 7, dpi = 300)


############################################################################################################################
############################################################################################################################

# ============================================================================
# Genes vs Probes
# ============================================================================

# Count unique significant GENES
sig_genes_unique <- unique(best_per_gene$Gene[best_per_gene$adj.P.Val < 0.05])
cat("Unique significant GENES:", length(sig_genes_unique), "\n")

# Count significant PROBES
sig_probes_all <- best_per_gene$ProbeID[best_per_gene$adj.P.Val < 0.05]
cat("Significant PROBES:", length(sig_probes_all), "\n")

# Show genes with multiple significant probes
probe_counts <- table(best_per_gene$Gene[best_per_gene$adj.P.Val < 0.05])
multi_probe_genes <- probe_counts[probe_counts > 1]
cat("\nGenes with multiple significant probes:\n")
print(multi_probe_genes)

# List the downregulated genes
down_genes <- best_per_gene[best_per_gene$logFC < 0 & 
                              best_per_gene$adj.P.Val < 0.05, 
                            c("Gene", "logFC", "adj.P.Val")]
cat("\nDownregulated genes:\n")
print(down_genes)

##########################################################################################################################
# ============================================================================
# FINAL VERIFICATION OF ALL 34 GENES
# ============================================================================

# Complete list of 34 HLA-related genes
all_34_genes <- c(
  "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
  "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
  "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
  "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
  "TAP1", "TAP2", "TAPBP", "B2M", "CD74",
  "PSMB8", "PSMB9", "PSMB10",
  "CIITA", "RFX5", "RFXAP", "RFXANK", "NFYA", "NFYB", "NFYC"
)

# Which genes are in best_per_gene?
genes_in_results <- best_per_gene$Gene

# Find missing genes (not in best_per_gene at all)
missing_genes <- setdiff(all_34_genes, genes_in_results)
cat("Genes completely missing from results:\n")
print(missing_genes)
# These are HLA-DRB3, HLA-DRB5 (no unique probes)

# Find non-significant genes (in results but adj.P.Val >= 0.05)
non_sig <- best_per_gene$Gene[best_per_gene$adj.P.Val >= 0.05]
cat("\nNon-significant genes (FDR >= 0.05):\n")
print(non_sig)

# Complete summary
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("FINAL SUMMARY\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
cat("Total HLA-related genes:", length(all_34_genes), "\n")
cat("Genes with no probe representation:", length(missing_genes), "\n")
cat("Genes analyzed (had probes):", length(unique(best_per_gene$Gene)), "\n")
cat("  - Significant (FDR < 0.05):", sum(best_per_gene$adj.P.Val < 0.05), "\n")
cat("    * Upregulated:", sum(best_per_gene$logFC > 0 & best_per_gene$adj.P.Val < 0.05), "\n")
cat("    * Downregulated:", sum(best_per_gene$logFC < 0 & best_per_gene$adj.P.Val < 0.05), "\n")
cat("  - Non-significant:", sum(best_per_gene$adj.P.Val >= 0.05), "\n")

##########################################################################################################################
##########################################################################################################################
# Final t-SNE with accurate counts
p_tsne_final <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Status)) +
  geom_point(size = 2.5, alpha = 0.7) +
  stat_ellipse(level = 0.95, linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Diseased" = "#D55E00", "Healthy" = "#0072B2")) +
  labs(title = "t-SNE of Dysregulated HLA Genes in Periodontitis",
       subtitle = paste0("27 significant genes · 24 upregulated, 3 downregulated · 247 samples"),
       x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2",
       caption = "Ellipses: 95% confidence regions") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray30"))

ggsave("HLA_Analysis/Plots/TSNE_HLA_final.png", p_tsne_final, width = 8, height = 6, dpi = 300)

###########################################################################################################################
###########################################################################################################################

# Check the format of row names in exprs_matrix
cat("Row names in exprs_matrix (first 5):\n")
print(head(rownames(exprs_matrix), 5))

cat("\nRow names in exprs_matrix (examples containing HLA):\n")
print(head(rownames(exprs_matrix)[grep("HLA", rownames(exprs_matrix))], 10))

# Check the format of your probe IDs from best_per_gene
cat("\nProbe IDs from best_per_gene (first 5):\n")
print(head(best_per_gene$ProbeID, 10))

# The issue: best_per_gene$ProbeID contains "HLA-DOB_205671_s_at" 
# But rownames(exprs_matrix) are just "205671_s_at" or "HLA-DOB_205671_s_at"? 
# Let's find out

# Check if probe IDs are just the numeric part
test_probe <- "205671_s_at"
if(test_probe %in% rownames(exprs_matrix)) {
  cat("\n✓ Probe IDs are stored as numeric IDs (e.g., 205671_s_at)\n")
} else {
  cat("\n✗ Probe IDs are not numeric IDs\n")
}

# Try to match by extracting the numeric part
# Extract the probe ID from the formatted name (remove gene prefix)
clean_probe_ids <- gsub("^.*_", "", best_per_gene$ProbeID)
cat("\nCleaned probe IDs (first 5):\n")
print(head(clean_probe_ids, 10))

# Check if cleaned IDs match
match_count <- sum(clean_probe_ids %in% rownames(exprs_matrix))
cat("\nCleaned probe IDs matching rownames:", match_count, "/", length(clean_probe_ids), "\n")

###########################################################################################################################

# ============================================================================
# FIGURE 3: CIITA Correlation with HLA Class II Genes
# Using Your Correct Table 3 Data
# ============================================================================

library(reshape2)
library(RColorBrewer)

# Your correct correlation data from Table 3
table3_data <- data.frame(
  Gene = c("HLA-DRA", "HLA-DPA1", "HLA-DRB1", "HLA-DPB1", "HLA-DQA1", 
           "HLA-DMA", "HLA-DMB", "HLA-DOB", "HLA-DQB1"),
  Correlation = c(0.72, 0.70, 0.68, 0.66, 0.65, 0.58, 0.55, 0.52, 0.48),
  CI_lower = c(0.66, 0.64, 0.62, 0.60, 0.58, 0.51, 0.48, 0.44, 0.40),
  CI_upper = c(0.77, 0.75, 0.73, 0.72, 0.71, 0.64, 0.62, 0.59, 0.55),
  P_value = rep("< 0.001", 9)
)

# Add category
table3_data$Category <- ifelse(table3_data$Gene %in% c("HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1"),
                               "Classical (DR/DQ/DP)", "Non-classical (DM/DO)")

# Order by correlation (descending)
table3_data <- table3_data[order(-table3_data$Correlation), ]
table3_data$Gene <- factor(table3_data$Gene, levels = table3_data$Gene)

# ============================================================================
# HEATMAP
# ============================================================================

# Create matrix for heatmap (1 row, multiple columns)
heatmap_matrix <- matrix(table3_data$Correlation, 
                         nrow = 1, 
                         ncol = nrow(table3_data),
                         dimnames = list("CIITA", as.character(table3_data$Gene)))

# Melt for ggplot
heatmap_melt <- melt(heatmap_matrix, varnames = c("Regulator", "Gene"), value.name = "Correlation")
heatmap_melt$Category <- table3_data$Category[match(heatmap_melt$Gene, table3_data$Gene)]

# Create heatmap
p_heatmap <- ggplot(heatmap_melt, aes(x = Gene, y = Regulator, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_gradient2(low = "#0571b0", mid = "white", high = "#ca0020",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Pearson\nCorrelation (r)") +
  geom_text(aes(label = sprintf("%.2f", Correlation)), 
            color = "black", size = 5, fontface = "bold") +
  labs(title = "CIITA Correlation with HLA Class II Gene Expression",
       subtitle = paste0("CIITA is the master regulator of MHC class II transcription\n",
                         "All correlations significant at p < 0.001"),
       x = "HLA Class II Gene",
       y = "Regulator") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "gray30", size = 10),
        legend.position = "right",
        legend.title = element_text(size = 10),
        panel.grid = element_blank()) +
  coord_fixed(ratio = 1)

# Save
ggsave("HLA_Analysis/Plots/Figure3_Heatmap.png", p_heatmap, 
       width = 12, height = 6, dpi = 300)
ggsave("HLA_Analysis/Plots/Figure3_Heatmap.pdf", p_heatmap, 
       width = 12, height = 6, dpi = 300)

print(p_heatmap)
cat("\n✓ Figure 3 (Heatmap) saved\n")


##########################################################################################################################
##########################################################################################################################

# ============================================================================
# Chronic vs Aggressive Periodontitis Comparison
# Only within DISEASED samples
# ============================================================================

# Create subtype assignment (based on original study: patients 1-63 = chronic, 64-90 = aggressive)
pheno_valid$Subtype <- NA
pheno_valid$Subtype[as.numeric(pheno_valid$PatientID) <= 63] <- "Chronic"
pheno_valid$Subtype[as.numeric(pheno_valid$PatientID) >= 64] <- "Aggressive"

# Verify distribution
cat("Subtype distribution in all samples:\n")
print(table(pheno_valid$Subtype, pheno_valid$Status))

# Filter to ONLY DISEASED samples
diseased_idx <- pheno_valid$Status == "Diseased" & !is.na(pheno_valid$Subtype)
pheno_diseased <- pheno_valid[diseased_idx, ]
hla_exprs_diseased <- hla_exprs_valid[, diseased_idx]

cat("\nDiseased samples by subtype:\n")
print(table(pheno_diseased$Subtype))

# Now compare chronic vs aggressive within diseased samples
library(limma)

# Create design matrix with patient blocking
subtype_factor <- factor(pheno_diseased$Subtype, levels = c("Chronic", "Aggressive"))
patient_factor <- factor(pheno_diseased$PatientID)

# Estimate correlation within patients
dupcor_subtype <- duplicateCorrelation(hla_exprs_diseased, block = patient_factor)
cat("\nCorrelation within patients:", round(dupcor_subtype$consensus.correlation, 3), "\n")

# Fit model
design_subtype <- model.matrix(~ subtype_factor)
fit_subtype <- lmFit(hla_exprs_diseased, design_subtype, 
                     block = patient_factor, 
                     correlation = dupcor_subtype$consensus.correlation)
fit_subtype <- eBayes(fit_subtype)

# Get results (Aggressive vs Chronic)
subtype_results <- topTable(fit_subtype, coef = "subtype_factorAggressive", 
                            number = Inf, adjust.method = "fdr")

# Add gene symbols
subtype_results$ProbeID <- rownames(subtype_results)
subtype_results$Gene <- gsub("_.*", "", subtype_results$ProbeID)
subtype_results <- subtype_results[, c("Gene", "ProbeID", "logFC", "P.Value", "adj.P.Val")]

# Sort by adjusted p-value
subtype_results <- subtype_results[order(subtype_results$adj.P.Val), ]

# Identify significant genes
sig_subtype <- subtype_results[subtype_results$adj.P.Val < 0.05, ]
cat("\nSignificant genes between chronic and aggressive (FDR < 0.05):", nrow(sig_subtype), "\n")

if(nrow(sig_subtype) == 0) {
  cat("✓ No significant differences between chronic and aggressive periodontitis\n")
} else {
  cat("⚠️ Significant differences found:\n")
  print(sig_subtype[, c("Gene", "logFC", "adj.P.Val")])
}

# Save results
write.csv(subtype_results, "HLA_Analysis/Results/Chronic_vs_Aggressive_DiseasedOnly.csv", row.names = FALSE)

########################################################################################################################
########################################################################################################################

# VALIDATION ANALYSIS
# ANALYSIS OF GSE16134 DATASET AND COMPARE IT WITH THE DISCOVERY DATASET
# Create the project structure

setwd("D:/MHPE7/periodontitis_validated")

# Create analysis folders
dir.create("HLA_Analysis", showWarnings = FALSE)
dir.create("HLA_Analysis/Results", showWarnings = FALSE)
dir.create("HLA_Analysis/Plots", showWarnings = FALSE)
dir.create("HLA_Analysis/Checkpoints", showWarnings = FALSE)

# Set analysis path
analysis_path <- file.path(getwd(), "HLA_Analysis")

###################################################################################################################

# ============================================================================
# LOAD REQUIRED PACKAGES
# ============================================================================

# Install if needed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages_needed <- c(
  # Bioconductor
  "oligo", "limma", "GEOquery",
  "hgu133plus2.db", "org.Hs.eg.db", "AnnotationDbi",
  # CRAN
  "tidyverse", "pheatmap", "RColorBrewer",
  "ggplot2", "reshape2", "corrplot"
)

for(pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE)) {
    if(pkg %in% c("oligo", "limma", "GEOquery", "hgu133plus2.db", "org.Hs.eg.db", "AnnotationDbi")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

cat("✓ All packages loaded\n")

##########################################################################################################################

# ============================================================================
# LOAD AND PROCESS GPL570 ANNOTATION
# ============================================================================

cat("\n========================================\n")
cat("Loading GPL570 annotation\n")
cat("========================================\n")

# Check if annotation file exists
if(file.exists("GPL570.soft.gz")) {
  cat("Found GPL570.soft.gz\n")
  
  # Read the SOFT file
  gpl <- getGEO(filename = "GPL570.soft.gz")
  
  # Extract annotation table
  gpl_table <- Table(gpl)
  cat("Annotation table dimensions:", dim(gpl_table), "\n")
  cat("Column names:", paste(colnames(gpl_table)[1:10], collapse=", "), "...\n")
  
  # Save for later use
  saveRDS(gpl_table, "HLA_Analysis/Checkpoints/gpl570_annotation.rds")
  cat("✓ Saved GPL570 annotation\n")
  
} else {
  cat("GPL570.soft.gz not found. Downloading from GEO...\n")
  gpl <- getGEO("GPL570", destdir = ".")
  gpl_table <- Table(gpl)
  saveRDS(gpl_table, "HLA_Analysis/Checkpoints/gpl570_annotation.rds")
}

##########################################################################################################################

# ============================================================================
# COMPILE COMPREHENSIVE HLA GENE LIST
# ============================================================================

cat("\n========================================\n")
cat("Creating HLA gene list\n")
cat("========================================\n")

# Comprehensive HLA and related genes
hla_genes <- c(
  # HLA Class I classical
  "HLA-A", "HLA-B", "HLA-C",
  # HLA Class I non-classical
  "HLA-E", "HLA-F", "HLA-G",
  # HLA Class II classical
  "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
  "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
  # HLA Class II non-classical (antigen loading)
  "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
  # Antigen processing machinery
  "TAP1", "TAP2", "TAPBP", "B2M", "CD74",
  # Immunoproteasome subunits
  "PSMB8", "PSMB9", "PSMB10",
  # MHC class II transcriptional regulators
  "CIITA", "RFX5", "RFXAP", "RFXANK", "NFYA", "NFYB", "NFYC"
)

cat("Total HLA-related genes:", length(hla_genes), "\n")
print(hla_genes)

# Save gene list
saveRDS(hla_genes, "HLA_Analysis/Checkpoints/hla_gene_list.rds")
cat("✓ Saved HLA gene list\n")

#########################################################################################################################

# ============================================================================
# MAP HLA GENES TO AFFYMETRIX PROBE IDs
# ============================================================================

cat("\n========================================\n")
cat("Mapping HLA genes to probe IDs\n")
cat("========================================\n")

# Method 1: Using annotation package (preferred) #######################################################################
library(hgu133plus2.db)

# Get all probe IDs for HLA genes
hla_probe_map <- mapIds(hgu133plus2.db,
                        keys = hla_genes,
                        column = "PROBEID",
                        keytype = "SYMBOL",
                        multiVals = "list")

# Convert to data frame for easy viewing
hla_probes_df <- data.frame(
  Gene = rep(names(hla_probe_map), lengths(hla_probe_map)),
  ProbeID = unlist(hla_probe_map),
  row.names = NULL
)

cat("\nProbe mapping summary:\n")
print(table(hla_probes_df$Gene))

# Get unique probe IDs
hla_probe_ids <- unique(unlist(hla_probe_map))
cat("\nTotal unique probe IDs:", length(hla_probe_ids), "\n")

# Method 2: Using GPL570 annotation (backup) #############################################################################
# Get gene symbols from annotation
gpl_table <- readRDS("HLA_Analysis/Checkpoints/gpl570_annotation.rds")

# Find which columns contain gene symbols
gene_cols <- grep("Gene.Symbol|GeneSymbol|gene_symbol|Gene Symbol", 
                  colnames(gpl_table), ignore.case = TRUE, value = TRUE)
cat("\nGene symbol columns in GPL570:", paste(gene_cols, collapse=", "), "\n")

# Use the first gene symbol column
gene_col <- gene_cols[1]

# Extract probe-gene mappings
probe_gene_map <- data.frame(
  ProbeID = gpl_table$ID,
  GeneSymbol = gpl_table[[gene_col]],
  stringsAsFactors = FALSE
)

# Clean up multiple gene symbols (some probes map to multiple genes)
probe_gene_map$GeneSymbol <- gsub(" /// .*", "", probe_gene_map$GeneSymbol)

# Filter for HLA genes
hla_from_gpl <- probe_gene_map[probe_gene_map$GeneSymbol %in% hla_genes, ]

cat("\nHLA probes found in GPL570:", nrow(hla_from_gpl), "\n")

# Save both mappings
saveRDS(hla_probes_df, "HLA_Analysis/Checkpoints/hla_probe_mapping.rds")
saveRDS(hla_probe_ids, "HLA_Analysis/Checkpoints/hla_probe_ids.rds")
saveRDS(hla_from_gpl, "HLA_Analysis/Checkpoints/hla_from_gpl.rds")

cat("\n✓ Saved HLA probe mappings\n")

##########################################################################################################################
##########################################################################################################################
# ============================================================================
# SEARCH FOR ALL RDS FILES
# ============================================================================

# Search current directory and subdirectories
all_rds <- list.files(pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
cat("All RDS files found:\n")
for(f in all_rds) {
  cat("  -", f, "\n")
}

# Check size of each file
file_info <- data.frame(
  File = all_rds,
  Size_MB = file.size(all_rds) / 1024^2
)
print(file_info)

# Look specifically for phenotype files
pheno_files <- all_rds[grepl("pheno", all_rds, ignore.case = TRUE)]
cat("\nPhenotype files found:\n")
print(pheno_files)

###########################################################################################################################
# ============================================================================
# Extract Phenotype for GSE16134 with Status Classification
# ============================================================================

if(file.exists("GSE16134_series_matrix.txt.gz")) {
  cat("Found GSE16134_series_matrix.txt.gz\n")
  
  # Read the compressed file
  con <- gzfile("GSE16134_series_matrix.txt.gz", "rt")
  lines <- readLines(con)
  close(con)
  
  # Find sample IDs
  data_start <- which(lines == "!series_matrix_table_begin") + 1
  header_line <- lines[data_start]
  sample_ids <- strsplit(header_line, "\t")[[1]][-1]
  cat("Found", length(sample_ids), "sample IDs\n")
  
  # Find sample titles (tab-separated, all in one line)
  title_line <- lines[grep("^!Sample_title", lines)]
  all_titles <- strsplit(title_line, "\t")[[1]][-1]
  cat("Found", length(all_titles), "sample titles\n")
  
  # Create phenotype data frame
  pheno_validation <- data.frame(
    SampleID = sample_ids,
    Title = all_titles,
    Cohort = "GSE16134",
    stringsAsFactors = FALSE
  )
  
  # ============================================================
  # CRITICAL FIX: Check for "Unaffected" FIRST, then "Affected"
  # ============================================================
  pheno_validation$Status <- ifelse(
    grepl("Unaffected site|Unaffected|Healthy", pheno_validation$Title, ignore.case = TRUE), 
    "Healthy",
    ifelse(grepl("Affected site|Affected", pheno_validation$Title, ignore.case = TRUE), 
           "Diseased", NA)
  )
  
  # Extract patient ID
  pheno_validation$PatientID <- sapply(strsplit(pheno_validation$Title, " "), function(words) {
    patient_pos <- which(words == "patient")
    if(length(patient_pos) > 0 && patient_pos < length(words)) {
      return(gsub(",", "", words[patient_pos + 1]))
    }
    return(NA)
  })
  
  # Extract sample number
  pheno_validation$SampleNum <- sapply(strsplit(pheno_validation$Title, " "), function(words) {
    sample_pos <- which(words == "sample")
    if(length(sample_pos) > 0 && sample_pos < length(words)) {
      return(gsub(",", "", words[sample_pos + 1]))
    }
    return(NA)
  })
  
  # Verify the fix - show first 10 samples
  cat("\n=== VERIFICATION: First 10 samples ===\n")
  print(head(pheno_validation[, c("SampleID", "Title", "Status", "PatientID")], 10))
  
  # Summary
  cat("\n=== GSE16134 PHENOTYPE SUMMARY (CORRECTED) ===\n")
  cat("Total samples:", nrow(pheno_validation), "\n")
  cat("Diseased:", sum(pheno_validation$Status == "Diseased", na.rm = TRUE), "\n")
  cat("Healthy:", sum(pheno_validation$Status == "Healthy", na.rm = TRUE), "\n")
  cat("Unknown:", sum(is.na(pheno_validation$Status)), "\n")
  cat("Unique patients:", length(unique(na.omit(pheno_validation$PatientID))), "\n")
  
  # Save corrected phenotype (overwrite the incorrect one)
  saveRDS(pheno_validation, "HLA_Analysis/Checkpoints/pheno_validation_corrected.rds")
  write.csv(pheno_validation, "HLA_Analysis/Checkpoints/pheno_validation_corrected.csv", row.names = FALSE)
  cat("\n✓ phenotype saved\n")
  
} else {
  cat("GSE16134_series_matrix.txt.gz not found\n")
}

###########################################################################################################################
###########################################################################################################################
# ============================================================================
# LOAD GSE10334 DISCOVERY DATA FROM YOUR OTHER DIRECTORY
# ============================================================================

# Option A: If you know the path to your GSE10334 data
# Replace with your actual path
gse10334_path <- "D:/MHPE7/periodontitis"  # Update this to your GSE10334 directory

# Check if the discovery expression matrix exists
if(file.exists(file.path(gse10334_path, "HLA_Analysis/Checkpoints/exprs_matrix_all.rds"))) {
  exprs_discovery <- readRDS(file.path(gse10334_path, "HLA_Analysis/Checkpoints/exprs_matrix_all.rds"))
  cat("Loaded discovery expression matrix:", dim(exprs_discovery), "\n")
} else {
  cat("Discovery expression matrix not found at that path\n")
  cat("Please provide the correct path\n")
}

# Load discovery phenotype
if(file.exists(file.path(gse10334_path, "HLA_Analysis/Checkpoints/pheno_data_corrected.rds"))) {
  pheno_discovery <- readRDS(file.path(gse10334_path, "HLA_Analysis/Checkpoints/pheno_data_corrected.rds"))
  cat("Loaded discovery phenotype:", nrow(pheno_discovery), "samples\n")
} else {
  cat("Discovery phenotype not found\n")
}

###########################################################################################################################
# ============================================================================
# LOAD DISCOVERY RESULTS (GSE10334)
# ============================================================================

# Set path for discovery data
discovery_path <- "D:/MHPE7/periodontitis"

# Load discovery results
if(file.exists(file.path(discovery_path, "HLA_Analysis/Results/HLA_best_per_gene_final.csv"))) {
  discovery_results <- read.csv(file.path(discovery_path, "HLA_Analysis/Results/HLA_best_per_gene_final.csv"))
  cat("Loaded discovery results:", nrow(discovery_results), "genes\n")
} else if(file.exists(file.path(discovery_path, "HLA_Analysis/Results/HLA_best_per_gene_corrected.csv"))) {
  discovery_results <- read.csv(file.path(discovery_path, "HLA_Analysis/Results/HLA_best_per_gene_corrected.csv"))
  cat("Loaded discovery results:", nrow(discovery_results), "genes\n")
} else {
  cat("Looking for discovery results in alternative location...\n")
  # Check other possible locations
  alt_files <- list.files(file.path(discovery_path, "HLA_Analysis/Results"), pattern = "best_per_gene", full.names = TRUE)
  print(alt_files)
}

# If still not found, use the discovery_results from your current environment
if(!exists("discovery_results") && exists("best_per_gene")) {
  discovery_results <- best_per_gene
  cat("Using best_per_gene from environment:", nrow(discovery_results), "genes\n")
}

# View first few rows
cat("\nFirst few discovery results:\n")
print(head(discovery_results[, c("Gene", "logFC", "adj.P.Val")]))

############################################################################################################################
############################################################################################################################
# ============================================================================
# RUN VALIDATION ANALYSIS ON GSE16134
# ============================================================================

# Set working directory to validation path
setwd("D:/MHPE7/periodontitis_validated")

# Load validation expression matrix
if(file.exists("HLA_Analysis/Checkpoints/exprs_matrix_all.rds")) {
  exprs_validation <- readRDS("HLA_Analysis/Checkpoints/exprs_matrix_all.rds")
  cat("Loaded validation expression matrix:", dim(exprs_validation), "\n")
} else {
  cat("ERROR: exprs_matrix_all.rds not found\n")
  cat("Available files in Checkpoints:\n")
  print(list.files("HLA_Analysis/Checkpoints/", pattern = "\\.rds$"))
}

# Load validation phenotype (corrected)
if(file.exists("HLA_Analysis/Checkpoints/pheno_validation_corrected.rds")) {
  pheno_validation <- readRDS("HLA_Analysis/Checkpoints/pheno_validation_corrected.rds")
  cat("Loaded validation phenotype:", nrow(pheno_validation), "samples\n")
} else {
  cat("ERROR: pheno_validation_corrected.rds not found\n")
}

# Load HLA probe mapping
if(file.exists("HLA_Analysis/Checkpoints/hla_probe_ids.rds")) {
  hla_probe_ids <- readRDS("HLA_Analysis/Checkpoints/hla_probe_ids.rds")
  cat("Loaded HLA probe IDs:", length(hla_probe_ids), "\n")
} else {
  cat("ERROR: hla_probe_ids.rds not found\n")
}

if(file.exists("HLA_Analysis/Checkpoints/hla_probe_mapping.rds")) {
  hla_probes_df <- readRDS("HLA_Analysis/Checkpoints/hla_probe_mapping.rds")
  cat("Loaded HLA probe mapping:", nrow(hla_probes_df), "probes\n")
} else {
  cat("ERROR: hla_probe_mapping.rds not found\n")
}

###########################################################################################################################
# ============================================================================
# EXTRACT PHENOTYPE FROM GSE16134 SERIES MATRIX
# ============================================================================

# Read the series matrix file
cat("Reading GSE16134_series_matrix.txt.gz...\n")
con <- gzfile("GSE16134_series_matrix.txt.gz", "rt")
lines <- readLines(con)
close(con)

# Get sample IDs from data table header
data_start <- which(lines == "!series_matrix_table_begin") + 1
header_line <- lines[data_start]
sample_ids <- strsplit(header_line, "\t")[[1]][-1]  # Remove "ID_REF"
cat("Sample IDs from header:", length(sample_ids), "\n")

# Clean sample IDs (remove .CEL suffix if present)
sample_ids <- gsub("\\.CEL$", "", sample_ids)
sample_ids <- gsub('"', '', sample_ids)
sample_ids <- trimws(sample_ids)

# Get sample titles
title_line <- lines[grep("^!Sample_title", lines)]
all_titles <- strsplit(title_line, "\t")[[1]][-1]
cat("Sample titles:", length(all_titles), "\n")

# Clean titles (remove quotes)
all_titles <- gsub('"', '', all_titles)

# Create NEW phenotype data frame
pheno_new <- data.frame(
  SampleID = sample_ids,
  Title = all_titles,
  stringsAsFactors = FALSE
)

# Extract status (CRITICAL: Check "Unaffected" FIRST)
pheno_new$Status <- ifelse(
  grepl("Unaffected site|Unaffected|Healthy", pheno_new$Title, ignore.case = TRUE), 
  "Healthy",
  ifelse(grepl("Affected site|Affected", pheno_new$Title, ignore.case = TRUE), 
         "Diseased", NA)
)

# Extract patient ID
pheno_new$PatientID <- sapply(strsplit(pheno_new$Title, " "), function(words) {
  patient_pos <- which(words == "patient")
  if(length(patient_pos) > 0 && patient_pos < length(words)) {
    return(gsub(",", "", words[patient_pos + 1]))
  }
  return(NA)
})

# Extract sample number
pheno_new$SampleNum <- sapply(strsplit(pheno_new$Title, " "), function(words) {
  sample_pos <- which(words == "sample")
  if(length(sample_pos) > 0 && sample_pos < length(words)) {
    return(gsub(",", "", words[sample_pos + 1]))
  }
  return(NA)
})

pheno_new$Cohort <- "GSE16134"

# Verify
cat("\n=== VERIFICATION ===\n")
cat("SampleIDs (first 5):", paste(head(pheno_new$SampleID, 5), collapse = ", "), "\n")
cat("Status distribution:\n")
print(table(pheno_new$Status, useNA = "always"))

# Save this clean phenotype
saveRDS(pheno_new, "HLA_Analysis/Checkpoints/pheno_validation_clean.rds")
write.csv(pheno_new, "HLA_Analysis/Checkpoints/pheno_validation_clean.csv", row.names = FALSE)
cat("\n✓ Clean phenotype saved\n")

###########################################################################################################################
# ============================================================================
# CLEAN EXPRESSION MATRIX COLUMN NAMES
# ============================================================================

# Clean column names of expression matrix
colnames(exprs_validation) <- gsub("\\.CEL$", "", colnames(exprs_validation))
colnames(exprs_validation) <- gsub('"', '', colnames(exprs_validation))
colnames(exprs_validation) <- trimws(colnames(exprs_validation))

cat("Cleaned column names (first 5):\n")
print(head(colnames(exprs_validation), 10))

###########################################################################################################################
# ============================================================================
# ALIGN PHENOTYPE AND EXPRESSION
# ============================================================================

# Reorder phenotype to match expression matrix
pheno_aligned <- pheno_new[match(colnames(exprs_validation), pheno_new$SampleID), ]

# Check alignment
cat("\nAlignment check:", all(colnames(exprs_validation) == pheno_aligned$SampleID), "\n")

# Filter to valid samples (with known status)
valid_idx <- !is.na(pheno_aligned$Status)
exprs_val_valid <- exprs_validation[, valid_idx]
pheno_val_valid <- pheno_aligned[valid_idx, ]

cat("\n=== VALIDATION SAMPLE SUMMARY ===\n")
cat("Total validation samples:", ncol(exprs_val_valid), "\n")
cat("Diseased:", sum(pheno_val_valid$Status == "Diseased"), "\n")
cat("Healthy:", sum(pheno_val_valid$Status == "Healthy"), "\n")
cat("Unique patients:", length(unique(pheno_val_valid$PatientID)), "\n")

# Save aligned data for future use
saveRDS(exprs_val_valid, "HLA_Analysis/Checkpoints/exprs_validation_aligned.rds")
saveRDS(pheno_val_valid, "HLA_Analysis/Checkpoints/pheno_validation_aligned.rds")

############################################################################################################################
###########################################################################################################################
# ============================================================================
# CHECK VALIDATION RESULTS
# ============================================================================

# Read the saved validation results
results_val <- read.csv("HLA_Analysis/Results/Validation_GSE16134_results.csv")

# Check structure
cat("=== VALIDATION RESULTS STRUCTURE ===\n")
cat("Dimensions:", dim(results_val), "\n")
cat("Column names:", paste(colnames(results_val), collapse = ", "), "\n")

# Check first few rows
cat("\nFirst 10 rows:\n")
print(head(results_val, 10))

# Check for NA genes
cat("\nNumber of rows with NA genes:", sum(is.na(results_val$Gene)), "\n")

# Check significance
cat("\nNumber of significant genes (adj.P.Val < 0.05):", 
    sum(results_val$adj.P.Val < 0.05, na.rm = TRUE), "\n")

# Display top significant genes
if(sum(results_val$adj.P.Val < 0.05) > 0) {
  top_sig <- results_val[order(results_val$adj.P.Val), ]
  cat("\nTop 10 significant genes:\n")
  print(head(top_sig[, c("Gene", "logFC", "adj.P.Val")], 10))
} else {
  cat("\nNo significant genes found. Check results.\n")
}

###########################################################################################################################
# ============================================================================
# RE-RUN VALIDATION WITH CORRECT GENE ASSIGNMENT
# ============================================================================

# Load the aligned validation data
exprs_val_valid <- readRDS("HLA_Analysis/Checkpoints/exprs_validation_aligned.rds")
pheno_val_valid <- readRDS("HLA_Analysis/Checkpoints/pheno_validation_aligned.rds")

# Load HLA probe mapping
hla_probes_df <- readRDS("HLA_Analysis/Checkpoints/hla_probe_mapping.rds")
hla_probe_ids <- readRDS("HLA_Analysis/Checkpoints/hla_probe_ids.rds")

# Extract HLA probes present
probes_present <- intersect(hla_probe_ids, rownames(exprs_val_valid))
cat("HLA probes present:", length(probes_present), "\n")

# Subset expression matrix
hla_exprs_val <- exprs_val_valid[probes_present, ]

# Create a clean mapping for these probes
probe_mapping <- data.frame(
  ProbeID = probes_present,
  Gene = hla_probes_df$Gene[match(probes_present, hla_probes_df$ProbeID)],
  stringsAsFactors = FALSE
)

# Remove any probes without gene assignment
probe_mapping <- probe_mapping[!is.na(probe_mapping$Gene), ]
hla_exprs_val <- hla_exprs_val[probe_mapping$ProbeID, ]

cat("After filtering:", nrow(hla_exprs_val), "probes with gene assignments\n")

# Run differential expression
library(limma)

patient_factor <- factor(pheno_val_valid$PatientID)
status_factor <- factor(pheno_val_valid$Status, levels = c("Healthy", "Diseased"))

# Estimate correlation
dupcor <- duplicateCorrelation(hla_exprs_val, block = patient_factor)
cat("Within-patient correlation:", round(dupcor$consensus.correlation, 3), "\n")

# Fit model
design <- model.matrix(~ status_factor)
fit <- lmFit(hla_exprs_val, design, 
             block = patient_factor, 
             correlation = dupcor$consensus.correlation)
fit <- eBayes(fit)

# Get results
results_val <- topTable(fit, coef = "status_factorDiseased", 
                        number = Inf, adjust.method = "fdr")

# Add gene symbols using the mapping
results_val$ProbeID <- rownames(results_val)
results_val$Gene <- probe_mapping$Gene[match(results_val$ProbeID, probe_mapping$ProbeID)]
results_val <- results_val[, c("Gene", "ProbeID", "logFC", "AveExpr", "P.Value", "adj.P.Val")]

# Remove rows with NA genes
results_val <- results_val[!is.na(results_val$Gene), ]

# Save results
write.csv(results_val, "HLA_Analysis/Results/Validation_GSE16134_results.csv", row.names = FALSE)
cat("\n✓ Validation results saved with correct gene assignments\n")

# Display top results
cat("\n=== TOP 10 VALIDATION RESULTS ===\n")
top_val <- head(results_val[order(results_val$adj.P.Val), c("Gene", "logFC", "adj.P.Val")], 10)
print(top_val)

###########################################################################################################################
# ============================================================================
# COMPARE DISCOVERY AND VALIDATION
# ============================================================================

# Set path for discovery results
discovery_path <- "D:/MHPE7/periodontitis"

# Load discovery results
if(file.exists(file.path(discovery_path, "HLA_Analysis/Results/HLA_best_per_gene_final.csv"))) {
  discovery_results <- read.csv(file.path(discovery_path, "HLA_Analysis/Results/HLA_best_per_gene_final.csv"))
} else {
  discovery_results <- best_per_gene
}
cat("Discovery results:", nrow(discovery_results), "genes\n")

# Load validation results (just saved)
results_val <- read.csv("HLA_Analysis/Results/Validation_GSE16134_results.csv")
cat("Validation results:", nrow(results_val), "genes\n")

# Merge
comparison <- merge(discovery_results[, c("Gene", "logFC", "adj.P.Val")],
                    results_val[, c("Gene", "logFC", "adj.P.Val")],
                    by = "Gene", suffixes = c("_disc", "_val"))

cat("\n=== VALIDATION SUMMARY ===\n")
cat("Genes compared:", nrow(comparison), "\n")

# Correlation
cor_val <- cor(comparison$logFC_disc, comparison$logFC_val, method = "pearson", use = "complete.obs")
cat("Pearson correlation:", round(cor_val, 3), "\n")

# Concordance
comparison$same_direction <- sign(comparison$logFC_disc) == sign(comparison$logFC_val)
concordance <- mean(comparison$same_direction, na.rm = TRUE) * 100
cat("Direction concordance:", round(concordance, 1), "%\n")

# Replication rate
comparison$sig_both <- comparison$adj.P.Val_disc < 0.05 & comparison$adj.P.Val_val < 0.05
replication_rate <- mean(comparison$sig_both[comparison$adj.P.Val_disc < 0.05], na.rm = TRUE) * 100
cat("Replication rate:", round(replication_rate, 1), "%\n")

# Key genes
key_genes <- c("HLA-DOB", "CD74", "TAPBP", "HLA-DMA", "CIITA", "HLA-G", "NFYB", "NFYC")
key_comparison <- comparison[comparison$Gene %in% key_genes, ]
cat("\n=== KEY GENES ===\n")
print(key_comparison[, c("Gene", "logFC_disc", "logFC_val", "adj.P.Val_disc", "adj.P.Val_val")])

# Create validation plot
library(ggplot2)

p_val <- ggplot(comparison, aes(x = logFC_disc, y = logFC_val, color = sig_both)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                     labels = c("TRUE" = "Significant in both", "FALSE" = "Not significant")) +
  labs(title = "Validation of HLA Expression Changes",
       subtitle = paste0("Pearson r = ", round(cor_val, 3), 
                         " · Concordance = ", round(concordance, 1), "%"),
       x = "log2 FC (Discovery: GSE10334)",
       y = "log2 FC (Validation: GSE16134)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("HLA_Analysis/Plots/Validation_Scatter.png", p_val, width = 8, height = 7, dpi = 300)
print(p_val)

# Save comparison
write.csv(comparison, "HLA_Analysis/Results/Discovery_vs_Validation.csv", row.names = FALSE)

###########################################################################################################################
# ============================================================================
# CREATE CLEAN VALIDATION TABLE FOR MANUSCRIPT
# ============================================================================

# Select the most significant probe per gene for validation
library(dplyr)

validation_best <- results_val %>%
  group_by(Gene) %>%
  slice_min(order_by = adj.P.Val, n = 1) %>%
  ungroup()

# Merge with discovery best
discovery_best <- discovery_results %>%
  group_by(Gene) %>%
  slice_min(order_by = adj.P.Val, n = 1) %>%
  ungroup()

# Create comparison table
validation_table <- merge(
  discovery_best[, c("Gene", "logFC", "adj.P.Val")],
  validation_best[, c("Gene", "logFC", "adj.P.Val")],
  by = "Gene", suffixes = c("_disc", "_val")
)

# Calculate replication status
validation_table$Replicated <- validation_table$adj.P.Val_disc < 0.05 & 
  validation_table$adj.P.Val_val < 0.05 &
  sign(validation_table$logFC_disc) == sign(validation_table$logFC_val)

# Sort by discovery logFC
validation_table <- validation_table[order(-abs(validation_table$logFC_disc)), ]

# Display top replicated genes
cat("\n=== TOP REPLICATED GENES (Best Probe Per Gene) ===\n")
top_replicated <- head(validation_table[validation_table$Replicated, 
                                        c("Gene", "logFC_disc", "logFC_val", "adj.P.Val_disc", "adj.P.Val_val")], 15)
print(top_replicated)

# Save
write.csv(validation_table, "HLA_Analysis/Results/Validation_Summary_Table.csv", row.names = FALSE)
