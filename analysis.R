library(tidyverse)
library(pheatmap)
library(ggrepel)

set.seed(123) # Ensures reproducible random imputation

# Load data
raw_data <- read.csv("proteins.csv", stringsAsFactors = FALSE)

# Strip column names so that they are "C3", "D1", etc.
names(raw_data) <- gsub("LFQ.intensity.", "", names(raw_data))

id_col <- grep("Majority.protein.IDs", names(raw_data), value = TRUE)
gene_col <- grep("Fasta.headers", names(raw_data), value = TRUE)

# Rename columns
raw_data <- raw_data %>% rename(ProteinID = all_of(id_col), FastaHeader = all_of(gene_col))

# Remove E11 and L6 
raw_data <- raw_data %>%
  select(-matches("E11")) %>%
  select(-matches("L6"))

# Identify current data columns
initial_data_cols <- names(raw_data)[!names(raw_data) %in% c("ProteinID", "FastaHeader", gene_col, id_col)]

# Calculate the median intensity for each sample
sample_medians <- raw_data %>%
  select(all_of(initial_data_cols)) %>%
  # Calculate median for each column (sample)
  summarise(across(everything(), ~ median(., na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "MedianIntensity")

# Determine which samples to keep (Top 3 per Group, based on median intensity)
samples_to_keep <- sample_medians %>%
  mutate(Group = substr(Sample, 1, 1)) %>%
  group_by(Group) %>%
  # Arrange by median intensity descending and take the top 3
  arrange(desc(MedianIntensity)) %>%
  slice_head(n = 3) %>%
  pull(Sample)

# Filter the raw_data to keep only the selected samples plus ID columns
id_cols_all <- c("ProteinID", "FastaHeader") # Gene column is created later
raw_data <- raw_data %>% select(all_of(id_cols_all), all_of(samples_to_keep))

# Update the data_cols variable for subsequent steps
data_cols <- samples_to_keep
unique_groups <- unique(substr(data_cols, 1, 1))

print(paste("Selected samples (Top 3 per group):", paste(data_cols, collapse = ", ")))
print(paste("Groups analyzed:", paste(unique_groups, collapse = ", ")))

# -----------------------------------------------------------------------------
# Preprocessing!
# -----------------------------------------------------------------------------

# Extract Gene Names for clarity
raw_data$Gene <- str_extract(raw_data$FastaHeader, "GN=[^ ]+")
raw_data$Gene <- gsub("GN=", "", raw_data$Gene)
# If no gene name found, default back to ProteinID
raw_data$Gene[is.na(raw_data$Gene)] <- raw_data$ProteinID[is.na(raw_data$Gene)]

# Use 'apply' to check each row for quantification completeness
keep_rows <- apply(raw_data[, data_cols], 1, function(row_vals) {
  # Check each group separately
  valid_counts <- sapply(unique_groups, function(g) {
    # Find columns belonging to this group (e.g., start with "C")
    group_cols <- grep(paste0("^", g), names(row_vals))
    # Count how many numbers are NOT NA
    sum(!is.na(row_vals[group_cols]))
  })
  # Keep row if any group has at least 3 valid values (since we only have 3 per group now)
  # This means the protein must be quantified in all 3 samples of at least one group
  any(valid_counts >= 3)
})

# Apply the filter mask
data_filtered <- raw_data[keep_rows, ]

print(paste("Proteins remaining after filter:", nrow(data_filtered)))

# Normalize (Median Centering)
data_norm <- data_filtered
numeric_mat <- as.matrix(data_norm[, data_cols])
global_median <- median(numeric_mat, na.rm = TRUE)

# Loop through columns to center them
for(col in data_cols) {
  col_median <- median(data_norm[[col]], na.rm = TRUE)
  data_norm[[col]] <- data_norm[[col]] - col_median + global_median
}

# Imputation
# Function to impute missing values with "noise"
impute_downshift <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  n_missing <- sum(is.na(x))
  
  if(n_missing > 0) {
    # Shift mean down by 1.8 SDs, width 0.3 SDs
    x[is.na(x)] <- rnorm(n_missing, mean = m - 1.8 * s, sd = 0.3 * s)
  }
  return(x)
}

# Apply imputation
data_imputed <- data_norm %>%
  mutate(across(all_of(data_cols), impute_downshift))

# ----------------------------------------------------------------------------
# PCA Plot!
pca_input <- t(data_imputed[, data_cols])
pca_res <- prcomp(pca_input, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$Group <- substr(rownames(pca_df), 1, 1)
pca_df$Sample <- rownames(pca_df)

# Plot
print(ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
        geom_point(size = 4) +
        geom_text_repel() +
        theme_minimal() +
        labs(title = "PCA Plot (Top 3 Samples Selected)"))

# ----------------------------------------------------------------------------
# ANOVA

long_data <- data_imputed %>%
  pivot_longer(cols = all_of(data_cols), names_to = "Sample", values_to = "Intensity") %>%
  mutate(Group = substr(Sample, 1, 1))

# Run ANOVA
anova_results <- long_data %>%
  group_by(ProteinID, Gene) %>%
  # Ensure there is variation within groups for ANOVA to run
  summarise(p_val = tryCatch({
    summary(aov(Intensity ~ Group, data = cur_data()))[[1]][["Pr(>F)"]][1]
  }, error = function(e) {
    # Return 1 if ANOVA fails (e.g., if no variation within a group)
    return(1)
  }), .groups = "drop")

# Adjust P-values
anova_results$p_adj <- p.adjust(anova_results$p_val, method = "BH")
sig_proteins <- anova_results %>% filter(p_adj < 0.05)

print(paste("Significant proteins found (p.adj < 0.05):", nrow(sig_proteins)))

# Calculate Fold Changes for Volcano Plot (L vs C)
fold_changes <- data_imputed %>%
  rowwise() %>%
  mutate(
    # Use selected data columns for mean calculation
    Mean_C = mean(c_across(all_of(data_cols[substr(data_cols, 1, 1) == "C"]))),
    Mean_L = mean(c_across(all_of(data_cols[substr(data_cols, 1, 1) == "L"]))),
    Log2FC_L_vs_C = Mean_L - Mean_C
  ) %>%
  ungroup() %>%
  select(ProteinID, Gene, Log2FC_L_vs_C)

# Combine results
final_results <- left_join(anova_results, fold_changes, by = c("ProteinID", "Gene"))

# -----------------------------------------------------------------------------
# Heat Map!
if(nrow(sig_proteins) > 0) {
  sig_matrix <- data_imputed %>%
    filter(ProteinID %in% sig_proteins$ProteinID) %>%
    select(all_of(data_cols)) %>%
    as.matrix()
  
  # Add Gene names as row labels
  rownames(sig_matrix) <- data_imputed$Gene[data_imputed$ProteinID %in% sig_proteins$ProteinID]
  
  # Scale by row (Z-score) for visualization
  print(pheatmap(sig_matrix, scale = "row", 
                 show_rownames = TRUE, 
                 fontsize_row = 4,
                 fontsize_col = 8,
                 main = "Heatmap of ANOVA Significant Proteins"))
} else {
  print("No proteins were significant (p.adj < 0.05), skipping heatmap generation.")
}

# ----------------------------------------------------------------------------
# Volcano Plot (Light vs Control)
final_results <- final_results %>%
  mutate(IsSig = case_when(
    p_adj < 0.05 & abs(Log2FC_L_vs_C) > 1 ~ "Significant",
    p_adj < 0.05 & abs(Log2FC_L_vs_C) <= 1 ~ "p < 0.05",
    TRUE ~ "Not Sig"
  ))

print(ggplot(final_results, aes(x = Log2FC_L_vs_C, y = -log10(p_adj), color = IsSig)) +
        geom_point(alpha = 0.6) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        scale_color_manual(values = c("p < 0.05" = "blue", "Significant" = "red", "Not Sig" = "grey")) +
        geom_text_repel(data = subset(final_results, IsSig == "Significant"), 
                        aes(label = Gene), max.overlaps = 15, size = 3) +
        theme_minimal() +
        labs(title = "Volcano Plot: Light vs Control (Top 3 Samples)", 
             x = "Log2 Fold Change (L - C)", 
             y = "-Log10 FDR",
             color = "Significance"))

# ---------------------------------------------------------------------
write.csv(final_results, "results.csv", row.names = FALSE)