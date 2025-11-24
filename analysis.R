library(tidyverse)
library(pheatmap)
library(ggrepel)

set.seed(123) # Ensures reproducible random imputation

# Load data
raw_data <- read.csv("proteins.csv", stringsAsFactors = FALSE)

# Strip column names so that they are "C3", "D1", et
names(raw_data) <- gsub("LFQ.intensity.", "", names(raw_data))

id_col <- grep("Majority.protein.IDs", names(raw_data), value = TRUE)
gene_col <- grep("Fasta.headers", names(raw_data), value = TRUE)

# Rename column
raw_data <- raw_data %>% rename(ProteinID = all_of(id_col), FastaHeader = all_of(gene_col))

# Remove E11 as outlier
raw_data <- raw_data %>% select(-matches("E11"))

# -----------------------------------------------------------------------------
# Preprocessing!
# -----------------------------------------------------------------------------

# Extract Gene Names for clarity
raw_data$Gene <- str_extract(raw_data$FastaHeader, "GN=[^ ]+")
raw_data$Gene <- gsub("GN=", "", raw_data$Gene)
# If no gene name found, default back to ProteinID
raw_data$Gene[is.na(raw_data$Gene)] <- raw_data$ProteinID[is.na(raw_data$Gene)]

# Identify data columns (excluding ID/Gene columns)
data_cols <- names(raw_data)[!names(raw_data) %in% c("ProteinID", "FastaHeader", "Gene")]

# Define unique groups (C, D, L, E)
unique_groups <- unique(substr(data_cols, 1, 1))

# Use 'apply' to check each row. It preserves column names reliably.
keep_rows <- apply(raw_data[, data_cols], 1, function(row_vals) {
  # Check each group separately
  valid_counts <- sapply(unique_groups, function(g) {
    # Find columns belonging to this group (e.g., start with "C")
    group_cols <- grep(paste0("^", g), names(row_vals))
    # Count how many numbers are NOT NA
    sum(!is.na(row_vals[group_cols]))
  })
  # Keep row if any group has at least 3 valid values
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
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel() +
  theme_minimal() +
  labs(title = "PCA Plot no E11")

# ----------------------------------------------------------------------------
# ANOVA

long_data <- data_imputed %>%
  pivot_longer(cols = all_of(data_cols), names_to = "Sample", values_to = "Intensity") %>%
  mutate(Group = substr(Sample, 1, 1))

# Run ANOVA
anova_results <- long_data %>%
  group_by(ProteinID, Gene) %>%
  summarise(p_val = summary(aov(Intensity ~ Group))[[1]][["Pr(>F)"]][1], .groups = "drop")

# Adjust P-values
anova_results$p_adj <- p.adjust(anova_results$p_val, method = "BH")
sig_proteins <- anova_results %>% filter(p_adj < 0.05)

print(paste("Significant proteins found:", nrow(sig_proteins)))

# Calculate Fold Changes for Volcano Plot
fold_changes <- data_imputed %>%
  rowwise() %>%
  mutate(
    Mean_C = mean(c_across(starts_with("C"))),
    Mean_L = mean(c_across(starts_with("L"))),
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
  
  pheatmap(sig_matrix, scale = "row", 
           show_rownames = TRUE, fontsize_row = 6,
           main = "Heatmap of ANOVA Significant Proteins")
}

# ----------------------------------------------------------------------------
# Volcano Plot (Light vs Control)
final_results <- final_results %>%
  mutate(IsSig = if_else(p_adj < 0.05 & abs(Log2FC_L_vs_C) > 1, "Significant", "Not Sig"))

ggplot(final_results, aes(x = Log2FC_L_vs_C, y = -log10(p_adj), color = IsSig)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(data = subset(final_results, IsSig == "Significant"), 
                  aes(label = Gene), max.overlaps = 10) +
  theme_minimal() +
  labs(title = "Volcano Plot: Light vs Control", x = "Log2 Fold Change (L - C)", y = "-Log10 FDR")

# ---------------------------------------------------------------------
write.csv(final_results, "results.csv", row.names = FALSE)