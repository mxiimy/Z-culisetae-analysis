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
  # CRITICAL FILTERING STEP: Keep row if any group has at least 3 valid values
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
# ANOVA (For Heatmap Generation)

long_data <- data_imputed %>%
  pivot_longer(cols = all_of(data_cols), names_to = "Sample", values_to = "Intensity") %>%
  mutate(Group = substr(Sample, 1, 1))

# Run ANOVA
anova_results <- long_data %>%
  group_by(ProteinID, Gene) %>%
  summarise(p_val = tryCatch({
    summary(aov(Intensity ~ Group, data = cur_data()))[[1]][["Pr(>F)"]][1]
  }, error = function(e) {
    return(1)
  }), .groups = "drop")

# Adjust P-values
anova_results$p_adj <- p.adjust(anova_results$p_val, method = "BH")
sig_proteins <- anova_results %>% filter(p_adj < 0.05)

print(paste("Proteins with significant overall change (ANOVA p.adj < 0.05):", nrow(sig_proteins)))

# -----------------------------------------------------------------------------
# Heat Map! (Uses overall significant proteins from ANOVA)
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
                 fontsize = 8,  # Reduced base font size to 8 to shrink the title
                 fontsize_row = 4, 
                 fontsize_col = 8, 
                 treeheight_col = 30, # Reduced column dendrogram height
                 main = "Heatmap of ANOVA Significant Proteins"))
} else {
  print("No proteins were significant (p.adj < 0.05), skipping heatmap generation.")
}

# ----------------------------------------------------------------------------
# Pairwise T-tests and Volcano Plots (All Combinations)
# ----------------------------------------------------------------------------

# Calculate Mean Intensity for all groups (needed for Log2FC)
group_means <- data_imputed %>%
  rowwise() %>%
  mutate(
    Mean_C = mean(c_across(all_of(data_cols[substr(data_cols, 1, 1) == "C"]))),
    Mean_D = mean(c_across(all_of(data_cols[substr(data_cols, 1, 1) == "D"]))),
    Mean_L = mean(c_across(all_of(data_cols[substr(data_cols, 1, 1) == "L"]))),
    Mean_E = mean(c_across(all_of(data_cols[substr(data_cols, 1, 1) == "E"])))
  ) %>%
  ungroup() %>%
  select(ProteinID, Gene, starts_with("Mean_"))


# Function to run pairwise T-test and generate volcano plot
generate_volcano <- function(g1, g2, data, group_means, data_cols) {
  
  comparison_name <- paste0(g2, " vs ", g1) # g2 is numerator, g1 is denominator
  
  # Calculate Log2FC
  fc_col_name <- paste0("Log2FC_", g2, "_vs_", g1)
  
  # Select the mean columns dynamically and calculate Log2FC
  fc_data <- group_means %>%
    mutate(!!fc_col_name := .data[[paste0("Mean_", g2)]] - .data[[paste0("Mean_", g1)]]) %>%
    select(ProteinID, Gene, !!fc_col_name)
  
  # Run T-tests for p-value (STATISTICALLY CORRECT APPROACH FOR PAIRWISE)
  cols_g1 <- data_cols[substr(data_cols, 1, 1) == g1]
  cols_g2 <- data_cols[substr(data_cols, 1, 1) == g2]
  
  pairwise_long_data <- data %>%
    select(ProteinID, Gene, all_of(c(cols_g1, cols_g2))) %>%
    pivot_longer(cols = all_of(c(cols_g1, cols_g2)), names_to = "Sample", values_to = "Intensity") %>%
    mutate(Group = substr(Sample, 1, 1))
  
  t_test_results <- pairwise_long_data %>%
    group_by(ProteinID, Gene) %>%
    summarise(
      p_val = tryCatch({
        # Two-sample T-test
        t.test(Intensity ~ Group, data = cur_data())$p.value
      }, error = function(e) {
        return(1) # Return 1 if t-test fails
      }),
      .groups = "drop"
    )
  
  # Combine results and adjust p-values
  final_df <- left_join(t_test_results, fc_data, by = c("ProteinID", "Gene"))
  final_df$p_adj <- p.adjust(final_df$p_val, method = "BH")
  
  log2fc_col <- sym(fc_col_name)
  
  # Define significance (FDR < 0.05 AND |Log2FC| > 1)
  final_df <- final_df %>%
    mutate(IsSig = case_when(
      p_adj < 0.05 & abs(!!log2fc_col) > 1 ~ "Significant",
      p_adj < 0.05 & abs(!!log2fc_col) <= 1 ~ "p < 0.05",
      TRUE ~ "Not Sig"
    ))
  
  # Generate Plot (Using T-test results)
  plot_title <- paste0("Volcano Plot: ", comparison_name, " (Top 3 Samples)")
  
  plot_out <- ggplot(final_df, aes(x = !!log2fc_col, y = -log10(p_adj), color = IsSig)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = c("p < 0.05" = "blue", "Significant" = "red", "Not Sig" = "grey")) +
    geom_text_repel(data = subset(final_df, IsSig %in% c("Significant", "p < 0.05")), 
                    aes(label = Gene), max.overlaps = 20, size = 3,
                    order = -log10(final_df$p_adj)) +
    theme_minimal() +
    labs(title = plot_title, 
         x = paste0("Log2 Fold Change (", g2, " - ", g1, ")"), 
         y = "-Log10 FDR",
         color = "Significance")
  
  print(plot_out)
  
  return(final_df)
}

# Main Plotting Loop
all_groups <- unique(substr(data_cols, 1, 1))
group_combinations <- combn(all_groups, 2) 

all_results_list <- list()

for(i in 1:ncol(group_combinations)) {
  # g1 is the denominator (baseline), g2 is the numerator (comparison)
  g1 <- group_combinations[1, i] 
  g2 <- group_combinations[2, i] 
  
  # Generate Volcano Plot and get results
  results_df <- generate_volcano(g1, g2, data_imputed, group_means, data_cols)
  
  # Store results
  results_df$Comparison <- paste0(g2, "_vs_", g1)
  all_results_list[[i]] <- results_df
}

# Combine all results into one large data frame and save
all_final_results <- bind_rows(all_results_list)

# ---------------------------------------------------------------------
# Saving all results
write.csv(all_final_results, "all_pairwise_results.csv", row.names = FALSE)