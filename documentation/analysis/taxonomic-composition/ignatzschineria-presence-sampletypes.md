```R
# ==================================================================================
# Relative Abundance of Ignatzschineria in Swab Samples by Estimated Diversity
# ==================================================================================
library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggthemes)
library(reshape2)
library(viridis)
library(tibble)
library(patchwork)  # For combining plots

# ==================================================================================
# Load custom ggplot2 theme
# ==================================================================================
theme_pub <- function(base_size = 14, base_family = "CMU Sans Serif", legend_pos = "bottom") {
  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.border = element_blank(),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_blank(),
      legend.position = legend_pos,
      legend.direction = "horizontal",
      legend.key.size = grid::unit(0.4, "cm"),
      legend.title = element_text(face = "bold"),
      plot.margin = grid::unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(fill = "#f0f0f0", colour = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
} 

# ==================================================================================
# Import data from QIIME2 artifacts
# ==================================================================================
ps <- qiime2R::qza_to_phyloseq(
  features = "merged-table-filtered-all.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

# ==================================================================================
# Filter to include only swab samples
# ==================================================================================
ps_swabs <- subset_samples(ps, SampleType == "Swab")

# Print filtered phyloseq object summary
cat("\nFiltered phyloseq object (Swab samples only):\n")
print(ps_swabs)

# Check if EstimatedDiversity_y exists in the metadata
sample_data_df <- as.data.frame(sample_data(ps_swabs))
if (!"EstimatedDiversity_y" %in% colnames(sample_data_df)) {
  stop("EstimatedDiversity_y column not found in metadata. Please check the column name.")
}

# Check if IntendedDiversity exists in the metadata
if (!"IntendedDiversity" %in% colnames(sample_data_df)) {
  cat("\nWARNING: IntendedDiversity column not found in metadata. Will only use EstimatedDiversity_y.\n")
  has_intended_diversity <- FALSE
} else {
  has_intended_diversity <- TRUE
  cat("\nIntendedDiversity column found in metadata.\n")
  cat("Unique values for IntendedDiversity:", paste(unique(sample_data_df$IntendedDiversity), collapse = ", "), "\n")
  cat("Number of samples with NA IntendedDiversity:", sum(is.na(sample_data_df$IntendedDiversity)), "\n")
}

# Check unique values for EstimatedDiversity_y
cat("\nUnique values for EstimatedDiversity_y:", paste(unique(sample_data_df$EstimatedDiversity_y), collapse = ", "), "\n")

# Check for NA values in EstimatedDiversity_y
cat("\nNumber of samples with NA EstimatedDiversity_y:", sum(is.na(sample_data_df$EstimatedDiversity_y)), "\n")

# Remove samples with NA EstimatedDiversity_y if any
ps_swabs <- subset_samples(ps_swabs, !is.na(EstimatedDiversity_y))

# Convert EstimatedDiversity_y to factor for grouping
sample_data(ps_swabs)$EstimatedDiversity_y <- as.factor(sample_data(ps_swabs)$EstimatedDiversity_y)

# If IntendedDiversity exists, convert it to factor as well
if (has_intended_diversity) {
  sample_data(ps_swabs)$IntendedDiversity <- as.factor(sample_data(ps_swabs)$IntendedDiversity)
}

# Check which diversity levels have samples
diversity_counts <- table(sample_data(ps_swabs)$EstimatedDiversity_y)
cat("\nNumber of samples per EstimatedDiversity_y level:\n")
print(diversity_counts)

# Check if we have at least one sample for each diversity level
if (any(diversity_counts == 0)) {
  cat("\nWARNING: Some EstimatedDiversity_y levels have no samples. This may cause plotting issues.\n")
}

# If IntendedDiversity exists, check which levels have samples
if (has_intended_diversity) {
  intended_diversity_counts <- table(sample_data(ps_swabs)$IntendedDiversity)
  cat("\nNumber of samples per IntendedDiversity level:\n")
  print(intended_diversity_counts)
  
  # Check if we have at least one sample for each intended diversity level
  if (any(intended_diversity_counts == 0)) {
    cat("\nWARNING: Some IntendedDiversity levels have no samples. This may cause plotting issues.\n")
  }
}

# Remove taxa that are not present in these samples (have a sum of 0)
ps_swabs <- prune_taxa(taxa_sums(ps_swabs) > 0, ps_swabs)

# ==================================================================================
# Transform to relative abundance
# ==================================================================================
ps_rel <- transform_sample_counts(ps_swabs, function(x) x / sum(x))

# Verify relative abundance transformation
sample_sums <- sample_sums(ps_rel)
cat("\nSample sums after relative abundance transformation:\n")
cat("Min:", min(sample_sums), "Max:", max(sample_sums), "Mean:", mean(sample_sums), "\n")
cat("All sample sums should be approximately 1.0\n")

# ==================================================================================
# Extract Ignatzschineria abundance
# ==================================================================================
extract_ignatzschineria <- function(ps_obj) {
  # Extract OTU table and taxonomy table
  otu_table <- as.data.frame(otu_table(ps_obj))
  tax_table <- as.data.frame(tax_table(ps_obj))
  sample_data <- as.data.frame(sample_data(ps_obj))
  
  # Ensure OTUs are rows
  if (!taxa_are_rows(ps_obj)) {
    otu_table <- t(otu_table)
  }
  
  # Print dimensions of tables for debugging
  cat("\nDimensions of OTU table:", dim(otu_table)[1], "rows,", dim(otu_table)[2], "columns\n")
  cat("Dimensions of taxonomy table:", dim(tax_table)[1], "rows,", dim(tax_table)[2], "columns\n")
  cat("Dimensions of sample data:", dim(sample_data)[1], "rows,", dim(sample_data)[2], "columns\n")
  
  # Check if Genus column exists in taxonomy table
  if (!"Genus" %in% colnames(tax_table)) {
    stop("ERROR: Genus column not found in taxonomy table. Cannot extract Ignatzschineria abundance.")
  }
  
  # Print unique genera in taxonomy table
  unique_genera <- unique(tax_table$Genus)
  cat("Number of unique genera:", length(unique_genera), "\n")
  cat("First 10 genera:", paste(head(unique_genera, 10), collapse=", "), "\n")
  
  # Check if Ignatzschineria exists in taxonomy table
  if (!"Ignatzschineria" %in% unique_genera) {
    cat("\nWARNING: Ignatzschineria not found in taxonomy table.\n")
    # Return a data frame with zeros
    result <- data.frame(
      SampleID = rownames(sample_data),
      Ignatzschineria = rep(0, nrow(sample_data)),
      Other = rep(1, nrow(sample_data))
    )
    return(result)
  }
  
  # Filter taxonomy table for Ignatzschineria
  ignatzschineria_taxa <- rownames(tax_table)[tax_table$Genus == "Ignatzschineria"]
  cat("Number of ASVs for Ignatzschineria:", length(ignatzschineria_taxa), "\n")
  
  # Print the first few Ignatzschineria ASVs
  if (length(ignatzschineria_taxa) > 0) {
    cat("First few Ignatzschineria ASVs:", paste(head(ignatzschineria_taxa), collapse=", "), "\n")
    
    # Check if these ASVs exist in the OTU table
    existing_taxa <- ignatzschineria_taxa[ignatzschineria_taxa %in% rownames(otu_table)]
    cat("Number of Ignatzschineria ASVs found in OTU table:", length(existing_taxa), "\n")
    
    if (length(existing_taxa) == 0) {
      cat("\nWARNING: None of the Ignatzschineria ASVs found in taxonomy table exist in the OTU table.\n")
      cat("This suggests a mismatch between taxonomy and OTU tables.\n")
    }
  }
  
  # Create a data frame to store the results
  result <- data.frame(
    SampleID = colnames(otu_table),
    Ignatzschineria = rep(0, ncol(otu_table)),
    Other = rep(0, ncol(otu_table))
  )
  
  # Calculate Ignatzschineria abundance
  if (length(ignatzschineria_taxa) > 0) {
    # Check if the ASVs exist in the OTU table
    existing_taxa <- ignatzschineria_taxa[ignatzschineria_taxa %in% rownames(otu_table)]
    
    if (length(existing_taxa) > 0) {
      if (length(existing_taxa) == 1) {
        # If only one taxon, get its abundance
        result$Ignatzschineria <- as.numeric(otu_table[existing_taxa, ])
        cat("Using single ASV for Ignatzschineria abundance\n")
      } else {
        # If multiple taxa, sum their abundances
        result$Ignatzschineria <- colSums(otu_table[existing_taxa, , drop = FALSE])
        cat("Summing", length(existing_taxa), "ASVs for Ignatzschineria abundance\n")
      }
      
      # Print summary statistics of Ignatzschineria abundance
      cat("Ignatzschineria abundance summary:\n")
      cat("Min:", min(result$Ignatzschineria), "Max:", max(result$Ignatzschineria), 
          "Mean:", mean(result$Ignatzschineria), "\n")
      cat("Number of samples with Ignatzschineria:", sum(result$Ignatzschineria > 0), 
          "out of", ncol(otu_table), "\n")
    } else {
      cat("\nWARNING: No matching Ignatzschineria ASVs found in OTU table. Using zeros.\n")
    }
  }
  
  # Calculate abundance for all other taxa
  result$Other <- 1 - result$Ignatzschineria
  
  # Add metadata
  result <- merge(result, sample_data, by.x = "SampleID", by.y = "row.names")
  
  # Print the first few rows of the result
  cat("\nFirst few rows of result data frame:\n")
  print(head(result[, c("SampleID", "Ignatzschineria", "Other", "EstimatedDiversity_y")]))
  
  return(result)
}

# Extract Ignatzschineria abundance
plot_data <- extract_ignatzschineria(ps_rel)

# Print debug information
cat("\nSummary of Ignatzschineria abundance in swab samples:\n")

# Check if we have any non-NA values
if (all(is.na(plot_data$Ignatzschineria))) {
  cat("All Ignatzschineria values are NA\n")
} else {
  # Only calculate stats if we have non-NA values
  cat("Min:", min(plot_data$Ignatzschineria, na.rm = TRUE), 
      "Max:", max(plot_data$Ignatzschineria, na.rm = TRUE), 
      "Mean:", mean(plot_data$Ignatzschineria, na.rm = TRUE), "\n")
}

# Check if we have any Ignatzschineria
if (all(is.na(plot_data$Ignatzschineria)) || all(plot_data$Ignatzschineria == 0, na.rm = TRUE)) {
  cat("\nWARNING: No Ignatzschineria found in any swab samples. Creating placeholder data for visualization.\n")
  
  # Create placeholder data with small values for visualization
  plot_data$Ignatzschineria <- 0.001  # Small placeholder value
  plot_data$Other <- 0.999
}

# Print the first few rows of the data
cat("\nFirst few rows of plot_data:\n")
print(head(plot_data))

# Convert to long format for plotting
plot_data_long <- plot_data %>%
  pivot_longer(
    cols = c("Ignatzschineria", "Other"),
    names_to = "Genus",
    values_to = "Abundance"
  )

# Print debug information for long format data
cat("\nSummary of long format data:\n")
cat("Number of rows:", nrow(plot_data_long), "\n")
cat("Number of unique EstimatedDiversity_y values:", length(unique(plot_data_long$EstimatedDiversity_y)), "\n")
cat("Number of unique Genus values:", length(unique(plot_data_long$Genus)), "\n")

# ==================================================================================
# Calculate mean abundance by diversity level
# ==================================================================================
bar_data <- plot_data_long %>%
  group_by(EstimatedDiversity_y, Genus) %>%
  summarise(
    MeanAbundance = mean(Abundance, na.rm = TRUE),
    SampleCount = n(),
    SE = sd(Abundance, na.rm = TRUE) / sqrt(n()),  # Calculate standard error
    .groups = "drop"
  )

# Print debug information
cat("\nSummary of grouped data before normalization:\n")
print(bar_data)

# Check for NaN values in MeanAbundance
if (any(is.nan(bar_data$MeanAbundance))) {
  cat("\nWARNING: NaN values found in MeanAbundance. Replacing with zeros.\n")
  bar_data$MeanAbundance[is.nan(bar_data$MeanAbundance)] <- 0
}

# Now normalize
bar_data <- bar_data %>%
  group_by(EstimatedDiversity_y) %>%
  mutate(
    GroupSum = sum(MeanAbundance, na.rm = TRUE),
    NormalizedAbundance = ifelse(GroupSum > 0, MeanAbundance / GroupSum, 0.5)  # Use 0.5 as default for empty groups
  ) %>%
  ungroup()

# Print debug information after normalization
cat("\nSummary of grouped data after normalization:\n")
print(bar_data)

# Check for NA or NaN values
if (any(is.na(bar_data$NormalizedAbundance)) || any(is.nan(bar_data$NormalizedAbundance))) {
  cat("\nWARNING: NA or NaN values found in normalized abundance. This may cause plotting issues.\n")
  # Replace NA/NaN with 0
  bar_data$NormalizedAbundance[is.na(bar_data$NormalizedAbundance) | is.nan(bar_data$NormalizedAbundance)] <- 0
}

# Verify that each group now sums to 1.0
verification <- bar_data %>%
  group_by(EstimatedDiversity_y) %>%
  summarise(GroupSum = sum(NormalizedAbundance), .groups = "drop")

cat("Group sums after normalization (should all be 1.0):\n")
cat("Min:", min(verification$GroupSum), "Max:", max(verification$GroupSum), "\n")

# ==================================================================================
# Create vertical bar plot
# ==================================================================================
create_vertical_barplot <- function(data) {
  # Set colors for Ignatzschineria and Other
  colors <- c("Ignatzschineria" = "#0ea5e9", "Other" = "#e5e5e5")
  
  # Add sample count information
  sample_counts <- data %>%
    group_by(EstimatedDiversity_y) %>%
    summarise(SampleCount = first(SampleCount), .groups = "drop") %>%
    mutate(Label = paste0("n=", SampleCount))
  
  # Create bar plot
  p <- ggplot(data, aes(x = EstimatedDiversity_y, y = NormalizedAbundance, fill = Genus)) +
    geom_col(position = "stack", color = "white", size = 0.1) +
    geom_text(data = sample_counts, aes(label = Label, y = -0.05), 
              position = position_dodge(width = 0.9), 
              vjust = 1, size = 3) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0.1), oob = scales::squish) +
    labs(
      title = "Relative Abundance of Ignatzschineria in Swab Samples by Estimated Diversity",
      subtitle = "Showing proportion of total bacterial abundance",
      y = "Mean Relative Abundance",
      x = "Estimated Diversity"
    ) +
    theme_pub(legend_pos = "right") +
    theme(
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      legend.direction = "vertical"
    )
  
  return(p)
}

# ==================================================================================
# Create horizontal bar plot
# ==================================================================================
create_horizontal_barplot <- function(data) {
  # Set colors for Ignatzschineria and Other
  colors <- c("Ignatzschineria" = "#0ea5e9", "Other" = "#e5e5e5")
  
  # Add sample count information
  sample_counts <- data %>%
    group_by(EstimatedDiversity_y) %>%
    summarise(SampleCount = first(SampleCount), .groups = "drop") %>%
    mutate(Label = paste0("n=", SampleCount))
  
  # Create horizontal bar plot
  p <- ggplot(data, aes(x = EstimatedDiversity_y, y = NormalizedAbundance, fill = Genus)) +
    geom_col(position = "stack", color = "white", size = 0.1) +
    geom_text(data = sample_counts, aes(label = Label, y = -0.05), 
              position = position_dodge(width = 0.9), 
              hjust = 1, size = 3) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0.1), oob = scales::squish) +
    coord_flip() +  # Flip coordinates to make horizontal bars
    labs(
      title = "Relative Abundance of Ignatzschineria in Swab Samples by Estimated Diversity",
      subtitle = "Showing proportion of total bacterial abundance",
      y = "Mean Relative Abundance",
      x = "Estimated Diversity"
    ) +
    theme_pub(legend_pos = "bottom") +
    theme(
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      legend.direction = "horizontal"
    )
  
  return(p)
}

# ==================================================================================
# Create Ignatzschineria-only bar plot by EstimatedDiversity_y
# ==================================================================================
create_ignatzschineria_only_barplot <- function(data) {
  # Filter for Ignatzschineria only
  ignatz_data <- data %>%
    filter(Genus == "Ignatzschineria")
  
  # Add sample count as a label
  ignatz_data <- ignatz_data %>%
    mutate(Label = paste0("n=", SampleCount))
  
  # Create bar plot
  p <- ggplot(ignatz_data, aes(x = EstimatedDiversity_y, y = NormalizedAbundance)) +
    geom_col(fill = "#808080", color = "white", size = 0.1) +
    # Add error bars
    geom_errorbar(
      aes(ymin = NormalizedAbundance - SE, 
          ymax = NormalizedAbundance + SE),
      width = 0.2,
      color = "black",
      size = 0.5
    ) +
    geom_text(aes(label = Label, y = NormalizedAbundance + SE + 0.02), 
              position = position_dodge(width = 0.9), 
              vjust = 0, size = 3) +
    scale_y_continuous(
      labels = scales::percent_format(), 
      limits = c(0, 0.4),  # Set y-axis range to 0-40%
      expand = c(0, 0), 
      oob = scales::squish
    ) +
    labs(
      # Removed title and subtitle as requested
      y = "Mean Relative Abundance",
      x = "Estimated Diversity"
    ) +
    theme_pub() +
    theme(legend.position = "none")
  
  return(p)
}

# ==================================================================================
# Create Ignatzschineria-only bar plot by IntendedDiversity
# ==================================================================================
create_ignatzschineria_intended_barplot <- function(data) {
  # Filter for Ignatzschineria only
  ignatz_data <- data %>%
    filter(Genus == "Ignatzschineria")
  
  # Add sample count as a label
  ignatz_data <- ignatz_data %>%
    mutate(Label = paste0("n=", SampleCount))
  
  # Create bar plot
  p <- ggplot(ignatz_data, aes(x = IntendedDiversity, y = NormalizedAbundance)) +
    geom_col(fill = "#808080", color = "white", size = 0.1) +
    # Add error bars
    geom_errorbar(
      aes(ymin = NormalizedAbundance - SE, 
          ymax = NormalizedAbundance + SE),
      width = 0.2,
      color = "black",
      size = 0.5
    ) +
    geom_text(aes(label = Label, y = NormalizedAbundance + SE + 0.02), 
              position = position_dodge(width = 0.9), 
              vjust = 0, size = 3) +
    scale_y_continuous(
      labels = scales::percent_format(), 
      limits = c(0, 0.4),  # Set y-axis range to 0-40%
      expand = c(0, 0), 
      oob = scales::squish
    ) +
    labs(
      # Removed title and subtitle as requested
      y = "Mean Relative Abundance",
      x = "Intended Diversity"
    ) +
    theme_pub() +
    theme(legend.position = "none")
  
  return(p)
}

# ==================================================================================
# Process data for IntendedDiversity if available
# ==================================================================================
if (has_intended_diversity) {
  # Create a new long format data frame for IntendedDiversity
  plot_data_intended_long <- plot_data %>%
    pivot_longer(
      cols = c("Ignatzschineria", "Other"),
      names_to = "Genus",
      values_to = "Abundance"
    )
  
  # Calculate mean abundance by intended diversity level
  bar_data_intended <- plot_data_intended_long %>%
    group_by(IntendedDiversity, Genus) %>%
    summarise(
      MeanAbundance = mean(Abundance, na.rm = TRUE),
      SampleCount = n(),
      SE = sd(Abundance, na.rm = TRUE) / sqrt(n()),  # Calculate standard error
      .groups = "drop"
    )
  
  # Check for NaN values in MeanAbundance
  if (any(is.nan(bar_data_intended$MeanAbundance))) {
    cat("\nWARNING: NaN values found in MeanAbundance for IntendedDiversity. Replacing with zeros.\n")
    bar_data_intended$MeanAbundance[is.nan(bar_data_intended$MeanAbundance)] <- 0
  }
  
  # Normalize
  bar_data_intended <- bar_data_intended %>%
    group_by(IntendedDiversity) %>%
    mutate(
      GroupSum = sum(MeanAbundance, na.rm = TRUE),
      NormalizedAbundance = ifelse(GroupSum > 0, MeanAbundance / GroupSum, 0.5)  # Use 0.5 as default for empty groups
    ) %>%
    ungroup()
  
  # Check for NA or NaN values
  if (any(is.na(bar_data_intended$NormalizedAbundance)) || any(is.nan(bar_data_intended$NormalizedAbundance))) {
    cat("\nWARNING: NA or NaN values found in normalized abundance for IntendedDiversity. This may cause plotting issues.\n")
    # Replace NA/NaN with 0
    bar_data_intended$NormalizedAbundance[is.na(bar_data_intended$NormalizedAbundance) | is.nan(bar_data_intended$NormalizedAbundance)] <- 0
  }
  
  # Verify that each group now sums to 1.0
  verification_intended <- bar_data_intended %>%
    group_by(IntendedDiversity) %>%
    summarise(GroupSum = sum(NormalizedAbundance), .groups = "drop")
  
  cat("Group sums after normalization for IntendedDiversity (should all be 1.0):\n")
  cat("Min:", min(verification_intended$GroupSum), "Max:", max(verification_intended$GroupSum), "\n")
}

# ==================================================================================
# Function to combine Ignatzschineria plots
# ==================================================================================
combine_ignatzschineria_plots <- function(estimated_plot, intended_plot) {
  # Use patchwork to combine plots side by side without a title
  combined_plot <- estimated_plot + intended_plot + 
    plot_layout(ncol = 2, widths = c(1, 1))
  
  return(combined_plot)
}

# ==================================================================================
# Generate and save visualizations
# ==================================================================================

# Create vertical bar plot
vertical_plot <- create_vertical_barplot(bar_data)

# Create horizontal bar plot
horizontal_plot <- create_horizontal_barplot(bar_data)

# Create Ignatzschineria-only bar plot by EstimatedDiversity_y
# First check if we have any non-zero Ignatzschineria data
ignatz_data <- bar_data %>% filter(Genus == "Ignatzschineria")
if (all(ignatz_data$NormalizedAbundance == 0)) {
  cat("\nWARNING: All Ignatzschineria abundance values are zero. Creating placeholder visualization.\n")
  # Create a placeholder plot with a message
  ignatz_only_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = "No Ignatzschineria detected in swab samples", 
             size = 6) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.margin = margin(50, 50, 50, 50)
    )
} else {
  ignatz_only_plot <- create_ignatzschineria_only_barplot(bar_data)
}

# Variable to store combined plot
combined_ignatz_plot <- NULL

# Create Ignatzschineria-only bar plot by IntendedDiversity if available
if (has_intended_diversity) {
  # Check if we have any non-zero Ignatzschineria data for IntendedDiversity
  ignatz_intended_data <- bar_data_intended %>% filter(Genus == "Ignatzschineria")
  if (all(ignatz_intended_data$NormalizedAbundance == 0)) {
    cat("\nWARNING: All Ignatzschineria abundance values are zero for IntendedDiversity. Creating placeholder visualization.\n")
    # Create a placeholder plot with a message
    ignatz_intended_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = "No Ignatzschineria detected in swab samples", 
               size = 6) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(50, 50, 50, 50)
      )
  } else {
    ignatz_intended_plot <- create_ignatzschineria_intended_barplot(bar_data_intended)
  }
  
  # Combine the two Ignatzschineria plots
  combined_ignatz_plot <- combine_ignatzschineria_plots(ignatz_only_plot, ignatz_intended_plot)
  
  # Save the combined plot
  ggsave("swab_ignatzschineria_combined_plots.png", combined_ignatz_plot, width = 12, height = 8, dpi = 300)
  
  # Display the combined plot
  print(combined_ignatz_plot)
}

# Save individual plots
 ggsave("swab_ignatzschineria_by_diversity_vertical.png", vertical_plot, width = 10, height = 8, dpi = 300)
ggsave("swab_ignatzschineria_by_diversity_horizontal.png", horizontal_plot, width = 10, height = 6, dpi = 300)
ggsave("swab_ignatzschineria_only_by_diversity.png", ignatz_only_plot, width = 8, height = 6, dpi = 300)
if (has_intended_diversity) {
  ggsave("swab_ignatzschineria_by_intended_diversity.png", ignatz_intended_plot, width = 8, height = 6, dpi = 300)
}

# Display plots
print(vertical_plot)
print(horizontal_plot)
print(ignatz_only_plot)

# Print summary message
cat("\nRelative abundance visualizations created:\n")
cat("1. Vertical bar plot (swab_ignatzschineria_by_diversity_vertical.png)\n")
cat("2. Horizontal bar plot (swab_ignatzschineria_by_diversity_horizontal.png)\n")
cat("3. Ignatzschineria-only plot by EstimatedDiversity_y (swab_ignatzschineria_only_by_diversity.png)\n")
if (has_intended_diversity) {
  cat("4. Ignatzschineria-only plot by IntendedDiversity (swab_ignatzschineria_by_intended_diversity.png)\n")
  cat("5. Combined Ignatzschineria plots (swab_ignatzschineria_combined_plots.png)\n")
}

```

```R
# ==================================================================================
# Relative Abundance of Ignatzschineria in Soil Samples by Estimated Diversity
# ==================================================================================
library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggthemes)
library(reshape2)
library(viridis)
library(tibble)
library(patchwork)  # For combining plots

# ==================================================================================
# Load custom ggplot2 theme
# ==================================================================================
theme_pub <- function(base_size = 14, base_family = "CMU Sans Serif", legend_pos = "bottom") {
  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.border = element_blank(),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_blank(),
      legend.position = legend_pos,
      legend.direction = "horizontal",
      legend.key.size = grid::unit(0.4, "cm"),
      legend.title = element_text(face = "bold"),
      plot.margin = grid::unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(fill = "#f0f0f0", colour = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
} 

# ==================================================================================
# Import data from QIIME2 artifacts
# ==================================================================================
ps <- qiime2R::qza_to_phyloseq(
  features = "merged-table-filtered-all.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

# ==================================================================================
# Filter to include only soil samples
# ==================================================================================
ps_soil <- subset_samples(ps, SampleType == "Soil")

# Print filtered phyloseq object summary
cat("\nFiltered phyloseq object (Soil samples only):\n")
print(ps_soil)

# Check if EstimatedDiversity_y exists in the metadata
sample_data_df <- as.data.frame(sample_data(ps_soil))
if (!"EstimatedDiversity_y" %in% colnames(sample_data_df)) {
  stop("EstimatedDiversity_y column not found in metadata. Please check the column name.")
}

# Check if IntendedDiversity exists in the metadata
if (!"IntendedDiversity" %in% colnames(sample_data_df)) {
  cat("\nWARNING: IntendedDiversity column not found in metadata. Will only use EstimatedDiversity_y.\n")
  has_intended_diversity <- FALSE
} else {
  has_intended_diversity <- TRUE
  cat("\nIntendedDiversity column found in metadata.\n")
  cat("Unique values for IntendedDiversity:", paste(unique(sample_data_df$IntendedDiversity), collapse = ", "), "\n")
  cat("Number of samples with NA IntendedDiversity:", sum(is.na(sample_data_df$IntendedDiversity)), "\n")
}

# Check unique values for EstimatedDiversity_y
cat("\nUnique values for EstimatedDiversity_y:", paste(unique(sample_data_df$EstimatedDiversity_y), collapse = ", "), "\n")

# Check for NA values in EstimatedDiversity_y
cat("\nNumber of samples with NA EstimatedDiversity_y:", sum(is.na(sample_data_df$EstimatedDiversity_y)), "\n")

# Remove samples with NA EstimatedDiversity_y if any
ps_soil <- subset_samples(ps_soil, !is.na(EstimatedDiversity_y))

# Convert EstimatedDiversity_y to factor for grouping
sample_data(ps_soil)$EstimatedDiversity_y <- as.factor(sample_data(ps_soil)$EstimatedDiversity_y)

# If IntendedDiversity exists, convert it to factor as well
if (has_intended_diversity) {
  sample_data(ps_soil)$IntendedDiversity <- as.factor(sample_data(ps_soil)$IntendedDiversity)
}

# Check which diversity levels have samples
diversity_counts <- table(sample_data(ps_soil)$EstimatedDiversity_y)
cat("\nNumber of samples per EstimatedDiversity_y level:\n")
print(diversity_counts)

# Check if we have at least one sample for each diversity level
if (any(diversity_counts == 0)) {
  cat("\nWARNING: Some EstimatedDiversity_y levels have no samples. This may cause plotting issues.\n")
}

# If IntendedDiversity exists, check which levels have samples
if (has_intended_diversity) {
  intended_diversity_counts <- table(sample_data(ps_soil)$IntendedDiversity)
  cat("\nNumber of samples per IntendedDiversity level:\n")
  print(intended_diversity_counts)
  
  # Check if we have at least one sample for each intended diversity level
  if (any(intended_diversity_counts == 0)) {
    cat("\nWARNING: Some IntendedDiversity levels have no samples. This may cause plotting issues.\n")
  }
}

# Remove taxa that are not present in these samples (have a sum of 0)
ps_soil <- prune_taxa(taxa_sums(ps_soil) > 0, ps_soil)

# ==================================================================================
# Transform to relative abundance
# ==================================================================================
ps_rel <- transform_sample_counts(ps_soil, function(x) x / sum(x))

# Verify relative abundance transformation
sample_sums <- sample_sums(ps_rel)
cat("\nSample sums after relative abundance transformation:\n")
cat("Min:", min(sample_sums), "Max:", max(sample_sums), "Mean:", mean(sample_sums), "\n")
cat("All sample sums should be approximately 1.0\n")

# ==================================================================================
# Extract Ignatzschineria abundance
# ==================================================================================
extract_ignatzschineria <- function(ps_obj) {
  # Extract OTU table and taxonomy table
  otu_table <- as.data.frame(otu_table(ps_obj))
  tax_table <- as.data.frame(tax_table(ps_obj))
  sample_data <- as.data.frame(sample_data(ps_obj))
  
  # Ensure OTUs are rows
  if (!taxa_are_rows(ps_obj)) {
    otu_table <- t(otu_table)
  }
  
  # Print dimensions of tables for debugging
  cat("\nDimensions of OTU table:", dim(otu_table)[1], "rows,", dim(otu_table)[2], "columns\n")
  cat("Dimensions of taxonomy table:", dim(tax_table)[1], "rows,", dim(tax_table)[2], "columns\n")
  cat("Dimensions of sample data:", dim(sample_data)[1], "rows,", dim(sample_data)[2], "columns\n")
  
  # Check if Genus column exists in taxonomy table
  if (!"Genus" %in% colnames(tax_table)) {
    stop("ERROR: Genus column not found in taxonomy table. Cannot extract Ignatzschineria abundance.")
  }
  
  # Print unique genera in taxonomy table
  unique_genera <- unique(tax_table$Genus)
  cat("Number of unique genera:", length(unique_genera), "\n")
  cat("First 10 genera:", paste(head(unique_genera, 10), collapse=", "), "\n")
  
  # Check if Ignatzschineria exists in taxonomy table
  if (!"Ignatzschineria" %in% unique_genera) {
    cat("\nWARNING: Ignatzschineria not found in taxonomy table.\n")
    # Return a data frame with zeros
    result <- data.frame(
      SampleID = rownames(sample_data),
      Ignatzschineria = rep(0, nrow(sample_data)),
      Other = rep(1, nrow(sample_data))
    )
    return(result)
  }
  
  # Filter taxonomy table for Ignatzschineria
  ignatzschineria_taxa <- rownames(tax_table)[tax_table$Genus == "Ignatzschineria"]
  cat("Number of ASVs for Ignatzschineria:", length(ignatzschineria_taxa), "\n")
  
  # Print the first few Ignatzschineria ASVs
  if (length(ignatzschineria_taxa) > 0) {
    cat("First few Ignatzschineria ASVs:", paste(head(ignatzschineria_taxa), collapse=", "), "\n")
    
    # Check if these ASVs exist in the OTU table
    existing_taxa <- ignatzschineria_taxa[ignatzschineria_taxa %in% rownames(otu_table)]
    cat("Number of Ignatzschineria ASVs found in OTU table:", length(existing_taxa), "\n")
    
    if (length(existing_taxa) == 0) {
      cat("\nWARNING: None of the Ignatzschineria ASVs found in taxonomy table exist in the OTU table.\n")
      cat("This suggests a mismatch between taxonomy and OTU tables.\n")
    }
  }
  
  # Create a data frame to store the results
  result <- data.frame(
    SampleID = colnames(otu_table),
    Ignatzschineria = rep(0, ncol(otu_table)),
    Other = rep(0, ncol(otu_table))
  )
  
  # Calculate Ignatzschineria abundance
  if (length(ignatzschineria_taxa) > 0) {
    # Check if the ASVs exist in the OTU table
    existing_taxa <- ignatzschineria_taxa[ignatzschineria_taxa %in% rownames(otu_table)]
    
    if (length(existing_taxa) > 0) {
      if (length(existing_taxa) == 1) {
        # If only one taxon, get its abundance
        result$Ignatzschineria <- as.numeric(otu_table[existing_taxa, ])
        cat("Using single ASV for Ignatzschineria abundance\n")
      } else {
        # If multiple taxa, sum their abundances
        result$Ignatzschineria <- colSums(otu_table[existing_taxa, , drop = FALSE])
        cat("Summing", length(existing_taxa), "ASVs for Ignatzschineria abundance\n")
      }
      
      # Print summary statistics of Ignatzschineria abundance
      cat("Ignatzschineria abundance summary:\n")
      cat("Min:", min(result$Ignatzschineria), "Max:", max(result$Ignatzschineria), 
          "Mean:", mean(result$Ignatzschineria), "\n")
      cat("Number of samples with Ignatzschineria:", sum(result$Ignatzschineria > 0), 
          "out of", ncol(otu_table), "\n")
    } else {
      cat("\nWARNING: No matching Ignatzschineria ASVs found in OTU table. Using zeros.\n")
    }
  }
  
  # Calculate abundance for all other taxa
  result$Other <- 1 - result$Ignatzschineria
  
  # Add metadata
  result <- merge(result, sample_data, by.x = "SampleID", by.y = "row.names")
  
  # Print the first few rows of the result
  cat("\nFirst few rows of result data frame:\n")
  print(head(result[, c("SampleID", "Ignatzschineria", "Other", "EstimatedDiversity_y")]))
  
  return(result)
}

# Extract Ignatzschineria abundance
plot_data <- extract_ignatzschineria(ps_rel)

# Print debug information
cat("\nSummary of Ignatzschineria abundance in soil samples:\n")

# Check if we have any non-NA values
if (all(is.na(plot_data$Ignatzschineria))) {
  cat("All Ignatzschineria values are NA\n")
} else {
  # Only calculate stats if we have non-NA values
  cat("Min:", min(plot_data$Ignatzschineria, na.rm = TRUE), 
      "Max:", max(plot_data$Ignatzschineria, na.rm = TRUE), 
      "Mean:", mean(plot_data$Ignatzschineria, na.rm = TRUE), "\n")
}

# Check if we have any Ignatzschineria
if (all(is.na(plot_data$Ignatzschineria)) || all(plot_data$Ignatzschineria == 0, na.rm = TRUE)) {
  cat("\nWARNING: No Ignatzschineria found in any soil samples. Creating placeholder data for visualization.\n")
  
  # Create placeholder data with small values for visualization
  plot_data$Ignatzschineria <- 0.001  # Small placeholder value
  plot_data$Other <- 0.999
}

# Print the first few rows of the data
cat("\nFirst few rows of plot_data:\n")
print(head(plot_data))

# Convert to long format for plotting
plot_data_long <- plot_data %>%
  pivot_longer(
    cols = c("Ignatzschineria", "Other"),
    names_to = "Genus",
    values_to = "Abundance"
  )

# Print debug information for long format data
cat("\nSummary of long format data:\n")
cat("Number of rows:", nrow(plot_data_long), "\n")
cat("Number of unique EstimatedDiversity_y values:", length(unique(plot_data_long$EstimatedDiversity_y)), "\n")
cat("Number of unique Genus values:", length(unique(plot_data_long$Genus)), "\n")

# ==================================================================================
# Calculate mean abundance by diversity level
# ==================================================================================
bar_data <- plot_data_long %>%
  group_by(EstimatedDiversity_y, Genus) %>%
  summarise(
    MeanAbundance = mean(Abundance, na.rm = TRUE),
    SampleCount = n(),
    SE = sd(Abundance, na.rm = TRUE) / sqrt(n()),  # Calculate standard error
    .groups = "drop"
  )

# Print debug information
cat("\nSummary of grouped data before normalization:\n")
print(bar_data)

# Check for NaN values in MeanAbundance
if (any(is.nan(bar_data$MeanAbundance))) {
  cat("\nWARNING: NaN values found in MeanAbundance. Replacing with zeros.\n")
  bar_data$MeanAbundance[is.nan(bar_data$MeanAbundance)] <- 0
}

# Now normalize
bar_data <- bar_data %>%
  group_by(EstimatedDiversity_y) %>%
  mutate(
    GroupSum = sum(MeanAbundance, na.rm = TRUE),
    NormalizedAbundance = ifelse(GroupSum > 0, MeanAbundance / GroupSum, 0.5)  # Use 0.5 as default for empty groups
  ) %>%
  ungroup()

# Print debug information after normalization
cat("\nSummary of grouped data after normalization:\n")
print(bar_data)

# Check for NA or NaN values
if (any(is.na(bar_data$NormalizedAbundance)) || any(is.nan(bar_data$NormalizedAbundance))) {
  cat("\nWARNING: NA or NaN values found in normalized abundance. This may cause plotting issues.\n")
  # Replace NA/NaN with 0
  bar_data$NormalizedAbundance[is.na(bar_data$NormalizedAbundance) | is.nan(bar_data$NormalizedAbundance)] <- 0
}

# Verify that each group now sums to 1.0
verification <- bar_data %>%
  group_by(EstimatedDiversity_y) %>%
  summarise(GroupSum = sum(NormalizedAbundance), .groups = "drop")

cat("Group sums after normalization (should all be 1.0):\n")
cat("Min:", min(verification$GroupSum), "Max:", max(verification$GroupSum), "\n")

# ==================================================================================
# Create vertical bar plot
# ==================================================================================
create_vertical_barplot <- function(data) {
  # Set colors for Ignatzschineria and Other
  colors <- c("Ignatzschineria" = "#0ea5e9", "Other" = "#e5e5e5")
  
  # Add sample count information
  sample_counts <- data %>%
    group_by(EstimatedDiversity_y) %>%
    summarise(SampleCount = first(SampleCount), .groups = "drop") %>%
    mutate(Label = paste0("n=", SampleCount))
  
  # Create bar plot
  p <- ggplot(data, aes(x = EstimatedDiversity_y, y = NormalizedAbundance, fill = Genus)) +
    geom_col(position = "stack", color = "white", size = 0.1) +
    geom_text(data = sample_counts, aes(label = Label, y = -0.05), 
              position = position_dodge(width = 0.9), 
              vjust = 1, size = 3) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0.1), oob = scales::squish) +
    labs(
      title = "Relative Abundance of Ignatzschineria in Soil Samples by Estimated Diversity",
      subtitle = "Showing proportion of total bacterial abundance",
      y = "Mean Relative Abundance",
      x = "Estimated Diversity"
    ) +
    theme_pub(legend_pos = "right") +
    theme(
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      legend.direction = "vertical"
    )
  
  return(p)
}

# ==================================================================================
# Create horizontal bar plot
# ==================================================================================
create_horizontal_barplot <- function(data) {
  # Set colors for Ignatzschineria and Other
  colors <- c("Ignatzschineria" = "#0ea5e9", "Other" = "#e5e5e5")
  
  # Add sample count information
  sample_counts <- data %>%
    group_by(EstimatedDiversity_y) %>%
    summarise(SampleCount = first(SampleCount), .groups = "drop") %>%
    mutate(Label = paste0("n=", SampleCount))
  
  # Create horizontal bar plot
  p <- ggplot(data, aes(x = EstimatedDiversity_y, y = NormalizedAbundance, fill = Genus)) +
    geom_col(position = "stack", color = "white", size = 0.1) +
    geom_text(data = sample_counts, aes(label = Label, y = -0.05), 
              position = position_dodge(width = 0.9), 
              hjust = 1, size = 3) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0.1), oob = scales::squish) +
    coord_flip() +  # Flip coordinates to make horizontal bars
    labs(
      title = "Relative Abundance of Ignatzschineria in Soil Samples by Estimated Diversity",
      subtitle = "Showing proportion of total bacterial abundance",
      y = "Mean Relative Abundance",
      x = "Estimated Diversity"
    ) +
    theme_pub(legend_pos = "bottom") +
    theme(
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      legend.direction = "horizontal"
    )
  
  return(p)
}

# ==================================================================================
# Create Ignatzschineria-only bar plot by EstimatedDiversity_y
# ==================================================================================
create_ignatzschineria_only_barplot <- function(data) {
  # Filter for Ignatzschineria only
  ignatz_data <- data %>%
    filter(Genus == "Ignatzschineria")
  
  # Add sample count as a label
  ignatz_data <- ignatz_data %>%
    mutate(Label = paste0("n=", SampleCount))
  
  # Create bar plot
  p <- ggplot(ignatz_data, aes(x = EstimatedDiversity_y, y = NormalizedAbundance)) +
    geom_col(fill = "#808080", color = "white", size = 0.1) +
    # Add error bars
    geom_errorbar(
      aes(ymin = NormalizedAbundance - SE, 
          ymax = NormalizedAbundance + SE),
      width = 0.2,
      color = "black",
      size = 0.5
    ) +
    geom_text(aes(label = Label, y = NormalizedAbundance + SE + 0.02), 
              position = position_dodge(width = 0.9), 
              vjust = 0, size = 3) +
    scale_y_continuous(
      labels = scales::percent_format(), 
      limits = c(0, 0.4),  # Set y-axis range to 0-40%
      expand = c(0, 0), 
      oob = scales::squish
    ) +
    labs(
      # Removed title and subtitle as requested
      y = "Mean Relative Abundance",
      x = "Estimated Diversity"
    ) +
    theme_pub() +
    theme(legend.position = "none")
  
  return(p)
}

# ==================================================================================
# Create Ignatzschineria-only bar plot by IntendedDiversity
# ==================================================================================
create_ignatzschineria_intended_barplot <- function(data) {
  # Filter for Ignatzschineria only
  ignatz_data <- data %>%
    filter(Genus == "Ignatzschineria")
  
  # Add sample count as a label
  ignatz_data <- ignatz_data %>%
    mutate(Label = paste0("n=", SampleCount))
  
  # Create bar plot
  p <- ggplot(ignatz_data, aes(x = IntendedDiversity, y = NormalizedAbundance)) +
    geom_col(fill = "#808080", color = "white", size = 0.1) +
    # Add error bars
    geom_errorbar(
      aes(ymin = NormalizedAbundance - SE, 
          ymax = NormalizedAbundance + SE),
      width = 0.2,
      color = "black",
      size = 0.5
    ) +
    geom_text(aes(label = Label, y = NormalizedAbundance + SE + 0.02), 
              position = position_dodge(width = 0.9), 
              vjust = 0, size = 3) +
    scale_y_continuous(
      labels = scales::percent_format(), 
      limits = c(0, 0.4),  # Set y-axis range to 0-40%
      expand = c(0, 0), 
      oob = scales::squish
    ) +
    labs(
      # Removed title and subtitle as requested
      y = "Mean Relative Abundance",
      x = "Intended Diversity"
    ) +
    theme_pub() +
    theme(legend.position = "none")
  
  return(p)
}

# ==================================================================================
# Process data for IntendedDiversity if available
# ==================================================================================
if (has_intended_diversity) {
  # Create a new long format data frame for IntendedDiversity
  plot_data_intended_long <- plot_data %>%
    pivot_longer(
      cols = c("Ignatzschineria", "Other"),
      names_to = "Genus",
      values_to = "Abundance"
    )
  
  # Calculate mean abundance by intended diversity level
  bar_data_intended <- plot_data_intended_long %>%
    group_by(IntendedDiversity, Genus) %>%
    summarise(
      MeanAbundance = mean(Abundance, na.rm = TRUE),
      SampleCount = n(),
      SE = sd(Abundance, na.rm = TRUE) / sqrt(n()),  # Calculate standard error
      .groups = "drop"
    )
  
  # Check for NaN values in MeanAbundance
  if (any(is.nan(bar_data_intended$MeanAbundance))) {
    cat("\nWARNING: NaN values found in MeanAbundance for IntendedDiversity. Replacing with zeros.\n")
    bar_data_intended$MeanAbundance[is.nan(bar_data_intended$MeanAbundance)] <- 0
  }
  
  # Normalize
  bar_data_intended <- bar_data_intended %>%
    group_by(IntendedDiversity) %>%
    mutate(
      GroupSum = sum(MeanAbundance, na.rm = TRUE),
      NormalizedAbundance = ifelse(GroupSum > 0, MeanAbundance / GroupSum, 0.5)  # Use 0.5 as default for empty groups
    ) %>%
    ungroup()
  
  # Check for NA or NaN values
  if (any(is.na(bar_data_intended$NormalizedAbundance)) || any(is.nan(bar_data_intended$NormalizedAbundance))) {
    cat("\nWARNING: NA or NaN values found in normalized abundance for IntendedDiversity. This may cause plotting issues.\n")
    # Replace NA/NaN with 0
    bar_data_intended$NormalizedAbundance[is.na(bar_data_intended$NormalizedAbundance) | is.nan(bar_data_intended$NormalizedAbundance)] <- 0
  }
  
  # Verify that each group now sums to 1.0
  verification_intended <- bar_data_intended %>%
    group_by(IntendedDiversity) %>%
    summarise(GroupSum = sum(NormalizedAbundance), .groups = "drop")
  
  cat("Group sums after normalization for IntendedDiversity (should all be 1.0):\n")
  cat("Min:", min(verification_intended$GroupSum), "Max:", max(verification_intended$GroupSum), "\n")
}

# ==================================================================================
# Function to combine Ignatzschineria plots
# ==================================================================================
combine_ignatzschineria_plots <- function(estimated_plot, intended_plot) {
  # Use patchwork to combine plots side by side without a title
  combined_plot <- estimated_plot + intended_plot + 
    plot_layout(ncol = 2, widths = c(1, 1))
  
  return(combined_plot)
}

# ==================================================================================
# Generate and save visualizations
# ==================================================================================

# Create vertical bar plot
vertical_plot <- create_vertical_barplot(bar_data)

# Create horizontal bar plot
horizontal_plot <- create_horizontal_barplot(bar_data)

# Create Ignatzschineria-only bar plot by EstimatedDiversity_y
# First check if we have any non-zero Ignatzschineria data
ignatz_data <- bar_data %>% filter(Genus == "Ignatzschineria")
if (all(ignatz_data$NormalizedAbundance == 0)) {
  cat("\nWARNING: All Ignatzschineria abundance values are zero. Creating placeholder visualization.\n")
  # Create a placeholder plot with a message
  ignatz_only_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = "No Ignatzschineria detected in soil samples", 
             size = 6) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.margin = margin(50, 50, 50, 50)
    )
} else {
  ignatz_only_plot <- create_ignatzschineria_only_barplot(bar_data)
}

# Variable to store combined plot
combined_ignatz_plot <- NULL

# Create Ignatzschineria-only bar plot by IntendedDiversity if available
if (has_intended_diversity) {
  # Check if we have any non-zero Ignatzschineria data for IntendedDiversity
  ignatz_intended_data <- bar_data_intended %>% filter(Genus == "Ignatzschineria")
  if (all(ignatz_intended_data$NormalizedAbundance == 0)) {
    cat("\nWARNING: All Ignatzschineria abundance values are zero for IntendedDiversity. Creating placeholder visualization.\n")
    # Create a placeholder plot with a message
    ignatz_intended_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = "No Ignatzschineria detected in soil samples", 
               size = 6) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white"),
        plot.margin = margin(50, 50, 50, 50)
      )
  } else {
    ignatz_intended_plot <- create_ignatzschineria_intended_barplot(bar_data_intended)
  }
  
  # Combine the two Ignatzschineria plots
  combined_ignatz_plot <- combine_ignatzschineria_plots(ignatz_only_plot, ignatz_intended_plot)
  
  # Save the combined plot
  ggsave("soil_ignatzschineria_combined_plots.png", combined_ignatz_plot, width = 12, height = 8, dpi = 300)
  
  # Display the combined plot
  print(combined_ignatz_plot)
}

# Save individual plots
ggsave("soil_ignatzschineria_by_diversity_vertical.png", vertical_plot, width = 10, height = 8, dpi = 300)
ggsave("soil_ignatzschineria_by_diversity_horizontal.png", horizontal_plot, width = 10, height = 6, dpi = 300)
ggsave("soil_ignatzschineria_only_by_diversity.png", ignatz_only_plot, width = 8, height = 6, dpi = 300)
if (has_intended_diversity) {
  ggsave("soil_ignatzschineria_by_intended_diversity.png", ignatz_intended_plot, width = 8, height = 6, dpi = 300)
}

# Display plots
print(vertical_plot)
print(horizontal_plot)
print(ignatz_only_plot)

# Print summary message
cat("\nRelative abundance visualizations created:\n")
cat("1. Vertical bar plot (soil_ignatzschineria_by_diversity_vertical.png)\n")
cat("2. Horizontal bar plot (soil_ignatzschineria_by_diversity_horizontal.png)\n")
cat("3. Ignatzschineria-only plot by EstimatedDiversity_y (soil_ignatzschineria_only_by_diversity.png)\n")
if (has_intended_diversity) {
  cat("4. Ignatzschineria-only plot by IntendedDiversity (soil_ignatzschineria_by_intended_diversity.png)\n")
  cat("5. Combined Ignatzschineria plots (soil_ignatzschineria_combined_plots.png)\n")
}

```

