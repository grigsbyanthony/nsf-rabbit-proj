## Load packages

```R
library(qiime2R)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(ape)
library(umap)
library(gridExtra)
library(patchwork)  
library(RColorBrewer)
library(pairwiseAdonis)  
```



## Custom theme

```R
theme_pub <- function(base_size=14, base_family="cmu sans serif") {
  library(grid)
  library(ggthemes)
  
  # Create a modified theme without using margin()
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size = unit(0.2, "cm"),
            # Avoid using margin() function
            legend.spacing = unit(0, "cm"),
            legend.box.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin = unit(c(10,5,5,5),"mm"),
            strip.background = element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}
```



## `phyloseq` object filtering

```R
# Import QIIME2 artifacts directly to phyloseq using qiime2R
ps <- qiime2R::qza_to_phyloseq(
  features = "merged-table-filtered-all.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

# Remove Phoridae sp. samples
ps <- subset_samples(ps, Source != "Phoridae sp.")

# Filter to include only Swab, Insect, and Soil samples
ps <- subset_samples(ps, SampleType %in% c("Swab", "Insect", "Soil"))


# Remove taxa that are not present in these samples (have a sum of 0)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Check what sample types are available
cat("\n=== AVAILABLE SAMPLE TYPES ===\n")
available_sample_types <- unique(sample_data(ps)$SampleType)
cat("Available sample types:", paste(available_sample_types, collapse = ", "), "\n")

# Check what microbiota values are available
cat("\n=== AVAILABLE MICROBIOTA VALUES ===\n")
available_microbiota <- unique(sample_data(ps)$Microbiota)
cat("Available microbiota:", paste(available_microbiota, collapse = ", "), "\n")

# Ensure SampleType is a factor (this is important for plotting)
sample_data(ps)$SampleType <- factor(sample_data(ps)$SampleType)
available_sample_types <- levels(sample_data(ps)$SampleType)

# Ensure Microbiota is a factor (this is important for shape aesthetics)
# For samples without Microbiota value, set a default value
sample_data(ps)$Microbiota <- as.character(sample_data(ps)$Microbiota)
sample_data(ps)$Microbiota[is.na(sample_data(ps)$Microbiota)] <- "Unknown"
sample_data(ps)$Microbiota <- factor(sample_data(ps)$Microbiota)
available_microbiota <- levels(sample_data(ps)$Microbiota)

# Ensure Stage is a factor with the desired order (Fr, Ac, Adv, NA)
sample_data(ps)$Stage <- factor(sample_data(ps)$Stage, levels = c("Fr", "Ac", "Adv", NA))

# Define a color palette for SampleType values
# Make sure we have enough colors for all sample types
n_colors <- max(3, length(available_sample_types))  # At least 3 colors
sample_type_colors <- brewer.pal(n_colors, "Set1")[1:length(available_sample_types)]
names(sample_type_colors) <- available_sample_types

# Print sample type and microbiota levels to verify
cat("\nSample type levels (as factors):", paste(levels(sample_data(ps)$SampleType), collapse = ", "), "\n")
cat("\nMicrobiota levels (as factors):", paste(levels(sample_data(ps)$Microbiota), collapse = ", "), "\n")

```



## Permanova function

```R
# ============================================================================
# BETA DIVERSITY CALCULATIONS
# ============================================================================

# Calculate distance matrices
dist_bray <- phyloseq::distance(ps, method = "bray")
dist_unifrac <- phyloseq::distance(ps, method = "unifrac")
dist_wunifrac <- phyloseq::distance(ps, method = "wunifrac")

# ============================================================================
# PERMANOVA TESTS
# ============================================================================

# Function to perform PERMANOVA and print results
perform_permanova <- function(dist_matrix, metadata) {
  # Convert metadata to data frame
  metadata_df <- as.data.frame(metadata)
  
  # Perform PERMANOVA for SampleType
  cat("PERMANOVA for SampleType:\n")
  permanova_sample_type <- adonis2(dist_matrix ~ metadata_df$SampleType, permutations = 999, by = "terms")
  print(permanova_sample_type)
  
  # Perform PERMANOVA for Microbiota
  cat("\nPERMANOVA for Microbiota:\n")
  permanova_microbiota <- adonis2(dist_matrix ~ metadata_df$Microbiota, permutations = 999, by = "terms")
  print(permanova_microbiota)
  
  # Perform PERMANOVA for SampleType and Microbiota interaction
  cat("\nPERMANOVA for SampleType * Microbiota interaction:\n")
  permanova_interaction <- adonis2(dist_matrix ~ metadata_df$SampleType * metadata_df$Microbiota, permutations = 999, by = "terms")
  print(permanova_interaction)
  
  # Return results (just the SampleType results for compatibility with existing code)
  return(permanova_sample_type)
}

# Get metadata
metadata <- as.data.frame(sample_data(ps))

# Perform PERMANOVA for each distance metric
cat("=== PERMANOVA RESULTS FOR BRAY-CURTIS DISTANCE ===\n")
permanova_bray <- perform_permanova(dist_bray, metadata)

cat("\n=== PERMANOVA RESULTS FOR UNWEIGHTED UNIFRAC DISTANCE ===\n")
permanova_unifrac <- perform_permanova(dist_unifrac, metadata)

cat("\n=== PERMANOVA RESULTS FOR WEIGHTED UNIFRAC DISTANCE ===\n")
permanova_wunifrac <- perform_permanova(dist_wunifrac, metadata)
```



## Pairwise permanova function

```R
# ============================================================================
# PAIRWISE PERMANOVA TESTS
# ============================================================================

# Function to perform pairwise PERMANOVA tests using pairwiseAdonis package
perform_pairwise_permanova <- function(ps_obj, distance_method) {
  # Calculate distance matrix
  dist_matrix <- phyloseq::distance(ps_obj, method = distance_method)
  
  # Get metadata
  metadata <- as.data.frame(sample_data(ps_obj))
  
  # Print sample types for debugging
  cat("Sample types for pairwise PERMANOVA:", paste(unique(metadata$SampleType), collapse = ", "), "\n")
  
  # Use pairwise.adonis function from pairwiseAdonis package
  cat("Performing pairwise PERMANOVA tests...\n")
  pairwise_result <- pairwiseAdonis::pairwise.adonis(
    dist_matrix, 
    metadata$SampleType, 
    p.adjust.m = "bonferroni",
    perm = 999
  )
  
  # Format results
  results <- data.frame(
    SampleType1 = character(),
    SampleType2 = character(),
    R2 = numeric(),
    p_value = numeric(),
    p_adjusted = numeric(),
    significance = character(),
    stringsAsFactors = FALSE
  )
  
  # Extract results from pairwise.adonis output
  for (i in 1:nrow(pairwise_result)) {
    # Parse the pairs from the row names
    pair <- strsplit(as.character(pairwise_result$pairs[i]), " vs ")[[1]]
    type1 <- pair[1]
    type2 <- pair[2]
    
    # Extract R2 and p-values
    r2 <- pairwise_result$R2[i]
    p_value <- pairwise_result$p.value[i]
    p_adjusted <- pairwise_result$p.adjusted[i]
    
    # Determine significance level based on adjusted p-value
    sig <- ""
    if (p_adjusted < 0.001) sig <- "***"
    else if (p_adjusted < 0.01) sig <- "**"
    else if (p_adjusted < 0.05) sig <- "*"
    else sig <- "ns"
    
    # Add to results
    results <- rbind(results, data.frame(
      SampleType1 = type1,
      SampleType2 = type2,
      R2 = r2,
      p_value = p_value,
      p_adjusted = p_adjusted,
      significance = sig,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

# Perform pairwise PERMANOVA tests for each distance metric
cat("\n=== PAIRWISE PERMANOVA RESULTS FOR BRAY-CURTIS DISTANCE ===\n")
pairwise_bray <- perform_pairwise_permanova(ps, "bray")
print(pairwise_bray)

cat("\n=== PAIRWISE PERMANOVA RESULTS FOR UNWEIGHTED UNIFRAC DISTANCE ===\n")
pairwise_unifrac <- perform_pairwise_permanova(ps, "unifrac")
print(pairwise_unifrac)

cat("\n=== PAIRWISE PERMANOVA RESULTS FOR WEIGHTED UNIFRAC DISTANCE ===\n")
pairwise_wunifrac <- perform_pairwise_permanova(ps, "wunifrac")
print(pairwise_wunifrac)
```



## Ordination functions

```R
# ============================================================================
# ORDINATION FUNCTIONS
# ============================================================================

# Function to perform PCoA ordination
perform_pcoa <- function(ps_obj, distance_method) {
  # Calculate distance matrix
  dist_matrix <- phyloseq::distance(ps_obj, method = distance_method)
  
  # Perform PCoA
  pcoa <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  
  # Calculate variance explained
  variance_explained <- round(100 * pcoa$eig / sum(pcoa$eig), 1)
  
  # Create data frame for plotting
  pcoa_df <- data.frame(
    PC1 = pcoa$points[, 1],
    PC2 = pcoa$points[, 2],
    Sample = rownames(pcoa$points)
  )
  
  # Add metadata
  metadata <- as.data.frame(sample_data(ps_obj))
  pcoa_df <- merge(pcoa_df, metadata, by.x = "Sample", by.y = "row.names")
  
  # Return results
  return(list(
    ordination = pcoa_df,
    variance_explained = variance_explained
  ))
}

# Function to perform UMAP ordination
perform_umap <- function(ps_obj, distance_method) {
  # Calculate distance matrix
  dist_matrix <- phyloseq::distance(ps_obj, method = distance_method)
  
  # Convert distance matrix to a matrix
  dist_matrix_mat <- as.matrix(dist_matrix)
  
  # Perform UMAP
  set.seed(42)  # For reproducibility
  umap_result <- umap(dist_matrix_mat, n_neighbors = min(15, nrow(dist_matrix_mat) - 1), min_dist = 0.1)
  
  # Create data frame for plotting
  umap_df <- data.frame(
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2],
    Sample = rownames(dist_matrix_mat)
  )
  
  # Add metadata
  metadata <- as.data.frame(sample_data(ps_obj))
  umap_df <- merge(umap_df, metadata, by.x = "Sample", by.y = "row.names")
  
  # Return results
  return(umap_df)
}
```



## Plotting functions

```R
# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

# Load ggExtra for marginal density plots
library(ggExtra)

# Function to create PCoA plot with density plots on margins
create_pcoa_plot <- function(pcoa_result, title, permanova_result) {
  # Extract data
  pcoa_df <- pcoa_result$ordination
  variance_explained <- pcoa_result$variance_explained
  
  # Create base plot
  p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = SampleType, shape = Microbiota)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = Stage, linetype = Stage), level = 0.95, size = 0.7, alpha = 0.5) +
    scale_color_manual(values = sample_type_colors) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
    labs(
      x = paste0("PC1 (", variance_explained[1], "%)"),
      y = paste0("PC2 (", variance_explained[2], "%)")
    ) +
    theme_pub() +
    theme(
      legend.position = "none",  # Remove legend
      plot.margin = unit(c(10, 10, 10, 10), "mm")  # top, right, bottom, left
    )
  
  # Add marginal density plots
  p_with_marginals <- ggExtra::ggMarginal(
    p, 
    type = "density", 
    margins = "both",
    groupColour = TRUE,
    groupFill = TRUE,
    alpha = 0.25,
    size = 5
  )
  
  return(p_with_marginals)
}

# Function to create UMAP plot with density plots on margins
create_umap_plot <- function(umap_df, title, permanova_result) {
  # Create base plot
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = SampleType, shape = Microbiota)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = Stage, linetype = Stage), level = 0.95, size = 0.7, alpha = 0.5) +
    scale_color_manual(values = sample_type_colors) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
    labs(
      x = "UMAP1",
      y = "UMAP2"
    ) +
    theme_pub() +
    theme(
      legend.position = "none",  # Remove legend
      plot.margin = unit(c(10, 10, 10, 10), "mm")  # top, right, bottom, left
    )
  
  # Add marginal density plots
  p_with_marginals <- ggExtra::ggMarginal(
    p, 
    type = "density", 
    margins = "both",
    groupColour = TRUE,
    groupFill = TRUE,
    alpha = 0.25,
    size = 5
  )
  
  return(p_with_marginals)
}

# ============================================================================
# PERFORM ORDINATIONS
# ============================================================================

# Bray-Curtis
pcoa_bray <- perform_pcoa(ps, "bray")
umap_bray <- perform_umap(ps, "bray")

# Unweighted UniFrac
pcoa_unifrac <- perform_pcoa(ps, "unifrac")
umap_unifrac <- perform_umap(ps, "unifrac")

# Weighted UniFrac
pcoa_wunifrac <- perform_pcoa(ps, "wunifrac")
umap_wunifrac <- perform_umap(ps, "wunifrac")

# ============================================================================
# CREATE PLOTS
# ============================================================================

# Bray-Curtis plots
pcoa_bray_plot <- create_pcoa_plot(pcoa_bray, "PCoA - Bray-Curtis", permanova_bray)
umap_bray_plot <- create_umap_plot(umap_bray, "UMAP - Bray-Curtis", permanova_bray)

# Unweighted UniFrac plots
pcoa_unifrac_plot <- create_pcoa_plot(pcoa_unifrac, "PCoA - Unweighted UniFrac", permanova_unifrac)
umap_unifrac_plot <- create_umap_plot(umap_unifrac, "UMAP - Unweighted UniFrac", permanova_unifrac)

# Weighted UniFrac plots
pcoa_wunifrac_plot <- create_pcoa_plot(pcoa_wunifrac, "PCoA - Weighted UniFrac", permanova_wunifrac)
umap_wunifrac_plot <- create_umap_plot(umap_wunifrac, "UMAP - Weighted UniFrac", permanova_wunifrac)

# ============================================================================
# COMBINE PLOTS SIDE BY SIDE
# ============================================================================

# Function to combine PCoA and UMAP plots side by side WITHOUT a legend
combine_plots <- function(pcoa_plot, umap_plot, distance_name) {
  # Combine plots without legend
  combined_plot <- cowplot::plot_grid(
    pcoa_plot, umap_plot,
    ncol = 2,
    align = "h"
  )
  
  return(combined_plot)
}

# Create a separate legend file
create_legend <- function() {
  # Create a simpler plot just for the legend, using the same theme as the plots
  legend_plot <- ggplot(pcoa_bray$ordination, aes(x = PC1, y = PC2, color = SampleType, shape = Microbiota)) +
    geom_point(size = 3) +
    stat_ellipse(aes(group = Stage, linetype = Stage), level = 0.95, size = 0.7) +
    scale_color_manual(values = sample_type_colors) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
    guides(
      linetype = guide_legend(title = "Stage", override.aes = list(linewidth = 1), keywidth = unit(2, "cm")),
      color = guide_legend(title = "Sample Type", override.aes = list(size = 3)),
      shape = guide_legend(title = "Microbiota", override.aes = list(size = 3))
    ) +
    theme_pub() +  # Apply the same theme as the plots
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      # Hide all plot elements except the legend
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank()
    )
  
  # Save the legend directly
  ggsave("simple_legend.png", legend_plot, width = 10, height = 3, dpi = 300, bg = "white")
  
  # Also return the legend for potential use in the script
  return(cowplot::get_legend(legend_plot))
}

# Combine plots
bray_combined <- combine_plots(pcoa_bray_plot, umap_bray_plot, "Bray-Curtis")
unifrac_combined <- combine_plots(pcoa_unifrac_plot, umap_unifrac_plot, "Unweighted UniFrac")
wunifrac_combined <- combine_plots(pcoa_wunifrac_plot, umap_wunifrac_plot, "Weighted UniFrac")

# Create a stacked plot of all distance metrics (without legends)
all_distances_combined <- cowplot::plot_grid(
  bray_combined,
  unifrac_combined,
  wunifrac_combined,
  ncol = 1,
  align = "v",
  axis = "lr"
)

# Create and save a separate legend file
legend_only <- create_legend()
# Note: The legend is already saved in the create_legend function

# Save plots without legends
ggsave("simple_bray_curtis.png", bray_combined, width = 16, height = 8, dpi = 300, bg = "white")
ggsave("simple_unweighted_unifrac.png", unifrac_combined, width = 16, height = 8, dpi = 300, bg = "white")
ggsave("simple_weighted_unifrac.png", wunifrac_combined, width = 16, height = 8, dpi = 300, bg = "white")
ggsave("simple_all_distances.png", all_distances_combined, width = 12, height = 24, dpi = 300, bg = "white")
```



