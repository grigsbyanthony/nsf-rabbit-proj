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



## Load custom theme

```R
# Define a custom theme function
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



## Load and filter `phyloseq` object

```R
ps <- qiime2R::qza_to_phyloseq(
  features = "merged-table-filtered-all.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

# Subset insect data
ps_insects <- subset_samples(ps, SampleType == "Insect")
ps_insects <- subset_samples(ps_insects, Source != "Phoridae sp.")


# Remove taxa that are not present in these samples (have a sum of 0)
ps_insects <- prune_taxa(taxa_sums(ps_insects) > 0, ps_insects)

# We don't need to check for stages as they're not relevant for this analysis

# Check what source values are available
cat("\n=== AVAILABLE SOURCE VALUES ===\n")
available_sources <- unique(sample_data(ps_insects)$Source)
cat("Available sources:", paste(available_sources, collapse = ", "), "\n")

# Check what microbiota values are available
cat("\n=== AVAILABLE MICROBIOTA VALUES ===\n")
available_microbiota <- unique(sample_data(ps_insects)$Microbiota)
cat("Available microbiota:", paste(available_microbiota, collapse = ", "), "\n")

# No need to process Stage variable as it's not used in this analysis

# Ensure Source is a factor (this is important for plotting)
sample_data(ps_insects)$Source <- factor(sample_data(ps_insects)$Source)
available_sources <- levels(sample_data(ps_insects)$Source)

# Ensure Microbiota is a factor (this is important for faceting)
sample_data(ps_insects)$Microbiota <- factor(sample_data(ps_insects)$Microbiota)
available_microbiota <- levels(sample_data(ps_insects)$Microbiota)

# Define a color palette for Source values
# Make sure we have enough colors for all source levels
n_colors <- max(3, length(available_sources))  # At least 3 colors
source_colors <- brewer.pal(n_colors, "Set1")[1:length(available_sources)]
names(source_colors) <- available_sources

# Print source and microbiota levels to verify
cat("\nSource levels (as factors):", paste(levels(sample_data(ps_insects)$Source), collapse = ", "), "\n")
cat("\nMicrobiota levels (as factors):", paste(levels(sample_data(ps_insects)$Microbiota), collapse = ", "), "\n")
```



## PERMANOVA

```R
# ============================================================================
# BETA DIVERSITY CALCULATIONS
# ============================================================================

# Calculate distance matrices
dist_bray <- phyloseq::distance(ps_insects, method = "bray")
dist_unifrac <- phyloseq::distance(ps_insects, method = "unifrac")
dist_wunifrac <- phyloseq::distance(ps_insects, method = "wunifrac")

# ============================================================================
# PERMANOVA TESTS
# ============================================================================

# Function to perform PERMANOVA and print results
perform_permanova <- function(dist_matrix, metadata) {
  # Convert metadata to data frame
  metadata_df <- as.data.frame(metadata)
  
  # Perform PERMANOVA for Source
  cat("PERMANOVA for Source:\n")
  permanova_source <- adonis2(dist_matrix ~ metadata_df$Source, permutations = 999, by = "terms")
  print(permanova_source)
  
  # Perform PERMANOVA for Microbiota
  cat("\nPERMANOVA for Microbiota:\n")
  permanova_microbiota <- adonis2(dist_matrix ~ metadata_df$Microbiota, permutations = 999, by = "terms")
  print(permanova_microbiota)
  
  # Perform PERMANOVA for Source and Microbiota interaction
  cat("\nPERMANOVA for Source * Microbiota interaction:\n")
  permanova_interaction <- adonis2(dist_matrix ~ metadata_df$Source * metadata_df$Microbiota, permutations = 999, by = "terms")
  print(permanova_interaction)
  
  # Return results (just the Source results for compatibility with existing code)
  return(permanova_source)
}

# Get metadata
metadata <- as.data.frame(sample_data(ps_insects))

# Perform PERMANOVA for each distance metric
cat("=== PERMANOVA RESULTS FOR BRAY-CURTIS DISTANCE ===\n")
permanova_bray <- perform_permanova(dist_bray, metadata)

cat("\n=== PERMANOVA RESULTS FOR UNWEIGHTED UNIFRAC DISTANCE ===\n")
permanova_unifrac <- perform_permanova(dist_unifrac, metadata)

cat("\n=== PERMANOVA RESULTS FOR WEIGHTED UNIFRAC DISTANCE ===\n")
permanova_wunifrac <- perform_permanova(dist_wunifrac, metadata)
```



## Pairwise PERMANOVA

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
  
  # Print sources for debugging
  cat("Sources for pairwise PERMANOVA:", paste(unique(metadata$Source), collapse = ", "), "\n")
  
  # Use pairwise.adonis function from pairwiseAdonis package
  cat("Performing pairwise PERMANOVA tests...\n")
  pairwise_result <- pairwiseAdonis::pairwise.adonis(
    dist_matrix, 
    metadata$Source, 
    p.adjust.m = "bonferroni",
    perm = 999
  )
  
  # Format results
  results <- data.frame(
    Source1 = character(),
    Source2 = character(),
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
    src1 <- pair[1]
    src2 <- pair[2]
    
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
      Source1 = src1,
      Source2 = src2,
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
pairwise_bray <- perform_pairwise_permanova(ps_insects, "bray")
print(pairwise_bray)

cat("\n=== PAIRWISE PERMANOVA RESULTS FOR UNWEIGHTED UNIFRAC DISTANCE ===\n")
pairwise_unifrac <- perform_pairwise_permanova(ps_insects, "unifrac")
print(pairwise_unifrac)

cat("\n=== PAIRWISE PERMANOVA RESULTS FOR WEIGHTED UNIFRAC DISTANCE ===\n")
pairwise_wunifrac <- perform_pairwise_permanova(ps_insects, "wunifrac")
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
  
  # Create a list to store plots for each microbiota
  plot_list <- list()
  
  # Get unique microbiota values
  microbiota_values <- unique(pcoa_df$Microbiota)
  
  # Create a plot for each microbiota with marginal density plots
  for (mb in microbiota_values) {
    # Subset data for this microbiota
    mb_data <- pcoa_df[pcoa_df$Microbiota == mb, ]
    
    # Create base plot
    p <- ggplot(mb_data, aes(x = PC1, y = PC2, color = Source)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(values = source_colors) +
      labs(
        # Remove title
        x = paste0("PC1 (", variance_explained[1], "%)"),
        y = paste0("PC2 (", variance_explained[2], "%)")
      ) +
      theme_pub() +
      theme(
        legend.position = "none",  # Remove legend
        plot.margin = unit(c(10, 10, 10, 10), "mm")  # Consistent margins
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
    
    # Add to list
    plot_list[[mb]] <- p_with_marginals
  }
  
  # Combine all plots into a grid
  combined_plot <- cowplot::plot_grid(
    plotlist = plot_list,
    ncol = 2,  # Adjust based on number of microbiota values
    align = "hv"
  )
  
  return(combined_plot)
}

# Function to create UMAP plot with density plots on margins
create_umap_plot <- function(umap_df, title, permanova_result) {
  # Create a list to store plots for each microbiota
  plot_list <- list()
  
  # Get unique microbiota values
  microbiota_values <- unique(umap_df$Microbiota)
  
  # Create a plot for each microbiota with marginal density plots
  for (mb in microbiota_values) {
    # Subset data for this microbiota
    mb_data <- umap_df[umap_df$Microbiota == mb, ]
    
    # Create base plot
    p <- ggplot(mb_data, aes(x = UMAP1, y = UMAP2, color = Source)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_manual(values = source_colors) +
      labs(
        # Remove title
        x = "UMAP1",
        y = "UMAP2"
      ) +
      theme_pub() +
      theme(
        legend.position = "none",  # Remove legend
        plot.margin = unit(c(10, 10, 10, 10), "mm")  # Consistent margins
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
    
    # Add to list
    plot_list[[mb]] <- p_with_marginals
  }
  
  # Combine all plots into a grid
  combined_plot <- cowplot::plot_grid(
    plotlist = plot_list,
    ncol = 2,  # Adjust based on number of microbiota values
    align = "hv"
  )
  
  return(combined_plot)
}

# ============================================================================
# PERFORM ORDINATIONS
# ============================================================================

# Bray-Curtis
pcoa_bray <- perform_pcoa(ps_insects, "bray")
umap_bray <- perform_umap(ps_insects, "bray")

# Unweighted UniFrac
pcoa_unifrac <- perform_pcoa(ps_insects, "unifrac")
umap_unifrac <- perform_umap(ps_insects, "unifrac")

# Weighted UniFrac
pcoa_wunifrac <- perform_pcoa(ps_insects, "wunifrac")
umap_wunifrac <- perform_umap(ps_insects, "wunifrac")

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
  # Combine plots without legend and without labels
  combined_plot <- cowplot::plot_grid(
    pcoa_plot, umap_plot,
    ncol = 1,
    align = "v"
    # Removed labels
  )
  
  return(combined_plot)
}

# Create a separate legend file
create_legend <- function() {
  # Create a simpler plot just for the legend, using the same theme as the plots
  legend_plot <- ggplot(pcoa_bray$ordination, aes(x = PC1, y = PC2, color = Source)) +
    geom_point(size = 3) +
    scale_color_manual(values = source_colors) +
    guides(
      color = guide_legend(title = "Source", override.aes = list(size = 3))
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
  ggsave("insects_legend.png", legend_plot, width = 10, height = 3, dpi = 300, bg = "white")
  
  # Also return the legend for potential use in the script
  return(cowplot::get_legend(legend_plot))
}

# Function to create a stacked plot of all distance metrics WITHOUT legends
create_stacked_plot <- function(bray_combined, unifrac_combined, wunifrac_combined) {
  # Combine all plots vertically without legends
  stacked_plot <- cowplot::plot_grid(
    bray_combined,
    unifrac_combined,
    wunifrac_combined,
    ncol = 1,
    align = "v",
    axis = "lr"
  )
  
  return(stacked_plot)
}

# Combine plots
bray_combined <- combine_plots(pcoa_bray_plot, umap_bray_plot, "Bray-Curtis")
unifrac_combined <- combine_plots(pcoa_unifrac_plot, umap_unifrac_plot, "Unweighted UniFrac")
wunifrac_combined <- combine_plots(pcoa_wunifrac_plot, umap_wunifrac_plot, "Weighted UniFrac")

# Create stacked plot with all distance metrics
all_distances_combined <- create_stacked_plot(bray_combined, unifrac_combined, wunifrac_combined)

# Create and save a separate legend file
legend_only <- create_legend()
# Note: The legend is already saved in the create_legend function

# Save plots without legends
ggsave("insects_bray_curtis.png", bray_combined, width = 14, height = 14, dpi = 300, bg = "white")
ggsave("insects_unweighted_unifrac.png", unifrac_combined, width = 14, height = 14, dpi = 300, bg = "white")
ggsave("insects_weighted_unifrac.png", wunifrac_combined, width = 14, height = 14, dpi = 300, bg = "white")
ggsave("insects_all_distances.png", all_distances_combined, width = 12, height = 18, dpi = 300, bg = "white")
```

