## Load necessary packages

```         R
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
library(patchwork)
library(ggpubr)  
library(rstatix)  
```

## Import and subset sequencing data using `qiime2R`

```         R
# Import data
ps <- qiime2R::qza_to_phyloseq(
  features = "merged-table-filtered-all.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

# Subset insect data
ps_insects <- subset_samples(ps, SampleType == "Insect")
r
# Further filter to exclude Phoridae sp.
ps_insects <- subset_samples(ps_insects, Source != "Phoridae sp.")

# Check if ps_insects has samples
print(paste("Number of samples in ps_insects after excluding Phoridae sp.:", nsamples(ps_insects)))
print(paste("Number of features in ps_insects:", ntaxa(ps_insects)))

# Print sample variables to check what's available
print("Sample variables:")
print(sample_variables(ps))
```

## Test `Microbiota`

```         R
# This section tests if the Microbiota (Endogenous vs. Exogenous) of specific insect species
# are significantly different in terms of alpha and beta diversity

# Alpha diversity comparison by Microbiota within each Source
cat("\n==== Alpha Diversity Comparison by Microbiota within each Insect Species ====\n")

# Function to perform t-tests for alpha diversity metrics by Microbiota within each Source
alpha_div_by_microbiota <- function(alpha_data, metric) {
  # Get unique Sources
  sources <- unique(alpha_data$Source)
  
  # Initialize results data frame
  results <- data.frame(
    Source = character(),
    Metric = character(),
    Endogenous_Mean = numeric(),
    Exogenous_Mean = numeric(),
    t_statistic = numeric(),
    p_value = numeric(),
    significant = character(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each Source
  for (i in 1:length(sources)) {
    source_name <- sources[i]
    
    # Subset data for this Source
    source_data <- alpha_data[alpha_data$Source == source_name, ]
    
    # Check if we have both Endogenous and Exogenous samples
    if (length(unique(source_data$Microbiota)) < 2) {
      cat(paste0("Skipping ", source_name, " - only one Microbiota type present\n"))
      next
    }
    
    # Perform t-test
    t_test_result <- t.test(
      source_data[source_data$Microbiota == "Endogenous", metric],
      source_data[source_data$Microbiota == "Exogenous", metric]
    )
    
    # Calculate means
    endo_mean <- mean(source_data[source_data$Microbiota == "Endogenous", metric], na.rm = TRUE)
    exo_mean <- mean(source_data[source_data$Microbiota == "Exogenous", metric], na.rm = TRUE)
    
    # Add to results
    results <- rbind(results, data.frame(
      Source = source_name,
      Metric = metric,
      Endogenous_Mean = endo_mean,
      Exogenous_Mean = exo_mean,
      t_statistic = t_test_result$statistic,
      p_value = t_test_result$p.value,
      significant = ifelse(t_test_result$p.value < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    ))
  }
  
  # Adjust p-values for multiple comparisons
  results$adj_p_value <- p.adjust(results$p_value, method = "bonferroni")
  results$significant <- ifelse(results$adj_p_value < 0.05, "Yes", "No")
  
  return(results)
}

# Perform t-tests for each alpha diversity metric
shannon_by_microbiota <- alpha_div_by_microbiota(alpha_div_meta, "Shannon")
chao1_by_microbiota <- alpha_div_by_microbiota(alpha_div_meta, "Chao1")
observed_by_microbiota <- alpha_div_by_microbiota(alpha_div_meta, "Observed")

# Print results
cat("\nShannon diversity comparison by Microbiota within each Source:\n")
print(shannon_by_microbiota)

cat("\nChao1 comparison by Microbiota within each Source:\n")
print(chao1_by_microbiota)

cat("\nObserved ASVs comparison by Microbiota within each Source:\n")
print(observed_by_microbiota)

# Save results to CSV
write.csv(shannon_by_microbiota, "alpha_diversity_shannon_by_microbiota.csv", row.names = FALSE)
write.csv(chao1_by_microbiota, "alpha_diversity_chao1_by_microbiota.csv", row.names = FALSE)
write.csv(observed_by_microbiota, "alpha_diversity_observed_by_microbiota.csv", row.names = FALSE)

# Create violin plots to visualize alpha diversity by Microbiota within each Source
# Shannon
p_shannon_microbiota <- ggplot(alpha_div_meta, aes(x = Microbiota, y = Shannon, fill = Microbiota)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~ Source, scales = "free_y") +
  labs(x = "Microbiota Type",
       y = "Shannon Diversity Index") +
  scale_fill_brewer(palette = "Set2") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"))

# Chao1
p_chao1_microbiota <- ggplot(alpha_div_meta, aes(x = Microbiota, y = Chao1, fill = Microbiota)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~ Source, scales = "free_y") +
  labs(x = "Microbiota Type",
       y = "Chao1 Richness") +
  scale_fill_brewer(palette = "Set2") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"))

# Observed ASVs
p_observed_microbiota <- ggplot(alpha_div_meta, aes(x = Microbiota, y = Observed, fill = Microbiota)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~ Source, scales = "free_y") +
  labs(x = "Microbiota Type",
       y = "Observed ASVs") +
  scale_fill_brewer(palette = "Set2") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"))

# Combine plots
alpha_microbiota_plots <- (p_shannon_microbiota / p_chao1_microbiota / p_observed_microbiota) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Alpha Diversity Comparison by Microbiota within each Insect Species",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))

# Display the combined plot
alpha_microbiota_plots

# Save the combined plot
ggsave("alpha_diversity_by_microbiota.pdf", alpha_microbiota_plots, width = 14, height = 18, units = "in")
ggsave("alpha_diversity_by_microbiota.png", alpha_microbiota_plots, width = 14, height = 18, units = "in", dpi = 300)

# Beta diversity comparison by Microbiota within each Source
cat("\n==== Beta Diversity Comparison by Microbiota within each Source ====\n")

# Function to perform PERMANOVA for beta diversity by Microbiota within each Source
beta_div_by_microbiota <- function(dist_matrix, metadata, distance_name) {
  # Get unique Sources
  sources <- unique(metadata$Source)
  
  # Initialize results data frame
  results <- data.frame(
    Source = character(),
    Distance_Metric = character(),
    F_value = numeric(),
    R2 = numeric(),
    p_value = numeric(),
    significant = character(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each Source
  for (i in 1:length(sources)) {
    source_name <- sources[i]
    
    # Subset metadata for this Source
    source_idx <- metadata$Source == source_name
    source_metadata <- metadata[source_idx, ]
    
    # Check if we have both Endogenous and Exogenous samples
    if (length(unique(source_metadata$Microbiota)) < 2) {
      cat(paste0("Skipping ", source_name, " - only one Microbiota type present\n"))
      next
    }
    
    # Subset distance matrix
    source_dist <- as.dist(as.matrix(dist_matrix)[source_idx, source_idx])
    
    # Run PERMANOVA
    perm_result <- adonis2(source_dist ~ Microbiota, data = source_metadata)
    
    # Add to results
    results <- rbind(results, data.frame(
      Source = source_name,
      Distance_Metric = distance_name,
      F_value = perm_result$F[1],
      R2 = perm_result$R2[1],
      p_value = perm_result$`Pr(>F)`[1],
      significant = ifelse(perm_result$`Pr(>F)`[1] < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    ))
  }
  
  # Adjust p-values for multiple comparisons
  results$adj_p_value <- p.adjust(results$p_value, method = "bonferroni")
  results$significant <- ifelse(results$adj_p_value < 0.05, "Yes", "No")
  
  return(results)
}

# Perform PERMANOVA for each distance metric
unifrac_by_microbiota <- beta_div_by_microbiota(unifrac_dist, metadata_df, "Unweighted UniFrac")
wunifrac_by_microbiota <- beta_div_by_microbiota(wunifrac_dist, metadata_df, "Weighted UniFrac")
bray_by_microbiota <- beta_div_by_microbiota(bray_dist, metadata_df, "Bray-Curtis")

# Print results
cat("\nUnweighted UniFrac comparison by Microbiota within each Source:\n")
print(unifrac_by_microbiota)

cat("\nWeighted UniFrac comparison by Microbiota within each Source:\n")
print(wunifrac_by_microbiota)

cat("\nBray-Curtis comparison by Microbiota within each Source:\n")
print(bray_by_microbiota)

# Save results to CSV
write.csv(unifrac_by_microbiota, "beta_diversity_unifrac_by_microbiota.csv", row.names = FALSE)
write.csv(wunifrac_by_microbiota, "beta_diversity_wunifrac_by_microbiota.csv", row.names = FALSE)
write.csv(bray_by_microbiota, "beta_diversity_bray_by_microbiota.csv", row.names = FALSE)

# Create PCoA plots for each Source, colored by Microbiota
# Function to create PCoA plot for a specific Source
create_source_pcoa <- function(ps, source_name, distance_method, weighted = FALSE) {
  # Subset phyloseq object for this Source
  ps_source <- subset_samples(ps, Source == source_name)
  
  # Check if we have enough samples
  if (nsamples(ps_source) < 3) {
    return(NULL)  # Not enough samples for meaningful ordination
  }
  
  # Calculate distance matrix
  if (distance_method == "unifrac") {
    dist_matrix <- phyloseq::distance(ps_source, method = "unifrac", weighted = weighted)
    method_name <- ifelse(weighted, "Weighted UniFrac", "Unweighted UniFrac")
  } else if (distance_method == "bray") {
    dist_matrix <- phyloseq::distance(ps_source, method = "bray")
    method_name <- "Bray-Curtis"
  } else {
    stop("Unsupported distance method")
  }
  
  # Perform PCoA
  pcoa_result <- ordinate(ps_source, method = "PCoA", distance = dist_matrix)
  
  # Extract variance explained
  var_explained <- pcoa_result$values$Eigenvalues / sum(pcoa_result$values$Eigenvalues) * 100
  
  # Create PCoA plot
  p <- plot_ordination(ps_source, pcoa_result, color = "Microbiota") +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = Microbiota), type = "t", level = 0.95, linetype = 2) +
    scale_color_brewer(palette = "Set2") +
    labs(title = paste(source_name, "-", method_name, "PCoA"),
         x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
         y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
    theme_pub() +
    theme(legend.position = "right",
          plot.title = element_text(size = 12))
  
  return(p)
}

# Create PCoA plots for each Source and distance metric
source_list <- unique(metadata_df$Source)
pcoa_plots <- list()

# Unweighted UniFrac
unifrac_plots <- list()
for (source in source_list) {
  p <- create_source_pcoa(ps_insects, source, "unifrac", weighted = FALSE)
  if (!is.null(p)) {
    unifrac_plots[[source]] <- p
  }
}

# Weighted UniFrac
wunifrac_plots <- list()
for (source in source_list) {
  p <- create_source_pcoa(ps_insects, source, "unifrac", weighted = TRUE)
  if (!is.null(p)) {
    wunifrac_plots[[source]] <- p
  }
}

# Bray-Curtis
bray_plots <- list()
for (source in source_list) {
  p <- create_source_pcoa(ps_insects, source, "bray")
  if (!is.null(p)) {
    bray_plots[[source]] <- p
  }
}

# Combine plots for each distance metric
if (length(unifrac_plots) > 0) {
  unifrac_combined <- wrap_plots(unifrac_plots, ncol = 2)
  ggsave("beta_diversity_unifrac_by_source.pdf", unifrac_combined, width = 14, height = 12, units = "in")
  ggsave("beta_diversity_unifrac_by_source.png", unifrac_combined, width = 14, height = 12, units = "in", dpi = 300)
}

if (length(wunifrac_plots) > 0) {
  wunifrac_combined <- wrap_plots(wunifrac_plots, ncol = 2)
  ggsave("beta_diversity_wunifrac_by_source.pdf", wunifrac_combined, width = 14, height = 12, units = "in")
  ggsave("beta_diversity_wunifrac_by_source.png", wunifrac_combined, width = 14, height = 12, units = "in", dpi = 300)
}

if (length(bray_plots) > 0) {
  bray_combined <- wrap_plots(bray_plots, ncol = 2)
  ggsave("beta_diversity_bray_by_source.pdf", bray_combined, width = 14, height = 12, units = "in")
  ggsave("beta_diversity_bray_by_source.png", bray_combined, width = 14, height = 12, units = "in", dpi = 300)
}

# Display one of the combined plots
if (length(unifrac_plots) > 0) {
  unifrac_combined
}
```
