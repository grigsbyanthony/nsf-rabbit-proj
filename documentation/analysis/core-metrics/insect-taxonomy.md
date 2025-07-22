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

## Create taxa barplot function

```         R
create_taxa_barplot <- function(ps, taxa_level, top_n = 10, custom_colors = NULL) {
  # Agglomerate taxa at the specified level
  ps_taxa <- tax_glom(ps, taxrank = taxa_level)
  
  # Convert to relative abundance
  ps_rel <- transform_sample_counts(ps_taxa, function(x) x / sum(x))
  
  # Extract abundance matrix
  otu_rel <- as.matrix(otu_table(ps_rel))
  
  # Extract taxonomy table
  tax_table_df <- as.data.frame(tax_table(ps_rel))
  
  # Get taxa names at the specified level
  taxa_names <- tax_table_df[[taxa_level]]
  
  # Replace NA with "Unknown"
  taxa_names[is.na(taxa_names)] <- "Unknown"
  
  # Replace underscores with spaces in taxa names
  taxa_names <- gsub("_", " ", taxa_names)
  
  # Create a data frame for plotting
  plot_data <- as.data.frame(t(otu_rel))
  colnames(plot_data) <- taxa_names
  
  # Add sample metadata directly without using sample_data function
  plot_data$SampleID <- rownames(plot_data)
  
  # Get metadata directly from the phyloseq object
  metadata <- data.frame(
    SampleID = sample_names(ps_rel),
    Source = sample_data(ps_rel)$Source,
    Microbiota = sample_data(ps_rel)$Microbiota
  )
  
  # Merge data
  plot_data <- merge(plot_data, metadata, by = "SampleID")
  
  # Ensure all abundance columns are numeric
  abundance_cols <- setdiff(colnames(plot_data), c("SampleID", "Source", "Microbiota"))
  plot_data[abundance_cols] <- lapply(plot_data[abundance_cols], as.numeric)
  
  # Melt the data for ggplot
  plot_data_melted <- reshape2::melt(plot_data, 
                                     id.vars = c("SampleID", "Source", "Microbiota"),
                                     variable.name = taxa_level,
                                     value.name = "Abundance")
  
  # Ensure Abundance is numeric
  plot_data_melted$Abundance <- as.numeric(plot_data_melted$Abundance)
  
  # First, sum abundances by sample and taxa
  sample_taxa_abundance <- plot_data_melted %>%
    group_by(SampleID, Source, !!sym(taxa_level)) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  # Get top N taxa based on overall abundance
  top_taxa <- sample_taxa_abundance %>%
    group_by(!!sym(taxa_level)) %>%
    summarise(total_abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_abundance)) %>%
    slice_head(n = top_n) %>%
    pull(!!sym(taxa_level))
  
  # Add "Other" category for remaining taxa
  sample_taxa_abundance <- sample_taxa_abundance %>%
    mutate(!!sym(taxa_level) := ifelse(!!sym(taxa_level) %in% top_taxa, 
                                      as.character(!!sym(taxa_level)), 
                                      "Other"))
  
  # Replace any remaining underscores with spaces in the taxa level column
  sample_taxa_abundance[[taxa_level]] <- gsub("_", " ", sample_taxa_abundance[[taxa_level]])
  
  # Aggregate abundances for the "Other" category
  sample_taxa_abundance <- sample_taxa_abundance %>%
    group_by(SampleID, Source, !!sym(taxa_level)) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  # Calculate mean relative abundance by Source and taxa
  # First, get total abundance per sample
  sample_totals <- sample_taxa_abundance %>%
    group_by(SampleID) %>%
    summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  # Join with sample_taxa_abundance
  sample_taxa_abundance <- left_join(sample_taxa_abundance, sample_totals, by = "SampleID")
  
  # Calculate relative abundance for each sample
  sample_taxa_abundance <- sample_taxa_abundance %>%
    mutate(rel_abundance = Abundance / total)
  
  # Now calculate mean relative abundance by Source and taxa
  mean_abundance <- sample_taxa_abundance %>%
    group_by(Source, !!sym(taxa_level)) %>%
    summarise(mean_abundance = mean(rel_abundance, na.rm = TRUE), .groups = "drop")
  
  # Order Source by the custom order
  source_order <- c("Chrysomya rufifacies", "Chrysomya megacephala", 
                    "Cochliomyia macellaria", "Lucilia coeruleiviridis", 
                    "Hydrotea aenescens", "Necrodes surinamensis")
  mean_abundance$Source <- factor(mean_abundance$Source, levels = source_order)
  
  # Create the barplot
  p <- ggplot(mean_abundance, aes(x = Source, y = mean_abundance, fill = !!sym(taxa_level))) +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(expand = c(0, 0)) +  # Remove padding at bottom of y-axis
    labs(x = "Insect Species",
         y = "Relative Abundance",
         fill = taxa_level) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          legend.direction = "vertical",
          legend.box = "vertical")
  
  # Use custom colors if provided, otherwise use default
  if (!is.null(custom_colors) && length(custom_colors) >= length(unique(mean_abundance[[taxa_level]]))) {
    p <- p + scale_fill_manual(values = custom_colors)
  }
  
  return(p)
}
```

## Create taxa barplots

```R
# Phylum level
phylum_plot <- create_taxa_barplot(ps_insects, "Phylum", top_n = 10)

# Class level
class_plot <- create_taxa_barplot(ps_insects, "Class", top_n = 10)

# Order level
order_plot <- create_taxa_barplot(ps_insects, "Order", top_n = 10)

# Family level
family_plot <- create_taxa_barplot(ps_insects, "Family", top_n = 10)

# Genus level
genus_plot <- create_taxa_barplot(ps_insects, "Genus", top_n = 15)

# Combine plots
taxa_plots <- (phylum_plot / class_plot / order_plot) | (family_plot / genus_plot / plot_spacer())
taxa_plots <- taxa_plots + 
  plot_annotation(title = "Taxonomic Composition Across Insect Species",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))

# Display the combined plot
taxa_plots

# Save the combined plot
ggsave("taxonomic_composition_barplots.pdf", taxa_plots, width = 18, height = 16, units = "in")
ggsave("taxonomic_composition_barplots.png", taxa_plots, width = 18, height = 16, units = "in", dpi = 300)

# Also save individual plots
ggsave("taxonomic_composition_phylum.pdf", phylum_plot, width = 10, height = 6, units = "in")
ggsave("taxonomic_composition_phylum.png", phylum_plot, width = 10, height = 6, units = "in", dpi = 300)

ggsave("taxonomic_composition_class.pdf", class_plot, width = 10, height = 6, units = "in")
ggsave("taxonomic_composition_class.png", class_plot, width = 10, height = 6, units = "in", dpi = 300)

ggsave("taxonomic_composition_order.pdf", order_plot, width = 10, height = 6, units = "in")
ggsave("taxonomic_composition_order.png", order_plot, width = 10, height = 6, units = "in", dpi = 300)

ggsave("taxonomic_composition_family.pdf", family_plot, width = 10, height = 6, units = "in")
ggsave("taxonomic_composition_family.png", family_plot, width = 10, height = 6, units = "in", dpi = 300)

ggsave("taxonomic_composition_genus.pdf", genus_plot, width = 10, height = 6, units = "in")
ggsave("taxonomic_composition_genus.png", genus_plot, width = 10, height = 6, units = "in", dpi = 300)

# Create faceted barplots by Microbiota
# Function to create faceted taxonomic barplots
create_faceted_taxa_barplot <- function(ps, taxa_level, top_n = 10, custom_colors = NULL) {
  # Agglomerate taxa at the specified level
  ps_taxa <- tax_glom(ps, taxrank = taxa_level)
  
  # Convert to relative abundance
  ps_rel <- transform_sample_counts(ps_taxa, function(x) x / sum(x))
  
  # Extract abundance matrix
  otu_rel <- as.matrix(otu_table(ps_rel))
  
  # Extract taxonomy table
  tax_table_df <- as.data.frame(tax_table(ps_rel))
  
  # Get taxa names at the specified level
  taxa_names <- tax_table_df[[taxa_level]]
  
  # Replace NA with "Unknown"
  taxa_names[is.na(taxa_names)] <- "Unknown"
  
  # Replace underscores with spaces in taxa names
  taxa_names <- gsub("_", " ", taxa_names)
  
  # Create a data frame for plotting
  plot_data <- as.data.frame(t(otu_rel))
  colnames(plot_data) <- taxa_names
  
  # Add sample metadata directly without using sample_data function
  plot_data$SampleID <- rownames(plot_data)
  
  # Get metadata directly from the phyloseq object
  metadata <- data.frame(
    SampleID = sample_names(ps_rel),
    Source = sample_data(ps_rel)$Source,
    Microbiota = sample_data(ps_rel)$Microbiota
  )
  
  # Merge data
  plot_data <- merge(plot_data, metadata, by = "SampleID")
  
  # Ensure all abundance columns are numeric
  abundance_cols <- setdiff(colnames(plot_data), c("SampleID", "Source", "Microbiota"))
  plot_data[abundance_cols] <- lapply(plot_data[abundance_cols], as.numeric)
  
  # Melt the data for ggplot
  plot_data_melted <- reshape2::melt(plot_data, 
                                     id.vars = c("SampleID", "Source", "Microbiota"),
                                     variable.name = taxa_level,
                                     value.name = "Abundance")
  
  # Ensure Abundance is numeric
  plot_data_melted$Abundance <- as.numeric(plot_data_melted$Abundance)
  
  # First, sum abundances by sample, microbiota, and taxa
  sample_taxa_abundance <- plot_data_melted %>%
    group_by(SampleID, Source, Microbiota, !!sym(taxa_level)) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  # Get top N taxa based on overall abundance
  top_taxa <- sample_taxa_abundance %>%
    group_by(!!sym(taxa_level)) %>%
    summarise(total_abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_abundance)) %>%
    slice_head(n = top_n) %>%
    pull(!!sym(taxa_level))
  
  # Add "Other" category for remaining taxa
  sample_taxa_abundance <- sample_taxa_abundance %>%
    mutate(!!sym(taxa_level) := ifelse(!!sym(taxa_level) %in% top_taxa, 
                                      as.character(!!sym(taxa_level)), 
                                      "Other"))
  
  # Replace any remaining underscores with spaces in the taxa level column
  sample_taxa_abundance[[taxa_level]] <- gsub("_", " ", sample_taxa_abundance[[taxa_level]])
  
  # Aggregate abundances for the "Other" category
  sample_taxa_abundance <- sample_taxa_abundance %>%
    group_by(SampleID, Source, Microbiota, !!sym(taxa_level)) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  # Calculate mean relative abundance by Source, Microbiota, and taxa
  # First, get total abundance per sample
  sample_totals <- sample_taxa_abundance %>%
    group_by(SampleID) %>%
    summarise(total = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  # Join with sample_taxa_abundance
  sample_taxa_abundance <- left_join(sample_taxa_abundance, sample_totals, by = "SampleID")
  
  # Calculate relative abundance for each sample
  sample_taxa_abundance <- sample_taxa_abundance %>%
    mutate(rel_abundance = Abundance / total)
  
  # Now calculate mean relative abundance by Source, Microbiota, and taxa
  mean_abundance <- sample_taxa_abundance %>%
    group_by(Source, Microbiota, !!sym(taxa_level)) %>%
    summarise(mean_abundance = mean(rel_abundance, na.rm = TRUE), .groups = "drop")
  
  # Order Source by the custom order
  source_order <- c("Chrysomya rufifacies", "Chrysomya megacephala", 
                    "Cochliomyia macellaria", "Lucilia coeruleiviridis", 
                    "Hydrotea aenescens", "Necrodes surinamensis")
  mean_abundance$Source <- factor(mean_abundance$Source, levels = source_order)
  
  # Create the faceted barplot
  p <- ggplot(mean_abundance, aes(x = Source, y = mean_abundance, fill = !!sym(taxa_level))) +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_continuous(expand = c(0, 0)) +  # Remove padding at bottom of y-axis
    facet_wrap(~ Microbiota) +
    labs(x = "Insect Species",
         y = "Relative Abundance",
         fill = taxa_level) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          legend.direction = "vertical",
          legend.box = "vertical")
  
  # Use custom colors if provided, otherwise use default
  if (!is.null(custom_colors) && length(custom_colors) >= length(unique(mean_abundance[[taxa_level]]))) {
    p <- p + scale_fill_manual(values = custom_colors)
  }
  
  return(p)
}

# Create faceted taxonomic barplots
phylum_facet_plot <- create_faceted_taxa_barplot(ps_insects, "Phylum", top_n = 10)
genus_facet_plot <- create_faceted_taxa_barplot(ps_insects, "Genus", top_n = 15)

# Save faceted plots
ggsave("taxonomic_composition_phylum_by_microbiota.pdf", phylum_facet_plot, width = 12, height = 8, units = "in")
ggsave("taxonomic_composition_phylum_by_microbiota.png", phylum_facet_plot, width = 12, height = 8, units = "in", dpi = 300)

ggsave("taxonomic_composition_genus_by_microbiota.pdf", genus_facet_plot, width = 12, height = 8, units = "in")
ggsave("taxonomic_composition_genus_by_microbiota.png", genus_facet_plot, width = 12, height = 8, units = "in", dpi = 300)
```

