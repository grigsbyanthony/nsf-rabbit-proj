```R
# ==================================================================================
# Relative Abundance Visualization for Soil Samples Faceted by Stage
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

# Check unique values for Stage
sample_data_df <- as.data.frame(sample_data(ps_soil))
cat("\nUnique values for Stage:", paste(unique(sample_data_df$Stage), collapse = ", "), "\n")

# Ensure Stage is a factor with the desired order (Fr, Ac, Adv, NA)
sample_data(ps_soil)$Stage <- factor(sample_data(ps_soil)$Stage, levels = c("Fr", "Ac", "Adv", NA))

# Check for NA values in Stage
sample_data_df <- as.data.frame(sample_data(ps_soil))
cat("\nNumber of samples with NA Stage:", sum(is.na(sample_data_df$Stage)), "\n")

# Remove samples with NA Stage if any
ps_soil <- subset_samples(ps_soil, !is.na(Stage))

# Check which stages have samples
stage_counts <- table(sample_data(ps_soil)$Stage)
cat("\nNumber of samples per stage:\n")
print(stage_counts)

# Check if we have at least one sample for each stage
if (any(stage_counts == 0)) {
  cat("\nWARNING: Some stages have no samples. This may cause plotting issues.\n")
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
# Function to prepare data for visualization
# ==================================================================================
prepare_abundance_data <- function(ps_obj, tax_level, top_n = 10) {
  # Extract data from phyloseq object
  otu_table <- as.data.frame(otu_table(ps_obj))
  tax_table <- as.data.frame(tax_table(ps_obj))
  sample_data <- as.data.frame(sample_data(ps_obj))
  
  # Ensure OTUs are rows
  if (!taxa_are_rows(ps_obj)) {
    otu_table <- t(otu_table)
  }
  
  # Merge taxonomy with OTU table
  otu_tax <- merge(otu_table, tax_table, by = "row.names")
  rownames(otu_tax) <- otu_tax$Row.names
  otu_tax$Row.names <- NULL
  
  # Aggregate by taxonomic level
  tax_level_data <- otu_tax %>%
    group_by(!!sym(tax_level)) %>%
    summarise(across(where(is.numeric), sum)) %>%
    filter(!is.na(!!sym(tax_level)) & !!sym(tax_level) != "")
  
  # Get top taxa based on total abundance across all samples
  total_abundance <- rowSums(select(tax_level_data, where(is.numeric)))
  tax_level_data$TotalAbundance <- total_abundance
  
  top_taxa <- tax_level_data %>%
    arrange(desc(TotalAbundance)) %>%
    head(top_n) %>%
    pull(!!sym(tax_level))
  
  # Create a copy for reshaping
  tax_data_for_reshape <- tax_level_data %>%
    select(-TotalAbundance) %>%
    # Convert to long format FIRST
    pivot_longer(cols = -!!sym(tax_level), names_to = "SampleID", values_to = "Abundance")
  
  # Now group low-abundance taxa into "Other" category PER SAMPLE
  # This ensures each sample still sums to 1.0
  plot_data <- tax_data_for_reshape %>%
    mutate(TaxaGroup = ifelse(!!sym(tax_level) %in% top_taxa, !!sym(tax_level), "Other")) %>%
    group_by(SampleID, TaxaGroup) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # Verify that each sample sums to 1.0
  sample_sums_check <- plot_data %>%
    group_by(SampleID) %>%
    summarise(SampleSum = sum(Abundance), .groups = "drop")
  
  cat("Sample sums after grouping (should all be ~1.0):\n")
  cat("Min:", min(sample_sums_check$SampleSum), "Max:", max(sample_sums_check$SampleSum), "\n")
  
  # Add sample metadata
  plot_data <- merge(plot_data, sample_data, by.x = "SampleID", by.y = "row.names")
  
  return(plot_data)
}

# ==================================================================================
# Create faceted bar plot visualization for phyla
# ==================================================================================
create_phylum_barplot <- function(ps_obj, top_n = 10) {
  tax_level <- "Phylum"
  
  # Prepare data
  plot_data <- prepare_abundance_data(ps_obj, tax_level, top_n)
  
  # Calculate mean abundance by taxa and stage
  bar_data <- plot_data %>%
    group_by(Stage, TaxaGroup) %>%
    summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    # Renormalize so each stage sums to 1.0
    group_by(Stage) %>%
    mutate(GroupSum = sum(MeanAbundance, na.rm = TRUE),
           NormalizedAbundance = ifelse(GroupSum > 0, MeanAbundance / GroupSum, 0)) %>%
    ungroup()
  
  # Check for NA or NaN values
  if (any(is.na(bar_data$NormalizedAbundance)) || any(is.nan(bar_data$NormalizedAbundance))) {
    cat("\nWARNING: NA or NaN values found in normalized abundance. This may cause plotting issues.\n")
    # Replace NA/NaN with 0
    bar_data$NormalizedAbundance[is.na(bar_data$NormalizedAbundance) | is.nan(bar_data$NormalizedAbundance)] <- 0
  }
  
  # Verify that each group now sums to 1.0
  verification <- bar_data %>%
    group_by(Stage) %>%
    summarise(GroupSum = sum(NormalizedAbundance), .groups = "drop")
  
  cat("Group sums after normalization (should all be 1.0):\n")
  cat("Min:", min(verification$GroupSum), "Max:", max(verification$GroupSum), "\n")
  
  # Order taxa by overall abundance
  taxa_order <- bar_data %>%
    group_by(TaxaGroup) %>%
    summarise(TotalAbundance = sum(NormalizedAbundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalAbundance)) %>%
    pull(TaxaGroup)
  
  # Replace underscores with spaces in phylum names
  bar_data$TaxaGroup <- gsub("_", " ", bar_data$TaxaGroup)
  taxa_order <- gsub("_", " ", taxa_order)
  
  # Make sure "Other" is at the end of the legend
  if("Other" %in% taxa_order) {
    taxa_order <- c(taxa_order[taxa_order != "Other"], "Other")
  }
  
  # Set factor levels for ordering
  bar_data$TaxaGroup <- factor(bar_data$TaxaGroup, levels = rev(taxa_order))
  
  # Use custom color palette with "Other" as light grey
  custom_colors <- c("#0ea5e9", "#06b6d4", "#14b8a6", "#10b981", "#22c55e", 
                    "#84cc16", "#eab308", "#f59e0b", "#f97316", "#ef4444")
  
  # Get unique taxa groups excluding "Other"
  taxa_groups <- levels(bar_data$TaxaGroup)
  other_index <- which(taxa_groups == "Other")
  
  # Create the final color vector with "Other" as light grey
  if (length(other_index) > 0) {
    # If "Other" is present, assign colors to all taxa except "Other"
    taxa_without_other <- taxa_groups[-other_index]
    n_taxa <- length(taxa_without_other)
    
    # Assign colors to taxa (excluding "Other")
    colors <- c()
    for (i in 1:length(taxa_groups)) {
      if (taxa_groups[i] == "Other") {
        colors <- c(colors, "#e5e5e5")  # Light grey for "Other"
      } else {
        # Find position in taxa_without_other
        pos <- which(taxa_without_other == taxa_groups[i])
        # Use custom color if available, otherwise use default
        if (pos <= length(custom_colors)) {
          colors <- c(colors, custom_colors[pos])
        } else {
          colors <- c(colors, "#999999")  # Fallback color
        }
      }
    }
  } else {
    # If "Other" is not present, just use the custom colors directly
    n_taxa <- length(taxa_groups)
    colors <- custom_colors[1:min(n_taxa, length(custom_colors))]
  }
  
  # Create bar plot with faceting by Stage
  p <- ggplot(bar_data, aes(x = Stage, y = NormalizedAbundance, fill = TaxaGroup)) +
    geom_col(position = "stack", color = "white", size = 0.1) +
    scale_fill_manual(values = colors, name = "Phylum") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0), oob = scales::squish) +
    labs(
      title = "Relative Abundance of Phyla in Soil Samples by Stage",
      y = "Mean Relative Abundance",
      x = "Stage"
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
# Create faceted bar plot visualization for genera
# ==================================================================================
create_genus_barplot <- function(ps_obj, top_n = 10) {
  tax_level <- "Genus"
  
  # Prepare data
  plot_data <- prepare_abundance_data(ps_obj, tax_level, top_n)
  
  # Calculate mean abundance by taxa and stage
  bar_data <- plot_data %>%
    group_by(Stage, TaxaGroup) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = "drop") %>%
    # Renormalize so each stage sums to 1.0
    group_by(Stage) %>%
    mutate(GroupSum = sum(MeanAbundance),
           NormalizedAbundance = MeanAbundance / GroupSum) %>%
    ungroup()
  
  # Verify that each group now sums to 1.0
  verification <- bar_data %>%
    group_by(Stage) %>%
    summarise(GroupSum = sum(NormalizedAbundance), .groups = "drop")
  
  cat("Group sums after normalization (should all be 1.0):\n")
  cat("Min:", min(verification$GroupSum), "Max:", max(verification$GroupSum), "\n")
  
  # Order taxa by overall abundance
  taxa_order <- bar_data %>%
    group_by(TaxaGroup) %>%
    summarise(TotalAbundance = sum(NormalizedAbundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalAbundance)) %>%
    pull(TaxaGroup)
  
  # Replace underscores with spaces in genus names
  bar_data$TaxaGroup <- gsub("_", " ", bar_data$TaxaGroup)
  taxa_order <- gsub("_", " ", taxa_order)
  
  # Make sure "Other" is at the end of the legend
  if("Other" %in% taxa_order) {
    taxa_order <- c(taxa_order[taxa_order != "Other"], "Other")
  }
  
  # Set factor levels for ordering
  bar_data$TaxaGroup <- factor(bar_data$TaxaGroup, levels = rev(taxa_order))
  
  # Use custom color palette with "Other" as light grey
  custom_colors <- c("#0ea5e9", "#06b6d4", "#14b8a6", "#10b981", "#22c55e", 
                    "#84cc16", "#eab308", "#f59e0b", "#f97316", "#ef4444")
  
  # Get unique taxa groups excluding "Other"
  taxa_groups <- levels(bar_data$TaxaGroup)
  other_index <- which(taxa_groups == "Other")
  
  # Create the final color vector with "Other" as light grey
  if (length(other_index) > 0) {
    # If "Other" is present, assign colors to all taxa except "Other"
    taxa_without_other <- taxa_groups[-other_index]
    n_taxa <- length(taxa_without_other)
    
    # Assign colors to taxa (excluding "Other")
    colors <- c()
    for (i in 1:length(taxa_groups)) {
      if (taxa_groups[i] == "Other") {
        colors <- c(colors, "#e5e5e5")  # Light grey for "Other"
      } else {
        # Find position in taxa_without_other
        pos <- which(taxa_without_other == taxa_groups[i])
        # Use custom color if available, otherwise use default
        if (pos <= length(custom_colors)) {
          colors <- c(colors, custom_colors[pos])
        } else {
          colors <- c(colors, "#999999")  # Fallback color
        }
      }
    }
  } else {
    # If "Other" is not present, just use the custom colors directly
    n_taxa <- length(taxa_groups)
    colors <- custom_colors[1:min(n_taxa, length(custom_colors))]
  }
  
  # Create bar plot with faceting by Stage
  p <- ggplot(bar_data, aes(x = Stage, y = NormalizedAbundance, fill = TaxaGroup)) +
    geom_col(position = "stack", color = "white", size = 0.1) +
    scale_fill_manual(values = colors, name = "Genus") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0), oob = scales::squish) +
    labs(
      title = "Relative Abundance of Genera in Soil Samples by Stage",
      y = "Mean Relative Abundance",
      x = "Stage"
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
# Create horizontal bar plots (with coord_flip)
# ==================================================================================
create_horizontal_phylum_barplot <- function(ps_obj, top_n = 10) {
  tax_level <- "Phylum"
  
  # Prepare data
  plot_data <- prepare_abundance_data(ps_obj, tax_level, top_n)
  
  # Calculate mean abundance by taxa and stage
  bar_data <- plot_data %>%
    group_by(Stage, TaxaGroup) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = "drop") %>%
    # Renormalize so each stage sums to 1.0
    group_by(Stage) %>%
    mutate(GroupSum = sum(MeanAbundance),
           NormalizedAbundance = MeanAbundance / GroupSum) %>%
    ungroup()
  
  # Order taxa by overall abundance
  taxa_order <- bar_data %>%
    group_by(TaxaGroup) %>%
    summarise(TotalAbundance = sum(NormalizedAbundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalAbundance)) %>%
    pull(TaxaGroup)
  
  # Replace underscores with spaces in phylum names
  bar_data$TaxaGroup <- gsub("_", " ", bar_data$TaxaGroup)
  taxa_order <- gsub("_", " ", taxa_order)
  
  # Make sure "Other" is at the end of the legend
  if("Other" %in% taxa_order) {
    taxa_order <- c(taxa_order[taxa_order != "Other"], "Other")
  }
  
  # Set factor levels for ordering
  bar_data$TaxaGroup <- factor(bar_data$TaxaGroup, levels = rev(taxa_order))
  
  # Use custom color palette with "Other" as light grey
  custom_colors <- c("#0ea5e9", "#06b6d4", "#14b8a6", "#10b981", "#22c55e", 
                    "#84cc16", "#eab308", "#f59e0b", "#f97316", "#ef4444")
  
  # Get unique taxa groups excluding "Other"
  taxa_groups <- levels(bar_data$TaxaGroup)
  other_index <- which(taxa_groups == "Other")
  
  # Create the final color vector with "Other" as light grey
  if (length(other_index) > 0) {
    # If "Other" is present, assign colors to all taxa except "Other"
    taxa_without_other <- taxa_groups[-other_index]
    n_taxa <- length(taxa_without_other)
    
    # Assign colors to taxa (excluding "Other")
    colors <- c()
    for (i in 1:length(taxa_groups)) {
      if (taxa_groups[i] == "Other") {
        colors <- c(colors, "#e5e5e5")  # Light grey for "Other"
      } else {
        # Find position in taxa_without_other
        pos <- which(taxa_without_other == taxa_groups[i])
        # Use custom color if available, otherwise use default
        if (pos <= length(custom_colors)) {
          colors <- c(colors, custom_colors[pos])
        } else {
          colors <- c(colors, "#999999")  # Fallback color
        }
      }
    }
  } else {
    # If "Other" is not present, just use the custom colors directly
    n_taxa <- length(taxa_groups)
    colors <- custom_colors[1:min(n_taxa, length(custom_colors))]
  }
  
  # Create horizontal bar plot
  p <- ggplot(bar_data, aes(x = Stage, y = NormalizedAbundance, fill = TaxaGroup)) +
    geom_col(position = "stack", color = "white", size = 0.1) +
    scale_fill_manual(values = colors, name = "Phylum") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0), oob = scales::squish) +
    coord_flip() +  # Flip coordinates to make horizontal bars
    labs(
      title = "Relative Abundance of Phyla in Soil Samples by Stage",
      y = "Mean Relative Abundance",
      x = "Stage"
    ) +
    theme_pub(legend_pos = "bottom") +
    theme(
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      legend.direction = "horizontal",
      legend.box = "horizontal"
    ) +
    guides(fill = guide_legend(nrow = 2))  # Arrange legend in 2 rows
  
  return(p)
}

create_horizontal_genus_barplot <- function(ps_obj, top_n = 10) {
  tax_level <- "Genus"
  
  # Prepare data
  plot_data <- prepare_abundance_data(ps_obj, tax_level, top_n)
  
  # Calculate mean abundance by taxa and stage
  bar_data <- plot_data %>%
    group_by(Stage, TaxaGroup) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = "drop") %>%
    # Renormalize so each stage sums to 1.0
    group_by(Stage) %>%
    mutate(GroupSum = sum(MeanAbundance),
           NormalizedAbundance = MeanAbundance / GroupSum) %>%
    ungroup()
  
  # Order taxa by overall abundance
  taxa_order <- bar_data %>%
    group_by(TaxaGroup) %>%
    summarise(TotalAbundance = sum(NormalizedAbundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalAbundance)) %>%
    pull(TaxaGroup)
  
  # Replace underscores with spaces in genus names
  bar_data$TaxaGroup <- gsub("_", " ", bar_data$TaxaGroup)
  taxa_order <- gsub("_", " ", taxa_order)
  
  # Make sure "Other" is at the end of the legend
  if("Other" %in% taxa_order) {
    taxa_order <- c(taxa_order[taxa_order != "Other"], "Other")
  }
  
  # Set factor levels for ordering
  bar_data$TaxaGroup <- factor(bar_data$TaxaGroup, levels = rev(taxa_order))
  
  # Use custom color palette with "Other" as light grey
  custom_colors <- c("#0ea5e9", "#06b6d4", "#14b8a6", "#10b981", "#22c55e", 
                    "#84cc16", "#eab308", "#f59e0b", "#f97316", "#ef4444")
  
  # Get unique taxa groups excluding "Other"
  taxa_groups <- levels(bar_data$TaxaGroup)
  other_index <- which(taxa_groups == "Other")
  
  # Create the final color vector with "Other" as light grey
  if (length(other_index) > 0) {
    # If "Other" is present, assign colors to all taxa except "Other"
    taxa_without_other <- taxa_groups[-other_index]
    n_taxa <- length(taxa_without_other)
    
    # Assign colors to taxa (excluding "Other")
    colors <- c()
    for (i in 1:length(taxa_groups)) {
      if (taxa_groups[i] == "Other") {
        colors <- c(colors, "#e5e5e5")  # Light grey for "Other"
      } else {
        # Find position in taxa_without_other
        pos <- which(taxa_without_other == taxa_groups[i])
        # Use custom color if available, otherwise use default
        if (pos <= length(custom_colors)) {
          colors <- c(colors, custom_colors[pos])
        } else {
          colors <- c(colors, "#999999")  # Fallback color
        }
      }
    }
  } else {
    # If "Other" is not present, just use the custom colors directly
    n_taxa <- length(taxa_groups)
    colors <- custom_colors[1:min(n_taxa, length(custom_colors))]
  }
  
  # Create horizontal bar plot
  p <- ggplot(bar_data, aes(x = Stage, y = NormalizedAbundance, fill = TaxaGroup)) +
    geom_col(position = "stack", color = "white", size = 0.1) +
    scale_fill_manual(values = colors, name = "Genus") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0), oob = scales::squish) +
    coord_flip() +  # Flip coordinates to make horizontal bars
    labs(
      title = "Relative Abundance of Genera in Soil Samples by Stage",
      y = "Mean Relative Abundance",
      x = "Stage"
    ) +
    theme_pub(legend_pos = "bottom") +
    theme(
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      legend.direction = "horizontal",
      legend.box = "horizontal"
    ) +
    guides(fill = guide_legend(nrow = 2))  # Arrange legend in 2 rows
  
  return(p)
}

# ==================================================================================
# Function to combine phylum and genus plots
# ==================================================================================

# Function to create a combined plot with phylum and genus
create_combined_plot <- function(phylum_plot, genus_plot, layout = "vertical") {
  # Remove titles from individual plots for better combined appearance
  phylum_plot <- phylum_plot + 
    labs(title = NULL) +
    theme(plot.margin = margin(5, 5, 10, 5))
  
  genus_plot <- genus_plot + 
    labs(title = NULL) +
    theme(plot.margin = margin(5, 5, 10, 5))
  
  # Combine plots based on layout without title
  if (layout == "vertical") {
    # Stack plots vertically
    combined_plot <- phylum_plot / genus_plot +
      plot_layout(heights = c(1, 1))
  } else {
    # Arrange plots side by side
    combined_plot <- phylum_plot + genus_plot +
      plot_layout(widths = c(1, 1))
  }
  
  return(combined_plot)
}

# ==================================================================================
# Generate and save visualizations
# ==================================================================================

# Create vertical bar plots
phylum_plot <- create_phylum_barplot(ps_rel, top_n = 10)
genus_plot <- create_genus_barplot(ps_rel, top_n = 10)

# Create horizontal bar plots
phylum_plot_horizontal <- create_horizontal_phylum_barplot(ps_rel, top_n = 10)
genus_plot_horizontal <- create_horizontal_genus_barplot(ps_rel, top_n = 10)

# Create combined plots
combined_vertical <- create_combined_plot(phylum_plot, genus_plot, "vertical")
combined_horizontal <- create_combined_plot(phylum_plot_horizontal, genus_plot_horizontal, "horizontal")

# Save individual plots
ggsave("soil_phylum_abundance_by_stage.png", phylum_plot, width = 10, height = 8, dpi = 300)
ggsave("soil_genus_abundance_by_stage.png", genus_plot, width = 10, height = 8, dpi = 300)
ggsave("soil_phylum_abundance_by_stage_horizontal.png", phylum_plot_horizontal, width = 10, height = 6, dpi = 300)
ggsave("soil_genus_abundance_by_stage_horizontal.png", genus_plot_horizontal, width = 10, height = 6, dpi = 300)

# Save combined plots
ggsave("soil_combined_vertical.png", combined_vertical, width = 10, height = 16, dpi = 300)
ggsave("soil_combined_horizontal.png", combined_horizontal, width = 20, height = 8, dpi = 300)

# Display the combined plots
print(combined_vertical)
print(combined_horizontal)

# Print summary message
cat("\nRelative abundance visualizations created:\n")
cat("1. Phylum-level bar plot by stage (soil_phylum_abundance_by_stage.png)\n")
cat("2. Genus-level bar plot by stage (soil_genus_abundance_by_stage.png)\n")
cat("3. Horizontal phylum-level bar plot by stage (soil_phylum_abundance_by_stage_horizontal.png)\n")
cat("4. Horizontal genus-level bar plot by stage (soil_genus_abundance_by_stage_horizontal.png)\n")
cat("5. Combined vertical plot with phylum and genus (soil_combined_vertical.png)\n")
cat("6. Combined horizontal plot with phylum and genus (soil_combined_horizontal.png)\n")
cat("\nAll bars in the stacked plots sum to 1.0 (100%) per stage!\n")

```

```R
# ==================================================================================
# Relative Abundance Visualization for Swab Samples Faceted by Stage
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

# Check unique values for Stage
sample_data_df <- as.data.frame(sample_data(ps_swabs))
cat("\nUnique values for Stage:", paste(unique(sample_data_df$Stage), collapse = ", "), "\n")

# Ensure Stage is a factor with the desired order (Fr, Ac, Adv, NA)
sample_data(ps_swabs)$Stage <- factor(sample_data(ps_swabs)$Stage, levels = c("Fr", "Ac", "Adv", NA))

# Remove samples with NA Stage if any
ps_swabs <- subset_samples(ps_swabs, !is.na(Stage))

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
# Function to prepare data for visualization
# ==================================================================================
prepare_abundance_data <- function(ps_obj, tax_level, top_n = 10) {
  # Extract data from phyloseq object
  otu_table <- as.data.frame(otu_table(ps_obj))
  tax_table <- as.data.frame(tax_table(ps_obj))
  sample_data <- as.data.frame(sample_data(ps_obj))
  
  # Ensure OTUs are rows
  if (!taxa_are_rows(ps_obj)) {
    otu_table <- t(otu_table)
  }
  
  # Merge taxonomy with OTU table
  otu_tax <- merge(otu_table, tax_table, by = "row.names")
  rownames(otu_tax) <- otu_tax$Row.names
  otu_tax$Row.names <- NULL
  
  # Aggregate by taxonomic level
  tax_level_data <- otu_tax %>%
    group_by(!!sym(tax_level)) %>%
    summarise(across(where(is.numeric), sum)) %>%
    filter(!is.na(!!sym(tax_level)) & !!sym(tax_level) != "")
  
  # Get top taxa based on total abundance across all samples
  total_abundance <- rowSums(select(tax_level_data, where(is.numeric)))
  tax_level_data$TotalAbundance <- total_abundance
  
  top_taxa <- tax_level_data %>%
    arrange(desc(TotalAbundance)) %>%
    head(top_n) %>%
    pull(!!sym(tax_level))
  
  # Create a copy for reshaping
  tax_data_for_reshape <- tax_level_data %>%
    select(-TotalAbundance) %>%
    # Convert to long format FIRST
    pivot_longer(cols = -!!sym(tax_level), names_to = "SampleID", values_to = "Abundance")
  
  # Now group low-abundance taxa into "Other" category PER SAMPLE
  # This ensures each sample still sums to 1.0
  plot_data <- tax_data_for_reshape %>%
    mutate(TaxaGroup = ifelse(!!sym(tax_level) %in% top_taxa, !!sym(tax_level), "Other")) %>%
    group_by(SampleID, TaxaGroup) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # Verify that each sample sums to 1.0
  sample_sums_check <- plot_data %>%
    group_by(SampleID) %>%
    summarise(SampleSum = sum(Abundance), .groups = "drop")
  
  cat("Sample sums after grouping (should all be ~1.0):\n")
  cat("Min:", min(sample_sums_check$SampleSum), "Max:", max(sample_sums_check$SampleSum), "\n")
  
  # Add sample metadata
  plot_data <- merge(plot_data, sample_data, by.x = "SampleID", by.y = "row.names")
  
  return(plot_data)
}

# ==================================================================================
# Create faceted bar plot visualization for phyla
# ==================================================================================
create_phylum_barplot <- function(ps_obj, top_n = 10) {
  tax_level <- "Phylum"
  
  # Prepare data
  plot_data <- prepare_abundance_data(ps_obj, tax_level, top_n)
  
  # Calculate mean abundance by taxa and stage
  bar_data <- plot_data %>%
    group_by(Stage, TaxaGroup) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = "drop") %>%
    # Renormalize so each stage sums to 1.0
    group_by(Stage) %>%
    mutate(GroupSum = sum(MeanAbundance),
           NormalizedAbundance = MeanAbundance / GroupSum) %>%
    ungroup()
  
  # Verify that each group now sums to 1.0
  verification <- bar_data %>%
    group_by(Stage) %>%
    summarise(GroupSum = sum(NormalizedAbundance), .groups = "drop")
  
  cat("Group sums after normalization (should all be 1.0):\n")
  cat("Min:", min(verification$GroupSum), "Max:", max(verification$GroupSum), "\n")
  
  # Order taxa by overall abundance
  taxa_order <- bar_data %>%
    group_by(TaxaGroup) %>%
    summarise(TotalAbundance = sum(NormalizedAbundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalAbundance)) %>%
    pull(TaxaGroup)
  
  # Replace underscores with spaces in phylum names
  bar_data$TaxaGroup <- gsub("_", " ", bar_data$TaxaGroup)
  taxa_order <- gsub("_", " ", taxa_order)
  
  # Make sure "Other" is at the end of the legend
  if("Other" %in% taxa_order) {
    taxa_order <- c(taxa_order[taxa_order != "Other"], "Other")
  }
  
  # Set factor levels for ordering
  bar_data$TaxaGroup <- factor(bar_data$TaxaGroup, levels = rev(taxa_order))
  
  # Use custom color palette with "Other" as light grey
  custom_colors <- c("#0ea5e9", "#06b6d4", "#14b8a6", "#10b981", "#22c55e", 
                    "#84cc16", "#eab308", "#f59e0b", "#f97316", "#ef4444")
  
  # Get unique taxa groups excluding "Other"
  taxa_groups <- levels(bar_data$TaxaGroup)
  other_index <- which(taxa_groups == "Other")
  
  # Create the final color vector with "Other" as light grey
  if (length(other_index) > 0) {
    # If "Other" is present, assign colors to all taxa except "Other"
    taxa_without_other <- taxa_groups[-other_index]
    n_taxa <- length(taxa_without_other)
    
    # Assign colors to taxa (excluding "Other")
    colors <- c()
    for (i in 1:length(taxa_groups)) {
      if (taxa_groups[i] == "Other") {
        colors <- c(colors, "#e5e5e5")  # Light grey for "Other"
      } else {
        # Find position in taxa_without_other
        pos <- which(taxa_without_other == taxa_groups[i])
        # Use custom color if available, otherwise use default
        if (pos <= length(custom_colors)) {
          colors <- c(colors, custom_colors[pos])
        } else {
          colors <- c(colors, "#999999")  # Fallback color
        }
      }
    }
  } else {
    # If "Other" is not present, just use the custom colors directly
    n_taxa <- length(taxa_groups)
    colors <- custom_colors[1:min(n_taxa, length(custom_colors))]
  }
  
  # Create bar plot with faceting by Stage
  p <- ggplot(bar_data, aes(x = Stage, y = NormalizedAbundance, fill = TaxaGroup)) +
    geom_col(position = "stack", color = "white", size = 0.1) +
    scale_fill_manual(values = colors, name = "Phylum") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = "Relative Abundance of Phyla in Swab Samples by Stage",
      y = "Mean Relative Abundance",
      x = "Stage"
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
# Create faceted bar plot visualization for genera
# ==================================================================================
create_genus_barplot <- function(ps_obj, top_n = 10) {
  tax_level <- "Genus"
  
  # Prepare data
  plot_data <- prepare_abundance_data(ps_obj, tax_level, top_n)
  
  # Calculate mean abundance by taxa and stage
  bar_data <- plot_data %>%
    group_by(Stage, TaxaGroup) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = "drop") %>%
    # Renormalize so each stage sums to 1.0
    group_by(Stage) %>%
    mutate(GroupSum = sum(MeanAbundance),
           NormalizedAbundance = MeanAbundance / GroupSum) %>%
    ungroup()
  
  # Verify that each group now sums to 1.0
  verification <- bar_data %>%
    group_by(Stage) %>%
    summarise(GroupSum = sum(NormalizedAbundance), .groups = "drop")
  
  cat("Group sums after normalization (should all be 1.0):\n")
  cat("Min:", min(verification$GroupSum), "Max:", max(verification$GroupSum), "\n")
  
  # Order taxa by overall abundance
  taxa_order <- bar_data %>%
    group_by(TaxaGroup) %>%
    summarise(TotalAbundance = sum(NormalizedAbundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalAbundance)) %>%
    pull(TaxaGroup)
  
  # Replace underscores with spaces in genus names
  bar_data$TaxaGroup <- gsub("_", " ", bar_data$TaxaGroup)
  taxa_order <- gsub("_", " ", taxa_order)
  
  # Set factor levels for ordering
  bar_data$TaxaGroup <- factor(bar_data$TaxaGroup, levels = rev(taxa_order))
  
  # Use custom color palette with "Other" as light grey
  custom_colors <- c("#0ea5e9", "#06b6d4", "#14b8a6", "#10b981", "#22c55e", 
                    "#84cc16", "#eab308", "#f59e0b", "#f97316", "#ef4444")
  
  # Get unique taxa groups excluding "Other"
  taxa_groups <- levels(bar_data$TaxaGroup)
  other_index <- which(taxa_groups == "Other")
  
  # Create the final color vector with "Other" as light grey
  if (length(other_index) > 0) {
    # If "Other" is present, assign colors to all taxa except "Other"
    taxa_without_other <- taxa_groups[-other_index]
    n_taxa <- length(taxa_without_other)
    
    # Assign colors to taxa (excluding "Other")
    colors <- c()
    for (i in 1:length(taxa_groups)) {
      if (taxa_groups[i] == "Other") {
        colors <- c(colors, "#e5e5e5")  # Light grey for "Other"
      } else {
        # Find position in taxa_without_other
        pos <- which(taxa_without_other == taxa_groups[i])
        # Use custom color if available, otherwise use default
        if (pos <= length(custom_colors)) {
          colors <- c(colors, custom_colors[pos])
        } else {
          colors <- c(colors, "#999999")  # Fallback color
        }
      }
    }
  } else {
    # If "Other" is not present, just use the custom colors directly
    n_taxa <- length(taxa_groups)
    colors <- custom_colors[1:min(n_taxa, length(custom_colors))]
  }
  
  # Create bar plot with faceting by Stage
  p <- ggplot(bar_data, aes(x = Stage, y = NormalizedAbundance, fill = TaxaGroup)) +
    geom_col(position = "stack", color = "white", size = 0.1) +
    scale_fill_manual(values = colors, name = "Genus") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = "Relative Abundance of Genera in Swab Samples by Stage",
      y = "Mean Relative Abundance",
      x = "Stage"
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
# Create horizontal bar plots (with coord_flip)
# ==================================================================================
create_horizontal_phylum_barplot <- function(ps_obj, top_n = 10) {
  tax_level <- "Phylum"
  
  # Prepare data
  plot_data <- prepare_abundance_data(ps_obj, tax_level, top_n)
  
  # Calculate mean abundance by taxa and stage
  bar_data <- plot_data %>%
    group_by(Stage, TaxaGroup) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = "drop") %>%
    # Renormalize so each stage sums to 1.0
    group_by(Stage) %>%
    mutate(GroupSum = sum(MeanAbundance),
           NormalizedAbundance = MeanAbundance / GroupSum) %>%
    ungroup()
  
  # Order taxa by overall abundance
  taxa_order <- bar_data %>%
    group_by(TaxaGroup) %>%
    summarise(TotalAbundance = sum(NormalizedAbundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalAbundance)) %>%
    pull(TaxaGroup)
  
  # Replace underscores with spaces in phylum names
  bar_data$TaxaGroup <- gsub("_", " ", bar_data$TaxaGroup)
  taxa_order <- gsub("_", " ", taxa_order)
  
  # Make sure "Other" is at the end of the legend
  if("Other" %in% taxa_order) {
    taxa_order <- c(taxa_order[taxa_order != "Other"], "Other")
  }
  
  # Set factor levels for ordering
  bar_data$TaxaGroup <- factor(bar_data$TaxaGroup, levels = rev(taxa_order))
  
  # Use custom color palette with "Other" as light grey
  custom_colors <- c("#0ea5e9", "#06b6d4", "#14b8a6", "#10b981", "#22c55e", 
                    "#84cc16", "#eab308", "#f59e0b", "#f97316", "#ef4444")
  
  # Get unique taxa groups excluding "Other"
  taxa_groups <- levels(bar_data$TaxaGroup)
  other_index <- which(taxa_groups == "Other")
  
  # Create the final color vector with "Other" as light grey
  if (length(other_index) > 0) {
    # If "Other" is present, assign colors to all taxa except "Other"
    taxa_without_other <- taxa_groups[-other_index]
    n_taxa <- length(taxa_without_other)
    
    # Assign colors to taxa (excluding "Other")
    colors <- c()
    for (i in 1:length(taxa_groups)) {
      if (taxa_groups[i] == "Other") {
        colors <- c(colors, "#e5e5e5")  # Light grey for "Other"
      } else {
        # Find position in taxa_without_other
        pos <- which(taxa_without_other == taxa_groups[i])
        # Use custom color if available, otherwise use default
        if (pos <= length(custom_colors)) {
          colors <- c(colors, custom_colors[pos])
        } else {
          colors <- c(colors, "#999999")  # Fallback color
        }
      }
    }
  } else {
    # If "Other" is not present, just use the custom colors directly
    n_taxa <- length(taxa_groups)
    colors <- custom_colors[1:min(n_taxa, length(custom_colors))]
  }
  
  # Create horizontal bar plot
  p <- ggplot(bar_data, aes(x = Stage, y = NormalizedAbundance, fill = TaxaGroup)) +
    geom_col(position = "stack", color = "white", size = 0.1) +
    scale_fill_manual(values = colors, name = "Phylum") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
    coord_flip() +  # Flip coordinates to make horizontal bars
    labs(
      title = "Relative Abundance of Phyla in Swab Samples by Stage",
      y = "Mean Relative Abundance",
      x = "Stage"
    ) +
    theme_pub(legend_pos = "bottom") +
    theme(
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      legend.direction = "horizontal",
      legend.box = "horizontal"
    ) +
    guides(fill = guide_legend(nrow = 2))  # Arrange legend in 2 rows
  
  return(p)
}

create_horizontal_genus_barplot <- function(ps_obj, top_n = 10) {
  tax_level <- "Genus"
  
  # Prepare data
  plot_data <- prepare_abundance_data(ps_obj, tax_level, top_n)
  
  # Calculate mean abundance by taxa and stage
  bar_data <- plot_data %>%
    group_by(Stage, TaxaGroup) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = "drop") %>%
    # Renormalize so each stage sums to 1.0
    group_by(Stage) %>%
    mutate(GroupSum = sum(MeanAbundance),
           NormalizedAbundance = MeanAbundance / GroupSum) %>%
    ungroup()
  
  # Order taxa by overall abundance
  taxa_order <- bar_data %>%
    group_by(TaxaGroup) %>%
    summarise(TotalAbundance = sum(NormalizedAbundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalAbundance)) %>%
    pull(TaxaGroup)
  
  # Replace underscores with spaces in genus names
  bar_data$TaxaGroup <- gsub("_", " ", bar_data$TaxaGroup)
  taxa_order <- gsub("_", " ", taxa_order)
  
  # Set factor levels for ordering
  bar_data$TaxaGroup <- factor(bar_data$TaxaGroup, levels = rev(taxa_order))
  
  # Use custom color palette with "Other" as light grey
  custom_colors <- c("#0ea5e9", "#06b6d4", "#14b8a6", "#10b981", "#22c55e", 
                    "#84cc16", "#eab308", "#f59e0b", "#f97316", "#ef4444")
  
  # Get unique taxa groups excluding "Other"
  taxa_groups <- levels(bar_data$TaxaGroup)
  other_index <- which(taxa_groups == "Other")
  
  # Create the final color vector with "Other" as light grey
  if (length(other_index) > 0) {
    # If "Other" is present, assign colors to all taxa except "Other"
    taxa_without_other <- taxa_groups[-other_index]
    n_taxa <- length(taxa_without_other)
    
    # Assign colors to taxa (excluding "Other")
    colors <- c()
    for (i in 1:length(taxa_groups)) {
      if (taxa_groups[i] == "Other") {
        colors <- c(colors, "#e5e5e5")  # Light grey for "Other"
      } else {
        # Find position in taxa_without_other
        pos <- which(taxa_without_other == taxa_groups[i])
        # Use custom color if available, otherwise use default
        if (pos <= length(custom_colors)) {
          colors <- c(colors, custom_colors[pos])
        } else {
          colors <- c(colors, "#999999")  # Fallback color
        }
      }
    }
  } else {
    # If "Other" is not present, just use the custom colors directly
    n_taxa <- length(taxa_groups)
    colors <- custom_colors[1:min(n_taxa, length(custom_colors))]
  }
  
  # Create horizontal bar plot
  p <- ggplot(bar_data, aes(x = Stage, y = NormalizedAbundance, fill = TaxaGroup)) +
    geom_col(position = "stack", color = "white", size = 0.1) +
    scale_fill_manual(values = colors, name = "Genus") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
    coord_flip() +  # Flip coordinates to make horizontal bars
    labs(
      title = "Relative Abundance of Genera in Swab Samples by Stage",
      y = "Mean Relative Abundance",
      x = "Stage"
    ) +
    theme_pub(legend_pos = "bottom") +
    theme(
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      legend.direction = "horizontal",
      legend.box = "horizontal"
    ) +
    guides(fill = guide_legend(nrow = 2))  # Arrange legend in 2 rows
  
  return(p)
}

# ==================================================================================
# Function to combine phylum and genus plots
# ==================================================================================

# Function to create a combined plot with phylum and genus
create_combined_plot <- function(phylum_plot, genus_plot, layout = "vertical") {
  # Remove titles from individual plots for better combined appearance
  phylum_plot <- phylum_plot + 
    labs(title = NULL) +
    theme(plot.margin = margin(5, 5, 10, 5))
  
  genus_plot <- genus_plot + 
    labs(title = NULL) +
    theme(plot.margin = margin(5, 5, 10, 5))
  
  # Combine plots based on layout without title
  if (layout == "vertical") {
    # Stack plots vertically
    combined_plot <- phylum_plot / genus_plot +
      plot_layout(heights = c(1, 1))
  } else {
    # Arrange plots side by side
    combined_plot <- phylum_plot + genus_plot +
      plot_layout(widths = c(1, 1))
  }
  
  return(combined_plot)
}

# ==================================================================================
# Generate and save visualizations
# ==================================================================================

# Create vertical bar plots
phylum_plot <- create_phylum_barplot(ps_rel, top_n = 10)
genus_plot <- create_genus_barplot(ps_rel, top_n = 10)

# Create horizontal bar plots
phylum_plot_horizontal <- create_horizontal_phylum_barplot(ps_rel, top_n = 10)
genus_plot_horizontal <- create_horizontal_genus_barplot(ps_rel, top_n = 10)

# Create combined plots
combined_vertical <- create_combined_plot(phylum_plot, genus_plot, "vertical")
combined_horizontal <- create_combined_plot(phylum_plot_horizontal, genus_plot_horizontal, "horizontal")

# Save individual plots
ggsave("swab_phylum_abundance_by_stage.png", phylum_plot, width = 10, height = 8, dpi = 300)
ggsave("swab_genus_abundance_by_stage.png", genus_plot, width = 10, height = 8, dpi = 300)
ggsave("swab_phylum_abundance_by_stage_horizontal.png", phylum_plot_horizontal, width = 10, height = 6, dpi = 300)
ggsave("swab_genus_abundance_by_stage_horizontal.png", genus_plot_horizontal, width = 10, height = 6, dpi = 300)

# Save combined plots
ggsave("swab_combined_vertical.png", combined_vertical, width = 10, height = 16, dpi = 300)
ggsave("swab_combined_horizontal.png", combined_horizontal, width = 20, height = 8, dpi = 300)

# Display the combined plots
print(combined_vertical)
print(combined_horizontal)

# Print summary message
cat("\nRelative abundance visualizations created:\n")
cat("1. Phylum-level bar plot by stage (swab_phylum_abundance_by_stage.png)\n")
cat("2. Genus-level bar plot by stage (swab_genus_abundance_by_stage.png)\n")
cat("3. Horizontal phylum-level bar plot by stage (swab_phylum_abundance_by_stage_horizontal.png)\n")
cat("4. Horizontal genus-level bar plot by stage (swab_genus_abundance_by_stage_horizontal.png)\n")
cat("5. Combined vertical plot with phylum and genus (swab_combined_vertical.png)\n")
cat("6. Combined horizontal plot with phylum and genus (swab_combined_horizontal.png)\n")
cat("\nAll bars in the stacked plots sum to 1.0 (100%) per stage!\n")

```

