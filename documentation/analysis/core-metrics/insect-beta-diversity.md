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

## Calculate diversity metrics and merge

```         R
# Shannon, Observed ASVs, and Chao1
alpha_div <- estimate_richness(ps_insects, measures = c("Shannon", "Observed", "Chao1"))

# Add sample IDs to alpha diversity data
alpha_div$SampleID <- rownames(alpha_div)

# Get sample metadata - using a safer approach
if (nsamples(ps_insects) > 0) {
  # Extract metadata directly from the phyloseq object
  sample_data_df <- data.frame(sample_data(ps_insects))
  
  # Check if sample_data_df has rows
  if (nrow(sample_data_df) > 0) {
    sample_data_df$SampleID <- rownames(sample_data_df)
    
    # Merge alpha diversity with metadata
    alpha_div_meta <- merge(alpha_div, sample_data_df, by = "SampleID")
    
    # Set custom order for Source variable
    source_order <- c("Chrysomya rufifacies", "Chrysomya megacephala", 
                      "Cochliomyia macellaria", "Lucilia coeruleiviridis", 
                      "Hydrotea aenescens", "Necrodes surinamensis")
    
    # Convert Source to factor with custom order
    alpha_div_meta$Source <- factor(alpha_div_meta$Source, levels = source_order)
    
    # Print the first few rows to verify the merge worked
    print("First few rows of merged data:")
    print(head(alpha_div_meta))
  } else {
    stop("Sample data frame has zero rows after extraction")
  }
} else {
  stop("No samples found in the filtered phyloseq object")
}
```

## Statistical testing

```         R
# Shannon diversity
shannon_stat <- alpha_div_meta %>%
  pairwise_t_test(Shannon ~ Source, p.adjust.method = "bonferroni") %>%
  filter(p.adj < 0.05)  # Keep only significant comparisons

# Shannon diversity by Microbiota
shannon_stat_by_microbiota <- alpha_div_meta %>%
  group_by(Microbiota) %>%
  pairwise_t_test(Shannon ~ Source, p.adjust.method = "bonferroni") %>%
  filter(p.adj < 0.05)  # Keep only significant comparisons

# Chao1
chao1_stat <- alpha_div_meta %>%
  pairwise_t_test(Chao1 ~ Source, p.adjust.method = "bonferroni") %>%
  filter(p.adj < 0.05)  # Keep only significant comparisons

# Chao1 by Microbiota
chao1_stat_by_microbiota <- alpha_div_meta %>%
  group_by(Microbiota) %>%
  pairwise_t_test(Chao1 ~ Source, p.adjust.method = "bonferroni") %>%
  filter(p.adj < 0.05)  # Keep only significant comparisons

# Observed ASVs
observed_stat <- alpha_div_meta %>%
  pairwise_t_test(Observed ~ Source, p.adjust.method = "bonferroni") %>%
  filter(p.adj < 0.05)  # Keep only significant comparisons

# Observed ASVs by Microbiota
observed_stat_by_microbiota <- alpha_div_meta %>%
  group_by(Microbiota) %>%
  pairwise_t_test(Observed ~ Source, p.adjust.method = "bonferroni") %>%
  filter(p.adj < 0.05)  # Keep only significant comparisons

# Print significant comparisons
cat("Significant pairwise comparisons for Shannon diversity:\n")
print(shannon_stat)

cat("\nSignificant pairwise comparisons for Shannon diversity by Microbiota:\n")
print(shannon_stat_by_microbiota)

cat("\nSignificant pairwise comparisons for Chao1:\n")
print(chao1_stat)

cat("\nSignificant pairwise comparisons for Chao1 by Microbiota:\n")
print(chao1_stat_by_microbiota)

cat("\nSignificant pairwise comparisons for Observed ASVs:\n")
print(observed_stat)

cat("\nSignificant pairwise comparisons for Observed ASVs by Microbiota:\n")
print(observed_stat_by_microbiota)

# Statistical tests for differences between sources
# ANOVA for each metric
shannon_anova <- aov(Shannon ~ Source, data = alpha_div_meta)
chao1_anova <- aov(Chao1 ~ Source, data = alpha_div_meta)
observed_anova <- aov(Observed ~ Source, data = alpha_div_meta)

# Print ANOVA summaries
cat("Shannon Diversity ANOVA:\n")
print(summary(shannon_anova))

cat("\nChao1 ANOVA:\n")
print(summary(chao1_anova))

cat("\nObserved ASVs ANOVA:\n")
print(summary(observed_anova))

# If ANOVA is significant, perform post-hoc Tukey test
if(summary(shannon_anova)[[1]][["Pr(>F)"]][1] < 0.05) {
  cat("\nTukey HSD for Shannon Diversity:\n")
  print(TukeyHSD(shannon_anova))
}

if(summary(chao1_anova)[[1]][["Pr(>F)"]][1] < 0.05) {
  cat("\nTukey HSD for Chao1:\n")
  print(TukeyHSD(chao1_anova))
}

if(summary(observed_anova)[[1]][["Pr(>F)"]][1] < 0.05) {
  cat("\nTukey HSD for Observed ASVs:\n")
  print(TukeyHSD(observed_anova))
}
```

## Plotting

```         R
# Shannon diversity
p1 <- ggplot(alpha_div_meta, aes(x = Source, y = Shannon, fill = Source)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(x = "Insect Species",
       y = "Shannon Diversity Index") +
  scale_fill_viridis_d() +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  # Remove legend from individual plot

# Add significance bars if there are significant comparisons
if(nrow(shannon_stat) > 0) {
  # Prepare the comparisons for ggpubr
  # We need to modify the format to work with stat_compare_means
  shannon_comparisons <- lapply(1:nrow(shannon_stat), function(i) {
    c(as.character(shannon_stat$group1[i]), as.character(shannon_stat$group2[i]))
  })
  
  # Add significance bars using stat_compare_means
  p1 <- p1 + stat_compare_means(
    comparisons = shannon_comparisons,
    method = "t.test",
    label = "p.signif"
  )
}

# Shannon diversity faceted by Microbiota
p1_facet <- ggplot(alpha_div_meta, aes(x = Source, y = Shannon, fill = Source)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~ Microbiota) +
  labs(x = "Insect Species",
       y = "Shannon Diversity Index") +
  scale_fill_viridis_d() +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  # Remove legend from individual plot

# Add significance bars if there are significant comparisons
if(nrow(shannon_stat_by_microbiota) > 0) {
  # For each Microbiota group, prepare the comparisons
  for (mb in unique(shannon_stat_by_microbiota$Microbiota)) {
    # Filter for this Microbiota
    mb_stats <- shannon_stat_by_microbiota %>% filter(Microbiota == mb)
    
    if(nrow(mb_stats) > 0) {
      # Prepare the comparisons for this Microbiota
      mb_comparisons <- lapply(1:nrow(mb_stats), function(i) {
        c(as.character(mb_stats$group1[i]), as.character(mb_stats$group2[i]))
      })
      
      # Add significance bars for this Microbiota
      p1_facet <- p1_facet + 
        stat_compare_means(
          comparisons = mb_comparisons,
          method = "t.test",
          label = "p.signif",
          data = alpha_div_meta %>% filter(Microbiota == mb)
        )
    }
  }
}

# Chao1
p3 <- ggplot(alpha_div_meta, aes(x = Source, y = Chao1, fill = Source)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(x = "Insect Species",
       y = "Chao1 Richness") +
  scale_fill_viridis_d() +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  # Remove legend from individual plot

# Add significance bars if there are significant comparisons
if(nrow(chao1_stat) > 0) {
  # Prepare the comparisons for ggpubr
  chao1_comparisons <- lapply(1:nrow(chao1_stat), function(i) {
    c(as.character(chao1_stat$group1[i]), as.character(chao1_stat$group2[i]))
  })
  
  # Add significance bars using stat_compare_means
  p3 <- p3 + stat_compare_means(
    comparisons = chao1_comparisons,
    method = "t.test",
    label = "p.signif"
  )
}

# Chao1 faceted by Microbiota
p3_facet <- ggplot(alpha_div_meta, aes(x = Source, y = Chao1, fill = Source)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~ Microbiota) +
  labs(x = "Insect Species",
       y = "Chao1 Richness") +
  scale_fill_viridis_d() +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  # Remove legend from individual plot

# Add significance bars if there are significant comparisons
if(nrow(chao1_stat_by_microbiota) > 0) {
  # For each Microbiota group, prepare the comparisons
  for (mb in unique(chao1_stat_by_microbiota$Microbiota)) {
    # Filter for this Microbiota
    mb_stats <- chao1_stat_by_microbiota %>% filter(Microbiota == mb)
    
    if(nrow(mb_stats) > 0) {
      # Prepare the comparisons for this Microbiota
      mb_comparisons <- lapply(1:nrow(mb_stats), function(i) {
        c(as.character(mb_stats$group1[i]), as.character(mb_stats$group2[i]))
      })
      
      # Add significance bars for this Microbiota
      p3_facet <- p3_facet + 
        stat_compare_means(
          comparisons = mb_comparisons,
          method = "t.test",
          label = "p.signif",
          data = alpha_div_meta %>% filter(Microbiota == mb)
        )
    }
  }
}

# Observed ASVs
p4 <- ggplot(alpha_div_meta, aes(x = Source, y = Observed, fill = Source)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(x = "Insect Species",
       y = "Observed ASVs") +
  scale_fill_viridis_d() +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  # Remove legend from individual plot

# Add significance bars if there are significant comparisons
if(nrow(observed_stat) > 0) {
  # Prepare the comparisons for ggpubr
  observed_comparisons <- lapply(1:nrow(observed_stat), function(i) {
    c(as.character(observed_stat$group1[i]), as.character(observed_stat$group2[i]))
  })
  
  # Add significance bars using stat_compare_means
  p4 <- p4 + stat_compare_means(
    comparisons = observed_comparisons,
    method = "t.test",
    label = "p.signif"
  )
}

# Observed ASVs faceted by Microbiota
p4_facet <- ggplot(alpha_div_meta, aes(x = Source, y = Observed, fill = Source)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~ Microbiota) +
  labs(x = "Insect Species",
       y = "Observed ASVs") +
  scale_fill_viridis_d() +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  # Remove legend from individual plot

# Add significance bars if there are significant comparisons
if(nrow(observed_stat_by_microbiota) > 0) {
  # For each Microbiota group, prepare the comparisons
  for (mb in unique(observed_stat_by_microbiota$Microbiota)) {
    # Filter for this Microbiota
    mb_stats <- observed_stat_by_microbiota %>% filter(Microbiota == mb)
    
    if(nrow(mb_stats) > 0) {
      # Prepare the comparisons for this Microbiota
      mb_comparisons <- lapply(1:nrow(mb_stats), function(i) {
        c(as.character(mb_stats$group1[i]), as.character(mb_stats$group2[i]))
      })
      
      # Add significance bars for this Microbiota
      p4_facet <- p4_facet + 
        stat_compare_means(
          comparisons = mb_comparisons,
          method = "t.test",
          label = "p.signif",
          data = alpha_div_meta %>% filter(Microbiota == mb)
        )
    }
  }
}
```

## Combining visualizations

```         R
# Shannon diversity combined plot (main + faceted)
shannon_combined <- p1 | p1_facet

# Chao1 combined plot (main + faceted)
chao1_combined <- p3 | p3_facet

# Observed ASVs combined plot (main + faceted)
observed_combined <- p4 | p4_facet

# Arrange all combined plots in a grid with a shared legend
all_combined_plots <- shannon_combined / chao1_combined / observed_combined +
  plot_layout(guides = "collect") +  # Collect all legends into one shared legend
  plot_annotation(title = "Alpha Diversity Metrics Across Insect Species and by Microbiota",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))) &
  theme(legend.direction = "vertical",
        legend.box = "vertical")

# Display the combined plot
all_combined_plots
```

## Exporting visualizations

```         R
# Save the combined plot
ggsave("alpha_diversity_combined_plots.pdf", all_combined_plots, width = 18, height = 12, units = "in")
ggsave("alpha_diversity_combined_plots.png", all_combined_plots, width = 18, height = 12, units = "in", dpi = 300)

# Also save individual metric plots for reference
# Shannon
ggsave("alpha_diversity_shannon.pdf", shannon_combined, width = 10, height = 12, units = "in")
ggsave("alpha_diversity_shannon.png", shannon_combined, width = 10, height = 12, units = "in", dpi = 300)

# Chao1
ggsave("alpha_diversity_chao1.pdf", chao1_combined, width = 10, height = 12, units = "in")
ggsave("alpha_diversity_chao1.png", chao1_combined, width = 10, height = 12, units = "in", dpi = 300)

# Observed ASVs
ggsave("alpha_diversity_observed.pdf", observed_combined, width = 10, height = 12, units = "in")
ggsave("alpha_diversity_observed.png", observed_combined, width = 10, height = 12, units = "in", dpi = 300)
```