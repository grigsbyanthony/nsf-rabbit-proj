## Agglomerate to genus level

```R
ps_genus <- tax_glom(ps_insects, taxrank = "Genus")
```



## Run ANCOM-BC for `Source` (insect species)

```R
ancombc_result <- ancombc2(
  data = ps_genus,
  assay_name = "counts",
  tax_level = "Genus",
  fix_formula = "Source",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.10,       
  lib_cut = 1000,       
  s0_perc = 0.05,
  group = "Source",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE
)
```



## Extract results and get signficant results

```R
results_df <- ancombc_results$res

print("ANCOM-BC Results Structure:")
print(paste("Number of genera tested:", nrow(results_df)))
print("Column names:")
print(colnames(results_df))

# Identify comparison groups (Source species)
comparison_cols <- colnames(results_df)[grepl("^diff_Source", colnames(results_df))]
print("Comparison groups found:")
print(comparison_cols)

# Create a summary of significant results for each comparison
all_sig_results <- data.frame()

for(comp_col in comparison_cols) {
  # Extract the source name from column name
  source_name <- gsub("diff_Source", "", comp_col)
  
  # Get corresponding columns for this comparison
  lfc_col <- paste0("lfc_Source", source_name)
  p_col <- paste0("p_Source", source_name)
  q_col <- paste0("q_Source", source_name)
  se_col <- paste0("se_Source", source_name)
  w_col <- paste0("W_Source", source_name)
  
  # Find significant results for this comparison
  sig_indices <- which(results_df[[comp_col]] == TRUE)
  
  if(length(sig_indices) > 0) {
    # Create summary for this comparison
    comp_results <- data.frame(
      taxon = results_df$taxon[sig_indices],
      comparison = source_name,
      lfc = results_df[[lfc_col]][sig_indices],
      se = results_df[[se_col]][sig_indices],
      W = results_df[[w_col]][sig_indices],
      p_val = results_df[[p_col]][sig_indices],
      q_val = results_df[[q_col]][sig_indices],
      diff_abn = results_df[[comp_col]][sig_indices]
    )
    
    all_sig_results <- rbind(all_sig_results, comp_results)
    
    print(paste("Significant genera for", source_name, "vs reference:", length(sig_indices)))
    print(comp_results)
  } else {
    print(paste("No significant genera found for", source_name, "vs reference"))
  }
}

print("Total significant results across all comparisons:")
print(nrow(all_sig_results))
```



## Visualize results

```R
if(nrow(all_sig_results) > 0) {
  print("All significant differentially abundant genera:")
  sig_summary <- all_sig_results %>%
    arrange(q_val)
  print(sig_summary)
  
  # Create visualizations
  
  # 1. Bar plot showing significant results by comparison
  if(nrow(all_sig_results) > 0) {
    comparison_plot <- all_sig_results %>%
      mutate(
        direction = ifelse(lfc > 0, "Enriched", "Depleted"),
        taxon_comparison = paste(taxon, "vs", comparison)
      ) %>%
      ggplot(aes(x = reorder(taxon_comparison, lfc), y = lfc, fill = direction)) +
      geom_col() +
      coord_flip() +
      facet_wrap(~comparison, scales = "free_y", ncol = 1) +
      scale_fill_manual(values = c("Enriched" = "steelblue", "Depleted" = "coral")) +
      labs(
        title = "Significant Differentially Abundant Genera",
        subtitle = "Each panel shows genera significantly different from reference",
        x = "Genus",
        y = "Log Fold Change",
        fill = "Direction"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        strip.text = element_text(size = 10, face = "bold")
      )
    
    print(comparison_plot)
  }
  
  # 2. Heatmap showing all comparisons for significant genera
  # Get all comparison data for significant genera
  sig_taxa <- unique(all_sig_results$taxon)
  
  if(length(sig_taxa) > 0) {
    heatmap_data <- data.frame()
    
    for(taxon in sig_taxa) {
      taxon_row <- results_df[results_df$taxon == taxon, ]
      
      for(comp_col in comparison_cols) {
        source_name <- gsub("diff_Source", "", comp_col)
        lfc_col <- paste0("lfc_Source", source_name)
        q_col <- paste0("q_Source", source_name)
        
        heatmap_data <- rbind(heatmap_data, data.frame(
          taxon = taxon,
          comparison = source_name,
          lfc = taxon_row[[lfc_col]],
          q_val = taxon_row[[q_col]],
          significant = taxon_row[[comp_col]]
        ))
      }
    }
    
    # Create heatmap
    heatmap_plot <- heatmap_data %>%
      mutate(
        significance = ifelse(significant, "Significant", "Not Significant"),
        lfc_capped = pmax(pmin(lfc, 3), -3)  # Cap extreme values for better visualization
      ) %>%
      ggplot(aes(x = comparison, y = taxon, fill = lfc_capped)) +
      geom_tile(color = "white", size = 0.5) +
      geom_point(aes(size = ifelse(significant, "Significant", "Not Significant")), 
                 color = "white", alpha = 0.8) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                          midpoint = 0, name = "Log FC") +
      scale_size_manual(values = c("Significant" = 2, "Not Significant" = 0), 
                       name = "Significance") +
      labs(
        title = "Log Fold Changes for Significant Genera",
        subtitle = "White dots indicate statistical significance (q < 0.05)",
        x = "Insect Species (vs Reference)",
        y = "Bacterial Genus"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10)
      )
    
    print(heatmap_plot)
  }
  
  print(sig_summary)
  
  # Create a volcano plot
  volcano_plot <- results_df %>%
    mutate(
      neg_log10_q = -log10(q_val),
      significant = ifelse(diff_abn == TRUE, "Significant", "Not Significant"),
      label = ifelse(diff_abn == TRUE, taxon, "")
    ) %>%
    ggplot(aes(x = lfc, y = neg_log10_q, color = significant)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_text(aes(label = label), hjust = 0, vjust = 0, size = 3, 
              check_overlap = TRUE, nudge_x = 0.02) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
    labs(
      title = "ANCOM-BC Volcano Plot",
      subtitle = "Differential abundance of bacterial genera between insect species",
      x = "Log Fold Change",
      y = "-log10(q-value)",
      color = "Status"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  print(volcano_plot)
  
  # Create a bar plot for significant genera
  if(nrow(sig_results) <= 20) {  # Only if not too many significant results
    bar_plot <- sig_results %>%
      arrange(desc(abs(lfc))) %>%
      slice_head(n = 15) %>%  # Top 15 by effect size
      mutate(
        taxon = reorder(taxon, lfc),
        direction = ifelse(lfc > 0, "Enriched", "Depleted")
      ) %>%
      ggplot(aes(x = taxon, y = lfc, fill = direction)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = c("Enriched" = "steelblue", "Depleted" = "coral")) +
      labs(
        title = "Top Differentially Abundant Genera",
        subtitle = "Based on log fold change from ANCOM-BC",
        x = "Genus",
        y = "Log Fold Change",
        fill = "Direction"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      )
    
    print(bar_plot)
  }
    
  } else {
    print("No significantly differentially abundant genera found.")
  }
} else {
  print("Could not find q-value column for significance testing")
  print("Available columns are:")
  print(colnames(actual_results))
}

# Save results to CSV
write.csv(actual_results, "ancombc_results_insect_microbiome.csv", row.names = FALSE)
```





## Alternatively running ANCOM-BC in QIIME2

```bash
qiime taxa collapse \
  --i-table merged-table-filtered-all.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table genus-level-table.qza
      
qiime composition ancombc \
  --i-table genus-level-table.qza \
  --m-metadata-file metadata.tsv \
  --p-formula "Source" \
  --p-reference-levels "Source::Necrodes surinamensis"\
  --p-alpha 0.05 \
  --p-p-adj-method "holm" \
  --p-prv-cut 0.1 \
  --p-lib-cut 1000 \
  --o-differentials ancombc-results.qza
  
qiime composition da-barplot \
	--i-data ancombc-results.qza \
	--p-significance-threshold 0.1 \
	--o-visualization ancombc-results.qzv

qiime tools export \
  --input-path ancombc-results.qza \
  --output-path exported-results
  
# WITHOUT reference group:
qiime composition ancombc \
  --i-table genus-level-table.qza \
  --m-metadata-file metadata.tsv \
  --p-formula "Source" \
  --p-alpha 0.05 \
  --p-p-adj-method "holm" \
  --p-prv-cut 0.1 \
  --p-lib-cut 1000 \
  --o-differentials ancombc-results.qza
  
qiime composition da-barplot \
	--i-data ancombc-results.qza \
	--p-significance-threshold 0.1 \
	--o-visualization ancombc-results.qzv
	
	
# By microbiota
qiime composition ancombc \
  --i-table genus-level-table.qza \
  --m-metadata-file metadata.tsv \
  --p-formula "Microbiota" \
  --p-alpha 0.05 \
  --p-p-adj-method "holm" \
  --p-prv-cut 0.1 \
  --p-lib-cut 1000 \
  --o-differentials ancombc-results.qza
  
qiime composition da-barplot \
	--i-data ancombc-results.qza \
	--p-significance-threshold 0.1 \
	--o-visualization ancombc-results.qzv
```

