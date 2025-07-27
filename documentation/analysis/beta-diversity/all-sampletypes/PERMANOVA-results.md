```R
=== PERMANOVA RESULTS FOR BRAY-CURTIS DISTANCE ===
> permanova_bray <- perform_permanova(dist_bray, metadata)
PERMANOVA for SampleType:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$SampleType, permutations = 999, by = "terms")
                        Df SumOfSqs      R2
metadata_df$SampleType   2   12.082 0.07101
Residual               415  158.056 0.92899
Total                  417  170.138 1.00000
                            F Pr(>F)    
metadata_df$SampleType 15.861  0.001 ***
Residual                                
Total                                   
---
Signif. codes:  
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1
  ‘ ’ 1

PERMANOVA for Microbiota:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$Microbiota, permutations = 999, by = "terms")
                        Df SumOfSqs      R2
metadata_df$Microbiota   2   11.138 0.06547
Residual               415  158.999 0.93453
Total                  417  170.138 1.00000
                            F Pr(>F)    
metadata_df$Microbiota 14.536  0.001 ***
Residual                                
Total                                   
---
Signif. codes:  
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1
  ‘ ’ 1

PERMANOVA for SampleType * Microbiota interaction:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$SampleType * metadata_df$Microbiota, permutations = 999, by = "terms")
                        Df SumOfSqs      R2
metadata_df$SampleType   2   12.082 0.07101
metadata_df$Microbiota   1    1.751 0.01029
Residual               414  156.305 0.91870
Total                  417  170.138 1.00000
                             F Pr(>F)    
metadata_df$SampleType 15.9999  0.001 ***
metadata_df$Microbiota  4.6367  0.001 ***
Residual                                 
Total                                    
---
Signif. codes:  
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1
  ‘ ’ 1
```

```R
=== PERMANOVA RESULTS FOR UNWEIGHTED UNIFRAC DISTANCE ===
> permanova_unifrac <- perform_permanova(dist_unifrac, metadata)
PERMANOVA for SampleType:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$SampleType, permutations = 999, by = "terms")
                        Df SumOfSqs      R2
metadata_df$SampleType   2    7.645 0.05649
Residual               415  127.702 0.94351
Total                  417  135.347 1.00000
                            F Pr(>F)    
metadata_df$SampleType 12.422  0.001 ***
Residual                                
Total                                   
---
Signif. codes:  
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1
  ‘ ’ 1

PERMANOVA for Microbiota:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$Microbiota, permutations = 999, by = "terms")
                        Df SumOfSqs      R2
metadata_df$Microbiota   2    5.911 0.04367
Residual               415  129.437 0.95633
Total                  417  135.347 1.00000
                            F Pr(>F)    
metadata_df$Microbiota 9.4753  0.001 ***
Residual                                
Total                                   
---
Signif. codes:  
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1
  ‘ ’ 1

PERMANOVA for SampleType * Microbiota interaction:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$SampleType * metadata_df$Microbiota, permutations = 999, by = "terms")
                        Df SumOfSqs      R2
metadata_df$SampleType   2    7.645 0.05649
metadata_df$Microbiota   1    0.748 0.00553
Residual               414  126.954 0.93799
Total                  417  135.347 1.00000
                             F Pr(>F)    
metadata_df$SampleType 12.4655  0.001 ***
metadata_df$Microbiota  2.4389  0.001 ***
Residual                                 
Total                                    
---
Signif. codes:  
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1
  ‘ ’ 1
```

```R
=== PERMANOVA RESULTS FOR WEIGHTED UNIFRAC DISTANCE ===
> permanova_wunifrac <- perform_permanova(dist_wunifrac, metadata)
PERMANOVA for SampleType:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$SampleType, permutations = 999, by = "terms")
                        Df SumOfSqs      R2
metadata_df$SampleType   2   0.2982 0.08049
Residual               415   3.4069 0.91951
Total                  417   3.7052 1.00000
                            F Pr(>F)    
metadata_df$SampleType 18.163  0.001 ***
Residual                                
Total                                   
---
Signif. codes:  
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1
  ‘ ’ 1

PERMANOVA for Microbiota:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$Microbiota, permutations = 999, by = "terms")
                        Df SumOfSqs      R2
metadata_df$Microbiota   2   0.4080 0.11011
Residual               415   3.2972 0.88989
Total                  417   3.7052 1.00000
                            F Pr(>F)    
metadata_df$Microbiota 25.676  0.001 ***
Residual                                
Total                                   
---
Signif. codes:  
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1
  ‘ ’ 1

PERMANOVA for SampleType * Microbiota interaction:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$SampleType * metadata_df$Microbiota, permutations = 999, by = "terms")
                        Df SumOfSqs      R2
metadata_df$SampleType   2   0.2982 0.08049
metadata_df$Microbiota   1   0.2011 0.05428
Residual               414   3.2058 0.86523
Total                  417   3.7052 1.00000
                            F Pr(>F)    
metadata_df$SampleType 19.256  0.001 ***
metadata_df$Microbiota 25.972  0.001 ***
Residual                                
Total                                   
---
Signif. codes:  
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1
  ‘ ’ 1
```

```R
=== PAIRWISE PERMANOVA RESULTS FOR BRAY-CURTIS DISTANCE ===
> pairwise_bray <- perform_pairwise_permanova(ps, "bray")
Sample types for pairwise PERMANOVA: Insect, Soil, Swab 
Performing pairwise PERMANOVA tests...
> print(pairwise_bray)
  SampleType1 SampleType2         R2 p_value
1      Insect        Soil 0.07417324   0.001
2      Insect        Swab 0.05918404   0.001
3        Soil        Swab 0.02624699   0.001
  p_adjusted significance
1      0.003           **
2      0.003           **
3      0.003           **
```

```R
=== PAIRWISE PERMANOVA RESULTS FOR UNWEIGHTED UNIFRAC DISTANCE ===
> pairwise_unifrac <- perform_pairwise_permanova(ps, "unifrac")
Sample types for pairwise PERMANOVA: Insect, Soil, Swab 
Performing pairwise PERMANOVA tests...
> print(pairwise_unifrac)
  SampleType1 SampleType2         R2 p_value
1      Insect        Soil 0.06573957   0.001
2      Insect        Swab 0.03322584   0.001
3        Soil        Swab 0.03042171   0.001
  p_adjusted significance
1      0.003           **
2      0.003           **
3      0.003           **
```

```R
=== PAIRWISE PERMANOVA RESULTS FOR WEIGHTED UNIFRAC DISTANCE ===
> pairwise_wunifrac <- perform_pairwise_permanova(ps, "wunifrac")
Sample types for pairwise PERMANOVA: Insect, Soil, Swab 
Performing pairwise PERMANOVA tests...
> print(pairwise_wunifrac)
  SampleType1 SampleType2         R2 p_value
1      Insect        Soil 0.04646537   0.001
2      Insect        Swab 0.06714453   0.001
3        Soil        Swab 0.09704759   0.001
  p_adjusted significance
1      0.003           **
2      0.003           **
3      0.003           **
```

