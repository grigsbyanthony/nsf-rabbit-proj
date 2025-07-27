```R
=== PERMANOVA RESULTS FOR BRAY-CURTIS DISTANCE ===
> permanova_bray <- perform_permanova(dist_bray, metadata)
PERMANOVA for Source:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$Source, permutations = 999, by = "terms")
                    Df SumOfSqs      R2
metadata_df$Source   5   14.895 0.25628
Residual           140   43.225 0.74372
Total              145   58.120 1.00000
                        F Pr(>F)    
metadata_df$Source 9.6487  0.001 ***
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
metadata_df$Microbiota   1    1.751 0.03012
Residual               144   56.369 0.96988
Total                  145   58.120 1.00000
                            F Pr(>F)    
metadata_df$Microbiota 4.4721  0.001 ***
Residual                                
Total                                   
---
Signif. codes:  
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1
  ‘ ’ 1

PERMANOVA for Source * Microbiota interaction:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$Source * metadata_df$Microbiota, permutations = 999, by = "terms")
                                           Df
metadata_df$Source                          5
metadata_df$Microbiota                      1
metadata_df$Source:metadata_df$Microbiota   5
Residual                                  134
Total                                     145
                                          SumOfSqs
metadata_df$Source                          14.895
metadata_df$Microbiota                       1.604
metadata_df$Source:metadata_df$Microbiota    3.225
Residual                                    38.396
Total                                       58.120
                                               R2
metadata_df$Source                        0.25628
metadata_df$Microbiota                    0.02760
metadata_df$Source:metadata_df$Microbiota 0.05548
Residual                                  0.66064
Total                                     1.00000
                                                F
metadata_df$Source                        10.3965
metadata_df$Microbiota                     5.5975
metadata_df$Source:metadata_df$Microbiota  2.2507
Residual                                         
Total                                            
                                          Pr(>F)
metadata_df$Source                         0.001
metadata_df$Microbiota                     0.001
metadata_df$Source:metadata_df$Microbiota  0.001
Residual                                        
Total                                           
                                             
metadata_df$Source                        ***
metadata_df$Microbiota                    ***
metadata_df$Source:metadata_df$Microbiota ***
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
PERMANOVA for Source:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$Source, permutations = 999, by = "terms")
                    Df SumOfSqs     R2
metadata_df$Source   5    8.082 0.1653
Residual           140   40.811 0.8347
Total              145   48.894 1.0000
                        F Pr(>F)    
metadata_df$Source 5.5452  0.001 ***
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
metadata_df$Microbiota   1    0.717 0.01466
Residual               144   48.177 0.98534
Total                  145   48.894 1.00000
                            F Pr(>F)   
metadata_df$Microbiota 2.1426  0.002 **
Residual                               
Total                                  
---
Signif. codes:  
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1
  ‘ ’ 1

PERMANOVA for Source * Microbiota interaction:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$Source * metadata_df$Microbiota, permutations = 999, by = "terms")
                                           Df
metadata_df$Source                          5
metadata_df$Microbiota                      1
metadata_df$Source:metadata_df$Microbiota   5
Residual                                  134
Total                                     145
                                          SumOfSqs
metadata_df$Source                           8.082
metadata_df$Microbiota                       0.729
metadata_df$Source:metadata_df$Microbiota    1.860
Residual                                    38.222
Total                                       48.894
                                               R2
metadata_df$Source                        0.16530
metadata_df$Microbiota                    0.01492
metadata_df$Source:metadata_df$Microbiota 0.03804
Residual                                  0.78174
Total                                     1.00000
                                               F
metadata_df$Source                        5.6670
metadata_df$Microbiota                    2.5571
metadata_df$Source:metadata_df$Microbiota 1.3041
Residual                                        
Total                                           
                                          Pr(>F)
metadata_df$Source                         0.001
metadata_df$Microbiota                     0.001
metadata_df$Source:metadata_df$Microbiota  0.017
Residual                                        
Total                                           
                                             
metadata_df$Source                        ***
metadata_df$Microbiota                    ***
metadata_df$Source:metadata_df$Microbiota *  
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
PERMANOVA for Source:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$Source, permutations = 999, by = "terms")
                    Df SumOfSqs     R2
metadata_df$Source   5  0.88034 0.3368
Residual           140  1.73353 0.6632
Total              145  2.61387 1.0000
                        F Pr(>F)    
metadata_df$Source 14.219  0.001 ***
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
metadata_df$Microbiota   1  0.19625 0.07508
Residual               144  2.41761 0.92492
Total                  145  2.61387 1.00000
                            F Pr(>F)   
metadata_df$Microbiota 11.689  0.002 **
Residual                               
Total                                  
---
Signif. codes:  
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1
  ‘ ’ 1

PERMANOVA for Source * Microbiota interaction:
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dist_matrix ~ metadata_df$Source * metadata_df$Microbiota, permutations = 999, by = "terms")
                                           Df
metadata_df$Source                          5
metadata_df$Microbiota                      1
metadata_df$Source:metadata_df$Microbiota   5
Residual                                  134
Total                                     145
                                          SumOfSqs
metadata_df$Source                         0.88034
metadata_df$Microbiota                     0.17979
metadata_df$Source:metadata_df$Microbiota  0.41621
Residual                                   1.13753
Total                                      2.61387
                                               R2
metadata_df$Source                        0.33680
metadata_df$Microbiota                    0.06878
metadata_df$Source:metadata_df$Microbiota 0.15923
Residual                                  0.43519
Total                                     1.00000
                                                F
metadata_df$Source                        20.7406
metadata_df$Microbiota                    21.1787
metadata_df$Source:metadata_df$Microbiota  9.8057
Residual                                         
Total                                            
                                          Pr(>F)
metadata_df$Source                         0.001
metadata_df$Microbiota                     0.001
metadata_df$Source:metadata_df$Microbiota  0.001
Residual                                        
Total                                           
                                             
metadata_df$Source                        ***
metadata_df$Microbiota                    ***
metadata_df$Source:metadata_df$Microbiota ***
Residual                                     
Total                                        
---
Signif. codes:  
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1
  ‘ ’ 1
```

```R
=== PAIRWISE PERMANOVA RESULTS FOR BRAY-CURTIS DISTANCE ===
> pairwise_bray <- perform_pairwise_permanova(ps_insects, "bray")
Sources for pairwise PERMANOVA: Chrysomya rufifacies, Cochliomyia macellaria, Hydrotea aenescens, Chrysomya megacephala, Lucilia coeruleiviridis, Necrodes surinamensis 
Performing pairwise PERMANOVA tests...
> print(pairwise_bray)
                   Source1
1     Chrysomya rufifacies
2     Chrysomya rufifacies
3     Chrysomya rufifacies
4     Chrysomya rufifacies
5     Chrysomya rufifacies
6   Cochliomyia macellaria
7   Cochliomyia macellaria
8   Cochliomyia macellaria
9   Cochliomyia macellaria
10      Hydrotea aenescens
11      Hydrotea aenescens
12      Hydrotea aenescens
13   Chrysomya megacephala
14   Chrysomya megacephala
15 Lucilia coeruleiviridis
                   Source2         R2
1   Cochliomyia macellaria 0.13283980
2       Hydrotea aenescens 0.24212766
3    Chrysomya megacephala 0.06202670
4  Lucilia coeruleiviridis 0.08242578
5    Necrodes surinamensis 0.20781697
6       Hydrotea aenescens 0.13638645
7    Chrysomya megacephala 0.09813407
8  Lucilia coeruleiviridis 0.07177935
9    Necrodes surinamensis 0.19541008
10   Chrysomya megacephala 0.19870689
11 Lucilia coeruleiviridis 0.24632265
12   Necrodes surinamensis 0.09787505
13 Lucilia coeruleiviridis 0.08305595
14   Necrodes surinamensis 0.10649895
15   Necrodes surinamensis 0.09958280
   p_value p_adjusted significance
1    0.001      0.015            *
2    0.001      0.015            *
3    0.064      0.960           ns
4    0.015      0.225           ns
5    0.001      0.015            *
6    0.001      0.015            *
7    0.001      0.015            *
8    0.001      0.015            *
9    0.001      0.015            *
10   0.016      0.240           ns
11   0.002      0.030            *
12   0.001      0.015            *
13   0.084      1.000           ns
14   0.001      0.015            *
15   0.001      0.015            *
```

```R
=== PAIRWISE PERMANOVA RESULTS FOR UNWEIGHTED UNIFRAC DISTANCE ===
> pairwise_unifrac <- perform_pairwise_permanova(ps_insects, "unifrac")
Sources for pairwise PERMANOVA: Chrysomya rufifacies, Cochliomyia macellaria, Hydrotea aenescens, Chrysomya megacephala, Lucilia coeruleiviridis, Necrodes surinamensis 
Performing pairwise PERMANOVA tests...
> print(pairwise_unifrac)
                   Source1
1     Chrysomya rufifacies
2     Chrysomya rufifacies
3     Chrysomya rufifacies
4     Chrysomya rufifacies
5     Chrysomya rufifacies
6   Cochliomyia macellaria
7   Cochliomyia macellaria
8   Cochliomyia macellaria
9   Cochliomyia macellaria
10      Hydrotea aenescens
11      Hydrotea aenescens
12      Hydrotea aenescens
13   Chrysomya megacephala
14   Chrysomya megacephala
15 Lucilia coeruleiviridis
                   Source2         R2
1   Cochliomyia macellaria 0.08896917
2       Hydrotea aenescens 0.19623348
3    Chrysomya megacephala 0.07917769
4  Lucilia coeruleiviridis 0.14903775
5    Necrodes surinamensis 0.13312355
6       Hydrotea aenescens 0.05930902
7    Chrysomya megacephala 0.05523098
8  Lucilia coeruleiviridis 0.07064945
9    Necrodes surinamensis 0.09749980
10   Chrysomya megacephala 0.15318868
11 Lucilia coeruleiviridis 0.13118184
12   Necrodes surinamensis 0.03072606
13 Lucilia coeruleiviridis 0.09839317
14   Necrodes surinamensis 0.05939381
15   Necrodes surinamensis 0.05596556
   p_value p_adjusted significance
1    0.001      0.015            *
2    0.001      0.015            *
3    0.003      0.045            *
4    0.001      0.015            *
5    0.001      0.015            *
6    0.001      0.015            *
7    0.002      0.030            *
8    0.001      0.015            *
9    0.001      0.015            *
10   0.001      0.015            *
11   0.007      0.105           ns
12   0.015      0.225           ns
13   0.014      0.210           ns
14   0.001      0.015            *
15   0.001      0.015            *
```

```R
=== PAIRWISE PERMANOVA RESULTS FOR WEIGHTED UNIFRAC DISTANCE ===
> pairwise_wunifrac <- perform_pairwise_permanova(ps_insects, "wunifrac")
Sources for pairwise PERMANOVA: Chrysomya rufifacies, Cochliomyia macellaria, Hydrotea aenescens, Chrysomya megacephala, Lucilia coeruleiviridis, Necrodes surinamensis 
Performing pairwise PERMANOVA tests...
> print(pairwise_wunifrac)
                   Source1
1     Chrysomya rufifacies
2     Chrysomya rufifacies
3     Chrysomya rufifacies
4     Chrysomya rufifacies
5     Chrysomya rufifacies
6   Cochliomyia macellaria
7   Cochliomyia macellaria
8   Cochliomyia macellaria
9   Cochliomyia macellaria
10      Hydrotea aenescens
11      Hydrotea aenescens
12      Hydrotea aenescens
13   Chrysomya megacephala
14   Chrysomya megacephala
15 Lucilia coeruleiviridis
                   Source2         R2
1   Cochliomyia macellaria 0.26887143
2       Hydrotea aenescens 0.22373107
3    Chrysomya megacephala 0.03375734
4  Lucilia coeruleiviridis 0.12885691
5    Necrodes surinamensis 0.47114944
6       Hydrotea aenescens 0.03371293
7    Chrysomya megacephala 0.11358659
8  Lucilia coeruleiviridis 0.02298911
9    Necrodes surinamensis 0.10371224
10   Chrysomya megacephala 0.20796529
11 Lucilia coeruleiviridis 0.11015605
12   Necrodes surinamensis 0.06575322
13 Lucilia coeruleiviridis 0.05371058
14   Necrodes surinamensis 0.33965686
15   Necrodes surinamensis 0.17574505
   p_value p_adjusted significance
1    0.001      0.015            *
2    0.010      0.150           ns
3    0.313      1.000           ns
4    0.039      0.585           ns
5    0.001      0.015            *
6    0.167      1.000           ns
7    0.010      0.150           ns
8    0.241      1.000           ns
9    0.001      0.015            *
10   0.034      0.510           ns
11   0.155      1.000           ns
12   0.047      0.705           ns
13   0.375      1.000           ns
14   0.001      0.015            *
15   0.001      0.015            *
```

