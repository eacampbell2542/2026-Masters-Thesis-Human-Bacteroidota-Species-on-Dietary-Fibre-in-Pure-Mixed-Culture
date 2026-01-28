# R version 4.4.3 (2025-02-28 ucrt) -- "Trophy Case"
# Date Created: # Tue Oct 14 16:28:17 2025 ------------------------------
#Date last modified:

#Project: Explore the relationship between Bacteroidota
#Author Ellie Campbell
#Email r01ec25@aberdeen.ac.uk

#Packages to load
library(ggplot2)
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)


####################
#Gene copies per mL
####################

Pair_Copies <- read.table(file = "~/Masters/Statistical Analysis/qPCR_Pair/sp1_sp2_fixed.txt", header = TRUE, sep = "\t",
                          stringsAsFactors = TRUE)
Pair_Copies$Substrate <- factor(Pair_Copies$Substrate, levels= c("NOCHO", "Glucose", "Apple_Pectin", 'B_D_glucan', 'Galactomannan', 'Potato_Starch', 'Xylan', 'Xyloglucan'))#Reorder by substrate
str(Pair_Copies$Substrate)
Pair_Copies$Substrate

Pair_Copies$Bacteria_Pair <- factor(Pair_Copies$Bacteria_Pair, levels= c("V975/DSM_19555", "V975/B54832", "V975/DSM_1447", 'DSM_19555/B5482', 'DSM_19555/DSM_1447', 'B5482/DSM_1447')) #reorder alphabetically by substrate
str(Pair_Copies$Bacteria_Pair)

str(Pair_Copies)
# 'data.frame':	144 obs. of  7 variables:
#   $ Bacteria_Pair: Factor w/ 6 levels "V975/DSM_19555",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Substrate    : Factor w/ 8 levels "NOCHO","Glucose",..: 1 1 1 2 2 2 3 3 3 4 ...
# $ S_Name       : Factor w/ 144 levels "B1_V9B5_NOCHO",..: 49 57 59 61 63 65 67 69 71 51 ...
# $ DNA_Number   : int  25 26 27 28 29 30 31 32 33 106 ...
# $ X16S_total   : num  8.50e+09 6.66e+09 4.34e+09 7.57e+09 1.08e+10 ...
# $ Species_1    : num  3.46e+08 2.78e+08 2.22e+08 1.10e+09 1.20e+09 ...
# $ Species_2    : num  1.63e+09 1.17e+09 6.84e+08 6.20e+08 1.04e+09 ...

boxplot(Species_1~Substrate, data=Pair_Copies, xlab="Substrate", ylab="Copier per mL", ylim=c(30000000, 4000000000))
boxplot(Species_2~Substrate, data=Pair_Copies, xlab="Substrate", ylab="Copier per mL", ylim=c(30000000, 4000000000))

#estimated marginal means was run to conduct a pairwise comparison between species pair, substrate, and the individual species in the pair (species 1 & species 2). 


#so each species has respective log of gene copies per ml
Pair_long <- Pair_Copies %>%
  pivot_longer(
    cols = c(Species_1, Species_2),
    names_to = "Species",
    values_to = "Copies"
  )


Pair_long <- Pair_long %>%
  mutate(logCopies = log10(Copies))

str(Pair_long)
# tibble [288 × 8] (S3: tbl_df/tbl/data.frame)
# $ Bacteria_Pair: Factor w/ 6 levels "V975/DSM_19555",..: 2 2 2 2 2 2 2 2 2 2 ...
# $ Substrate    : Factor w/ 8 levels "NOCHO","Glucose",..: 1 1 1 1 1 1 2 2 2 2 ...
# $ S_Name       : Factor w/ 144 levels "B1_V9B5_NOCHO",..: 1 1 9 9 11 11 13 13 15 15 ...
# $ DNA_Number   : int [1:288] 1 1 2 2 3 3 4 4 5 5 ...
# $ X16S_total   : num [1:288] 5.77e+09 5.77e+09 6.10e+09 6.10e+09 4.39e+09 4.39e+09 4.85e+09 4.85e+09 5.03e+09 5.03e+09 ...
# $ Species      : chr [1:288] "Species_1" "Species_2" "Species_1" "Species_2" ...
# $ Copies       : num [1:288] 8.20e+08 4.78e+08 8.78e+08 5.36e+08 6.42e+08 4.41e+08 5.07e+08 5.28e+08 5.04e+08 5.67e+08 ...
# $ logCopies    : num [1:288] 8.91 8.68 8.94 8.73 8.81 ...


#Fit a linear mixed effects model
model <- lmer(
  logCopies ~ Species * Bacteria_Pair * Substrate +
    (1 | S_Name),
  data = Pair_long
) #linear mixed model

anova(model, test = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#                                 Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
# Species                         1.3161 1.31610     1    96 511.3746 < 2.2e-16 ***
# Bacteria_Pair                   0.2759 0.05518     5    96  21.4421 2.387e-14 ***
# Substrate                       1.2354 0.17648     7    96  68.5723 < 2.2e-16 ***
# Species:Bacteria_Pair           2.4101 0.48201     5    96 187.2862 < 2.2e-16 ***
# Species:Substrate               1.7797 0.25425     7    96  98.7877 < 2.2e-16 ***
# Bacteria_Pair:Substrate         0.8357 0.02388    35    96   9.2774 < 2.2e-16 ***
# Species:Bacteria_Pair:Substrate 8.7850 0.25100    35    96  97.5271 < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(model)
#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: logCopies ~ Species * Bacteria_Pair * Substrate + (1 | S_Name)
#    Data: Pair_long
# 
# REML criterion at convergence: -307.8
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -2.88974 -0.30156  0.00531  0.26622  2.02771 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  S_Name   (Intercept) 0.007702 0.08776 
#  Residual             0.002574 0.05073 
# Number of obs: 288, groups:  S_Name, 144
# 
# Fixed effects:
#                                                                           Estimate Std. Error         df t value
# (Intercept)                                                              9.038e+00  5.853e-02  1.229e+02 154.437
# SpeciesSpecies_2                                                        -5.953e-01  4.142e-02  9.600e+01 -14.372
# Bacteria_PairV975/B54832                                                -1.502e-01  8.277e-02  1.229e+02  -1.815
# Bacteria_PairV975/DSM_1447                                              -4.253e-02  8.277e-02  1.229e+02  -0.514
# Bacteria_PairDSM_19555/B5482                                            -7.080e-01  8.277e-02  1.229e+02  -8.554
# Bacteria_PairDSM_19555/DSM_1447                                         -5.681e-01  8.277e-02  1.229e+02  -6.864
# Bacteria_PairB5482/DSM_1447                                             -1.796e-01  8.277e-02  1.229e+02  -2.170
# SubstrateGlucose                                                        -9.632e-02  8.277e-02  1.229e+02  -1.164
# SubstrateApple_Pectin                                                    2.426e-02  8.277e-02  1.229e+02   0.293
# SubstrateB_D_glucan                                                     -1.955e-01  8.277e-02  1.229e+02  -2.362
# SubstrateGalactomannan                                                   5.301e-02  8.277e-02  1.229e+02   0.641
# SubstratePotato_Starch                                                  -6.360e-01  8.277e-02  1.229e+02  -7.684
# SubstrateXylan                                                           2.196e-01  8.277e-02  1.229e+02   2.653
# SubstrateXyloglucan                                                      2.403e-01  8.277e-02  1.229e+02   2.903
# SpeciesSpecies_2:Bacteria_PairV975/B54832                                3.914e-01  5.858e-02  9.600e+01   6.681
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447                              1.112e-01  5.858e-02  9.600e+01   1.898
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482                            1.134e+00  5.858e-02  9.600e+01  19.357
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447                         7.352e-01  5.858e-02  9.600e+01  12.551
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447                             3.105e-01  5.858e-02  9.600e+01   5.300
# SpeciesSpecies_2:SubstrateGlucose                                        6.976e-01  5.858e-02  9.600e+01  11.909
# SpeciesSpecies_2:SubstrateApple_Pectin                                   8.331e-01  5.858e-02  9.600e+01  14.222
# SpeciesSpecies_2:SubstrateB_D_glucan                                     8.304e-02  5.858e-02  9.600e+01   1.418
# SpeciesSpecies_2:SubstrateGalactomannan                                 -1.801e-01  5.858e-02  9.600e+01  -3.074
# SpeciesSpecies_2:SubstratePotato_Starch                                  1.183e+00  5.858e-02  9.600e+01  20.197
# SpeciesSpecies_2:SubstrateXylan                                         -3.046e-02  5.858e-02  9.600e+01  -0.520
# SpeciesSpecies_2:SubstrateXyloglucan                                    -2.434e-01  5.858e-02  9.600e+01  -4.154
# Bacteria_PairV975/B54832:SubstrateGlucose                               -5.861e-03  1.171e-01  1.229e+02  -0.050
# Bacteria_PairV975/DSM_1447:SubstrateGlucose                             -2.380e-01  1.171e-01  1.229e+02  -2.034
# Bacteria_PairDSM_19555/B5482:SubstrateGlucose                            7.714e-01  1.171e-01  1.229e+02   6.591
# Bacteria_PairDSM_19555/DSM_1447:SubstrateGlucose                         7.326e-01  1.171e-01  1.229e+02   6.259
# Bacteria_PairB5482/DSM_1447:SubstrateGlucose                             3.229e-02  1.171e-01  1.229e+02   0.276
# Bacteria_PairV975/B54832:SubstrateApple_Pectin                          -7.334e-04  1.171e-01  1.229e+02  -0.006
# Bacteria_PairV975/DSM_1447:SubstrateApple_Pectin                         2.502e-01  1.171e-01  1.229e+02   2.137
# Bacteria_PairDSM_19555/B5482:SubstrateApple_Pectin                       1.005e+00  1.171e-01  1.229e+02   8.585
# Bacteria_PairDSM_19555/DSM_1447:SubstrateApple_Pectin                    9.716e-01  1.171e-01  1.229e+02   8.300
# Bacteria_PairB5482/DSM_1447:SubstrateApple_Pectin                        3.060e-01  1.171e-01  1.229e+02   2.614
# Bacteria_PairV975/B54832:SubstrateB_D_glucan                             3.542e-01  1.171e-01  1.229e+02   3.026
# Bacteria_PairV975/DSM_1447:SubstrateB_D_glucan                           4.046e-01  1.171e-01  1.229e+02   3.457
# Bacteria_PairDSM_19555/B5482:SubstrateB_D_glucan                         1.022e-01  1.171e-01  1.229e+02   0.874
# Bacteria_PairDSM_19555/DSM_1447:SubstrateB_D_glucan                     -1.248e-01  1.171e-01  1.229e+02  -1.066
# Bacteria_PairB5482/DSM_1447:SubstrateB_D_glucan                          1.053e-01  1.171e-01  1.229e+02   0.900
# Bacteria_PairV975/B54832:SubstrateGalactomannan                          3.816e-01  1.171e-01  1.229e+02   3.260
# Bacteria_PairV975/DSM_1447:SubstrateGalactomannan                        1.623e-01  1.171e-01  1.229e+02   1.386
# Bacteria_PairDSM_19555/B5482:SubstrateGalactomannan                     -1.164e-01  1.171e-01  1.229e+02  -0.994
# Bacteria_PairDSM_19555/DSM_1447:SubstrateGalactomannan                  -8.251e-02  1.171e-01  1.229e+02  -0.705
# Bacteria_PairB5482/DSM_1447:SubstrateGalactomannan                      -5.129e-01  1.171e-01  1.229e+02  -4.382
# Bacteria_PairV975/B54832:SubstratePotato_Starch                          2.936e-01  1.171e-01  1.229e+02   2.508
# Bacteria_PairV975/DSM_1447:SubstratePotato_Starch                       -2.222e-01  1.171e-01  1.229e+02  -1.898
# Bacteria_PairDSM_19555/B5482:SubstratePotato_Starch                      1.369e+00  1.171e-01  1.229e+02  11.692
# Bacteria_PairDSM_19555/DSM_1447:SubstratePotato_Starch                   1.193e+00  1.171e-01  1.229e+02  10.194
# Bacteria_PairB5482/DSM_1447:SubstratePotato_Starch                       2.953e-01  1.171e-01  1.229e+02   2.523
# Bacteria_PairV975/B54832:SubstrateXylan                                  1.762e-01  1.171e-01  1.229e+02   1.505
# Bacteria_PairV975/DSM_1447:SubstrateXylan                                1.050e-01  1.171e-01  1.229e+02   0.897
# Bacteria_PairDSM_19555/B5482:SubstrateXylan                             -7.874e-03  1.171e-01  1.229e+02  -0.067
# Bacteria_PairDSM_19555/DSM_1447:SubstrateXylan                          -6.596e-02  1.171e-01  1.229e+02  -0.564
# Bacteria_PairB5482/DSM_1447:SubstrateXylan                               1.328e-02  1.171e-01  1.229e+02   0.113
# Bacteria_PairV975/B54832:SubstrateXyloglucan                             1.783e-01  1.171e-01  1.229e+02   1.523
# Bacteria_PairV975/DSM_1447:SubstrateXyloglucan                           2.104e-01  1.171e-01  1.229e+02   1.798
# Bacteria_PairDSM_19555/B5482:SubstrateXyloglucan                        -2.165e-01  1.171e-01  1.229e+02  -1.850
# Bacteria_PairDSM_19555/DSM_1447:SubstrateXyloglucan                      8.802e-03  1.171e-01  1.229e+02   0.075
# Bacteria_PairB5482/DSM_1447:SubstrateXyloglucan                         -1.544e-01  1.171e-01  1.229e+02  -1.319
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstrateGlucose              -5.342e-01  8.284e-02  9.600e+01  -6.449
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstrateGlucose            -2.250e-01  8.284e-02  9.600e+01  -2.716
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstrateGlucose          -1.428e+00  8.284e-02  9.600e+01 -17.233
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstrateGlucose       -1.337e+00  8.284e-02  9.600e+01 -16.137
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstrateGlucose           -6.028e-01  8.284e-02  9.600e+01  -7.276
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstrateApple_Pectin         -3.741e-01  8.284e-02  9.600e+01  -4.516
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstrateApple_Pectin       -6.948e-01  8.284e-02  9.600e+01  -8.387
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstrateApple_Pectin     -1.666e+00  8.284e-02  9.600e+01 -20.116
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstrateApple_Pectin  -1.776e+00  8.284e-02  9.600e+01 -21.439
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstrateApple_Pectin      -6.174e-01  8.284e-02  9.600e+01  -7.452
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstrateB_D_glucan            5.699e-02  8.284e-02  9.600e+01   0.688
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstrateB_D_glucan         -2.366e-01  8.284e-02  9.600e+01  -2.856
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstrateB_D_glucan        1.970e-02  8.284e-02  9.600e+01   0.238
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstrateB_D_glucan     4.837e-02  8.284e-02  9.600e+01   0.584
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstrateB_D_glucan         1.452e-02  8.284e-02  9.600e+01   0.175
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstrateGalactomannan        -4.235e-01  8.284e-02  9.600e+01  -5.112
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstrateGalactomannan      -2.547e-01  8.284e-02  9.600e+01  -3.075
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstrateGalactomannan    -2.299e-01  8.284e-02  9.600e+01  -2.775
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstrateGalactomannan -8.163e-02  8.284e-02  9.600e+01  -0.985
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstrateGalactomannan      4.343e-01  8.284e-02  9.600e+01   5.243
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstratePotato_Starch        -6.667e-01  8.284e-02  9.600e+01  -8.048
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstratePotato_Starch      -1.308e-01  8.284e-02  9.600e+01  -1.579
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstratePotato_Starch    -2.042e+00  8.284e-02  9.600e+01 -24.649
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstratePotato_Starch -1.578e+00  8.284e-02  9.600e+01 -19.051
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstratePotato_Starch     -7.712e-01  8.284e-02  9.600e+01  -9.310
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstrateXylan                -2.977e-02  8.284e-02  9.600e+01  -0.359
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstrateXylan               2.485e-01  8.284e-02  9.600e+01   2.999
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstrateXylan             8.603e-03  8.284e-02  9.600e+01   0.104
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstrateXylan          2.998e-01  8.284e-02  9.600e+01   3.619
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstrateXylan              2.439e-01  8.284e-02  9.600e+01   2.944
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstrateXyloglucan            2.583e-01  8.284e-02  9.600e+01   3.118
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstrateXyloglucan         -9.023e-02  8.284e-02  9.600e+01  -1.089
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstrateXyloglucan        3.653e-01  8.284e-02  9.600e+01   4.410
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstrateXyloglucan     1.546e-01  8.284e-02  9.600e+01   1.866
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstrateXyloglucan         2.084e-01  8.284e-02  9.600e+01   2.515
#                                                                         Pr(>|t|)    
# (Intercept)                                                              < 2e-16 ***
# SpeciesSpecies_2                                                         < 2e-16 ***
# Bacteria_PairV975/B54832                                                0.072013 .  
# Bacteria_PairV975/DSM_1447                                              0.608315    
# Bacteria_PairDSM_19555/B5482                                            3.96e-14 ***
# Bacteria_PairDSM_19555/DSM_1447                                         2.92e-10 ***
# Bacteria_PairB5482/DSM_1447                                             0.031908 *  
# SubstrateGlucose                                                        0.246764    
# SubstrateApple_Pectin                                                   0.769953    
# SubstrateB_D_glucan                                                     0.019730 *  
# SubstrateGalactomannan                                                  0.523027    
# SubstratePotato_Starch                                                  4.19e-12 ***
# SubstrateXylan                                                          0.009032 ** 
# SubstrateXyloglucan                                                     0.004379 ** 
# SpeciesSpecies_2:Bacteria_PairV975/B54832                               1.55e-09 ***
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447                             0.060712 .  
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482                            < 2e-16 ***
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447                         < 2e-16 ***
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447                            7.37e-07 ***
# SpeciesSpecies_2:SubstrateGlucose                                        < 2e-16 ***
# SpeciesSpecies_2:SubstrateApple_Pectin                                   < 2e-16 ***
# SpeciesSpecies_2:SubstrateB_D_glucan                                    0.159550    
# SpeciesSpecies_2:SubstrateGalactomannan                                 0.002746 ** 
# SpeciesSpecies_2:SubstratePotato_Starch                                  < 2e-16 ***
# SpeciesSpecies_2:SubstrateXylan                                         0.604276    
# SpeciesSpecies_2:SubstrateXyloglucan                                    7.08e-05 ***
# Bacteria_PairV975/B54832:SubstrateGlucose                               0.960148    
# Bacteria_PairV975/DSM_1447:SubstrateGlucose                             0.044146 *  
# Bacteria_PairDSM_19555/B5482:SubstrateGlucose                           1.15e-09 ***
# Bacteria_PairDSM_19555/DSM_1447:SubstrateGlucose                        5.90e-09 ***
# Bacteria_PairB5482/DSM_1447:SubstrateGlucose                            0.783094    
# Bacteria_PairV975/B54832:SubstrateApple_Pectin                          0.995011    
# Bacteria_PairV975/DSM_1447:SubstrateApple_Pectin                        0.034543 *  
# Bacteria_PairDSM_19555/B5482:SubstrateApple_Pectin                      3.35e-14 ***
# Bacteria_PairDSM_19555/DSM_1447:SubstrateApple_Pectin                   1.57e-13 ***
# Bacteria_PairB5482/DSM_1447:SubstrateApple_Pectin                       0.010054 *  
# Bacteria_PairV975/B54832:SubstrateB_D_glucan                            0.003020 ** 
# Bacteria_PairV975/DSM_1447:SubstrateB_D_glucan                          0.000751 ***
# Bacteria_PairDSM_19555/B5482:SubstrateB_D_glucan                        0.384070    
# Bacteria_PairDSM_19555/DSM_1447:SubstrateB_D_glucan                     0.288503    
# Bacteria_PairB5482/DSM_1447:SubstrateB_D_glucan                         0.370038    
# Bacteria_PairV975/B54832:SubstrateGalactomannan                         0.001440 ** 
# Bacteria_PairV975/DSM_1447:SubstrateGalactomannan                       0.168196    
# Bacteria_PairDSM_19555/B5482:SubstrateGalactomannan                     0.321968    
# Bacteria_PairDSM_19555/DSM_1447:SubstrateGalactomannan                  0.482200    
# Bacteria_PairB5482/DSM_1447:SubstrateGalactomannan                      2.50e-05 ***
# Bacteria_PairV975/B54832:SubstratePotato_Starch                         0.013436 *  
# Bacteria_PairV975/DSM_1447:SubstratePotato_Starch                       0.060054 .  
# Bacteria_PairDSM_19555/B5482:SubstratePotato_Starch                      < 2e-16 ***
# Bacteria_PairDSM_19555/DSM_1447:SubstratePotato_Starch                   < 2e-16 ***
# Bacteria_PairB5482/DSM_1447:SubstratePotato_Starch                      0.012913 *  
# Bacteria_PairV975/B54832:SubstrateXylan                                 0.134782    
# Bacteria_PairV975/DSM_1447:SubstrateXylan                               0.371551    
# Bacteria_PairDSM_19555/B5482:SubstrateXylan                             0.946475    
# Bacteria_PairDSM_19555/DSM_1447:SubstrateXylan                          0.574101    
# Bacteria_PairB5482/DSM_1447:SubstrateXylan                              0.909883    
# Bacteria_PairV975/B54832:SubstrateXyloglucan                            0.130301    
# Bacteria_PairV975/DSM_1447:SubstrateXyloglucan                          0.074707 .  
# Bacteria_PairDSM_19555/B5482:SubstrateXyloglucan                        0.066756 .  
# Bacteria_PairDSM_19555/DSM_1447:SubstrateXyloglucan                     0.940183    
# Bacteria_PairB5482/DSM_1447:SubstrateXyloglucan                         0.189682    
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstrateGlucose              4.54e-09 ***
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstrateGlucose            0.007840 ** 
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstrateGlucose           < 2e-16 ***
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstrateGlucose        < 2e-16 ***
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstrateGlucose           9.35e-11 ***
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstrateApple_Pectin         1.79e-05 ***
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstrateApple_Pectin       4.28e-13 ***
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstrateApple_Pectin      < 2e-16 ***
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstrateApple_Pectin   < 2e-16 ***
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstrateApple_Pectin      4.02e-11 ***
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstrateB_D_glucan           0.493161    
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstrateB_D_glucan         0.005256 ** 
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstrateB_D_glucan       0.812533    
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstrateB_D_glucan    0.560643    
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstrateB_D_glucan        0.861242    
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstrateGalactomannan        1.62e-06 ***
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstrateGalactomannan      0.002743 ** 
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstrateGalactomannan    0.006632 ** 
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstrateGalactomannan 0.326956    
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstrateGalactomannan     9.40e-07 ***
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstratePotato_Starch        2.25e-12 ***
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstratePotato_Starch      0.117702    
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstratePotato_Starch     < 2e-16 ***
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstratePotato_Starch  < 2e-16 ***
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstratePotato_Starch     4.53e-15 ***
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstrateXylan                0.720082    
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstrateXylan              0.003445 ** 
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstrateXylan            0.917505    
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstrateXylan         0.000474 ***
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstrateXylan             0.004065 ** 
# SpeciesSpecies_2:Bacteria_PairV975/B54832:SubstrateXyloglucan           0.002400 ** 
# SpeciesSpecies_2:Bacteria_PairV975/DSM_1447:SubstrateXyloglucan         0.278823    
# SpeciesSpecies_2:Bacteria_PairDSM_19555/B5482:SubstrateXyloglucan       2.70e-05 ***
# SpeciesSpecies_2:Bacteria_PairDSM_19555/DSM_1447:SubstrateXyloglucan    0.065145 .  
# SpeciesSpecies_2:Bacteria_PairB5482/DSM_1447:SubstrateXyloglucan        0.013565 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Correlation matrix not shown by default, as p = 96 > 12.
# Use print(x, correlation=TRUE)  or
#     vcov(x)        if you need it


#Estimated means of the model
emm <- emmeans(model, ~ Species | Bacteria_Pair * Substrate) #estimated means
pairs(emm) #pairwise comparison from the marginal estimated means

# Bacteria_Pair = V975/DSM_19555, Substrate = NOCHO:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.5953 0.0414 96  14.372  <.0001
# 
# Bacteria_Pair = V975/B54832, Substrate = NOCHO:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.2039 0.0414 96   4.923  <.0001
# 
# Bacteria_Pair = V975/DSM_1447, Substrate = NOCHO:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.4841 0.0414 96  11.688  <.0001
# 
# Bacteria_Pair = DSM_19555/B5482, Substrate = NOCHO:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.5386 0.0414 96 -13.002  <.0001
# 
# Bacteria_Pair = DSM_19555/DSM_1447, Substrate = NOCHO:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.1399 0.0414 96  -3.378  0.0011
# 
# Bacteria_Pair = B5482/DSM_1447, Substrate = NOCHO:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.2848 0.0414 96   6.876  <.0001
# 
# Bacteria_Pair = V975/DSM_19555, Substrate = Glucose:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.1023 0.0414 96  -2.470  0.0153
# 
# Bacteria_Pair = V975/B54832, Substrate = Glucose:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.0406 0.0414 96   0.979  0.3301
# 
# Bacteria_Pair = V975/DSM_1447, Substrate = Glucose:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.0115 0.0414 96   0.278  0.7820
# 
# Bacteria_Pair = DSM_19555/B5482, Substrate = Glucose:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.1914 0.0414 96   4.622  <.0001
# 
# Bacteria_Pair = DSM_19555/DSM_1447, Substrate = Glucose:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.4993 0.0414 96  12.055  <.0001
# 
# Bacteria_Pair = B5482/DSM_1447, Substrate = Glucose:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.1900 0.0414 96   4.586  <.0001
# 
# Bacteria_Pair = V975/DSM_19555, Substrate = Apple_Pectin:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.2378 0.0414 96  -5.741  <.0001
# 
# Bacteria_Pair = V975/B54832, Substrate = Apple_Pectin:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.2551 0.0414 96  -6.158  <.0001
# 
# Bacteria_Pair = V975/DSM_1447, Substrate = Apple_Pectin:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.3458 0.0414 96   8.349  <.0001
# 
# Bacteria_Pair = DSM_19555/B5482, Substrate = Apple_Pectin:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.2948 0.0414 96   7.116  <.0001
# 
# Bacteria_Pair = DSM_19555/DSM_1447, Substrate = Apple_Pectin:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.8031 0.0414 96  19.388  <.0001
# 
# Bacteria_Pair = B5482/DSM_1447, Substrate = Apple_Pectin:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.0691 0.0414 96   1.668  0.0987
# 
# Bacteria_Pair = V975/DSM_19555, Substrate = B_D_glucan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.5123 0.0414 96  12.367  <.0001
# 
# Bacteria_Pair = V975/B54832, Substrate = B_D_glucan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.0639 0.0414 96   1.543  0.1262
# 
# Bacteria_Pair = V975/DSM_1447, Substrate = B_D_glucan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.6377 0.0414 96  15.396  <.0001
# 
# Bacteria_Pair = DSM_19555/B5482, Substrate = B_D_glucan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.6413 0.0414 96 -15.483  <.0001
# 
# Bacteria_Pair = DSM_19555/DSM_1447, Substrate = B_D_glucan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.2713 0.0414 96  -6.550  <.0001
# 
# Bacteria_Pair = B5482/DSM_1447, Substrate = B_D_glucan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.1873 0.0414 96   4.521  <.0001
# 
# Bacteria_Pair = V975/DSM_19555, Substrate = Galactomannan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.7754 0.0414 96  18.720  <.0001
# 
# Bacteria_Pair = V975/B54832, Substrate = Galactomannan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.8075 0.0414 96  19.496  <.0001
# 
# Bacteria_Pair = V975/DSM_1447, Substrate = Galactomannan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.9190 0.0414 96  22.186  <.0001
# 
# Bacteria_Pair = DSM_19555/B5482, Substrate = Galactomannan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.1286 0.0414 96  -3.104  0.0025
# 
# Bacteria_Pair = DSM_19555/DSM_1447, Substrate = Galactomannan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.1218 0.0414 96   2.941  0.0041
# 
# Bacteria_Pair = B5482/DSM_1447, Substrate = Galactomannan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.0306 0.0414 96   0.739  0.4616
# 
# Bacteria_Pair = V975/DSM_19555, Substrate = Potato_Starch:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.5878 0.0414 96 -14.190  <.0001
# 
# Bacteria_Pair = V975/B54832, Substrate = Potato_Starch:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.3124 0.0414 96  -7.543  <.0001
# 
# Bacteria_Pair = V975/DSM_1447, Substrate = Potato_Starch:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.5682 0.0414 96 -13.717  <.0001
# 
# Bacteria_Pair = DSM_19555/B5482, Substrate = Potato_Starch:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.3204 0.0414 96   7.734  <.0001
# 
# Bacteria_Pair = DSM_19555/DSM_1447, Substrate = Potato_Starch:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.2553 0.0414 96   6.163  <.0001
# 
# Bacteria_Pair = B5482/DSM_1447, Substrate = Potato_Starch:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.1270 0.0414 96  -3.067  0.0028
# 
# Bacteria_Pair = V975/DSM_19555, Substrate = Xylan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.6258 0.0414 96  15.107  <.0001
# 
# Bacteria_Pair = V975/B54832, Substrate = Xylan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.2642 0.0414 96   6.378  <.0001
# 
# Bacteria_Pair = V975/DSM_1447, Substrate = Xylan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.2661 0.0414 96   6.425  <.0001
# 
# Bacteria_Pair = DSM_19555/B5482, Substrate = Xylan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.5167 0.0414 96 -12.475  <.0001
# 
# Bacteria_Pair = DSM_19555/DSM_1447, Substrate = Xylan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.4093 0.0414 96  -9.880  <.0001
# 
# Bacteria_Pair = B5482/DSM_1447, Substrate = Xylan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.0714 0.0414 96   1.724  0.0879
# 
# Bacteria_Pair = V975/DSM_19555, Substrate = Xyloglucan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.8387 0.0414 96  20.247  <.0001
# 
# Bacteria_Pair = V975/B54832, Substrate = Xyloglucan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.1889 0.0414 96   4.562  <.0001
# 
# Bacteria_Pair = V975/DSM_1447, Substrate = Xyloglucan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.8177 0.0414 96  19.741  <.0001
# 
# Bacteria_Pair = DSM_19555/B5482, Substrate = Xyloglucan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.6605 0.0414 96 -15.946  <.0001
# 
# Bacteria_Pair = DSM_19555/DSM_1447, Substrate = Xyloglucan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2  -0.0511 0.0414 96  -1.234  0.2203
# 
# Bacteria_Pair = B5482/DSM_1447, Substrate = Xyloglucan:
#   contrast              estimate     SE df t.ratio p.value
# Species_1 - Species_2   0.3198 0.0414 96   7.722  <.0001
# 
# Degrees-of-freedom method: kenward-roger 

posthoc <- pairs(emm, adjust = "fdr") |>  # adjust optional
  as.data.frame()


#Labels for significance
posthoc$stars <- cut(
  posthoc$p.value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns")
)


#Group for use in graph
y_pos <- Pair_long |>
  group_by(Bacteria_Pair, Substrate) |>
  summarise(y = max(logCopies, na.rm = TRUE) + 0.25)

plot_df <- left_join(
  posthoc,
  y_pos,
  by = c("Bacteria_Pair", "Substrate")
)



#clean so only those with significance have bars
plot_df_clean <- plot_df |>
  dplyr::filter(!is.na(stars), stars != "ns")


#visualize the graph 
ggplot(Pair_long,
       aes(x = Species, y = logCopies, fill = Species)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
#appropriate labeling
  facet_grid(
    Bacteria_Pair ~ Substrate,
    labeller = labeller(
      Bacteria_Pair = c(
        "V975/DSM_19555" = "B. ova & B. ster",
        "V975/B54832" = "B. ova & B. theta",
        "V975/DSM_1447" = "B. ova & P. vul",
        "DSM_19555/B5482" = "B. ster & B. theta",
        "DSM_19555/DSM_1447" = "B. ster & P. vul",
        "B5482/DSM_1447" = "B. theta & P. vul"
      ),
      Substrate = c(
        "NOCHO" = "NO CHO",
        "Glucose" = "Glucose",
        "Apple_Pectin" = "Apple Pectin",
        "B_D_glucan" = "β-D-glucan",
        "Galactomannan" = "Galactomannan",
        "Potato_Starch" = "Potato Starch",
        "Xylan" = "Xylan",
        "Xyloglucan" = "Xyloglucan"
      )
    )
  ) +
  
  ## significance bars
  geom_segment(
    data = plot_df_clean,
    aes(x = 1, xend = 2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  
  ## significance stars
  geom_text(
    data = plot_df_clean,
    aes(x = 1.5, y = y + 0.05, label = stars),
    inherit.aes = FALSE,
    size = 4
  ) +
  
  scale_x_discrete(
    labels = c(
      "Species_1" = "Species 1",
      "Species_2" = "Species 2"
    )
  ) +
  scale_fill_discrete(
    labels = c(
      "Species_1" = "Species 1",
      "Species_2" = "Species 2"
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15)) #move bars up
  ) +
  
  labs(y = "log10(gene copies per ml)", fill = "Species") +
  
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    legend.position = "right"
  ) #position of the legend





