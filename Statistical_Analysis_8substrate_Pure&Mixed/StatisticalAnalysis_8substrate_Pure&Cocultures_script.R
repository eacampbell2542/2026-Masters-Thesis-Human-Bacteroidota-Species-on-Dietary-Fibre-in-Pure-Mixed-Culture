# R version 4.4.3 (2025-02-28 ucrt) -- "Trophy Case"
# Date Created: # Tue Oct 14 16:28:17 2025 ------------------------------


#Project: Explore the relationship between Bacteroidota
#Author Ellie Campbell
#Email r01ec25@aberdeen.ac.uk

#Packages to load
library(ggplot2)
library(emmeans)
library(patchwork)
library(psych)

#####################################
#Six species grown on 8 substrates#
####################################

#Maximum OD for pure cultures

plate_max_od <- read.table(file = "~/Masters/Statistical Analysis/qPCR_8substrate/Plate_ind_max.txt", header = TRUE, sep = "\t",
                                stringsAsFactors = TRUE)
plate_max_od$Substrate <- factor(plate_max_od$Substrate, levels= c("NO_CHO", "Glucose", "Apple_Pecitn", 'B_D_glucan', 'Galactomannan', 'Potato_Starch', 'Xylan', 'Xyloglucan'))
str(plate_max_od)
# 'data.frame':	144 obs. of  4 variables:
#   $ Well     : Factor w/ 72 levels "A1","A10","A11",..: 61 65 66 67 68 69 70 71 72 62 ...
# $ Bacteria : Factor w/ 6 levels "B5482","DSM_1447",..: 3 3 3 3 3 3 3 3 3 3 ...
# $ Substrate: Factor w/ 8 levels "NO_CHO","Glucose",..: 1 1 1 2 2 2 3 3 3 6 ...
# $ MaxOD    : num  0.005 0.004 0.059 0.231 0.228 0.239 0.051 0.052 0.055 0.235 ...
table(plate_max_od$Substrate)

boxplot(MaxOD~Bacteria, data=plate_max_od, xlab="Substrate", ylab="Maximum OD")
boxplot(MaxOD~Substrate, data=plate_max_od, xlab="Substrate", ylab="Maximum OD")

######
#Comparison of Maximum OD between each species. Analysis was okay to be done seperately are each well is an individual.
####

#Initally, done based on substrate for figure S2
#B. cellulosilyticus

max_od_14838 <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_dsm14838.txt", header = TRUE, sep = "\t",
                                          stringsAsFactors = TRUE)
str(max_od_14838)
max_od_14838$Substrate <- factor(max_od_14838$Substrate, levels= c("NO_CHO", "Glucose", "Apple_Pecitn", 'B_D_glucan', 'Galactomannan', 'Potato_Starch', 'Xylan', 'Xyloglucan'))
str(max_od_14838$Substrate)
table(max_od_14838$Substrate)

boxplot(MaxOD~Substrate, data=max_od_14838, xlab="Species", ylab="Substrate")

max_od_14838_lm<- lm(MaxOD~Substrate, data=max_od_14838)
drop1(max_od_14838_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Substrate
# Df Sum of Sq      RSS     AIC F value    Pr(>F)    
# <none>                 0.019056 -155.32                      
# Substrate  7   0.22316 0.242218 -108.30  26.768 1.087e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(max_od_14838_lm)
# Call:
#   lm(formula = MaxOD ~ Substrate, data = max_od_14838)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.100667 -0.001667 -0.000167  0.003250  0.089333 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            0.005667   0.019925   0.284  0.77975    
# SubstrateGlucose       0.227000   0.028178   8.056 5.07e-07 ***
#   SubstrateApple_Pecitn  0.047000   0.028178   1.668  0.11477    
# SubstrateB_D_glucan    0.012667   0.028178   0.450  0.65908    
# SubstrateGalactomannan 0.228000   0.028178   8.091 4.78e-07 ***
#   SubstratePotato_Starch 0.232333   0.028178   8.245 3.74e-07 ***
#   SubstrateXylan         0.090000   0.028178   3.194  0.00565 ** 
#   SubstrateXyloglucan    0.033333   0.028178   1.183  0.25411    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03451 on 16 degrees of freedom
# Multiple R-squared:  0.9213,	Adjusted R-squared:  0.8869 
# F-statistic: 26.77 on 7 and 16 DF,  p-value: 1.087e-07

anova(max_od_14838_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df   Sum Sq  Mean Sq F value    Pr(>F)    
# Substrate  7 0.223162 0.031880  26.768 1.087e-07 ***
#   Residuals 16 0.019056 0.001191                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#B. ovatus

max_od_v975 <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_v975.txt", header = TRUE, sep = "\t",
                           stringsAsFactors = TRUE)
str(max_od_v975)
max_od_v975$Substrate <- factor(max_od_v975$Substrate, levels= c("NO_CHO", "Glucose", "Apple_Pecitn", 'B_D_glucan', 'Galactomannan', 'Potato_Starch', 'Xylan', 'Xyloglucan'))
str(max_od_v975$Substrate)
table(max_od_v975$Substrate)
table(max_od_v975$MaxOD)


boxplot(MaxOD~Substrate, data=max_od_v975, xlab="Species", ylab="Substrate")

max_od_v975_lm<- lm(MaxOD~Substrate, data=max_od_v975)
drop1(max_od_v975_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Substrate
# Df Sum of Sq      RSS     AIC F value    Pr(>F)    
# <none>                 0.004899 -187.92                      
# Substrate  7   0.13739 0.142287 -121.07  64.096 1.624e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(max_od_v975_lm)
# Call:
#   lm(formula = MaxOD ~ Substrate, data = max_od_v975)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.020667 -0.009833 -0.002000  0.004917  0.037333 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.13200    0.01010   13.07 5.94e-10 ***
#   SubstrateGlucose        0.24133    0.01429   16.89 1.27e-11 ***
#   SubstrateApple_Pecitn   0.26000    0.01429   18.20 4.08e-12 ***
#   SubstrateB_D_glucan     0.17233    0.01429   12.06 1.91e-09 ***
#   SubstrateGalactomannan  0.16467    0.01429   11.53 3.69e-09 ***
#   SubstratePotato_Starch  0.21367    0.01429   14.96 8.00e-11 ***
#   SubstrateXylan          0.21800    0.01429   15.26 5.92e-11 ***
#   SubstrateXyloglucan     0.16067    0.01429   11.24 5.24e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0175 on 16 degrees of freedom
# Multiple R-squared:  0.9656,	Adjusted R-squared:  0.9505 
# F-statistic:  64.1 on 7 and 16 DF,  p-value: 1.624e-10

anova(max_od_v975_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Substrate  7 0.137388 0.0196269  64.096 1.624e-10 ***
#   Residuals 16 0.004899 0.0003062                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#B. stercoris

max_od_dsm19555 <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_dsm19555.txt", header = TRUE, sep = "\t",
                           stringsAsFactors = TRUE)
str(max_od_dsm19555)
max_od_dsm19555$Substrate <- factor(max_od_dsm19555$Substrate, levels= c("NO_CHO", "Glucose", "Apple_Pecitn", 'B_D_glucan', 'Galactomannan', 'Potato_Starch', 'Xylan', 'Xyloglucan'))
str(max_od_dsm19555$Substrate)
table(max_od_dsm19555$Substrate)

boxplot(MaxOD~Substrate, data=max_od_dsm19555, xlab="Species", ylab="Substrate")

max_od_dsm19555_lm<- lm(MaxOD~Substrate, data=max_od_dsm19555)
drop1(max_od_dsm19555_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Substrate
# Df Sum of Sq     RSS     AIC F value    Pr(>F)    
# <none>                 0.00225 -206.56                      
# Substrate  7   0.75142 0.75368  -81.06  762.22 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(max_od_dsm19555_lm)
# 
# Call:
#   lm(formula = MaxOD ~ Substrate, data = max_od_dsm19555)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0290000 -0.0041667  0.0003333  0.0020833  0.0240000 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             1.087e-01  6.852e-03  15.860 3.30e-11 ***
#   SubstrateGlucose        3.283e-01  9.690e-03  33.885 2.51e-16 ***
#   SubstrateApple_Pecitn   3.633e-01  9.690e-03  37.497  < 2e-16 ***
#   SubstrateB_D_glucan    -3.138e-16  9.690e-03   0.000  1.00000    
# SubstrateGalactomannan -3.633e-02  9.690e-03  -3.750  0.00175 ** 
#   SubstratePotato_Starch  3.873e-01  9.690e-03  39.974  < 2e-16 ***
#   SubstrateXylan          1.300e-02  9.690e-03   1.342  0.19845    
# SubstrateXyloglucan     6.000e-03  9.690e-03   0.619  0.54449    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01187 on 16 degrees of freedom
# Multiple R-squared:  0.997,	Adjusted R-squared:  0.9957 
# F-statistic: 762.2 on 7 and 16 DF,  p-value: < 2.2e-16

anova(max_od_dsm19555_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df  Sum Sq  Mean Sq F value    Pr(>F)    
# Substrate  7 0.75142 0.107346  762.22 < 2.2e-16 ***
#   Residuals 16 0.00225 0.000141                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#B. thetaiotaomicron

max_od_b5482 <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_b5482.txt", header = TRUE, sep = "\t",
                           stringsAsFactors = TRUE)
str(max_od_b5482)
max_od_b5482$Substrate <- factor(max_od_b5482$Substrate, levels= c("NO_CHO", "Glucose", "Apple_Pecitn", 'B_D_glucan', 'Galactomannan', 'Potato_Starch', 'Xylan', 'Xyloglucan'))
str(max_od_b5482$Substrate)
table(max_od_b5482$Substrate)

boxplot(MaxOD~Substrate, data=max_od_b5482, xlab="Species", ylab="Substrate")

max_od_b5482_lm<- lm(MaxOD~Substrate, data=max_od_b5482)
drop1(max_od_b5482_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Substrate
# Df Sum of Sq     RSS      AIC F value    Pr(>F)    
# <none>                 0.01234 -165.749                      
# Substrate  7   0.53139 0.54373  -88.897  98.417 5.945e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(max_od_b5482_lm)
# Call:
#   lm(formula = MaxOD ~ Substrate, data = max_od_b5482)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.072667 -0.005667 -0.001667  0.003250  0.078333 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.13900    0.01604   8.669 1.93e-07 ***
#   SubstrateGlucose        0.34867    0.02268  15.376 5.27e-11 ***
#   SubstrateApple_Pecitn   0.13133    0.02268   5.792 2.75e-05 ***
#   SubstrateB_D_glucan    -0.00700    0.02268  -0.309  0.76154    
# SubstrateGalactomannan -0.07667    0.02268  -3.381  0.00381 ** 
#   SubstratePotato_Starch  0.31367    0.02268  13.832 2.56e-10 ***
#   SubstrateXylan          0.00500    0.02268   0.220  0.82828    
# SubstrateXyloglucan     0.01267    0.02268   0.559  0.58418    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02777 on 16 degrees of freedom
# Multiple R-squared:  0.9773,	Adjusted R-squared:  0.9674 
# F-statistic: 98.42 on 7 and 16 DF,  p-value: 5.945e-12


anova(max_od_b5482_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df  Sum Sq  Mean Sq F value    Pr(>F)    
# Substrate  7 0.53139 0.075912  98.417 5.945e-12 ***
#   Residuals 16 0.01234 0.000771                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#B. xylanisolvens

max_od_18836 <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_dsm18836.txt", header = TRUE, sep = "\t",
                           stringsAsFactors = TRUE)
str(max_od_18836)
max_od_18836$Substrate <- factor(max_od_18836$Substrate, levels= c("NO_CHO", "Glucose", "Apple_Pecitn", 'B_D_glucan', 'Galactomannan', 'Potato_Starch', 'Xylan', 'Xyloglucan'))
str(max_od_18836$Substrate)
table(max_od_18836$Substrate)

boxplot(MaxOD~Substrate, data=max_od_18836, xlab="Species", ylab="Substrate")

max_od_18836_lm<- lm(MaxOD~Substrate, data=max_od_18836)
drop1(max_od_18836_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Substrate
# Df Sum of Sq      RSS     AIC F value    Pr(>F)    
# <none>                 0.002835 -201.05                      
# Substrate  7   0.28703 0.289868 -103.99  231.45 7.267e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(max_od_18836_lm)

# Call:
#   lm(formula = MaxOD ~ Substrate, data = max_od_18836)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.025000 -0.004750 -0.000333  0.001750  0.033000 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.160667   0.007685  20.907 4.82e-13 ***
#   SubstrateGlucose        0.228667   0.010868  21.041 4.37e-13 ***
#   SubstrateApple_Pecitn   0.099333   0.010868   9.140 9.45e-08 ***
#   SubstrateB_D_glucan     0.221333   0.010868  20.366 7.23e-13 ***
#   SubstrateGalactomannan -0.078000   0.010868  -7.177 2.20e-06 ***
#   SubstratePotato_Starch  0.191667   0.010868  17.636 6.59e-12 ***
#   SubstrateXylan          0.194333   0.010868  17.881 5.34e-12 ***
#   SubstrateXyloglucan     0.021333   0.010868   1.963   0.0673 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01331 on 16 degrees of freedom
# Multiple R-squared:  0.9902,	Adjusted R-squared:  0.9859 
# F-statistic: 231.4 on 7 and 16 DF,  p-value: 7.267e-15

anova(max_od_18836_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df   Sum Sq  Mean Sq F value    Pr(>F)    
# Substrate  7 0.287033 0.041005  231.45 7.267e-15 ***
#   Residuals 16 0.002835 0.000177                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#P. vulgatus

max_od_1447 <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_dsm1447.txt", header = TRUE, sep = "\t",
                           stringsAsFactors = TRUE)
str(max_od_1447)
max_od_1447$Substrate <- factor(max_od_1447$Substrate, levels= c("NO_CHO", "Glucose", "Apple_Pecitn", 'B_D_glucan', 'Galactomannan', 'Potato_Starch', 'Xylan', 'Xyloglucan'))
str(max_od_1447$Substrate)
table(max_od_1447$Substrate)

boxplot(MaxOD~Substrate, data=max_od_1447, xlab="Species", ylab="Substrate")

max_od_1447_lm<- lm(MaxOD~Substrate, data=max_od_1447)
drop1(max_od_1447_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Substrate
# Df Sum of Sq      RSS     AIC F value    Pr(>F)    
# <none>                 0.002248 -206.62                      
# Substrate  7   0.20408 0.206330 -112.15  207.51 1.721e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(max_od_1447_lm)

# Call:
#   lm(formula = MaxOD ~ Substrate, data = max_od_1447)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.016000 -0.005500 -0.001333  0.003083  0.024333 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.142000   0.006843  20.750 5.42e-13 ***
#   SubstrateGlucose        0.199667   0.009678  20.631 5.93e-13 ***
#   SubstrateApple_Pecitn   0.110667   0.009678  11.435 4.13e-09 ***
#   SubstrateB_D_glucan    -0.050000   0.009678  -5.166 9.37e-05 ***
#   SubstrateGalactomannan -0.032000   0.009678  -3.306  0.00446 ** 
#   SubstratePotato_Starch  0.154000   0.009678  15.912 3.14e-11 ***
#   SubstrateXylan          0.081000   0.009678   8.369 3.07e-07 ***
#   SubstrateXyloglucan    -0.054667   0.009678  -5.648 3.63e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01185 on 16 degrees of freedom
# Multiple R-squared:  0.9891,	Adjusted R-squared:  0.9843 
# F-statistic: 207.5 on 7 and 16 DF,  p-value: 1.721e-14

anova(max_od_1447_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Substrate  7 0.204082 0.0291545  207.51 1.721e-14 ***
#   Residuals 16 0.002248 0.0001405                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#####individual based on substrate

#NO CHO
max_od_nc <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_nc.txt", header = TRUE, sep = "\t",
                          stringsAsFactors = TRUE)
str(max_od_nc)
max_od_nc$Bacteria <- factor(max_od_nc$Bacteria, levels= c('V975', 'DSM_14838', 'DSM_19555', 'B5482', 'DSM_18836', 'DSM_1447'))
str(max_od_nc$Bacteria)
table(max_od_nc$Bacteria)

boxplot(MaxOD~Bacteria, data=max_od_nc, xlab="Species", ylab="Substrate")

max_od_nc_lm<- lm(MaxOD~Bacteria, data=max_od_nc)
drop1(max_od_nc_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Bacteria
# Df Sum of Sq      RSS     AIC F value   Pr(>F)    
# <none>                0.000428 -179.64                     
# Bacteria  5  0.047018 0.047446 -104.89  263.65 7.81e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(max_od_nc_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df   Sum Sq   Mean Sq F value   Pr(>F)    
# Bacteria   5 0.047018 0.0094036  263.65 7.81e-12 ***
#   Residuals 12 0.000428 0.0000357                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(max_od_nc_lm)
# Call:
#   lm(formula = MaxOD ~ Bacteria, data = max_od_nc)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0120000 -0.0020000 -0.0001667  0.0013333  0.0140000 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.132000   0.003448  38.283 6.50e-14 ***
#   BacteriaDSM_14838 -0.126333   0.004876 -25.908 6.67e-12 ***
#   BacteriaDSM_19555 -0.023333   0.004876  -4.785 0.000445 ***
#   BacteriaB5482      0.007000   0.004876   1.436 0.176686    
# BacteriaDSM_18836  0.028667   0.004876   5.879 7.49e-05 ***
#   BacteriaDSM_1447   0.010000   0.004876   2.051 0.062788 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.005972 on 12 degrees of freedom
# Multiple R-squared:  0.991,	Adjusted R-squared:  0.9872 
# F-statistic: 263.7 on 5 and 12 DF,  p-value: 7.81e-12

nctuk<- aov(max_od_nc_lm)
TukeyHSD(nctuk)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = max_od_nc_lm)
# 
# $Bacteria
# diff          lwr          upr     p adj
# DSM_14838-V975      -0.12633333 -0.142712259 -0.109954408 0.0000000
# DSM_19555-V975      -0.02333333 -0.039712259 -0.006954408 0.0045867
# B5482-V975           0.00700000 -0.009378925  0.023378925 0.7070857
# DSM_18836-V975       0.02866667  0.012287741  0.045045592 0.0008179
# DSM_1447-V975        0.01000000 -0.006378925  0.026378925 0.3710288
# DSM_19555-DSM_14838  0.10300000  0.086621075  0.119378925 0.0000000
# B5482-DSM_14838      0.13333333  0.116954408  0.149712259 0.0000000
# DSM_18836-DSM_14838  0.15500000  0.138621075  0.171378925 0.0000000
# DSM_1447-DSM_14838   0.13633333  0.119954408  0.152712259 0.0000000
# B5482-DSM_19555      0.03033333  0.013954408  0.046712259 0.0004911
# DSM_18836-DSM_19555  0.05200000  0.035621075  0.068378925 0.0000021
# DSM_1447-DSM_19555   0.03333333  0.016954408  0.049712259 0.0002032
# DSM_18836-B5482      0.02166667  0.005287741  0.038045592 0.0080616
# DSM_1447-B5482       0.00300000 -0.013378925  0.019378925 0.9877152
# DSM_1447-DSM_18836  -0.01866667 -0.035045592 -0.002287741 0.0226696


#Tukey HSD plot where any relationship with a p<0.01 was marked in red for significance
par(mar=c(5,6,4,1)+.1)
plot(TukeyHSD(nctuk, conf.level=.95), las = 2)

par(mar=c(5,6,4,1)+.1)
species_map <- c(
  "DSM_18836" = "B. xyl",
  "V975"     = "B. ova",
  "B5482"    = "B. theta",
  "DSM_19555" = "B. ster",
  "DSM_1447"  = "P. vul",
  "DSM_14838" = "B. cel"
)

tuk_plate_nc <- TukeyHSD(nctuk, conf.level = 0.95)$Bacteria

tuk_df <- as.data.frame(tuk_plate_nc)
tuk_df$Comparison <- rownames(tuk_plate_nc)

# Replace species codes with names
for (code in names(species_map)) {
  tuk_df$Comparison <- gsub(code, species_map[code], tuk_df$Comparison)
}

# Identify significance
tuk_df$Significant <- tuk_df$`p adj` < 0.01


#Makes it an order comparison
tuk_df$Comparison <- factor(
  tuk_df$Comparison,
  levels = rev(tuk_df$Comparison)
)

# Colors
cols <- ifelse(tuk_df$Significant, "red", "black")

#make the yplot
y_pos <- seq_len(nrow(tuk_df))

# Plot points
plot(
  tuk_df$diff,
  y_pos,
  xlim = range(c(tuk_df$lwr, tuk_df$upr)),
  yaxt = "n",
  ylab = "",
  xlab = "Differences in Mean Levels between Species Grown on NO CHO in 96-well Plate",
  pch = 16,
  col = cols,
  cex = 0.9
)

# Add confidence intervals
segments(
  tuk_df$lwr,
  y_pos,
  tuk_df$upr,
  y_pos,
  col = cols,
  lwd = 2
)

# Add vertical zero line
abline(v = 0, lty = 2)

# Add species comparison labels
axis(
  2,
  at = y_pos,
  labels = tuk_df$Comparison,
  las = 2,
  cex.axis = 0.8
)




#Glucose
max_od_glucose <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_glucose.txt", header = TRUE, sep = "\t",
                        stringsAsFactors = TRUE)
str(max_od_glucose)
max_od_glucose$Bacteria <- factor(max_od_glucose$Bacteria, levels= c('V975', 'DSM_14838', 'DSM_19555', 'B5482', 'DSM_18836', 'DSM_1447'))
str(max_od_glucose$Bacteria)
table(max_od_glucose$Bacteria)

boxplot(MaxOD~Bacteria, data=max_od_glucose, xlab="Species", ylab="Substrate")

max_od_glucose_lm<- lm(MaxOD~Bacteria, data=max_od_glucose)
drop1(max_od_glucose_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Bacteria
# Df Sum of Sq      RSS      AIC F value    Pr(>F)    
# <none>                0.001159 -161.705                      
# Bacteria  5   0.11428 0.115439  -88.889  236.58 1.485e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(max_od_glucose_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Bacteria   5 0.114280 0.0228559  236.58 1.485e-11 ***
#   Residuals 12 0.001159 0.0000966                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(max_od_glucose_lm)

# Call:
#   lm(formula = MaxOD ~ Bacteria, data = max_od_glucose)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.010000 -0.005917 -0.002167  0.005583  0.016000 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.373333   0.005675  65.788  < 2e-16 ***
#   BacteriaDSM_14838 -0.140667   0.008025 -17.528 6.47e-10 ***
#   BacteriaDSM_19555  0.063667   0.008025   7.933 4.10e-06 ***
#   BacteriaB5482      0.114333   0.008025  14.246 7.00e-09 ***
#   BacteriaDSM_18836  0.016000   0.008025   1.994  0.06942 .  
# BacteriaDSM_1447  -0.031667   0.008025  -3.946  0.00194 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.009829 on 12 degrees of freedom
# Multiple R-squared:   0.99,	Adjusted R-squared:  0.9858 
# F-statistic: 236.6 on 5 and 12 DF,  p-value: 1.485e-11

glucosetuk<- aov(max_od_glucose_lm)
TukeyHSD(glucosetuk)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = max_od_glucose_lm)
# 
# $Bacteria
# diff         lwr          upr     p adj
# DSM_14838-V975      -0.14066667 -0.16762343 -0.113709907 0.0000000
# DSM_19555-V975       0.06366667  0.03670991  0.090623426 0.0000471
# B5482-V975           0.11433333  0.08737657  0.141290093 0.0000001
# DSM_18836-V975       0.01600000 -0.01095676  0.042956759 0.3986082
# DSM_1447-V975       -0.03166667 -0.05862343 -0.004709907 0.0185805
# DSM_19555-DSM_14838  0.20433333  0.17737657  0.231290093 0.0000000
# B5482-DSM_14838      0.25500000  0.22804324  0.281956759 0.0000000
# DSM_18836-DSM_14838  0.15666667  0.12970991  0.183623426 0.0000000
# DSM_1447-DSM_14838   0.10900000  0.08204324  0.135956759 0.0000001
# B5482-DSM_19555      0.05066667  0.02370991  0.077623426 0.0004288
# DSM_18836-DSM_19555 -0.04766667 -0.07462343 -0.020709907 0.0007464
# DSM_1447-DSM_19555  -0.09533333 -0.12229009 -0.068376574 0.0000006
# DSM_18836-B5482     -0.09833333 -0.12529009 -0.071376574 0.0000005
# DSM_1447-B5482      -0.14600000 -0.17295676 -0.119043241 0.0000000
# DSM_1447-DSM_18836  -0.04766667 -0.07462343 -0.020709907 0.0007464


#Tukey HSD plot where any relationship with a p<0.01 was marked in red for significance
par(mar=c(5,6,4,1)+.1)
plot(TukeyHSD(glucosetuk, conf.level=.95), las = 2)

par(mar=c(5,6,4,1)+.1)
species_map <- c(
  "DSM_18836" = "B. xyl",
  "V975"     = "B. ova",
  "B5482"    = "B. theta",
  "DSM_19555" = "B. ster",
  "DSM_1447"  = "P. vul",
  "DSM_14838" = "B. cel"
)

tuk_plate_glucose <- TukeyHSD(glucosetuk, conf.level = 0.95)$Bacteria

tuk_df <- as.data.frame(tuk_plate_glucose)
tuk_df$Comparison <- rownames(tuk_plate_glucose)

# Replace species codes with names
for (code in names(species_map)) {
  tuk_df$Comparison <- gsub(code, species_map[code], tuk_df$Comparison)
}

# Identify significance
tuk_df$Significant <- tuk_df$`p adj` < 0.01


#Makes it an order comparison
tuk_df$Comparison <- factor(
  tuk_df$Comparison,
  levels = rev(tuk_df$Comparison)
)

# Colors
cols <- ifelse(tuk_df$Significant, "red", "black")

#make the yplot
y_pos <- seq_len(nrow(tuk_df))

# Plot points
plot(
  tuk_df$diff,
  y_pos,
  xlim = range(c(tuk_df$lwr, tuk_df$upr)),
  yaxt = "n",
  ylab = "",
  xlab = "Differences in Mean Levels between Species Grown on Glucose in 96-well Plate",
  pch = 16,
  col = cols,
  cex = 0.9
)

# Add confidence intervals
segments(
  tuk_df$lwr,
  y_pos,
  tuk_df$upr,
  y_pos,
  col = cols,
  lwd = 2
)

# Add vertical zero line
abline(v = 0, lty = 2)

# Add species comparison labels
axis(
  2,
  at = y_pos,
  labels = tuk_df$Comparison,
  las = 2,
  cex.axis = 0.8
)


#Apple pectin
max_od_ap <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_ap.txt", header = TRUE, sep = "\t",
                             stringsAsFactors = TRUE)
str(max_od_ap)
max_od_ap$Bacteria <- factor(max_od_ap$Bacteria, levels= c('V975', 'DSM_14838', 'DSM_19555', 'B5482', 'DSM_18836', 'DSM_1447'))
str(max_od_ap$Bacteria)
table(max_od_ap$Bacteria)

boxplot(MaxOD~Bacteria, data=max_od_ap, xlab="Species", ylab="Substrate")

max_od_ap_lm<- lm(MaxOD~Bacteria, data=max_od_ap)
drop1(max_od_ap_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Bacteria
# Df Sum of Sq      RSS     AIC F value    Pr(>F)    
# <none>                0.001632 -155.55                      
# Bacteria  5   0.30679 0.308426  -71.20  451.17 3.197e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(max_od_ap_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df   Sum Sq  Mean Sq F value    Pr(>F)    
# Bacteria   5 0.306794 0.061359  451.17 3.197e-13 ***
#   Residuals 12 0.001632 0.000136                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(max_od_ap_lm)
# 
# Call:
#   lm(formula = MaxOD ~ Bacteria, data = max_od_ap)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.013667 -0.004250 -0.001833  0.002583  0.024333 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.392000   0.006733  58.221 4.35e-16 ***
#   BacteriaDSM_14838 -0.339333   0.009522 -35.637 1.52e-13 ***
#   BacteriaDSM_19555  0.080000   0.009522   8.402 2.27e-06 ***
#   BacteriaB5482     -0.121667   0.009522 -12.778 2.40e-08 ***
#   BacteriaDSM_18836 -0.132000   0.009522 -13.863 9.54e-09 ***
#   BacteriaDSM_1447  -0.139333   0.009522 -14.633 5.16e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01166 on 12 degrees of freedom
# Multiple R-squared:  0.9947,	Adjusted R-squared:  0.9925 
# F-statistic: 451.2 on 5 and 12 DF,  p-value: 3.197e-13

applepectintuk<- aov(max_od_ap_lm)
TukeyHSD(applepectintuk)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = max_od_ap_lm)
# 
# $Bacteria
# diff         lwr         upr     p adj
# DSM_14838-V975      -0.339333333 -0.37131666 -0.30735001 0.0000000
# DSM_19555-V975       0.080000000  0.04801668  0.11198332 0.0000262
# B5482-V975          -0.121666667 -0.15364999 -0.08968334 0.0000003
# DSM_18836-V975      -0.132000000 -0.16398332 -0.10001668 0.0000001
# DSM_1447-V975       -0.139333333 -0.17131666 -0.10735001 0.0000001
# DSM_19555-DSM_14838  0.419333333  0.38735001  0.45131666 0.0000000
# B5482-DSM_14838      0.217666667  0.18568334  0.24964999 0.0000000
# DSM_18836-DSM_14838  0.207333333  0.17535001  0.23931666 0.0000000
# DSM_1447-DSM_14838   0.200000000  0.16801668  0.23198332 0.0000000
# B5482-DSM_19555     -0.201666667 -0.23364999 -0.16968334 0.0000000
# DSM_18836-DSM_19555 -0.212000000 -0.24398332 -0.18001668 0.0000000
# DSM_1447-DSM_19555  -0.219333333 -0.25131666 -0.18735001 0.0000000
# DSM_18836-B5482     -0.010333333 -0.04231666  0.02164999 0.8781857
# DSM_1447-B5482      -0.017666667 -0.04964999  0.01431666 0.4699814
# DSM_1447-DSM_18836  -0.007333333 -0.03931666  0.02464999 0.9676177


#Tukey HSD plot where any relationship with a p<0.01 was marked in red for significance
par(mar=c(5,6,4,1)+.1)
plot(TukeyHSD(applepectintuk, conf.level=.95), las = 2)

par(mar=c(5,6,4,1)+.1)
species_map <- c(
  "DSM_18836" = "B. xyl",
  "V975"     = "B. ova",
  "B5482"    = "B. theta",
  "DSM_19555" = "B. ster",
  "DSM_1447"  = "P. vul",
  "DSM_14838" = "B. cel"
)

tuk_plate_apple <- TukeyHSD(applepectintuk, conf.level = 0.95)$Bacteria

tuk_df <- as.data.frame(tuk_plate_apple)
tuk_df$Comparison <- rownames(tuk_plate_apple)

# Replace species codes with names
for (code in names(species_map)) {
  tuk_df$Comparison <- gsub(code, species_map[code], tuk_df$Comparison)
}

# Identify significance
tuk_df$Significant <- tuk_df$`p adj` < 0.01


#Makes it an order comparison
tuk_df$Comparison <- factor(
  tuk_df$Comparison,
  levels = rev(tuk_df$Comparison)
)

# Colors
cols <- ifelse(tuk_df$Significant, "red", "black")

#make the yplot
y_pos <- seq_len(nrow(tuk_df))

# Plot points
plot(
  tuk_df$diff,
  y_pos,
  xlim = range(c(tuk_df$lwr, tuk_df$upr)),
  yaxt = "n",
  ylab = "",
  xlab = "Differences in Mean Levels between Species Grown on Apple Pectin in 96-well Plate",
  pch = 16,
  col = cols,
  cex = 0.9
)

# Add confidence intervals
segments(
  tuk_df$lwr,
  y_pos,
  tuk_df$upr,
  y_pos,
  col = cols,
  lwd = 2
)

# Add vertical zero line
abline(v = 0, lty = 2)

# Add species comparison labels
axis(
  2,
  at = y_pos,
  labels = tuk_df$Comparison,
  las = 2,
  cex.axis = 0.8
)


#beta-D-glucan
max_od_bg <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_bg.txt", header = TRUE, sep = "\t",
                        stringsAsFactors = TRUE)
str(max_od_bg)
max_od_bg$Bacteria <- factor(max_od_bg$Bacteria, levels= c('V975', 'DSM_14838', 'DSM_19555', 'B5482', 'DSM_18836', 'DSM_1447'))
str(max_od_bg$Bacteria)
table(max_od_bg$Bacteria)

boxplot(MaxOD~Bacteria, data=max_od_bg, xlab="Species", ylab="Substrate")

max_od_bg_lm<- lm(MaxOD~Bacteria, data=max_od_bg)
drop1(max_od_bg_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Bacteria
# Df Sum of Sq    RSS      AIC F value    Pr(>F)    
# <none>                0.0040 -139.413                      
# Bacteria  5    0.2917 0.2957  -71.958  175.02 8.829e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(max_od_bg_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df Sum Sq  Mean Sq F value    Pr(>F)    
# Bacteria   5 0.2917 0.058339  175.02 8.829e-11 ***
#   Residuals 12 0.0040 0.000333                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(max_od_bg_lm)

# Call:
#   lm(formula = MaxOD ~ Bacteria, data = max_od_bg)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.02500 -0.00800 -0.00050  0.00425  0.03300 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.30433    0.01054   28.87 1.85e-12 ***
#   BacteriaDSM_14838 -0.28600    0.01491  -19.18 2.27e-10 ***
#   BacteriaDSM_19555 -0.19567    0.01491  -13.13 1.77e-08 ***
#   BacteriaB5482     -0.17233    0.01491  -11.56 7.32e-08 ***
#   BacteriaDSM_18836  0.07767    0.01491    5.21 0.000218 ***
#   BacteriaDSM_1447  -0.21233    0.01491  -14.24 7.01e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01826 on 12 degrees of freedom
# Multiple R-squared:  0.9865,	Adjusted R-squared:  0.9808 
# F-statistic:   175 on 5 and 12 DF,  p-value: 8.829e-11

betatuk<- aov(max_od_bg_lm)
TukeyHSD(betatuk)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = max_od_bg_lm)
# 
# $Bacteria
# diff         lwr         upr     p adj
# DSM_14838-V975      -0.28600000 -0.33607184 -0.23592816 0.0000000 ***
# DSM_19555-V975      -0.19566667 -0.24573850 -0.14559483 0.0000002 ***
# B5482-V975          -0.17233333 -0.22240517 -0.12226150 0.0000009 ***
# DSM_18836-V975       0.07766667  0.02759483  0.12773850 0.0023105 **
# DSM_1447-V975       -0.21233333 -0.26240517 -0.16226150 0.0000001 ***
# DSM_19555-DSM_14838  0.09033333  0.04026150  0.14040517 0.0006233 ***
# B5482-DSM_14838      0.11366667  0.06359483  0.16373850 0.0000700 ***
# DSM_18836-DSM_14838  0.36366667  0.31359483  0.41373850 0.0000000 ***
# DSM_1447-DSM_14838   0.07366667  0.02359483  0.12373850 0.0035548 ***
# B5482-DSM_19555      0.02333333 -0.02673850  0.07340517 0.6334212
# DSM_18836-DSM_19555  0.27333333  0.22326150  0.32340517 0.0000000 ***
# DSM_1447-DSM_19555  -0.01666667 -0.06673850  0.03340517 0.8649141
# DSM_18836-B5482      0.25000000  0.19992816  0.30007184 0.0000000 ***
# DSM_1447-B5482      -0.04000000 -0.09007184  0.01007184 0.1500657
# DSM_1447-DSM_18836  -0.29000000 -0.34007184 -0.23992816 0.0000000 ***


#Tukey HSD plot where any relationship with a p<0.01 was marked in red for significance
par(mar=c(5,6,4,1)+.1)
plot(TukeyHSD(betatuk, conf.level=.95), las = 2)

par(mar=c(5,6,4,1)+.1)
species_map <- c(
  "DSM_18836" = "B. xyl",
  "V975"     = "B. ova",
  "B5482"    = "B. theta",
  "DSM_19555" = "B. ster",
  "DSM_1447"  = "P. vul",
  "DSM_14838" = "B. cel"
)

tuk_plate_beta <- TukeyHSD(betatuk, conf.level = 0.95)$Bacteria

tuk_df <- as.data.frame(tuk_plate_beta)
tuk_df$Comparison <- rownames(tuk_plate_beta)

# Replace species codes with names
for (code in names(species_map)) {
  tuk_df$Comparison <- gsub(code, species_map[code], tuk_df$Comparison)
}

# Identify significance
tuk_df$Significant <- tuk_df$`p adj` < 0.01


#Makes it an order comparison
tuk_df$Comparison <- factor(
  tuk_df$Comparison,
  levels = rev(tuk_df$Comparison)
)

# Colors
cols <- ifelse(tuk_df$Significant, "red", "black")

#make the yplot
y_pos <- seq_len(nrow(tuk_df))

# Plot points
plot(
  tuk_df$diff,
  y_pos,
  xlim = range(c(tuk_df$lwr, tuk_df$upr)),
  yaxt = "n",
  ylab = "",
  xlab = "Differences in Mean Levels between Species Grown on BetaDglucan in 96-well Plate",
  pch = 16,
  col = cols,
  cex = 0.9
)

# Add confidence intervals
segments(
  tuk_df$lwr,
  y_pos,
  tuk_df$upr,
  y_pos,
  col = cols,
  lwd = 2
)

# Add vertical zero line
abline(v = 0, lty = 2)

# Add species comparison labels
axis(
  2,
  at = y_pos,
  labels = tuk_df$Comparison,
  las = 2,
  cex.axis = 0.8
)


#galactomannan
max_od_gala <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_gala.txt", header = TRUE, sep = "\t",
                        stringsAsFactors = TRUE)
str(max_od_gala)
max_od_gala$Bacteria <- factor(max_od_gala$Bacteria, levels= c('V975', 'DSM_14838', 'DSM_19555', 'B5482', 'DSM_18836', 'DSM_1447'))
str(max_od_gala$Bacteria)
table(max_od_gala$Bacteria)

boxplot(MaxOD~Bacteria, data=max_od_gala, xlab="Species", ylab="Substrate")

max_od_gala_lm<- lm(MaxOD~Bacteria, data=max_od_gala)
drop1(max_od_gala_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Bacteria
# Df Sum of Sq      RSS      AIC F value    Pr(>F)    
# <none>                0.018937 -111.426                      
# Bacteria  5   0.14419 0.163129  -82.665  18.274 3.067e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(max_od_gala_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Bacteria   5 0.144192 0.0288383  18.274 3.067e-05 ***
#   Residuals 12 0.018937 0.0015781                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(max_od_gala_lm)
# Call:
#   lm(formula = MaxOD ~ Bacteria, data = max_od_gala)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.100667 -0.001583  0.000167  0.001583  0.089333 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.29667    0.02294  12.935 2.09e-08 ***
#   BacteriaDSM_14838 -0.06300    0.03244  -1.942   0.0759 .  
# BacteriaDSM_19555 -0.22433    0.03244  -6.916 1.61e-05 ***
#   BacteriaB5482     -0.23433    0.03244  -7.225 1.05e-05 ***
#   BacteriaDSM_18836 -0.21400    0.03244  -6.598 2.55e-05 ***
#   BacteriaDSM_1447  -0.18667    0.03244  -5.755 9.09e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03973 on 12 degrees of freedom
# Multiple R-squared:  0.8839,	Adjusted R-squared:  0.8355 
# F-statistic: 18.27 on 5 and 12 DF,  p-value: 3.067e-05

galatuk<- aov(max_od_gala_lm)
TukeyHSD(galatuk)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = max_od_gala_lm)
# 
# $Bacteria
# diff         lwr         upr     p adj
# DSM_14838-V975      -0.06300000 -0.17194892  0.04594892 0.4244011
# DSM_19555-V975      -0.22433333 -0.33328225 -0.11538441 0.0001817 ***
# B5482-V975          -0.23433333 -0.34328225 -0.12538441 0.0001191 ***
# DSM_18836-V975      -0.21400000 -0.32294892 -0.10505108 0.0002845 ***
# DSM_1447-V975       -0.18666667 -0.29561559 -0.07771775 0.0009874 ***
# DSM_19555-DSM_14838 -0.16133333 -0.27028225 -0.05238441 0.0033741 ***
# B5482-DSM_14838     -0.17133333 -0.28028225 -0.06238441 0.0020603 ***
# DSM_18836-DSM_14838 -0.15100000 -0.25994892 -0.04205108 0.0056749 ***
# DSM_1447-DSM_14838  -0.12366667 -0.23261559 -0.01471775 0.0232678
# B5482-DSM_19555     -0.01000000 -0.11894892  0.09894892 0.9995197
# DSM_18836-DSM_19555  0.01033333 -0.09861559  0.11928225 0.9994369
# DSM_1447-DSM_19555   0.03766667 -0.07128225  0.14661559 0.8463881
# DSM_18836-B5482      0.02033333 -0.08861559  0.12928225 0.9866485
# DSM_1447-B5482       0.04766667 -0.06128225  0.15661559 0.6879937
# DSM_1447-DSM_18836   0.02733333 -0.08161559  0.13628225 0.9531762

#Tukey HSD plot where any relationship with a p<0.01 was marked in red for significance
par(mar=c(5,6,4,1)+.1)
plot(TukeyHSD(galatuk, conf.level=.95), las = 2)

par(mar=c(5,6,4,1)+.1)
species_map <- c(
  "DSM_18836" = "B. xyl",
  "V975"     = "B. ova",
  "B5482"    = "B. theta",
  "DSM_19555" = "B. ster",
  "DSM_1447"  = "P. vul",
  "DSM_14838" = "B. cel"
)

tuk_plate_gala <- TukeyHSD(galatuk, conf.level = 0.95)$Bacteria

tuk_df <- as.data.frame(tuk_plate_gala)
tuk_df$Comparison <- rownames(tuk_plate_gala)

# Replace species codes with names
for (code in names(species_map)) {
  tuk_df$Comparison <- gsub(code, species_map[code], tuk_df$Comparison)
}

# Identify significance
tuk_df$Significant <- tuk_df$`p adj` < 0.01


#Makes it an order comparison
tuk_df$Comparison <- factor(
  tuk_df$Comparison,
  levels = rev(tuk_df$Comparison)
)

# Colors
cols <- ifelse(tuk_df$Significant, "red", "black")

#make the yplot
y_pos <- seq_len(nrow(tuk_df))

# Plot points
plot(
  tuk_df$diff,
  y_pos,
  xlim = range(c(tuk_df$lwr, tuk_df$upr)),
  yaxt = "n",
  ylab = "",
  xlab = "Differences in Mean Levels between Species Grown on Galactomannan in 96-well Plate",
  pch = 16,
  col = cols,
  cex = 0.9
)

# Add confidence intervals
segments(
  tuk_df$lwr,
  y_pos,
  tuk_df$upr,
  y_pos,
  col = cols,
  lwd = 2
)

# Add vertical zero line
abline(v = 0, lty = 2)

# Add species comparison labels
axis(
  2,
  at = y_pos,
  labels = tuk_df$Comparison,
  las = 2,
  cex.axis = 0.8
)



#potato starch
max_od_ps <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_ps.txt", header = TRUE, sep = "\t",
                          stringsAsFactors = TRUE)
str(max_od_ps)
max_od_ps$Bacteria <- factor(max_od_ps$Bacteria, levels= c('V975', 'DSM_14838', 'DSM_19555', 'B5482', 'DSM_18836', 'DSM_1447'))
str(max_od_ps$Bacteria)
table(max_od_ps$Bacteria)

boxplot(MaxOD~Bacteria, data=max_od_ps, xlab="Species", ylab="Substrate")

max_od_ps_lm<- lm(MaxOD~Bacteria, data=max_od_ps)
drop1(max_od_ps_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Bacteria
# Df Sum of Sq      RSS      AIC F value    Pr(>F)    
# <none>                0.013532 -117.475                      
# Bacteria  5   0.13877 0.152300  -83.901  24.612 6.406e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(max_od_ps_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Bacteria   5 0.138768 0.0277537  24.612 6.406e-06 ***
#   Residuals 12 0.013532 0.0011277                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(max_od_ps_lm)
# Call:
#   lm(formula = MaxOD ~ Bacteria, data = max_od_ps)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.072667 -0.006667 -0.001000  0.008750  0.078333 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.345667   0.019388  17.829 5.31e-10 ***
#   BacteriaDSM_14838 -0.107667   0.027419  -3.927  0.00201 ** 
#   BacteriaDSM_19555  0.150333   0.027419   5.483  0.00014 ***
#   BacteriaB5482      0.107000   0.027419   3.902  0.00210 ** 
#   BacteriaDSM_18836  0.006667   0.027419   0.243  0.81200    
# BacteriaDSM_1447  -0.049667   0.027419  -1.811  0.09516 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03358 on 12 degrees of freedom
# Multiple R-squared:  0.9111,	Adjusted R-squared:  0.8741 
# F-statistic: 24.61 on 5 and 12 DF,  p-value: 6.406e-06

pstuk<- aov(max_od_ps_lm)
TukeyHSD(pstuk)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = max_od_ps_lm)
# 
# $Bacteria
# diff         lwr          upr     p adj
# DSM_14838-V975      -0.107666667 -0.19976346 -0.015569872 0.0191869
# DSM_19555-V975       0.150333333  0.05823654  0.242430128 0.0015030 ***
# B5482-V975           0.107000000  0.01490321  0.199096795 0.0199914
# DSM_18836-V975       0.006666667 -0.08543013  0.098763461 0.9998493
# DSM_1447-V975       -0.049666667 -0.14176346  0.042430128 0.4938298
# DSM_19555-DSM_14838  0.258000000  0.16590321  0.350096795 0.0000081 ***
# B5482-DSM_14838      0.214666667  0.12256987  0.306763461 0.0000537 ***
# DSM_18836-DSM_14838  0.114333333  0.02223654  0.206430128 0.0127358
# DSM_1447-DSM_14838   0.058000000 -0.03409679  0.150096795 0.3412963
# B5482-DSM_19555     -0.043333333 -0.13543013  0.048763461 0.6246899
# DSM_18836-DSM_19555 -0.143666667 -0.23576346 -0.051569872 0.0022039 ***
# DSM_1447-DSM_19555  -0.200000000 -0.29209679 -0.107903205 0.0001085
# DSM_18836-B5482     -0.100333333 -0.19243013 -0.008236539 0.0301531
# DSM_1447-B5482      -0.156666667 -0.24876346 -0.064569872 0.0010515 ***
# DSM_1447-DSM_18836  -0.056333333 -0.14843013  0.035763461 0.3692305

#Tukey HSD plot where any relationship with a p<0.01 was marked in red for significance
par(mar=c(5,6,4,1)+.1)
plot(TukeyHSD(pstuk, conf.level=.95), las = 2)

par(mar=c(5,6,4,1)+.1)
species_map <- c(
  "DSM_18836" = "B. xyl",
  "V975"     = "B. ova",
  "B5482"    = "B. theta",
  "DSM_19555" = "B. ster",
  "DSM_1447"  = "P. vul",
  "DSM_14838" = "B. cel"
)

tuk_plate_ps <- TukeyHSD(pstuk, conf.level = 0.95)$Bacteria

tuk_df <- as.data.frame(tuk_plate_ps)
tuk_df$Comparison <- rownames(tuk_plate_ps)

# Replace species codes with names
for (code in names(species_map)) {
  tuk_df$Comparison <- gsub(code, species_map[code], tuk_df$Comparison)
}

# Identify significance
tuk_df$Significant <- tuk_df$`p adj` < 0.01


#Makes it an order comparison
tuk_df$Comparison <- factor(
  tuk_df$Comparison,
  levels = rev(tuk_df$Comparison)
)

# Colors
cols <- ifelse(tuk_df$Significant, "red", "black")

#make the yplot
y_pos <- seq_len(nrow(tuk_df))

# Plot points
plot(
  tuk_df$diff,
  y_pos,
  xlim = range(c(tuk_df$lwr, tuk_df$upr)),
  yaxt = "n",
  ylab = "",
  xlab = "Differences in Mean Levels between Species Grown on Potato Starch in 96-well Plate",
  pch = 16,
  col = cols,
  cex = 0.9
)

# Add confidence intervals
segments(
  tuk_df$lwr,
  y_pos,
  tuk_df$upr,
  y_pos,
  col = cols,
  lwd = 2
)

# Add vertical zero line
abline(v = 0, lty = 2)

# Add species comparison labels
axis(
  2,
  at = y_pos,
  labels = tuk_df$Comparison,
  las = 2,
  cex.axis = 0.8
)


#xylan
max_od_xylan <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_xylan.txt", header = TRUE, sep = "\t",
                        stringsAsFactors = TRUE)
str(max_od_xylan)
max_od_xylan$Bacteria <- factor(max_od_xylan$Bacteria, levels= c('V975', 'DSM_14838', 'DSM_19555', 'B5482', 'DSM_18836', 'DSM_1447'))
str(max_od_xylan$Bacteria)
table(max_od_xylan$Bacteria)

boxplot(MaxOD~Bacteria, data=max_od_xylan, xlab="Species", ylab="Substrate")

max_od_xylan_lm<- lm(MaxOD~Bacteria, data=max_od_xylan)
drop1(max_od_xylan_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Bacteria
# Df Sum of Sq      RSS      AIC F value    Pr(>F)    
# <none>                0.000091 -207.445                      
# Bacteria  5   0.19764 0.197736  -79.202  5193.6 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(max_od_xylan_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df   Sum Sq  Mean Sq F value    Pr(>F)    
# Bacteria   5 0.197644 0.039529  5193.6 < 2.2e-16 ***
#   Residuals 12 0.000091 0.000008                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(max_od_xylan_lm)
# Call:
#   lm(formula = MaxOD ~ Bacteria, data = max_od_xylan)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0040000 -0.0019167  0.0006667  0.0012500  0.0050000 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.350000   0.001593  219.74  < 2e-16 ***
#   BacteriaDSM_14838 -0.254333   0.002253 -112.91  < 2e-16 ***
#   BacteriaDSM_19555 -0.228333   0.002253 -101.37  < 2e-16 ***
#   BacteriaB5482     -0.206000   0.002253  -91.45  < 2e-16 ***
#   BacteriaDSM_18836  0.005000   0.002253    2.22   0.0465 *  
#   BacteriaDSM_1447  -0.127000   0.002253  -56.38 6.39e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.002759 on 12 degrees of freedom
# Multiple R-squared:  0.9995,	Adjusted R-squared:  0.9993 
# F-statistic:  5194 on 5 and 12 DF,  p-value: < 2.2e-16


xylantuk<- aov(max_od_xylan_lm)
TukeyHSD(xylantuk)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = max_od_xylan_lm)
# 
# $Bacteria
# diff          lwr         upr     p adj
# DSM_14838-V975      -0.25433333 -0.261899540 -0.24676713 0.0000000
# DSM_19555-V975      -0.22833333 -0.235899540 -0.22076713 0.0000000
# B5482-V975          -0.20600000 -0.213566206 -0.19843379 0.0000000
# DSM_18836-V975       0.00500000 -0.002566206  0.01256621 0.2967326
# DSM_1447-V975       -0.12700000 -0.134566206 -0.11943379 0.0000000
# DSM_19555-DSM_14838  0.02600000  0.018433794  0.03356621 0.0000009
# B5482-DSM_14838      0.04833333  0.040767127  0.05589954 0.0000000
# DSM_18836-DSM_14838  0.25933333  0.251767127  0.26689954 0.0000000
# DSM_1447-DSM_14838   0.12733333  0.119767127  0.13489954 0.0000000
# B5482-DSM_19555      0.02233333  0.014767127  0.02989954 0.0000046
# DSM_18836-DSM_19555  0.23333333  0.225767127  0.24089954 0.0000000
# DSM_1447-DSM_19555   0.10133333  0.093767127  0.10889954 0.0000000
# DSM_18836-B5482      0.21100000  0.203433794  0.21856621 0.0000000
# DSM_1447-B5482       0.07900000  0.071433794  0.08656621 0.0000000
# DSM_1447-DSM_18836  -0.13200000 -0.139566206 -0.12443379 0.0000000

#Tukey HSD plot where any relationship with a p<0.01 was marked in red for significance
par(mar=c(5,6,4,1)+.1)
plot(TukeyHSD(xylantuk, conf.level=.95), las = 2)

par(mar=c(5,6,4,1)+.1)
species_map <- c(
  "DSM_18836" = "B. xyl",
  "V975"     = "B. ova",
  "B5482"    = "B. theta",
  "DSM_19555" = "B. ster",
  "DSM_1447"  = "P. vul",
  "DSM_14838" = "B. cel"
)

tuk_plate_xylan <- TukeyHSD(xylantuk, conf.level = 0.95)$Bacteria

tuk_df <- as.data.frame(tuk_plate_xylan)
tuk_df$Comparison <- rownames(tuk_plate_xylan)

# Replace species codes with names
for (code in names(species_map)) {
  tuk_df$Comparison <- gsub(code, species_map[code], tuk_df$Comparison)
}

# Identify significance
tuk_df$Significant <- tuk_df$`p adj` < 0.01


#Makes it an order comparison
tuk_df$Comparison <- factor(
  tuk_df$Comparison,
  levels = rev(tuk_df$Comparison)
)

# Colors
cols <- ifelse(tuk_df$Significant, "red", "black")

#make the yplot
y_pos <- seq_len(nrow(tuk_df))

# Plot points
plot(
  tuk_df$diff,
  y_pos,
  xlim = range(c(tuk_df$lwr, tuk_df$upr)),
  yaxt = "n",
  ylab = "",
  xlab = "Differences in Mean Levels between Species Grown on Xylan in 96-well Plate",
  pch = 16,
  col = cols,
  cex = 0.9
)

# Add confidence intervals
segments(
  tuk_df$lwr,
  y_pos,
  tuk_df$upr,
  y_pos,
  col = cols,
  lwd = 2
)

# Add vertical zero line
abline(v = 0, lty = 2)

# Add species comparison labels
axis(
  2,
  at = y_pos,
  labels = tuk_df$Comparison,
  las = 2,
  cex.axis = 0.8
)


#xyloglucan
max_od_xylo <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/max_xyloglucan.txt", header = TRUE, sep = "\t",
                           stringsAsFactors = TRUE)
str(max_od_xylo)
max_od_xylo$Bacteria <- factor(max_od_xylo$Bacteria, levels= c('V975', 'DSM_14838', 'DSM_19555', 'B5482', 'DSM_18836', 'DSM_1447'))
str(max_od_xylo$Bacteria)
table(max_od_xylo$Bacteria)

boxplot(MaxOD~Bacteria, data=max_od_xylo, xlab="Species", ylab="Substrate")

max_od_xylo_lm<- lm(MaxOD~Bacteria, data=max_od_xylo)
drop1(max_od_xylo_lm, t='F')
# Single term deletions
# 
# Model:
#   MaxOD ~ Bacteria
# Df Sum of Sq      RSS      AIC F value    Pr(>F)    
# <none>                0.003853 -140.089                      
# Bacteria  5    0.1161 0.119950  -88.199  72.323 1.544e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(max_od_xylo_lm)
# Analysis of Variance Table
# 
# Response: MaxOD
# Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Bacteria   5 0.116098 0.0232196  72.323 1.544e-08 ***
#   Residuals 12 0.003853 0.0003211                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(max_od_xylo_lm)
# Call:
#   lm(formula = MaxOD ~ Bacteria, data = max_od_xylo)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.020667 -0.008167 -0.005000  0.011250  0.037333 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.29267    0.01034  28.291 2.36e-12 ***
#   BacteriaDSM_14838 -0.25367    0.01463 -17.339 7.34e-10 ***
#   BacteriaDSM_19555 -0.17800    0.01463 -12.167 4.15e-08 ***
#   BacteriaB5482     -0.14100    0.01463  -9.638 5.33e-07 ***
#   BacteriaDSM_18836 -0.11067    0.01463  -7.564 6.64e-06 ***
#   BacteriaDSM_1447  -0.20533    0.01463 -14.035 8.29e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01792 on 12 degrees of freedom
# Multiple R-squared:  0.9679,	Adjusted R-squared:  0.9545 
# F-statistic: 72.32 on 5 and 12 DF,  p-value: 1.544e-08

xylotuk<- aov(max_od_xylo_lm)
TukeyHSD(xylotuk)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = max_od_xylo_lm)
# 
# $Bacteria
# diff           lwr         upr     p adj
# DSM_14838-V975      -0.25366667 -0.3028076947 -0.20452564 0.0000000
# DSM_19555-V975      -0.17800000 -0.2271410280 -0.12885897 0.0000005
# B5482-V975          -0.14100000 -0.1901410280 -0.09185897 0.0000062
# DSM_18836-V975      -0.11066667 -0.1598076947 -0.06152564 0.0000758
# DSM_1447-V975       -0.20533333 -0.2544743614 -0.15619231 0.0000001
# DSM_19555-DSM_14838  0.07566667  0.0265256386  0.12480769 0.0024548
# B5482-DSM_14838      0.11266667  0.0635256386  0.16180769 0.0000634
# DSM_18836-DSM_14838  0.14300000  0.0938589720  0.19214103 0.0000054
# DSM_1447-DSM_14838   0.04833333 -0.0008076947  0.09747436 0.0548360
# B5482-DSM_19555      0.03700000 -0.0121410280  0.08614103 0.1900171
# DSM_18836-DSM_19555  0.06733333  0.0181923053  0.11647436 0.0061927
# DSM_1447-DSM_19555  -0.02733333 -0.0764743614  0.02180769 0.4630584
# DSM_18836-B5482      0.03033333 -0.0188076947  0.07947436 0.3604418
# DSM_1447-B5482      -0.06433333 -0.1134743614 -0.01519231 0.0087028
# DSM_1447-DSM_18836  -0.09466667 -0.1438076947 -0.04552564 0.0003412

#Tukey HSD plot where any relationship with a p<0.01 was marked in red for significance
par(mar=c(5,6,4,1)+.1)
plot(TukeyHSD(xylotuk, conf.level=.95), las = 2)

par(mar=c(5,6,4,1)+.1)
species_map <- c(
  "DSM_18836" = "B. xyl",
  "V975"     = "B. ova",
  "B5482"    = "B. theta",
  "DSM_19555" = "B. ster",
  "DSM_1447"  = "P. vul",
  "DSM_14838" = "B. cel"
)

tuk_plate_xylo <- TukeyHSD(xylotuk, conf.level = 0.95)$Bacteria

tuk_df <- as.data.frame(tuk_plate_xylo)
tuk_df$Comparison <- rownames(tuk_plate_xylo)

# Replace species codes with names
for (code in names(species_map)) {
  tuk_df$Comparison <- gsub(code, species_map[code], tuk_df$Comparison)
}

# Identify significance
tuk_df$Significant <- tuk_df$`p adj` < 0.01


#Makes it an order comparison
tuk_df$Comparison <- factor(
  tuk_df$Comparison,
  levels = rev(tuk_df$Comparison)
)

# Colors
cols <- ifelse(tuk_df$Significant, "red", "black")

#make the yplot
y_pos <- seq_len(nrow(tuk_df))

# Plot points
plot(
  tuk_df$diff,
  y_pos,
  xlim = range(c(tuk_df$lwr, tuk_df$upr)),
  yaxt = "n",
  ylab = "",
  xlab = "Differences in Mean Levels between Species Grown on Xyloglucan in 96-well Plate",
  pch = 16,
  col = cols,
  cex = 0.9
)

# Add confidence intervals
segments(
  tuk_df$lwr,
  y_pos,
  tuk_df$upr,
  y_pos,
  col = cols,
  lwd = 2
)

# Add vertical zero line
abline(v = 0, lty = 2)

# Add species comparison labels
axis(
  2,
  at = y_pos,
  labels = tuk_df$Comparison,
  las = 2,
  cex.axis = 0.8
)






###############
#growth rate for individual pure culture
###############

Plate_8Sub_Growth <- read.table(file = "~/Masters/Statistical Analysis/qPCR_8substrate/plate8_growthrate_corrected.txt", header = TRUE, sep = "\t",
                   stringsAsFactors = TRUE)
Plate_8Sub_Growth$Substrate <- factor(Plate_8Sub_Growth$Substrate, levels= c("NOCHO", "Glucose", "Apple_Pectin", 'BDglucan', 'Galactomannans', 'Potato_Starch', 'Xylan', 'Xyloglucan')) #order so comparison is to NOCHO
str(Plate_8Sub_Growth$Substrate)
table(Plate_8Sub_Growth$Substrate)


Plate_8Sub_Growth$Bacterium <- factor(Plate_8Sub_Growth$Bacterium, levels= c( 'SynCom', "DSM14838", "V975", "DSM19555", "DSM18836", 'DSM1447'))
str(Plate_8Sub_Growth$Bacterium)
table(Plate_8Sub_Growth$Bacterium)


str(Plate_8Sub_Growth)
# 'data.frame':	216 obs. of  6 variables:
#   $ Number    : int  25 26 27 55 56 57 73 74 75 88 ...
# $ Well      : Factor w/ 120 levels "A1","A10","A11",..: 1 5 6 32 34 36 50 52 54 70 ...
# $ Bacterium : Factor w/ 6 levels "DSM14838","V975",..: 2 2 2 5 5 5 3 3 3 4 ...
# $ Substrate : Factor w/ 8 levels "NOCHO","Glucose",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Replicate : Factor w/ 3 levels "R1","R2","R3": 1 2 3 1 2 3 1 2 3 1 ...
# $ GrowthRate: num  0.1618 0.1557 0.1317 0.1184 0.0937 ...

boxplot(GrowthRate~Substrate, data=Plate_8Sub_Growth, xlab="Substrate", ylab="Growth Rate")

#####
#individual growth rates: all species are compared to B. ovatus (V975) as it had growth on all species
#####

#NO CHO
NOCHO_GrowthRate <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/plate8_growthrate_NOCHO.txt", header = TRUE, sep = "\t",
                           stringsAsFactors = TRUE)
str(NOCHO_GrowthRate)
NOCHO_GrowthRate$Bacterium <- factor(NOCHO_GrowthRate$Bacterium, levels= c('V975', 'DSM14838', 'DSM19555', 'B5482', 'DSM18836', 'DSM1447')) 
str(NOCHO_GrowthRate$Bacterium)
table(NOCHO_GrowthRate$Bacterium)

boxplot(GrowthRate~Bacterium, data=NOCHO_GrowthRate, xlab="Species", ylab="Substrate")

NOCHO_GrowthRate_lm<- lm(GrowthRate~Bacterium, data=NOCHO_GrowthRate)
drop1(NOCHO_GrowthRate_lm, t='F')
# Single term deletions
# 
# Model:
#   GrowthRate ~ Bacterium
# Df Sum of Sq      RSS      AIC F value    Pr(>F)    
# <none>                 0.004356 -112.165                      
# Bacterium  4   0.02738 0.031736  -90.375  15.715 0.0002588 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 


anova(NOCHO_GrowthRate_lm)
# Analysis of Variance Table
# 
# Response: GrowthRate
# Df    Sum Sq   Mean Sq F value    Pr(>F)    
# Bacterium  4 0.0273802 0.0068451  15.715 0.0002588 ***
#   Residuals 10 0.0043557 0.0004356                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(NOCHO_GrowthRate_lm)
# Call:
#   lm(formula = GrowthRate ~ Bacterium, data = NOCHO_GrowthRate)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0280667 -0.0124500  0.0005333  0.0099667  0.0291333 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.14973    0.01205  12.427  2.1e-07 ***
#   BacteriumDSM19555 -0.06187    0.01704  -3.631  0.00461 ** 
#   BacteriumB5482    -0.02490    0.01704  -1.461  0.17465    
# BacteriumDSM18836  0.05853    0.01704   3.435  0.00639 ** 
#   BacteriumDSM1447  -0.04967    0.01704  -2.915  0.01544 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02087 on 10 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.8628,	Adjusted R-squared:  0.8079 
# F-statistic: 15.72 on 4 and 10 DF,  p-value: 0.0002588


#GLUCOSE
Glucose_GrowthRate <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/plate8_growthrate_Glucose.txt", header = TRUE, sep = "\t",
                               stringsAsFactors = TRUE)
str(Glucose_GrowthRate)
Glucose_GrowthRate$Bacterium <- factor(Glucose_GrowthRate$Bacterium, levels= c('V975', 'DSM14838', 'DSM19555', 'B5482', 'DSM18836', 'DSM1447'))
str(Glucose_GrowthRate$Bacterium)
table(Glucose_GrowthRate$Bacterium)

boxplot(GrowthRate~Bacterium, data=Glucose_GrowthRate, xlab="Species", ylab="Growth Rate")

Glucose_GrowthRate_lm<- lm(GrowthRate~Bacterium, data=Glucose_GrowthRate)
drop1(Glucose_GrowthRate_lm, t='F')
# Single term deletions
# 
# Model:
#   GrowthRate ~ Bacterium
# Df Sum of Sq      RSS     AIC F value Pr(>F)
# <none>                 0.045324 -95.717               
# Bacterium  5  0.041779 0.087103 -93.959  2.2123 0.1206


anova(Glucose_GrowthRate_lm)
# Analysis of Variance Table
# 
# Response: GrowthRate
# Df   Sum Sq   Mean Sq F value Pr(>F)
# Bacterium  5 0.041779 0.0083558  2.2123 0.1206
# Residuals 12 0.045324 0.0037770               

summary(Glucose_GrowthRate_lm)
# Call:
#   lm(formula = GrowthRate ~ Bacterium, data = Glucose_GrowthRate)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.108200 -0.015242  0.004133  0.009850  0.114767 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.17657    0.03548   4.976 0.000322 ***
#   BacteriumDSM14838  0.01673    0.05018   0.333 0.744535    
# BacteriumDSM19555  0.09273    0.05018   1.848 0.089382 .  
# BacteriumB5482     0.09997    0.05018   1.992 0.069603 .  
# BacteriumDSM18836  0.10453    0.05018   2.083 0.059287 .  
# BacteriumDSM1447  -0.00510    0.05018  -0.102 0.920725    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.06146 on 12 degrees of freedom
# Multiple R-squared:  0.4797,	Adjusted R-squared:  0.2628 
# F-statistic: 2.212 on 5 and 12 DF,  p-value: 0.1206

#APPLE PECTIN
AP_GrowthRate <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/plate8_growthrate_ap.txt", header = TRUE, sep = "\t",
                                 stringsAsFactors = TRUE)
str(AP_GrowthRate)
AP_GrowthRate$Bacterium <- factor(AP_GrowthRate$Bacterium, levels= c('V975', 'DSM14838', 'DSM19555', 'B5482', 'DSM18836', 'DSM1447'))
str(AP_GrowthRate$Bacterium)
table(AP_GrowthRate$Bacterium)

boxplot(GrowthRate~Bacterium, data=AP_GrowthRate, xlab="Species", ylab="Substrate")

AP_GrowthRate_lm<- lm(GrowthRate~Bacterium, data=AP_GrowthRate)
drop1(Glucose_GrowthRate_lm, t='F')
# Single term deletions
# 
# Model:
#   GrowthRate ~ Bacterium
# Df Sum of Sq      RSS     AIC F value Pr(>F)
# <none>                 0.045324 -95.717               
# Bacterium  5  0.041779 0.087103 -93.959  2.2123 0.1206

anova(AP_GrowthRate_lm)
# Analysis of Variance Table
# 
# Response: GrowthRate
# Df  Sum Sq  Mean Sq F value    Pr(>F)    
# Bacterium  4 0.96762 0.241905  243.47 6.453e-10 ***
#   Residuals 10 0.00994 0.000994                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(AP_GrowthRate_lm)
# Call:
#   lm(formula = GrowthRate ~ Bacterium, data = AP_GrowthRate)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.04133 -0.01718 -0.00760  0.01373  0.06007 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.39550    0.01820  21.732 9.51e-10 ***
#   BacteriumDSM19555  0.48207    0.02574  18.731 4.07e-09 ***
#   BacteriumB5482    -0.20900    0.02574  -8.121 1.03e-05 ***
#   BacteriumDSM18836 -0.15457    0.02574  -6.006 0.000131 ***
#   BacteriumDSM1447  -0.15097    0.02574  -5.866 0.000158 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03152 on 10 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.9898,	Adjusted R-squared:  0.9858 
# F-statistic: 243.5 on 4 and 10 DF,  p-value: 6.453e-10


#Beta-D-glucan
bg_GrowthRate <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/plate8_growthrate_bg.txt", header = TRUE, sep = "\t",
                                 stringsAsFactors = TRUE)
str(bg_GrowthRate)
bg_GrowthRate$Bacterium <- factor(bg_GrowthRate$Bacterium, levels= c('V975','DSM14838',  'DSM19555', 'B5482', 'DSM18836', 'DSM1447'))
str(bg_GrowthRate$Bacterium)
table(bg_GrowthRate$Bacterium)

boxplot(GrowthRate~Bacterium, data=bg_GrowthRate, xlab="Species", ylab="GrowthRate")

bg_GrowthRate_lm<- lm(GrowthRate~Bacterium, data=bg_GrowthRate)
drop1(bg_GrowthRate_lm, t='F')
# Single term deletions
# 
# Model:
#   GrowthRate ~ Bacterium
# Df Sum of Sq      RSS     AIC F value   Pr(>F)   
# <none>                 0.034511 -67.109                    
# Bacterium  4   0.24025 0.274765 -48.138  13.924 0.001119 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


anova(bg_GrowthRate_lm)
# Analysis of Variance Table
# 
# Response: GrowthRate
# Df   Sum Sq  Mean Sq F value   Pr(>F)   
# Bacterium  4 0.240254 0.060064  13.924 0.001119 **
#   Residuals  8 0.034511 0.004314                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(bg_GrowthRate_lm)
# Call:
#   lm(formula = GrowthRate ~ Bacterium, data = bg_GrowthRate)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.10177 -0.02563  0.00000  0.01860  0.10307 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.39193    0.03792  10.336 6.63e-06 ***
#   BacteriumDSM19555 -0.34357    0.05363  -6.407 0.000208 ***
#   BacteriumB5482    -0.13853    0.05363  -2.583 0.032451 *  
#   BacteriumDSM18836 -0.05487    0.05363  -1.023 0.336194    
# BacteriumDSM1447  -0.33113    0.07584  -4.366 0.002393 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.06568 on 8 degrees of freedom
# (14 observations deleted due to missingness)
# Multiple R-squared:  0.8744,	Adjusted R-squared:  0.8116 
# F-statistic: 13.92 on 4 and 8 DF,  p-value: 0.001119

bg_tukM<- aov(bg_GrowthRate_lm)
TukeyHSD(bg_tukM, conf.level=.95)

par(mar=c(5,6,4,1)+.1)
plot(TukeyHSD(bg_tukM, conf.level=.95), las = 2)


#GALACTO
Gala_GrowthRate <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/plate8_growthrate_galactomannan.txt", header = TRUE, sep = "\t",
                                 stringsAsFactors = TRUE)
str(Gala_GrowthRate)
Gala_GrowthRate$Bacterium <- factor(Gala_GrowthRate$Bacterium, levels= c('V975', 'DSM14838', 'DSM19555', 'B5482', 'DSM18836', 'DSM1447'))
str(Gala_GrowthRate$Bacterium)
table(Gala_GrowthRate$Bacterium)

boxplot(GrowthRate~Bacterium, data=Gala_GrowthRate, xlab="Species", ylab="Growth Rate")

Gala_GrowthRate_lm<- lm(GrowthRate~Bacterium, data=Gala_GrowthRate)
drop1(Gala_GrowthRate_lm, t='F')
# Single term deletions
# 
# Model:
#   GrowthRate ~ Bacterium
# Df Sum of Sq      RSS     AIC F value    Pr(>F)    
# <none>                 0.001166 -74.566                      
# Bacterium  2  0.095287 0.096453 -38.823  245.27 1.764e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


anova(Gala_GrowthRate_lm)
# Analysis of Variance Table
# 
# Response: GrowthRate
# Df   Sum Sq  Mean Sq F value    Pr(>F)    
# Bacterium  2 0.095287 0.047644  245.27 1.764e-06 ***
#   Residuals  6 0.001166 0.000194                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(Gala_GrowthRate_lm)

# 
# Call:
#   lm(formula = GrowthRate ~ Bacterium, data = Gala_GrowthRate)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0218667 -0.0047000 -0.0002667  0.0037000  0.0221333 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.281867   0.008047   35.03 3.61e-08 ***
#   BacteriumDSM14838 -0.194867   0.011380  -17.12 2.54e-06 ***
#   BacteriumDSM1447  -0.235867   0.011380  -20.73 8.21e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01394 on 6 degrees of freedom
# (9 observations deleted due to missingness)
# Multiple R-squared:  0.9879,	Adjusted R-squared:  0.9839 
# F-statistic: 245.3 on 2 and 6 DF,  p-value: 1.764e-06


#Potato Starch
ps_GrowthRate <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/plate8_growthrate_ps.txt", header = TRUE, sep = "\t",
                                 stringsAsFactors = TRUE)
str(ps_GrowthRate)
ps_GrowthRate$Bacterium <- factor(Glucose_GrowthRate$Bacterium, levels= c('V975', 'DSM14838', 'DSM19555', 'B5482', 'DSM18836', 'DSM1447'))
str(ps_GrowthRate$Bacterium)
table(ps_GrowthRate$Bacterium)

boxplot(GrowthRate~Bacterium, data=ps_GrowthRate, xlab="Species", ylab="Growth Rate")

ps_GrowthRate_lm<- lm(GrowthRate~Bacterium, data=ps_GrowthRate)
drop1(ps_GrowthRate_lm, t='F')
# Single term deletions
# 
# Model:
#   GrowthRate ~ Bacterium
# Df Sum of Sq      RSS     AIC F value    Pr(>F)    
# <none>                 0.046204 -95.371                      
# Bacterium  5   0.25155 0.297751 -71.834  13.066 0.0001653 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


anova(ps_GrowthRate_lm)
# Analysis of Variance Table
# 
# Response: GrowthRate
# Df   Sum Sq  Mean Sq F value    Pr(>F)    
# Bacterium  5 0.251547 0.050309  13.066 0.0001653 ***
#   Residuals 12 0.046204 0.003850                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(ps_GrowthRate_lm)
# Call:
#   lm(formula = GrowthRate ~ Bacterium, data = ps_GrowthRate)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.09837 -0.03421  0.02213  0.03707  0.07123 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.39073    0.03583  10.907 1.39e-07 ***
#   BacteriumDSM14838 -0.20133    0.05066  -3.974  0.00185 ** 
#   BacteriumDSM19555  0.07040    0.05066   1.390  0.18991    
# BacteriumB5482     0.17943    0.05066   3.542  0.00406 ** 
#   BacteriumDSM18836 -0.01490    0.05066  -0.294  0.77371    
# BacteriumDSM1447   0.09273    0.05066   1.830  0.09213 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.06205 on 12 degrees of freedom
# Multiple R-squared:  0.8448,	Adjusted R-squared:  0.7802 
# F-statistic: 13.07 on 5 and 12 DF,  p-value: 0.0001653

#XYLAN
xylan_GrowthRate <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/plate8_growthrate_xylan.txt", header = TRUE, sep = "\t",
                                 stringsAsFactors = TRUE)
str(xylan_GrowthRate)
xylan_GrowthRate$Bacterium <- factor(xylan_GrowthRate$Bacterium, levels= c('V975', 'DSM14838', 'DSM19555', 'B5482', 'DSM18836', 'DSM1447'))
str(xylan_GrowthRate$Bacterium)
table(xylan_GrowthRate$Bacterium)

boxplot(GrowthRate~Bacterium, data=xylan_GrowthRate, xlab="Species", ylab="Growth Rate")

xylan_GrowthRate_lm<- lm(GrowthRate~Bacterium, data=xylan_GrowthRate)
drop1(xylan_GrowthRate_lm, t='F')
# Single term deletions
# 
# Model:
#   GrowthRate ~ Bacterium
# Df Sum of Sq     RSS     AIC F value    Pr(>F)    
# <none>                 0.04284 -96.732                      
# Bacterium  5    0.3627 0.40554 -66.272   20.32 1.767e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


anova(xylan_GrowthRate_lm)
# Analysis of Variance Table
# 
# Response: GrowthRate
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# Bacterium  5 0.36270 0.07254   20.32 1.767e-05 ***
#   Residuals 12 0.04284 0.00357                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1  

summary(xylan_GrowthRate_lm)
# Call:
#   lm(formula = GrowthRate ~ Bacterium, data = xylan_GrowthRate)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.093400 -0.017450 -0.000783  0.021600  0.088567 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.49260    0.03450  14.280 6.81e-09 ***
#   BacteriumDSM14838 -0.45470    0.04878  -9.321 7.62e-07 ***
#   BacteriumDSM19555 -0.38770    0.04878  -7.947 4.02e-06 ***
#   BacteriumB5482    -0.29133    0.04878  -5.972 6.49e-05 ***
#   BacteriumDSM18836 -0.26937    0.04878  -5.522 0.000132 ***
#   BacteriumDSM1447  -0.29360    0.04878  -6.018 6.04e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05975 on 12 degrees of freedom
# Multiple R-squared:  0.8944,	Adjusted R-squared:  0.8504 
# F-statistic: 20.32 on 5 and 12 DF,  p-value: 1.767e-05


#XYLOGLUCAN

xylo_GrowthRate <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/plate8_growthrate_xyloglucan.txt", header = TRUE, sep = "\t",
                                 stringsAsFactors = TRUE)
str(xylo_GrowthRate)
xylo_GrowthRate$Bacterium <- factor(Glucose_GrowthRate$Bacterium, levels= c('V975', 'DSM14838', 'DSM19555', 'B5482', 'DSM18836', 'DSM1447'))
str(xylo_GrowthRate$Bacterium)
table(xylo_GrowthRate$Bacterium)

boxplot(GrowthRate~Bacterium, data=xylo_GrowthRate, xlab="Species", ylab="Substrate")

xylo_GrowthRate_lm<- lm(GrowthRate~Bacterium, data=xylo_GrowthRate)
drop1(xylo_GrowthRate_lm, t='F')
# Single term deletions
# 
# Model:
#   GrowthRate ~ Bacterium
# Df Sum of Sq      RSS     AIC F value   Pr(>F)   
# <none>                 0.024439 -71.595                    
# Bacterium  4   0.11622 0.140657 -56.843  9.5109 0.003923 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


anova(xylo_GrowthRate_lm)
# Analysis of Variance Table
# 
# Response: GrowthRate
# Df   Sum Sq   Mean Sq F value   Pr(>F)   
# Bacterium  4 0.116218 0.0290545  9.5109 0.003923 **
#   Residuals  8 0.024439 0.0030549                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(xylo_GrowthRate_lm)
# Call:
#   lm(formula = GrowthRate ~ Bacterium, data = xylo_GrowthRate)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.096933 -0.020167  0.000000  0.008067  0.088867 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.30527    0.03191   9.566 1.18e-05 ***
#   BacteriumDSM19555 -0.23313    0.04513  -5.166 0.000857 ***
#   BacteriumB5482    -0.02523    0.04513  -0.559 0.591367    
# BacteriumDSM18836 -0.12610    0.04513  -2.794 0.023404 *  
#   BacteriumDSM1447  -0.22107    0.06382  -3.464 0.008520 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05527 on 8 degrees of freedom
# (5 observations deleted due to missingness)
# Multiple R-squared:  0.8263,	Adjusted R-squared:  0.7394 
# F-statistic: 9.511 on 4 and 8 DF,  p-value: 0.003923





####################################################
#SYN COM
####################################################

#Correlation pH, Max OD, 16S

correlation_three <- read.table(file = "~/Masters/Statistical Analysis/qPCR_8substrate/threevariablecorrelation.txt", header = TRUE, sep = "\t",
                                stringsAsFactors = TRUE)
correlation_three$Substrate <- factor(correlation_three$Substrate, levels= c("NO_CHO", "Glucose", "Apple_Pectin", 'B_D_glucan', 'Galactomannan', 'Potato_Starch', 'Xylan', 'Xyloglucan')) #reorder for the substrate to compare to NO CHO
table(correlation_three$Substrate)

structure(correlation_three)

cor(correlation_three[,c('MAXOD', 'copperml', 'ph_check')], method='pearson')

cor(correlation_three[,c('MAXOD', 'copperml', 'ph_check')], method='pearson')


pairs(
  correlation_three[, c("MAXOD", "copperml", "ph_check")],
  pch = 19,
  main = "Scatterplot Matrix of MAXOD, Copper, and pH"
)




plot1<-ggplot(correlation_three, aes(MAXOD, ph_check)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(x='Maximum OD', y='pH')
plot1
plot2<-ggplot(correlation_three, aes(MAXOD, copperml)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()+
  labs(x='Maximum OD', y='16S rRNA Gene Copies per ml')

plot3<-ggplot(correlation_three, aes(copperml, ph_check)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()+
  labs(y='pH', x='16S rRNA Gene Copies per ml')
plot3
# Combine them in one row
combined_plot <- plot1 | plot2 | plot3

# Display
combined_plot

#Correlation pH, Max OD, 16S, growth rate

correlation_four <- read.table(file = "~/Masters/Statistical Analysis/qPCR_8substrate/fourvariablecorrelation.txt", header = TRUE, sep = "\t",
                                stringsAsFactors = TRUE)
correlation_four$Substrate <- factor(correlation_four$Substrate, levels= c("NO_CHO", "Glucose", "Apple_Pectin", 'B_D_glucan', 'Galactomannan', 'Potato_Starch', 'Xylan', 'Xyloglucan'))
table(correlation_four$Substrate)

structure(correlation_four)

cor(correlation_four[,c('MAXOD', 'copperml', 'ph_check', 'growth')], method='pearson')

cor(correlation_four[,c('MAXOD', 'copperml', 'ph_check', 'growth')], method='pearson')
#           MAXOD    copperml    ph_check     growth
# MAXOD     1.0000000  0.47899144 -0.71419438  0.4021868
# copperml  0.4789914  1.00000000 -0.02912269  0.3434600
# ph_check -0.7141944 -0.02912269  1.00000000 -0.3552570
# growth    0.4021868  0.34346001 -0.35525697  1.0000000

pairs(
  correlation_four[, c("MAXOD", "copperml", "ph_check", 'growth')],
  pch = 19,
  main = "Scatterplot Matrix of MAXOD, growth rate, Copper, and pH"
)



#plot all seperately so can then combine plots at the end
plot1<-ggplot(correlation_four, aes(MAXOD, ph_check)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(x='Maximum OD', y='pH')
plot1
plot2<-ggplot(correlation_four, aes(MAXOD, copperml)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()+
  labs(x='Maximum OD', y='16S rRNA Gene Copies per ml')

plot3<-ggplot(correlation_four, aes(copperml, ph_check)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()+
  labs(y='pH', x='16S rRNA Gene Copies per ml')
plot3

plot4<-ggplot(correlation_four, aes(growth, ph_check)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(x='Growth Rate', y='pH')
plot4
plot5<-ggplot(correlation_four, aes(growth, MAXOD)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()+
  labs(x='Growth Rate', y='Maximum OD')
plot5
plot6<-ggplot(correlation_four, aes(growth, copperml)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal()+
  labs(x='Growth Rate', y='16S rRNA Gene Copies per ml')
plot6


# Combine them in one row
combined_plot <- plot1 | plot2 | plot3 
combined_plot2<-  plot4 | plot5 | plot6
# Display
combined_plot
combined_plot2

#linear model & anova to determine if substrate and dietary fibre was significant

#16S for substrate only

relative_16S<- lm(log(copperml)~Substrate, data=correlation_three)
summary(relative_16S)
# Call:
#   lm(formula = log(copperml) ~ Substrate, data = correlation_three)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.40865 -0.11610  0.00132  0.07569  0.52962 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             20.8962     0.1269 164.637  < 2e-16 ***
#   SubstrateGlucose         0.2310     0.1795   1.287 0.216331    
# SubstrateApple_Pectin    1.0449     0.1795   5.821 2.60e-05 ***
#   SubstrateB_D_glucan      0.1304     0.1795   0.726 0.478176    
# SubstrateGalactomannan   0.1924     0.1795   1.072 0.299713    
# SubstratePotato_Starch   0.1255     0.1795   0.699 0.494504    
# SubstrateXylan           1.0497     0.1795   5.848 2.47e-05 ***
#   SubstrateXyloglucan      0.8667     0.1795   4.828 0.000185 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2198 on 16 degrees of freedom
# Multiple R-squared:  0.8456,	Adjusted R-squared:  0.778 
# F-statistic: 12.51 on 7 and 16 DF,  p-value: 1.988e-05
anova(relative_16S)
# Analysis of Variance Table
# 
# Response: log(copperml)
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Substrate  7 4.2337 0.60482  12.515 1.988e-05 ***
#   Residuals 16 0.7732 0.04833                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#maximmum OD for the synthetic community
relative_maxod<- lm(MAXOD~Substrate, data=correlation_four)
anova(relative_maxod)
# Analysis of Variance Table
# 
# Response: MAXOD
# Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Substrate  7 0.149384 0.0213405  16.201 3.627e-06 ***
#   Residuals 16 0.021076 0.0013172                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(relative_maxod)
# Call:
#   lm(formula = MAXOD ~ Substrate, data = correlation_four)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.08617 -0.01121 -0.00225  0.01158  0.06933 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.15500    0.02095   7.397 1.51e-06 ***
#   SubstrateGlucose        0.25500    0.02963   8.605 2.13e-07 ***
#   SubstrateApple_Pectin   0.27083    0.02963   9.139 9.46e-08 ***
#   SubstrateB_D_glucan     0.18611    0.02963   6.280 1.10e-05 ***
#   SubstrateGalactomannan  0.17561    0.02963   5.926 2.13e-05 ***
#   SubstratePotato_Starch  0.18978    0.02963   6.404 8.73e-06 ***
#   SubstrateXylan          0.21306    0.02963   7.190 2.15e-06 ***
#   SubstrateXyloglucan     0.23678    0.02963   7.990 5.64e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03629 on 16 degrees of freedom
# Multiple R-squared:  0.8764,	Adjusted R-squared:  0.8223 
# F-statistic:  16.2 on 7 and 16 DF,  p-value: 3.627e-06

maxOD_plateTukey <- aov(relative_maxod)

par(mar=c(5,6,8,1)+.1)
plot(TukeyHSD(maxOD_plateTukey, conf.level=.95), las = 2)



#TukeyHSD plot for When symbol bar is red it has a p<0.01
par(mar=c(6,4,4,2)+.1)


tuk_plate <- TukeyHSD(maxOD_plateTukey, conf.level = 0.95)$Substrate

tuk_df <- as.data.frame(tuk_plate)
tuk_df$Comparison <- rownames(tuk_plate)

# Replace species codes with names
for (code in names(species_map)) {
  tuk_df$Comparison <- gsub(code, species_map[code], tuk_df$Comparison)
}

# Identify significance
tuk_df$Significant <- tuk_df$`p adj` < 0.05


#Makes it an order comparison
tuk_df$Comparison <- factor(
  tuk_df$Comparison,
  levels = rev(tuk_df$Comparison)
)

# Colors
cols <- ifelse(tuk_df$Significant, "red", "black")

#make the yplot
y_pos <- seq_len(nrow(tuk_df))

# Plot points
plot(
  tuk_df$diff,
  y_pos,
  xlim = range(c(tuk_df$lwr, tuk_df$upr)),
  yaxt = "n",
  ylab = "",
  xlab = "Differences in Mean Levels between Species in 96-well Plate",
  pch = 16,
  col = cols,
  cex = 0.9
)

# Add confidence intervals
segments(
  tuk_df$lwr,
  y_pos,
  tuk_df$upr,
  y_pos,
  col = cols,
  lwd = 2
)

# Add vertical zero line
abline(v = 0, lty = 2)

# Add species comparison labels
axis(
  2,
  at = y_pos,
  labels = tuk_df$Comparison,
  las = 2,
  cex.axis = 0.8
)


#########
#For relative abundance
#########
relative_abundance <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/qPCR_8substrate/relative_abundance.txt", header = TRUE, sep = "\t",
                          stringsAsFactors = TRUE)
str(relative_abundance)
# 'data.frame':	144 obs. of  4 variables:
#   $ DNACodes  : int  6 6 6 6 6 6 8 8 8 8 ...
# $ Substrate : Factor w/ 9 levels "Apple_Pectin",..: 5 5 5 5 5 5 5 5 5 5 ...
# $ Species   : Factor w/ 6 levels "B5482","DSM_1447",..: 2 6 4 1 3 5 2 6 4 1 ...
# $ Percentage: num  15.25 43.32 8.6 17.08 1.94 ...



#Reorder to compare to NO CHO
relative_abundance$Substrate <- factor(relative_abundance$Substrate, levels= c("NOCHO", "Glucose", "Apple_Pectin", 'B_Glucan', 'Galactomannan', 'Potato_Starch', 'Xylan', 'Xyloglucan'))
str(relative_abundance$Substrate)
table(relative_abundance$Substrate)

#Reorder species alphabetically
relative_abundance$Species <- factor(relative_abundance$Species, levels= c('DSM_14838', 'V975', 'B5482','DSM_19555','DSM_18836','DSM_1447'))
str(relative_abundance$Species)
table(relative_abundance$Species)

boxplot(Percentage~Species, data=relative_abundance, xlab="Species", ylab="Percentage")
boxplot(Percentage~Substrate, data=relative_abundance, xlab="Species", ylab="Substrate")

relative_ab_add<- lm(Percentage~Species+Substrate, data=relative_abundance)
drop1(relative_ab_add, t='F')
# Single term deletions
# 
# Model:
#   Percentage ~ Species + Substrate
# Df Sum of Sq   RSS    AIC F value Pr(>F)    
# <none>                 10161 638.93                   
# Species    5     19281 29442 782.13  49.716 <2e-16 ***
#   Substrate  7         0 10161 624.93   0.000      1    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

relative_substrate<- lm(Percentage~Substrate, data=relative_abundance)
drop1(relative_substrate, t='F')
# Single term deletions
# 
# Model:
#   Percentage ~ Substrate
# Df Sum of Sq   RSS    AIC F value Pr(>F)
# <none>                 29442 782.13               
# Substrate  7 9.375e-05 29442 768.13       0      1

#substrate has non-significant impact

#Impact of species on relative abundance
relative_species<- lm(Percentage~Species, data=relative_abundance)
drop1(relative_species, t='F')
# Single term deletions
# 
# Model:
#   Percentage ~ Species
# Df Sum of Sq   RSS    AIC F value    Pr(>F)    
# <none>               10161 624.93                      
# Species  5     19281 29442 768.13  52.373 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(relative_species)
# Analysis of Variance Table
# 
# Response: Percentage
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Species     5  19281  3856.2  52.373 < 2.2e-16 ***
#   Residuals 138  10161    73.6                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(relative_species) #all have significant interacions
# Call:
#   lm(formula = Percentage ~ Species, data = relative_abundance)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -17.4878  -5.1094  -0.8081   3.4587  30.7193 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         2.579      1.752   1.472 0.143173    
# SpeciesV975        37.229      2.477  15.029  < 2e-16 ***
#   SpeciesB5482       16.289      2.477   6.576 9.22e-10 ***
#   SpeciesDSM_19555   14.185      2.477   5.726 6.15e-08 ***
#   SpeciesDSM_18836    8.168      2.477   3.297 0.001241 ** 
#   SpeciesDSM_1447     8.656      2.477   3.494 0.000639 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 8.581 on 138 degrees of freedom
# Multiple R-squared:  0.6549,	Adjusted R-squared:  0.6424 
# F-statistic: 52.37 on 5 and 138 DF,  p-value: < 2.2e-16

relative_species_aov_plate<- aov(relative_species) #anova format needed for TukeyHSD analysis
summary(relative_species_aov_plate)
# Df Sum Sq Mean Sq F value Pr(>F)    
# Species       5  19281    3856   52.37 <2e-16 ***
#   Residuals   138  10161      74                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(relative_species_aov_plate, conf.level=.95)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = relative_species)
# 
# $Species
#                         diff        lwr         upr     p adj
# V975-DSM_14838       37.2287099  30.069586  44.3878338 0.0000000 *
# B5482-DSM_14838      16.2890267   9.129903  23.4481506 0.0000000 *
# DSM_19555-DSM_14838  14.1845885   7.025465  21.3437124 0.0000009 *
# DSM_18836-DSM_14838   8.1679594   1.008836  15.3270833 0.0153995 .
# DSM_1447-DSM_14838    8.6556222   1.496498  15.8147461 0.0082348 .
# B5482-V975          -20.9396832 -28.098807 -13.7805593 0.0000000 *
# DSM_19555-V975      -23.0441214 -30.203245 -15.8849975 0.0000000 *
# DSM_18836-V975      -29.0607504 -36.219874 -21.9016266 0.0000000 *
# DSM_1447-V975       -28.5730877 -35.732212 -21.4139638 0.0000000 *
# DSM_19555-B5482      -2.1044381  -9.263562   5.0546857 0.9574745 
# DSM_18836-B5482      -8.1210672 -15.280191  -0.9619434 0.0163242 .
# DSM_1447-B5482       -7.6334045 -14.792528  -0.4742806 0.0293366 .
# DSM_18836-DSM_19555  -6.0166291 -13.175753   1.1424948 0.1535952 
# DSM_1447-DSM_19555   -5.5289663 -12.688090   1.6301575 0.2299652
# DSM_1447-DSM_18836    0.4876628  -6.671461   7.6467866 0.9999584

par(mar=c(5,6,4,1)+.1)
plot(TukeyHSD(relative_species_aov_plate, conf.level=.95), las = 2)



#TukeyHSD plot for When symbol bar is red it has a p<0.01
par(mar=c(5,6,4,1)+.1)
species_map <- c(
  "DSM_18836" = "B. xyl",
  "V975"     = "B. ova",
  "B5482"    = "B. theta",
  "DSM_19555" = "B. ster",
  "DSM_1447"  = "P. vul",
  "DSM_14838" = "B. cel"
)

tuk_plate <- TukeyHSD(relative_species_aov_plate, conf.level = 0.95)$Species

tuk_df <- as.data.frame(tuk_plate)
tuk_df$Comparison <- rownames(tuk_plate)

# Replace species codes with names
for (code in names(species_map)) {
  tuk_df$Comparison <- gsub(code, species_map[code], tuk_df$Comparison)
}

# Identify significance
tuk_df$Significant <- tuk_df$`p adj` < 0.05


#Makes it an order comparison
tuk_df$Comparison <- factor(
  tuk_df$Comparison,
  levels = rev(tuk_df$Comparison)
)

# Colors
cols <- ifelse(tuk_df$Significant, "red", "black")

#make the yplot
y_pos <- seq_len(nrow(tuk_df))

# Plot points
plot(
  tuk_df$diff,
  y_pos,
  xlim = range(c(tuk_df$lwr, tuk_df$upr)),
  yaxt = "n",
  ylab = "",
  xlab = "Differences in Mean Levels between Species in 96-well Plate",
  pch = 16,
  col = cols,
  cex = 0.9
)

# Add confidence intervals
segments(
  tuk_df$lwr,
  y_pos,
  tuk_df$upr,
  y_pos,
  col = cols,
  lwd = 2
)

# Add vertical zero line
abline(v = 0, lty = 2)

# Add species comparison labels
axis(
  2,
  at = y_pos,
  labels = tuk_df$Comparison,
  las = 2,
  cex.axis = 0.8
)


