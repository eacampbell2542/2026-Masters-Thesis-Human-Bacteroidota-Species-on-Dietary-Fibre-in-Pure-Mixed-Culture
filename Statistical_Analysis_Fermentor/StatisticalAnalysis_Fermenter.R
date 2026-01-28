# R version 4.4.3 (2025-02-28 ucrt) -- "Trophy Case"
# Date Created: # Tue Oct 14 16:28:17 2025 ------------------------------


#Project: Explore the relationship between Bacteroidota
#Author Ellie Campbell
#Email r01ec25@aberdeen.ac.uk

#Packages to load
library(ggplot2)


###########################################
#Total 16S rRNA gene copies per ml analysis 
###########################################

Fermenter_16S <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/Fermenter/16s_Fermenter_3.txt", header = TRUE, sep = "\t",
                            stringsAsFactors = TRUE)

str(Fermenter_16S)
Fermenter_16S$Substrate <- factor(Fermenter_16S$Substrate, levels= c('Mix', 'Apple_Pectin', 'Potato_Starch', 'Xylan', 'Glucose', 'Basal'))Fermenter_16S$Substrate <- factor(Fermenter_16S$Substrate, levels= c("Mix","Apple_Pectin",'Potato_Starch','Xylan','Glucose', 'Basal'))
str(Fermenter_16S$Substrate)

Fermenter_16S$Phase <- factor(Fermenter_16S$Phase, levels= c('Batch', 'Continous')) #sorted into phase
str(Fermenter_16S$Phase)
table(Fermenter_16S$Phase)

boxplot(Copies~Substrate, data=Fermenter_16S, xlab="Substrate", ylab="Copies per ml")
boxplot(Copies~Phase, data=Fermenter_16S, xlab="Phase", ylab="Copies per ml")

#simplified
Fermenter_Copies_phase<-lm(Copies~Phase, data=Fermenter_16S)
anova(Fermenter_Copies_phase)
summary(Fermenter_Copies_phase)

Fermenter_Copies_substrate<-lm(Copies~Substrate, data=Fermenter_16S)
anova(Fermenter_Copies_substrate)
summary(Fermenter_Copies_substrate)

Fermenter_Copies_both<-lm(log10(Copies)~Substrate*Phase, data=Fermenter_16S)
anova(Fermenter_Copies_both)
summary(Fermenter_Copies_both)


#Model selection
Fermenter_Copies_MAX<-lm(Copies~Substrate+Phase+Substrate:Phase, data=Fermenter_16S)
drop1(Fermenter_Copies_MAX, test='F')

# Single term deletions
# 
# Model:
#   Copies ~ Substrate + Phase + Substrate:Phase
# Df Sum of Sq        RSS    AIC F value    Pr(>F)    
# <none>                       3.3645e+19 1513.6                      
# Substrate:Phase  4 3.883e+19 7.2475e+19 1533.3  6.9248 0.0007442 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Both substrate and the phase have a significant effect
Fermenter_Copies_add<-lm(Copies~Substrate+Phase, data=Fermenter_16S)
drop1(Fermenter_Copies_add, test='F')
# Single term deletions
# 
# Model:
#   Copies ~ Substrate + Phase
# Df  Sum of Sq        RSS    AIC F value  Pr(>F)   
# <none>                  7.2475e+19 1533.3                   
# Substrate  5 5.7074e+19 1.2955e+20 1544.2  4.4100 0.00436 **
#   Phase      1 1.8685e+19 9.1160e+19 1539.5  7.2188 0.01200 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#took log of copies to make data more manageable
Fermenter_Copies_both<-lm(log(Copies)~Substrate*Phase, data=Fermenter_16S)
drop1(Fermenter_Copies_both, test='F')
# Single term deletions
# 
# Model:
#   log(Copies) ~ Substrate * Phase
# Df Sum of Sq    RSS     AIC F value    Pr(>F)    
# <none>                        5.567 -43.200                      
# Substrate:Phase  4    20.428 25.995   4.278  22.017 9.706e-08 ***


anova(Fermenter_Copies_both)
# Analysis of Variance Table
# 
# Response: log(Copies)
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# Substrate        5 22.0427  4.4085  40.164 2.430e-10 ***
#   Phase            1  6.4390  6.4390  58.663 1.207e-07 ***
#   Substrate:Phase  4 20.4282  5.1071  46.528 1.997e-10 ***
#   Residuals       22  2.4148  0.1098                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


summary(Fermenter_Copies_both)
# 
# Call:
#   lm(formula = log(Copies) ~ Substrate * Phase, data = Fermenter_16S)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.94691 -0.09825  0.00497  0.10812  0.71969 
# 
# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            18.2888     0.1913  95.613  < 2e-16 ***
#   SubstrateApple_Pectin                   4.1798     0.2705  15.451 2.70e-13 ***
#   SubstratePotato_Starch                  4.3890     0.2705  16.225 1.00e-13 ***
#   SubstrateXylan                          4.1903     0.2705  15.491 2.57e-13 ***
#   SubstrateGlucose                        3.9573     0.2705  14.629 8.12e-13 ***
#   SubstrateBasal                         -0.1427     0.2705  -0.527    0.603    
# PhaseContinous                          4.2208     0.2705  15.603 2.22e-13 ***
#   SubstrateApple_Pectin:PhaseContinous   -4.0439     0.3826 -10.571 4.36e-10 ***
#   SubstratePotato_Starch:PhaseContinous  -4.2553     0.3826 -11.123 1.68e-10 ***
#   SubstrateXylan:PhaseContinous          -4.2027     0.3826 -10.986 2.13e-10 ***
#   SubstrateGlucose:PhaseContinous        -3.9693     0.3826 -10.376 6.14e-10 ***
#   SubstrateBasal:PhaseContinous               NA         NA      NA       NA    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3313 on 22 degrees of freedom
# Multiple R-squared:  0.953,	Adjusted R-squared:  0.9316 
# F-statistic: 44.56 on 10 and 22 DF,  p-value: 2.862e-12
# 
# 
#Tukey HSD test was run but not used
fermenter_copies_posthoc<-aov(Fermenter_Copies_both)
fermenter_copies_posthoc
# Call:
#   aov(formula = Fermenter_Copies_multi_good)
# 
# Terms:
#   Substrate     Phase Substrate:Phase Residuals
# Sum of Squares  22.042741  6.439044       20.428213  2.414789
# Deg. of Freedom         5         1               4        22
# 
# Residual standard error: 0.3313052
# 1 out of 12 effects not estimable
# Estimated effects may be unbalanced

TukeyHSD(fermenter_copies_posthoc, conf.level=.95)

par(mar=c(5,6,4,1)+.1)
plot(TukeyHSD(fermenter_copies_posthoc, conf.level=.95), las = 2)






##########
#Relative abundance
##########


Fermenter_Percentage <- read.table(file = "C:/Users/r01ec25/Documents/Masters/Statistical Analysis/Fermenter/Species_Fermenter.txt", header = TRUE, sep = "\t",
                   stringsAsFactors = TRUE)

str(Fermenter_Percentage)
Fermenter_Percentage$Fibre <- factor(Fermenter_Percentage$Fibre, levels= c("Initial", "Mix", "Apple_Pectin", 'Potato_Starch', 'Xylan', 'Basal')) #reorder by substrate
str(Fermenter_Percentage$Fibre)

Fermenter_Percentage$Species <- factor(Fermenter_Percentage$Species, levels= c("DSM14838", "V975", "DSM19555", 'B5482', 'DSM18836', 'DSM1447')) #reoder by species
str(Fermenter_Percentage$Species)

str(Fermenter_Percentage)
# 'data.frame':	216 obs. of  7 variables:
#   $ DNA_Sample     : int  37 37 37 37 37 37 4 4 4 4 ...
# $ Fermenter      : Factor w/ 3 levels "F1","F2","F3": 1 1 1 1 1 1 2 2 2 2 ...
# $ Time           : int  0 0 0 0 0 0 0 0 0 0 ...
# $ Fibre          : Factor w/ 6 levels "Initial","Mix",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Incubation.Type: Factor w/ 4 levels "Batch","Continous",..: 4 4 4 4 4 4 4 4 4 4 ...
# $ Species        : Factor w/ 6 levels "DSM14838","V975",..: 1 2 3 4 5 6 1 2 3 4 ...
# $ percentage     : num  16.4 20.4 24 15.2 11.1 ...

boxplot(percentage~Species, data=Fermenter_Percentage, xlab="Species", ylab="Copies per ml")
boxplot(percentage~Fibre, data=Fermenter_Percentage, xlab="Substrate", ylab="Copies per ml")

#Confirmed that there is no statisitcal difference of the relative abundances between fermentors- so it means they can be used as replicates
Fermentertype<-lm(percentage~Fermenter, data=Fermenter_Percentage)
summary(Fermentertype)
# Call:
#   lm(formula = percentage ~ Fermenter, data = Fermenter_Percentage)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -16.667 -11.949  -4.131   7.224  51.660 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.667e+01  1.861e+00   8.956   <2e-16 ***
#   FermenterF2 -1.528e-10  2.632e+00   0.000        1    
# FermenterF3 -2.778e-11  2.632e+00   0.000        1    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 15.79 on 213 degrees of freedom
# Multiple R-squared:  1.796e-23,	Adjusted R-squared:  -0.00939 
# F-statistic: 1.912e-21 on 2 and 213 DF,  p-value: 1
anova(Fermentertype)
# Analysis of Variance Table
# 
# Response: percentage
# Df Sum Sq Mean Sq F value Pr(>F)
# Fermenter   2      0    0.00       0      1
# Residuals 213  53113  249.35   

#all possible factors
Fermenter_Both_MAX<-lm(percentage~Species+Fibre+Species*Fibre, data=Fermenter_Percentage)
drop1(Fermenter_Both_MAX, test='F')
# Single term deletions
# 
# Model:
#   percentage ~ Species + Fibre + Species * Fibre
# Df Sum of Sq     RSS    AIC F value    Pr(>F)    
# <none>                      4031.6 605.90                      
# Species:Fibre 25    7856.5 11888.1 737.57  10.289 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#interaction significant
Fermenter_Both_multi<-lm(percentage~Species*Fibre, data=Fermenter_Percentage)
drop1(Fermenter_Both_multi, test='F')
# Single term deletions
# 
# Model:
#   percentage ~ Species * Fibre
# Df Sum of Sq     RSS    AIC F value    Pr(>F)    
# <none>                      4031.6 605.90                      
# Species:Fibre 25    7856.5 11888.1 737.57  10.289 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#fibre didn't have an effect, so significance due to species
Fermenter_Both_add<-lm(percentage~Species+Fibre, data=Fermenter_Percentage)
drop1(Fermenter_Both_add, test='F')
# Single term deletions
# 
# Model:
#   percentage ~ Species + Fibre
# Df Sum of Sq   RSS    AIC F value Pr(>F)    
# <none>               11888 737.57                   
# Species  5     22958 34846 908.24   60.64 <2e-16 ***
#   Fibre    5         0 11888 727.57    0.00      1    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#time aka another way of saying fibre did not have an effect
Fermenter_Both_time<-lm(percentage~Species+Time, data=Fermenter_Percentage)
drop1(Fermenter_Both_time, test='F')
# Single term deletions
# 
# Model:
#   percentage ~ Species + Time
# Df Sum of Sq   RSS     AIC F value Pr(>F)    
# <none>               18592  976.33                   
# Species  5     34520 53113 1193.06   77.61 <2e-16 ***
#   Time     1         0 18592  974.33    0.00      1    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Only species had a significant effect

#IF FOLLOW THE DROP ON THE SPECIES MATTER
Fermenter_Species<-lm(percentage~Species, data=Fermenter_Percentage)
drop1(Fermenter_Species, test='F')
# Single term deletions
# 
# Model:
#   percentage ~ Species
# Df Sum of Sq   RSS     AIC F value    Pr(>F)    
# <none>               18592  974.33                      
# Species  5     34520 53113 1191.06  77.981 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(Fermenter_Species)
# 
# Call:
#   lm(formula = percentage ~ Species, data = Fermenter_Percentage)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -23.873  -5.418  -2.215   5.446  29.228 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       7.2423     1.5682   4.618 6.74e-06 ***
#   SpeciesV975      31.8560     2.2178  14.364  < 2e-16 ***
#   SpeciesDSM19555   4.2062     2.2178   1.897   0.0593 .  
# SpeciesB5482     21.3449     2.2178   9.624  < 2e-16 ***
#   SpeciesDSM18836  -1.7229     2.2178  -0.777   0.4381    
# SpeciesDSM1447    0.8623     2.2178   0.389   0.6978    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 9.409 on 210 degrees of freedom
# Multiple R-squared:  0.6499,	Adjusted R-squared:  0.6416 
# F-statistic: 77.98 on 5 and 210 DF,  p-value: < 2.2e-16

anova(Fermenter_Species) #partitions the total variance observed in a dataset into different components attributed to different sources of variation
# Analysis of Variance Table
# 
# Response: percentage
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Species     5  34520  6904.1  77.981 < 2.2e-16 ***
#   Residuals 210  18592    88.5                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Tukey Honest Significane Difference Test
anova_Fermenter_Species<-aov(percentage~Species, data=Fermenter_Percentage)
summary(anova_Fermenter_Species)
# Df Sum Sq Mean Sq F value Pr(>F)    
# Species       5  34520    6904   77.98 <2e-16 ***
#   Residuals   210  18592      89                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(anova_Fermenter_Species, conf.level=.95)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = percentage ~ Species, data = Fermenter_Percentage)
# 
# $Species
# diff        lwr        upr     p adj
# V975-DSM14838      31.8560015  25.477149  38.234854 0.0000000 
# DSM19555-DSM14838   4.2061897  -2.172662  10.585042 0.4071425 #none
# B5482-DSM14838     21.3449100  14.966058  27.723762 0.0000000
# DSM18836-DSM14838  -1.7228684  -8.101721   4.655984 0.9711639 #none
# DSM1447-DSM14838    0.8622595  -5.516593   7.241112 0.9988376 #none
# DSM19555-V975     -27.6498118 -34.028664 -21.270960 0.0000000
# B5482-V975        -10.5110915 -16.889944  -4.132239 0.0000576
# DSM18836-V975     -33.5788700 -39.957722 -27.200018 0.0000000
# DSM1447-V975      -30.9937420 -37.372594 -24.614890 0.0000000
# B5482-DSM19555     17.1387203  10.759868  23.517572 0.0000000
# DSM18836-DSM19555  -5.9290582 -12.307910   0.449794 0.0849930 #none
# DSM1447-DSM19555   -3.3439302  -9.722782   3.034922 0.6596713 #none
# DSM18836-B5482    -23.0677784 -29.446631 -16.688926 0.0000000
# DSM1447-B5482     -20.4826505 -26.861503 -14.103798 0.0000000
# DSM1447-DSM18836    2.5851279  -3.793724   8.963980 0.8526376 #none



plot(TukeyHSD(anova_Fermenter_Species, conf.level=.95), las = 2)


#Tukey HSD plot where if p<0.01 is shown in red
par(mar=c(5,6,4,1)+.1)
species_map <- c(
  "DSM18836" = "B. xyl",
  "V975"     = "B. ova",
  "B5482"    = "B. theta",
  "DSM19555" = "B. ster",
  "DSM1447"  = "P. vul",
  "DSM14838" = "B. cel"
)

tuk <- TukeyHSD(anova_Fermenter_Species, conf.level = 0.95)$Species

tuk_df <- as.data.frame(tuk)
tuk_df$Comparison <- rownames(tuk)

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
  xlab = "Differences in Mean Levels between Species in Fermenter",
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



