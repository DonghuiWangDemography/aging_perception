#multiple imputation 

install.packages("mice")

library("mice")
library(tidyverse)
library(ggplot2)
library(lme4)
library(haven)
library(foreign)
library(lmerTest)
library(sjPlot)
library(broom.mixed) #for tidy, glance, and augment functions for lme4 models
library(dplyr)
library(MASS)
library(table1)
library(Hmisc)  #label 

rm(list = ls())

#setwd("C:/Users/Donghui/SynologyDrive/Aging/CHLS/cleaned") #home desktop
#setwd("/Users/priscillawang/SynologyDrive/Aging/CHLS/cleaned") #mac
setwd ("C:/Users/donghuiwang/SynologyDrive/Aging/CHLS/cleaned") #um office desktop 

aging <-read_dta("aging_long_forR_mi.dta")%>%
  arrange(ID, year)%>%
  mutate(cohort=ch+1,
         byrsq= byr*byr,
         ch2sq=ch2*ch2,
         ch2quar=ch2sq*ch2,
         age653=agesq65*age65)

aging_df <- as.data.frame(aging)

vars_to_impute <- c( "highed", 
                    "female", "han", "married", "alone", "rurallive", "ghlth", 
                      "ID")

pred_matrix <- make.predictorMatrix(aging_df[, vars_to_impute])
pred_matrix[,"ID"] <- 0 

imputed_data <- mice(aging, m = 5, method = 'rf', maxit = 50, 
                     predictorMatrix = pred_matrix, seed = 500)




fits <- with(imputed_data, lmer(std_miniscale ~ 1 + age65*ch2*factor(highed) + ch2sq + agesq65 +
                         female + han + factor(highed) +
                         married + alone + rurallive + factor(ghlth) +
                         factor(status) + factor(surveyear) +
                         (age65 || ID),
                       REML = FALSE, 
                       control = lmerControl(optimizer = "bobyqa")))

pooled_results <- pool(fits)


pooled_estimates <- tidy(pooled_results, conf.int = TRUE)

# Now use tab_model on the tidied data
sjPlot::tab_model(
  pooled_estimates,
  show.ci = FALSE, 
  show.r2 = TRUE,
  p.style = "stars",
  digits = 3
)


BIC(fits)


