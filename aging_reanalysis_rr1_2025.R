#aging_reanalysis RR1 EJA
#created on 2025 March 17 

install.packages("Hmisc")
install.packages("broom.mixed")
install.packages("table1")

install.packages("performance")  # 如果未安装
library(performance)

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
library(merTools)

rm(list = ls())

#setwd("C:/Users/Donghui/SynologyDrive/Aging/CHLS/cleaned") #home desktop
#setwd("/Users/priscillawang/SynologyDrive/Aging/CHLS/cleaned") #mac
setwd ("C:/Users/donghuiwang/SynologyDrive/Aging/CHLS/cleaned") #um office desktop 

aging <-read_dta("aging_long_forR.dta")%>%
  arrange(ID, year)%>%
  mutate(cohort=ch+1,
         byrsq= byr*byr,
         ch2sq=ch2*ch2,
         ch2quar=ch2sq*ch2,
         age653=agesq65*age65)

table(aging$pattern2)

label(aging$agesq) <- "Age squared"
label(aging$useful)<- "Feel useful with age"


aging_baseline<-aging%>%
  dplyr::select(ID, TRUEAGE, highed,ch2,female, han,status )%>%
  unique()

table1(~std_miniscale+happy+useful+
        factor(married)+factor(alone)+factor(rurallive)+factor(ghlth), data = aging,decimal=2)

table1(~factor(ch2)+factor(highed)+factor(female)+factor(han)+factor(status), data=aging_baseline, decimal = 4)


#anova test 

dependent_vars <- c( "std_miniscale", "happy", "useful")
anova_results <- list()

for (var in dependent_vars) {
  formula <- as.formula(paste(var, "~ factor(ch2)"))
  anova_results[[var]] <- aov(formula, data = aging)
}


lapply(anova_results, summary)

#anova for age 
summary(aov(TRUEAGE ~ factor(ch2), data = aging_baseline))

aging_forcat<-aging%>%
  mutate(
    highed = highed+1,
    female = female +1,
    status = status +1,
     married = married +1,
    alone = alone +1,
    rurallive= rurallive +1,
    ghlth= ghlth+1
  )

#chi-squared test for time-varying 
cat_vars <- c( "married", "alone", "rurallive", "ghlth")

chi_sq_results <- list()
for (var in cat_vars) {
  tbl <- table(aging_forcat[[var]], aging_forcat$ch2)
  chi_sq_results[[var]] <- chisq.test(tbl)
}
  
#show test statistics and p value 
for (var in names(chi_sq_results)) {
  cat("Variable:", var, "\n")
  cat("Chi-squared statistic:", chi_sq_results[[var]]$statistic, "\n")
  cat("Degrees of freedom:", chi_sq_results[[var]]$parameter, "\n")
  cat("P-value:", chi_sq_results[[var]]$p.value, "\n\n")
}



#chi-squared for time-constant (aka baseline)
cat_vars <- c("highed", "female", "han", "status")

chi_sq_results <- list()
for (var in cat_vars) {
  tbl <- table(aging_baseline[[var]], aging_baseline$ch2)
  chi_sq_results[[var]] <- chisq.test(tbl)
}

#show test statistics and p value 
for (var in names(chi_sq_results)) {
  cat("Variable:", var, "\n")
  cat("Chi-squared statistic:", chi_sq_results[[var]]$statistic, "\n")
  cat("Degrees of freedom:", chi_sq_results[[var]]$parameter, "\n")
  cat("P-value:", chi_sq_results[[var]]$p.value, "\n\n")
}




#by cohort 

table1(~factor(highed)+factor(female)+factor(han)+factor(status) |factor(ch2), 
       data=aging_baseline, overall= F, test = TRUE, decimal = 4)

table1(~std_miniscale+happy+useful+
         factor(married)+factor(alone)+factor(rurallive)+factor(ghlth)|ch2,
       overall= F, data = aging,decimal=2)


#graphing of the descriptive results 
  aging_ch2<-aging%>%
    mutate(age=age65+65)%>%
    group_by(ch2, age)%>%
    summarise(meanscale=weighted.mean(miniscale, weight))%>%
    filter(meanscale<8, meanscale>5)
  
  
ggplot(aging_ch2, aes(x=age, y=meanscale))+
  geom_line(aes(color=factor(ch2)))

#-------------------------------
#main results 
#cohort as continuous ; SES is mesured as education 
#null model 
m0_std_miniscale<-lmer(std_miniscale ~ 
                         (age65 || ID), data = aging, REML = FALSE, 
                       control = lmerControl(optimizer ="bobyqa"))

summary(m0_std_miniscale)

tab_model(m0_std_miniscale,
          collapse.se = TRUE, show.ci = FALSE, show.r2 = TRUE, show.re.var = TRUE,
          p.style =  "stars",  digits = 3,  show.icc = TRUE)


#null model cohort2 
  m1_std_miniscale<-lmer(std_miniscale ~ 1 + age65*ch2 +ch2sq+ agesq65+
                       (age65 || ID), data = aging, REML = FALSE, 
                     control = lmerControl(optimizer ="bobyqa"))
  
  summary(m1_std_miniscale)

#model 2: time-invariant controls ->cohort effect would reduce
  m2_std_miniscale<-lmer(std_miniscale ~ 1 + age65*ch2 +ch2sq+ agesq65+
                           #controls
                           +female+han+factor(highed)+
                           #mortality pattern
                           +factor(status)+factor(surveyear)+
                           (age65 || ID), data = aging, REML = FALSE, 
                         control = lmerControl(optimizer ="bobyqa"))
  
  summary(m2_std_miniscale)
  
# model 3 : time-varying controls -> mainly affect age effect 
  m3_std_miniscale<-lmer(std_miniscale ~ 1 + age65*ch2 +ch2sq+ agesq65+
                           #controls
                           +female+han+factor(highed)+
                           #time-varying controls
                           +married+alone+rurallive+factor(ghlth)+
                           #mortality pattern
                           +factor(status)+factor(surveyear)+
                           (age65 || ID), data = aging, REML = FALSE, 
                         control = lmerControl(optimizer ="bobyqa"))
  
  summary(m3_std_miniscale)
  
  
 
#model 4: ses interactions 
m4_std_miniscale<-lmer(std_miniscale ~ 1 + age65*ch2*factor(highed) +ch2sq+ agesq65+
                         +female+han+factor(highed)+
                         #time-varying controls
                         +married+alone+rurallive+factor(ghlth)+
                         #mortality pattern
                         +factor(status)+factor(surveyear)+
                         (age65 || ID), data = aging, REML = FALSE, 
                       control = lmerControl(optimizer ="bobyqa"))

summary(m4_std_miniscale)



# m5_std_miniscale<-lmer(std_miniscale ~ 1 + age65*factor(ch2)*factor(highed) + agesq65+
#                          +female+han+factor(highed)+
#                          #time-varying controls
#                          +married+alone+rurallive+factor(ghlth)+
#                          #mortality pattern
#                          +factor(status)+factor(surveyear)+
#                          (age65 || ID), data = aging, REML = FALSE, 
#                        control = lmerControl(optimizer ="bobyqa"))
# 
# summary(m5_std_miniscale)

v_fixed <- var(predict(m1_std_miniscale, re.form = NA))  # Fixed effects prediction
vc <- as.data.frame(VarCorr(m1_std_miniscale))
var_random <- sum(vc$vcov)  # Sum of all random effect variances


v_resid <- attr(VarCorr(m1_std_miniscale), "sc")^2  # Residual variance

 <- v_fixed / (v_fixed + var_random + v_resid)



BIC(m1_std_miniscale, m2_std_miniscale, m3_std_miniscale, m4_std_miniscale )
icc(m1_std_miniscale )

tab_model(m1_std_miniscale, m2_std_miniscale, m3_std_miniscale, m4_std_miniscale ,
          collapse.se = TRUE, show.ci = FALSE, show.r2 = TRUE, show.re.var = TRUE,
          p.style =  "stars",  digits = 3,  show.icc = TRUE)



tab_model(m1_std_miniscale, m2_std_miniscale, m3_std_miniscale, m4_std_miniscale ,
          collapse.se = TRUE, show.ci = FALSE, show.r2 = TRUE, show.re.var = TRUE,
          p.style =  "stars",  digits = 3,  show.icc = TRUE, 
          file = "main.rtf")

#==============================================================

#supplemental analysis 
#A. separate analysis of happy and useful 

#null model cohort2 
m1_happy<-lmer(happy ~ 1 + age65*ch2 +ch2sq+ agesq65+
                         (age65 || ID), data = aging, REML = FALSE, 
                       control = lmerControl(optimizer ="bobyqa"))

summary(m1_happy)

#model 2: time-invariant controls ->cohort effect would reduce
m2_happy<-lmer(happy ~ 1 + age65*ch2 +ch2sq+ agesq65+
                         #controls
                         +female+han+factor(highed)+
                         #mortality pattern
                         +factor(status)+factor(surveyear)+
                         (age65 || ID), data = aging, REML = FALSE, 
                       control = lmerControl(optimizer ="bobyqa"))

summary(m2_happy)

# model 3 : time-varying controls -> mainly affect age effect 
m3_happy<-lmer(happy ~ 1 + age65*ch2 +ch2sq+ agesq65+
                         #controls
                         +female+han+factor(highed)+
                         #time-varying controls
                         +married+alone+rurallive+factor(ghlth)+
                         #mortality pattern
                         +factor(status)+factor(surveyear)+
                         (age65 || ID), data = aging, REML = FALSE, 
                       control = lmerControl(optimizer ="bobyqa"))

summary(m3_happy)


#model 4: ses interactions 
m4_happy<-lmer(happy ~ 1 + age65*ch2*factor(highed) +ch2sq+ agesq65+
                         +female+han+factor(highed)+
                         #time-varying controls
                         +married+alone+rurallive+factor(ghlth)+
                         #mortality pattern
                         +factor(status)+factor(surveyear)+
                         (age65 || ID), data = aging, REML = FALSE, 
                       control = lmerControl(optimizer ="bobyqa"))

summary(m4_happy)

BIC(m1_happy, m2_happy, m3_happy, m4_happy )
tab_model(m1_happy, m2_happy, m3_happy, m4_happy ,
          collapse.se = TRUE, show.ci = FALSE, show.r2 = TRUE, show.re.var = TRUE,
          p.style =  "stars",  digits = 3,  show.icc = TRUE, 
          file= "happy.rtf")



#--------------------------
#B usefull 

#null model cohort2 
m1_useful<-lmer(useful ~ 1 + age65*ch2 +ch2sq+ agesq65+
                 (age65 || ID), data = aging, REML = FALSE, 
               control = lmerControl(optimizer ="bobyqa"))

summary(m1_useful)

#model 2: time-invariant controls ->cohort effect would reduce
m2_useful<-lmer(useful ~ 1 + age65*ch2 +ch2sq+ agesq65+
                 #controls
                 +female+han+factor(highed)+
                 #mortality pattern
                 +factor(status)+factor(surveyear)+
                 (age65 || ID), data = aging, REML = FALSE, 
               control = lmerControl(optimizer ="bobyqa"))

summary(m2_useful)

# model 3 : time-varying controls -> mainly affect age effect 
m3_useful<-lmer(useful ~ 1 + age65*ch2 +ch2sq+ agesq65+
                 #controls
                 +female+han+factor(highed)+
                 #time-varying controls
                 +married+alone+rurallive+
                 #mortality pattern
                 +factor(status)+factor(surveyear)+
                 (age65 || ID), data = aging, REML = FALSE, 
               control = lmerControl(optimizer ="bobyqa"))

summary(m3_useful)


#model 4: ses interactions 
m4_useful<-lmer(useful ~ 1 + age65*ch2*factor(highed) +ch2sq+ agesq65+
                 +female+han+factor(highed)+
                 #time-varying controls
                 +married+alone+rurallive+
                 #mortality pattern
                 +factor(status)+factor(surveyear)+
                 (age65 || ID), data = aging, REML = FALSE, 
               control = lmerControl(optimizer ="bobyqa"))

summary(m4_useful)

BIC(m1_useful, m2_useful, m3_useful, m4_useful )
tab_model(m1_useful, m2_useful, m3_useful, m4_useful ,
          collapse.se = TRUE, show.ci = FALSE, show.r2 = TRUE, show.re.var = TRUE,
          p.style =  "stars",  digits = 3,
          file= "useful.rtf")

#--------------------------
#graph results 
age_seq <- seq(0, 49, by = 1)


pred_df <- expand.grid(age65 = age_seq, ch2 = 1:10, highed = 0:1)%>%
  as_tibble()%>%
  mutate(ID = ch2+1 ,
         ch2sq=ch2*ch2,
         agesq65=age65*age65,
         age= age65+65,
         
         age65 = case_when(
           (age65 >= 28  & age65 <= 40 & ch2 == 1) ~ age65,
           (age65 >= 23 & age65 <= 35 & ch2 == 2) ~ age65,
           (age65 >= 18 & age65 <= 30 & ch2 == 3) ~ age65,
           (age65 >= 13 & age65 <= 25 & ch2 == 4) ~ age65,
           (age65 >= 12 & age65 <= 20 & ch2 == 5) ~ age65,
           (age65 >= 7 & age65 <= 15 & ch2 == 6) ~ age65,
           (age65 >= 2 & age65 <= 10 & ch2 == 7) ~ age65,
           (age65 >= 0 & age65 <= 10 & ch2 == 8) ~ age65,
           (age65 >= 0 & age65 <= 8 & ch2 == 9) ~ age65,
           (age65 >= 0 & age65 <= 9 & ch2 == 10) ~ age65,
           
           TRUE ~ NA_real_
         ),
         
         #others fixed at the largest share 
         female=1,
         han=1,
         married=1,
         alone=0,
         status = 0,
         surveyear= 1998,
         ghlth = 4,
         rurallive = 1
  )

#predictions 
pred_df <- pred_df %>% filter(!is.na(age65))

#pred_df$null <- predict(m1_std_miniscale, newdata = pred_df, 
#                        allow.new.levels = TRUE, re.form= NA )

pred_df$full <- predict(m4_std_miniscale, newdata = pred_df, 
                        allow.new.levels = TRUE, re.form= NA)

pred_df$happy <- predict(m4_happy, newdata = pred_df, 
                        allow.new.levels = TRUE, re.form= NA)

pred_df$useful <- predict(m4_useful, newdata = pred_df, 
                         allow.new.levels = TRUE, re.form= NA)


#Version 2: predict interval 
pred_df$full<-predictInterval(m4_std_miniscale, newdata = pred_df, n.sims = 999, level= 0.95 )






#graph results 

pred_df$highed<-factor(pred_df$highed, labels = c("Low education", "High education"))
pred_df$ch2<-factor(pred_df$ch2, labels = c("1898-1904", 
                                            "1905-1909" , 
                                            "1910-1914" , 
                                            "1915-1919" ,
                                            "1920-1924" ,
                                            "1925-1929" ,
                                            "1930-1934" ,
                                            "1935-1939" ,
                                            "1940-1944" ,
                                            "1945-1949" ))

#cohort effect: high as example


plot_std_miniscale <- ggplot(filter(pred_df, highed == "High education"), aes(x = age65, y = full)) +
  geom_line(linewidth = 1, aes(color = factor(ch2))) +
  geom_text(
    data = filter(pred_df, highed == "High education") %>% group_by(ch2) %>% filter(age65 == min(age65)),  # Label at the end of each line
    aes(x = age65, y = full, label = ch2, color = factor(ch2)),
    hjust = -0.2,  # Adjust horizontal position of the label
    vjust = 0.5,   # Adjust vertical position of the label
    size = 4.5       # Adjust text size
  ) +
  scale_x_continuous(breaks = seq(0, 49, by = 5),  # Adjust breaks every 5 integers
                     labels = seq(65, 114, by = 5),  # Adjust labels accordingly
                     expand = c(0.01, 0)) +  # Adjust expansion for better fit
  scale_y_continuous(breaks = seq(0.2, 0.8, by = 0.2), limits = c(0.2, 0.8)) +
  labs(color = "Birth cohort") + 
  xlab("Age") +
  ylab("Self-perception of aging") +
  theme_minimal()+
  theme(
    legend.position = "none",  # Remove the legend
    text = element_text(size = 13))

print(plot_std_miniscale)
ggsave("miniscale_high.png")


#SES difference 
ggplot(pred_df, aes(x=age, y=full))+
  geom_line(linewidth = 1, aes(color=factor(ch2), linetype=factor(highed)))+
  theme_minimal()




plot_ses<-ggplot(filter(pred_df), aes(x=age, y=full))+
  geom_line(linewidth = 1, aes(linetype=factor(highed) , color = factor(highed)))+
  scale_linetype_manual(values = c("dashed",  "solid"))+
  scale_color_manual(values = c("black",  "red"))+
  facet_wrap(~ch2,scales = "free_x")+
  theme_minimal()+
  theme(legend.position = c(0.95, 0.05),  # Move legend to lower-right corner
        legend.justification = c(1, 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        text = element_text(size = 14))+
  scale_x_continuous(breaks = seq(floor(min(pred_df$age)), ceiling(max(pred_df$age)), by = 2)) +
  xlab("Age") +
  ylab("Self-perception of aging") 


print(plot_ses)

ggsave("plot_ses.png")



#try age effect only 
pred_df_age<-pred_df%>%
  group_by(age, highed)%>%
  summarise(full_age= mean(full))

ggplot(pred_df_age, aes(x=age, y=full_age))+
  geom_point(aes(color=factor(highed)))+
  stat_smooth(aes(linetype = factor(highed),
                  color = factor(highed)),
              se = FALSE)+
  theme_bw()


plot_happy_ses<-ggplot(filter(pred_df), aes(x=age, y=happy))+
  geom_line(linewidth = 1, aes(linetype=factor(highed) , color = factor(highed)))+
  scale_linetype_manual(values = c("dashed",  "solid"))+
  scale_color_manual(values = c("black",  "red"))+
  facet_wrap(~ch2,scales = "free_x")+
  theme_minimal()+
  theme(legend.position = c(0.95, 0.05),  # Move legend to lower-right corner
        legend.justification = c(1, 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        text = element_text(size = 14))+
  scale_x_continuous(breaks = seq(floor(min(pred_df$age)), ceiling(max(pred_df$age)), by = 2)) +
  xlab("Age") +
  ylab("Feel happy with age") 


print(plot_happy_ses)
ggsave("plot_happy_ses.png")



plot_useful_ses<-ggplot(filter(pred_df), aes(x=age, y=useful))+
  geom_line(linewidth = 1, aes(linetype=factor(highed) , color = factor(highed)))+
  scale_linetype_manual(values = c("dashed",  "solid"))+
  scale_color_manual(values = c("black",  "red"))+
  facet_wrap(~ch2,scales = "free_x")+
  theme_minimal()+
  theme(legend.position = c(0.95, 0.05),  # Move legend to lower-right corner
        legend.justification = c(1, 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        text = element_text(size = 14))+
  scale_x_continuous(breaks = seq(floor(min(pred_df$age)), ceiling(max(pred_df$age)), by = 2)) +
  xlab("Age") +
  ylab("Feel useful with age") 


print(plot_useful_ses)
ggsave("plot_useful_ses.png")



#----------
#happy

plot_happy <- ggplot(filter(pred_df, highed == "High education"), aes(x = age65, y = happy)) +
  geom_line(linewidth = 1, aes(color = factor(ch2))) +
  geom_text(
    data = filter(pred_df, highed == "High education") %>% group_by(ch2) %>% filter(age65 == min(age65)),  # Label at the end of each line
    aes(x = age65, y = happy, label = ch2, color = factor(ch2)),
    hjust = -0.2,  # Adjust horizontal position of the label
    vjust = 0.5,   # Adjust vertical position of the label
    size = 4.5       # Adjust text size
  ) +
  scale_x_continuous(breaks = seq(0, 49, by = 5),  # Adjust breaks every 5 integers
                     labels = seq(65, 114, by = 5),  # Adjust labels accordingly
                     expand = c(0.01, 0)) +  # Adjust expansion for better fit
  labs(color = "Birth cohort") + 
  xlab("Age") +
  ylab("Feeling happy with age") +
  theme_minimal()+
  theme(
    legend.position = "none",  # Remove the legend
    text = element_text(size = 13))

print(plot_happy)
ggsave("happy.png")




# useful 
plot_useful <- ggplot(filter(pred_df, highed == "High education"), aes(x = age65, y = useful)) +
  geom_line(linewidth = 1, aes(color = factor(ch2))) +
  geom_text(
    data = filter(pred_df, highed == "High education") %>% group_by(ch2) %>% filter(age65 == min(age65)),  # Label at the end of each line
    aes(x = age65, y = useful, label = ch2, color = factor(ch2)),
    hjust = -0.2,  # Adjust horizontal position of the label
    vjust = 0.5,   # Adjust vertical position of the label
    size = 4.5       # Adjust text size
  ) +
  scale_x_continuous(breaks = seq(0, 49, by = 5),  # Adjust breaks every 5 integers
                     labels = seq(65, 114, by = 5),  # Adjust labels accordingly
                     expand = c(0.01, 0)) +  # Adjust expansion for better fit
  labs(color = "Birth cohort") + 
  xlab("Age") +
  ylab("Feel useful with age") +
  theme_minimal()+
  theme(
    legend.position = "none",  # Remove the legend
    text = element_text(size = 13))

print(plot_useful)
ggsave("useful.png")



#--------------------------
# sensitivity analysis 
#cohort as birth year  

std_miniscale_byr<-lmer(std_miniscale ~ 1 + age65*byr*factor(highed) +agesq65 + byrsq+
                         +female+han+factor(highed)+
                         #time-varying controls
                         +married+alone+rurallive+
                         #mortality pattern
                         +factor(status)+
                         (age65 || ID), data = aging, REML = FALSE, 
                       control = lmerControl(optimizer ="bobyqa"))

summary(std_miniscale_byr)



happy_byr<-lmer(happy ~ 1 + age65*byr*factor(highed) +agesq65 +byrsq+
                             +female+han+factor(highed)+
                             #time-varying controls
                             +married+alone+rurallive+
                             #mortality pattern
                             +factor(status)+
                             (age65 || ID), data = aging, REML = FALSE, 
                           control = lmerControl(optimizer ="bobyqa"))

summary(happy_byr)


useful_byr<-lmer(useful ~ 1 + age65*byr*factor(highed) +agesq65+ byrsq+
                  +female+han+factor(highed)+
                  #time-varying controls
                  +married+alone+rurallive+
                  #mortality pattern
                  +factor(status)+
                  (age65 || ID), data = aging, REML = FALSE, 
                control = lmerControl(optimizer ="bobyqa"))

summary(useful_byr)



tab_model(std_miniscale_byr, happy_byr, useful_byr, 
          collapse.se = TRUE, show.ci = FALSE, show.r2 = TRUE, show.re.var = TRUE,
          p.style =  "stars",  digits = 3,
          file= "byr.rtf")

BIC(std_miniscale_byr, happy_byr, useful_byr )

#--------------------------
#C. occupation instead of education 

std_miniscale_occ<-lmer(std_miniscale ~ 1 + age65*ch2*factor(prof) +ch2sq+ agesq65+
                  +female+han+factor(prof)+
                  #time-varying controls
                  +married+alone+rurallive+
                  #mortality pattern
                  +factor(status)++factor(surveyear)+
                  (age65 || ID)+(1 | pattern2), data = aging, REML = FALSE, 
                control = lmerControl(optimizer ="bobyqa"))

summary(std_miniscale_occ)




happy_occ<-lmer(happy ~ 1 + age65*ch2*factor(prof) +ch2sq+ agesq65+
                   +female+han+factor(prof)+
                   #time-varying controls
                   +married+alone+rurallive+
                   #mortality pattern
                   +factor(status)++factor(surveyear)+
                   (age65 || ID)+(1 | pattern2), data = aging, REML = FALSE, 
                 control = lmerControl(optimizer ="bobyqa"))

summary(happy_occ)


useful_occ<-lmer(useful ~ 1 + age65*ch2*factor(prof) +ch2sq+ agesq65+
                  +female+han+factor(prof)+
                  #time-varying controls
                  +married+alone+rurallive+
                  #mortality pattern
                  +factor(status)+factor(surveyear)+
                  (age65 || ID)+(1 | pattern2), data = aging, REML = FALSE, 
                control = lmerControl(optimizer ="bobyqa"))

summary(useful_occ)


tab_model(std_miniscale_occ, happy_occ, useful_occ, 
          collapse.se = TRUE, show.ci = FALSE, show.r2 = TRUE, show.re.var = TRUE,
          p.style =  "stars",  digits = 3,
          file= "occ.rtf")
BIC(std_miniscale_occ, happy_occ, useful_occ)

#--------------------------
#D. level2 interactions with other demo char: 

std_miniscale_inter<-lmer(std_miniscale ~ 1 + age65*ch2*factor(highed) +ch2sq+ agesq65+
                            age65*ch2*factor(female)+ 
                            age65*ch2*factor(han)+ 
                          +female+han+factor(highed)+
                          #time-varying controls
                          +married+alone+rurallive+
                          #mortality pattern
                          +factor(status)+factor(surveyear)+
                          (age65 || ID), data = aging, REML = FALSE, 
                        control = lmerControl(optimizer ="bobyqa"))

summary(std_miniscale_inter)



#E. quadratic term 

std_miniscale_quar<-lmer(std_miniscale ~ 1 + age65*ch2*factor(highed) +ch2sq+ agesq65*ch2*factor(highed)+
                         +female+han+factor(highed)+
                         #time-varying controls
                         +married+alone+rurallive+factor(ghlth)+
                         #mortality pattern
                         +factor(status)+factor(surveyear)+
                         (age65 || ID), data = aging, REML = FALSE, 
                       control = lmerControl(optimizer ="bobyqa"))

summary(std_miniscale_quar)

tab_model(std_miniscale_inter,std_miniscale_quar ,
          collapse.se = TRUE, show.ci = FALSE, show.r2 = TRUE, show.re.var = TRUE,
          p.style =  "stars",  digits = 3,
          file= "alter_model.rtf")




#--------------------------
#E. pattern mixture 
patternmixture<-lmer(std_miniscale ~ 1 + age65*ch2*factor(highed) +ch2sq+ agesq65+
                         +female+han+factor(highed)+
                         #time-varying controls
                         +married+alone+rurallive+factor(ghlth)+factor(pattern1)*age65+
                         #mortality pattern
                         (1 | pattern1) + (age65 || ID), data = aging, REML = FALSE, 
                       control = lmerControl(optimizer ="bobyqa"))

summary(patternmixture)


tab_model(patternmixture,
          collapse.se = TRUE, show.ci = FALSE, show.r2 = TRUE, show.re.var = TRUE,
          p.style =  "stars",  digits = 3,
          file= "patternmixture.rtf")

#----------------------------------
# F: focusing on later born only : only the lastest 5 cohorts 
age_small<-aging%>%
  filter(ch2>4 )

std_miniscale_small<-lmer(std_miniscale ~ 1 + age65*ch2*factor(highed) +ch2sq+ agesq65+
                           +female+han+factor(highed)+
                           #time-varying controls
                           +married+alone+rurallive+factor(ghlth)+factor(surveyear)+
                           #mortality pattern
                           +factor(status)+
                           (age65 || ID), data = age_small, REML = FALSE, 
                         control = lmerControl(optimizer ="bobyqa"))

summary(std_miniscale_small)

tab_model(std_miniscale_small,
          collapse.se = TRUE, show.ci = FALSE, show.r2 = TRUE, show.re.var = TRUE,
          p.style =  "stars",  digits = 3,
          file= "std_miniscale_small.rtf")



#-----------------------------------------------
# separate analysis 
age_1998<-aging%>%
  filter(surveyear== 1998 )%>%
  filter(TRUEAGE>=85)


#model 4: ses interactions 
std_miniscale_1998<-lmer(std_miniscale ~ 1 + age65*ch2*factor(highed) +ch2sq+ agesq65+
                         +female+han+factor(highed)+
                         #time-varying controls
                         +married+alone+rurallive+
                         #mortality pattern
                         +factor(status)+
                         (age65 || ID), data = age_1998, REML = FALSE, 
                       control = lmerControl(optimizer ="bobyqa"))

summary(std_miniscale_1998)


age_2000<-aging%>%
  filter(surveyear== 2000 )

std_miniscale_2000<-lmer(std_miniscale ~ 1 + age65*ch2*factor(highed) +ch2sq+ agesq65+
                           +female+han+factor(highed)+
                           #time-varying controls
                           +married+alone+rurallive+factor(ghlth)+
                           #mortality pattern
                           factor(status)+
                           (age65 || ID), data = age_2000, REML = FALSE, 
                         control = lmerControl(optimizer ="bobyqa"))

summary(std_miniscale_2000)


age_2002<-aging%>%
  filter(surveyear== 2002 )
std_miniscale_2002<-lmer(std_miniscale ~ 1 + age65*ch2*factor(highed) +ch2sq+ agesq65+
                           +female+han+factor(highed)+
                           #time-varying controls
                           +married+alone+rurallive+factor(ghlth)+
                           #mortality pattern
                           factor(status)+
                           (age65 || ID), data = age_2002, REML = FALSE, 
                         control = lmerControl(optimizer ="bobyqa"))

summary(std_miniscale_2002)


