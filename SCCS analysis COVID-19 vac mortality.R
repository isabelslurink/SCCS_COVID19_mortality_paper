# COVID-19 vaccination and short-term mortality risk: a nationwide SCCS in the Netherlands
# Author: IAL Slurink

set.seed(123)

##### Load packages ##### 
library(SCCS)
library(data.table)
library(eeptools)
library(dplyr)
library(gnm)
library(lubridate)
library(splines)

##### Preparations analyses ######
# Load in functions
source("H:/Annemarijn/oversterfte_covidvac/SCCS/Model_Functions.R")

#Load dataset
sccs_cohort <- fread(file = "H:/Datafiles/Oversterfte_covidvac/sccs_cohort_20241105.csv", header = TRUE, fill = TRUE)

##### Primary series vaccination ##### 

###### Set parameters - same for all analyses ###### 

## Parameters for stratification of days of the start of exposure-related risk

expogrp_1 = c(0,1) # Ignore first day after vaccination (no mortality to be expected related to vaccination on day 1)
expogrp_3 = c(0,1,7,14) # Ignore first day after vaccination, risk period is 3 separate weeks
expogrp_12 = c(0,1,7,14,21,28,35,42,49,56,63,70,77) # Ignore first day after vaccination, risk period is 12 separate weeks

## Parameters for 14 day interval calendar adjustment for seasonality

agegrp = seq(14,316,14) # Divide follow-up time into 14 day intervals starting from beginning to end of exposure period (called agegrp in original eventdepenexp function, but has nothing to do with age here)

## Parameters spline adjustment for seasonality

### Parameter agegrp is now the calender day of the start of the week included as continuous variable in the analysis

start_date <- as.Date("2021-01-06")
end_date <- as.Date("2021-11-18")

week_start <- start_date - wday(start_date, week_start = 1) + 1
week_start_days <- seq.Date(week_start, end_date, by = "week")

agegrp_spline = as.numeric(week_start_days - start_date)[-1]

####### All ages #####  
df <- sccs_cohort 
results_list1 <- run_models_primary(df, outcome = "alldeaths")
df <- sccs_cohort %>% filter(covid == 0)
results_list2 <- run_models_primary(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "allages"))
results <- export_outcome_primary(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_allages_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_allages_list_20241105")
rm(results, results_list, results_list1, results_list2)
gc()

####### Age group 12 and 29 #####
df <- sccs_cohort %>% filter(age_start >=12 & age_start <= 29)
results_list1 <- run_models_primary(df, outcome = "alldeaths")
df <- sccs_cohort %>% filter(age_start >=12 & age_start <= 29) %>% filter(covid == 0)
results_list2 <- run_models_primary(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "12_29"))
results <- export_outcome_primary(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_1229_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_1229_list_20241105")
rm(results, results_list, results_list1, results_list2)
gc()

####### Age group 30 and 39 #####  
df <- sccs_cohort %>% filter(age_start >=30 & age_start <= 39)
results_list1 <- run_models_primary(df, outcome = "alldeaths")
df <- sccs_cohort %>% filter(age_start >=30 & age_start <= 39) %>% filter(covid == 0)
results_list2 <- run_models_primary(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "30_39"))
results <- export_outcome_primary(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_3039_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_3039_list_20241105")
rm(results, results_list, results_list1, results_list2)
gc()

####### Age group 40 and 49 #####  
df <- sccs_cohort %>% filter(age_start >=40 & age_start <= 49)
results_list1 <- run_models_primary(df, outcome = "alldeaths")
df <- sccs_cohort %>% filter(age_start >=40 & age_start <= 49) %>% filter(covid == 0)
results_list2 <- run_models_primary(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "40_49"))
results <- export_outcome_primary(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_4049_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_4049_list_20241105")
rm(results, results_list, results_list1, results_list2)
gc()

####### Age group 50 and 59 ##### 
df <- sccs_cohort %>% filter(age_start >=50 & age_start <= 59)
results_list1 <- run_models_primary(df, outcome = "alldeaths")
df <- sccs_cohort %>% filter(age_start >=50 & age_start <= 59) %>% filter(covid == 0)
results_list2 <- run_models_primary(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "50_59"))
results <- export_outcome_primary(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_5059_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_5059_list_20241105")
rm(results, results_list, results_list1, results_list2)
gc()

####### Age group 60 and 69 #####  
df <- sccs_cohort %>% filter(age_start >=60 & age_start <= 69)
results_list1 <- run_models_primary(df, outcome = "alldeaths")
df <- sccs_cohort %>% filter(age_start >=60 & age_start <= 69) %>% filter(covid == 0)
results_list2 <- run_models_primary(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "6069"))
results <- export_outcome_primary(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_6069_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_6069_list_20241105")
rm(results, results_list, results_list1, results_list2)
gc()

####### Age group 70 and 79 #####  
df <- sccs_cohort %>% filter(age_start >=70 & age_start <= 79)
results_list1 <- run_models_primary(df, outcome = "alldeaths")
df <- sccs_cohort %>% filter(age_start >=70 & age_start <= 79) %>% filter(covid == 0)
results_list2 <- run_models_primary(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "70_79"))
results <- export_outcome_primary(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_7079_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_7079_list_20241105")
rm(results, results_list, results_list1, results_list2)
gc()

####### Age group 80 and 89 #####  
df <- sccs_cohort %>% filter(age_start >=80 & age_start <= 89)
results_list1 <- run_models_primary(df, outcome = "alldeaths")
df <- sccs_cohort %>% filter(age_start >=80 & age_start <= 89) %>% filter(covid == 0)
results_list2 <- run_models_primary(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "80_89"))
results <- export_outcome_primary(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_8089_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_8089_list_20241105")
rm(results, results_list, results_list1, results_list2)
gc()

####### Age group >90 #####  
df <- sccs_cohort %>% filter(age_start >=90)
results_list1 <- run_models_primary(df, outcome = "alldeaths")
df <- sccs_cohort %>% filter(age_start >=90) %>% filter(covid == 0)
results_list2 <- run_models_primary(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "90"))
results <- export_outcome_primary(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_90_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_90_list_20241105")
rm(results, results_list, results_list1, results_list2)
gc()

#### SARS-CoV-2 infection ####

#Load dataset
sccs_cohort_inf <- fread(file = "H:/Datafiles/Oversterfte_covidvac/sccs_cohort_inf_20240711.csv", header = TRUE, fill = TRUE)

## Rename the infection SCCS variables 

sccs_cohort_inf_sel <- sccs_cohort_inf %>% dplyr::select(!c(start, end, event)) %>%
  dplyr::rename(case = case_inf, start = start_inf, end = end_inf, event = event_inf)

###### Set parameters - same for all analyses ###### 

## Parameters for stratification of days of the start of exposure-related risk

expogrp_1 = 0 # First day of infection is not ignored
expogrp_3 = c(0,7,14) # Risk period is 3 separate weeks
expogrp_12 = c(0,7,14,21,28,35,42,49,56,63,70,77) # Risk period is 12 separate weeks

## Parameters for 14 day interval calendar adjustment for seasonality

agegrp = seq(14,578,14) # Divide age groups into 14 day intervals starting from beginning to end of exposure period, accounts for seasonality

## Parameters spline adjustment for seasonality

### Parameter agegrp is now the calender day of the start of the week included as continuous variable in the analysis

start_date <- as.Date("2020-06-01")
end_date <- as.Date("2021-12-31")

week_start <- start_date - wday(start_date, week_start = 1) + 1
week_start_days <- seq.Date(week_start, end_date, by = "week")

agegrp_spline = as.numeric(week_start_days - start_date)[-1]

######   All ages ######
df <- sccs_cohort_inf_sel
results_list1 <- run_models_inf(df, outcome = "alldeaths")
df <- sccs_cohort_inf_sel %>% filter(covid == 0)
results_list2 <- run_models_inf(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "allages"))
results <- export_outcome_inf(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_allages_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_allages_list_20241105")
gc()

####### Age group 12 and 29 #####
df <- sccs_cohort_inf_sel %>% filter(age_start >=12 & age_start <= 29)
results_list1 <- run_models_inf(df, outcome = "alldeaths")
df <- sccs_cohort_inf_sel %>% filter(age_start >=12 & age_start <= 29) %>% filter(covid == 0)
results_list2 <- run_models_inf(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "12_29"))
results <- export_outcome_inf(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_1229_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_1229_list_20241105")
gc()

####### Age group 30 and 39 #####
df <- sccs_cohort_inf_sel %>% filter(age_start >=30 & age_start <= 39)
results_list1 <- run_models_inf(df, outcome = "alldeaths")
df <- sccs_cohort_inf_sel %>% filter(age_start >=30 & age_start <= 39) %>% filter(covid == 0)
results_list2 <- run_models_inf(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "30_39"))
results <- export_outcome_inf(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_3039_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_3039_list_20241105")
gc()

######## Age group 40 and 49 #####
df <- sccs_cohort_inf_sel %>% filter(age_start >=40 & age_start <= 49)
results_list1 <- run_models_inf(df, outcome = "alldeaths")
df <- sccs_cohort_inf_sel %>% filter(age_start >=40 & age_start <= 49) %>% filter(covid == 0)
results_list2 <- run_models_inf(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "40_49"))
results <- export_outcome_inf(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_4049_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_4049_list_20241105")
gc()

####### Age group 50 and 59 #####
df <- sccs_cohort_inf_sel %>% filter(age_start >=50 & age_start <= 59)
results_list1 <- run_models_inf(df, outcome = "alldeaths")
df <- sccs_cohort_inf_sel %>% filter(age_start >=50 & age_start <= 59) %>% filter(covid == 0)
results_list2 <- run_models_inf(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "50_59"))
results <- export_outcome_inf(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_5059_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_5059_list_20241105")
gc()

####### Age group 60 and 69 #####
df <- sccs_cohort_inf_sel %>% filter(age_start >=60 & age_start <= 69)
results_list1 <- run_models_inf(df, outcome = "alldeaths")
df <- sccs_cohort_inf_sel %>% filter(age_start >=60 & age_start <= 69) %>% filter(covid == 0)
results_list2 <- run_models_inf(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "60_69"))
results <- export_outcome_inf(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_6069_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_6069_list_20241105")
gc()

####### Age group 70 and 79 #####
df <- sccs_cohort_inf_sel %>% filter(age_start >=70 & age_start <= 79)
results_list1 <- run_models_inf(df, outcome = "alldeaths")
df <- sccs_cohort_inf_sel %>% filter(age_start >=70 & age_start <= 79) %>% filter(covid == 0)
results_list2 <- run_models_inf(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "70_79"))
results <- export_outcome_inf(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_7079_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_7079_list_20241105")
gc()

####### Age group 80 and 89 #####
df <- sccs_cohort_inf_sel %>% filter(age_start >=80 & age_start <= 89)
results_list1 <- run_models_inf(df, outcome = "alldeaths")
df <- sccs_cohort_inf_sel %>% filter(age_start >=80 & age_start <= 89) %>% filter(covid == 0)
results_list2 <- run_models_inf(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "80_89"))
results <- export_outcome_inf(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_8089_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_8089_list_20241105")
gc()

####### Age group 90 #####
df <- sccs_cohort_inf_sel %>% filter(age_start >=90)
results_list1 <- run_models_inf(df, outcome = "alldeaths")
df <- sccs_cohort_inf_sel %>% filter(age_start >=90) %>% filter(covid == 0)
results_list2 <- run_models_inf(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "90"))
results <- export_outcome_inf(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_90_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_inf_90_list_20241105")
gc()

##### Booster series vaccination #####

###### Set parameters - same for all analyses ###### 

# Load dataset

sccs_cohort_booster <- fread(file = "H:/Datafiles/Oversterfte_covidvac/sccs_cohort_booster_20240711.csv", header = TRUE, fill = TRUE)

# Set parameters - same for all analyses

# Parameters for stratification of days of the start of exposure-related risk

expogrp_1 = c(0,1)  # Ignore first day after vaccination (no mortality to be expected related to vaccination on day 1)
expogrp_3 = c(0,1,7,14) # Ignore first day after vaccination, risk period is 3 separate weeks
expogrp_12 = c(0,1,7,14,21,28,35,42,49,56,63,70,77) # Ignore first day after vaccination, risk period is 12 separate weeks

## Parameters for 14 day interval calendar adjustment for seasonality

agegrp = seq(14,844,14) # Divide age groups into 14 day intervals starting from beginning to end of exposure period, accounts for seasonality

######   All ages #####  
df <- sccs_cohort_booster
results_list1 <- run_models_booster(df, outcome = "alldeaths")
df <- sccs_cohort_booster %>% filter(covid == 0)
results_list2 <- run_models_booster(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "allages"))
results <- export_outcome_booster(results_list)
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_allages_20241105.csv")
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_allages_list_20241105")
rm(results, results_list, results_list1, results_list2)
gc()

######## Per age group #####
results_agegroup <- list()
######## Age group 12 and 29 #####
df <- sccs_cohort_booster %>% filter(age_start >=12 & age_start <= 29)
results_list1 <- run_models_booster(df, outcome = "alldeaths")
df <- sccs_cohort_booster %>% filter(age_start >=12 & age_start <= 29) %>% filter(covid == 0)
results_list2 <- run_models_booster(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "12_29"))
results <- export_outcome_booster(results_list)
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_1229_20241105")
results_agegroup[[1]] <- results
save(results_agegroup, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_peragegroup_20241105.csv")
rm(results, results_list, results_list1, results_list2)
gc()

######## Age group 30 and 39 #####
df <- sccs_cohort_booster %>% filter(age_start >=30 & age_start <= 39)
results_list1 <- run_models_booster(df, outcome = "alldeaths")
df <- sccs_cohort_booster %>% filter(age_start >=30 & age_start <= 39) %>% filter(covid == 0)
results_list2 <- run_models_booster(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "30_39"))
results <- export_outcome_booster(results_list)
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_3039_20241105")
results_agegroup[[2]] <- results
save(results_agegroup, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_peragegroup_20241105.csv")
rm(results, results_list, results_list1, results_list2)
gc()

######## Age group 40 and 49 #####
df <- sccs_cohort_booster %>% filter(age_start >=40 & age_start <= 49)
results_list1 <- run_models_booster(df, outcome = "alldeaths")
df <- sccs_cohort_booster %>% filter(age_start >=40 & age_start <= 49) %>% filter(covid == 0)
results_list2 <- run_models_booster(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "40_49"))
results <- export_outcome_booster(results_list)
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_4049_20241105")
results_agegroup[[2]] <- results
save(results_agegroup, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_peragegroup_20241105.csv")
rm(results, results_list, results_list1, results_list2)
gc()

######## Age group 50 and 59 #####
df <- sccs_cohort_booster %>% filter(age_start >=50 & age_start <= 59)
results_list1 <- run_models_booster(df, outcome = "alldeaths")
df <- sccs_cohort_booster %>% filter(age_start >=50 & age_start <= 59) %>% filter(covid == 0)
results_list2 <- run_models_booster(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "50_59"))
results <- export_outcome_booster(results_list)
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_5059_20241105")
results_agegroup[[2]] <- results
save(results_agegroup, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_peragegroup_20241105.csv")
rm(results, results_list, results_list1, results_list2)
gc()

######## Age group 60 and 69 #####
df <- sccs_cohort_booster %>% filter(age_start >=60 & age_start <= 69)
results_list1 <- run_models_booster(df, outcome = "alldeaths")
df <- sccs_cohort_booster %>% filter(age_start >=60 & age_start <= 69) %>% filter(covid == 0)
results_list2 <- run_models_booster(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "60_69"))
results <- export_outcome_booster(results_list)
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_6069_20241105")
results_agegroup[[2]] <- results
save(results_agegroup, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_peragegroup_20241105.csv")
rm(results, results_list, results_list1, results_list2)
gc()

######## Age group 70 and 79 #####
df <- sccs_cohort_booster %>% filter(age_start >=70 & age_start <= 79)
results_list1 <- run_models_booster(df, outcome = "alldeaths")
df <- sccs_cohort_booster %>% filter(age_start >=70 & age_start <= 79) %>% filter(covid == 0)
results_list2 <- run_models_booster(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "70_79"))
results <- export_outcome_booster(results_list)
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_7079_20241105")
results_agegroup[[2]] <- results
save(results_agegroup, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_peragegroup_20241105.csv")
rm(results, results_list, results_list1, results_list2)
gc()

######## Age group 80 and 89 #####
df <- sccs_cohort_booster %>% filter(age_start >=80 & age_start <= 89)
results_list1 <- run_models_booster(df, outcome = "alldeaths")
df <- sccs_cohort_booster %>% filter(age_start >=80 & age_start <= 89) %>% filter(covid == 0)
results_list2 <- run_models_booster(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "80_89"))
results <- export_outcome_booster(results_list)
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_8089_20241105")
results_agegroup[[2]] <- results
save(results_agegroup, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_peragegroup_20241105.csv")
rm(results, results_list, results_list1, results_list2)
gc()

######## Age group >90 #####
df <- sccs_cohort_booster %>% filter(age_start >=90)
results_list1 <- run_models_booster(df, outcome = "alldeaths")
df <- sccs_cohort_booster %>% filter(age_start >=90) %>% filter(covid == 0)
results_list2 <- run_models_booster(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, age_group = "90"))
results <- export_outcome_booster(results_list)
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_90_20241105")
results_agegroup[[2]] <- results
save(results_agegroup, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_peragegroup_20241105.csv")
rm(results, results_list, results_list1, results_list2)
gc()

######## Per sex ######
results_sex <- list()

df <- sccs_cohort_booster %>% filter(sex == "Man")
results_list1 <- run_models_booster(df, outcome = "alldeaths")
df <- sccs_cohort_booster %>% filter(sex == "Man") %>% filter(covid == 0)
results_list2 <- run_models_booster(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, sex = "Man"))
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_men_20241105")
results_sex[[1]] <- export_outcome_booster(results_list)
save(results_sex, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_persex_20241105.csv")
gc()

df <- sccs_cohort_booster %>% filter(sex == "Woman")
results_list1 <- run_models_booster(df, outcome = "alldeaths")
df <- sccs_cohort_booster %>% filter(sex == "Woman") %>% filter(covid == 0)
results_list2 <- run_models_booster(df, outcome = "noncovdeaths")
results_list <- c(results_list1, results_list2)
results_list <- lapply(results_list, function(list) transform(list, sex = "Woman"))
save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_woman_20241105")
results_sex[[2]] <- export_outcome_booster(results_list)
save(results_sex, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_persex_20241105.csv")
gc()

results <- do.call(rbind, results_sex)
row.names(results) <- NULL
write.csv(results, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_booster_persex_20241105.csv")
