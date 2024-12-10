# Function to create custom parameters for SCCS analysis according to number of exposures
create_parameters <- function(df) { 
  # Rename dataset and put in global environment
  df <- df
  # Create adrug and aedrug for 3 and 12-week exposure intervals according to number of exposures.
  if (all(is.na(df$expo2)) || (all(is.na(df$expo3)))) {
    if (all(is.na(df$expo2)) & (all(is.na(df$expo3)))) {
      adrug <<- df$expo1
      aedrug_3 <<- cbind(df$expo1+21) # Three week exposure interval
      aedrug_12 <<- cbind(df$expo1+84) # Twelve week exposure interval
    }
    if (!all(is.na(df$expo2)) & (all(is.na(df$expo3)))) {
      adrug <<- cbind(df$expo1, df$expo2)
      aedrug_3 <<- cbind(df$expo1+21, df$expo2+21)
      aedrug_12 <<- cbind(df$expo1+84, df$expo2+84)
    }} else {
      adrug <<- cbind(df$expo1, df$expo2, df$expo3)
      aedrug_3 <<- cbind(df$expo1+21, df$expo2+21, df$expo3+21)
      aedrug_12 <<- cbind(df$expo1+84, df$expo2+84, df$expo3+84)
    }
  return(df)
  
}

# Function to create custom parameters for SCCS analysis according to number of exposures
create_parameters_booster <- function(df) { 
  # Rename dataset and put in global environment
  df <- df
  
  adrug <<- cbind(df$expo1, df$expo2, df$expo3, df$expo4, df$expo5, df$expo6, df$expo7)
  aedrug_3 <<- cbind(df$expo1+21, df$expo2+21, df$expo3+21, df$expo4+21, df$expo5+21, df$expo6+21, df$expo7+21)
  aedrug_12 <<- cbind(df$expo1+84, df$expo2+84, df$expo3+84, df$expo4+84, df$expo5+84, df$expo6+84, df$expo7+84)
  
  return(df)
}

# Function to run SCCS models for primary series vaccination
run_models_primary <- function(df, outcome) { 
  results_list <- list()
  indicator = 1
  
  df <- df %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_3", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = F, agegrp = agegrp, dataformat = "multi", verbose = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_3_dose", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_3, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_3_week", aedrug = "aedrug_3", expogrp = "expogrp_3", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_3, sameexpopar = F, agegrp = agegrp, dataformat = "multi", verbose = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_3_dose_week", aedrug = "aedrug_3", expogrp = "expogrp_3", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  # Stratified analyses
  
  ## Vaccin type
  
  df_sel <- df %>% filter(vaccin_type_mRNA == 1 | (is.na(vaccin_type_mRNA) & is.na(vaccin_type_nonmRNA))) %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_mRNA, expo2 = expo2_mRNA, expo3 = expo3_mRNA) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_mRNA", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = F, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel, itermax = 500)
  results_list[[indicator]] <- cbind(mod = "mod_3_dose_mRNA", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  df_sel <- df %>% filter(vaccin_type_nonmRNA == 1 | (is.na(vaccin_type_mRNA) & is.na(vaccin_type_nonmRNA))) %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_nonmRNA, expo2 = expo2_nonmRNA, expo3 = expo3_nonmRNA) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_nonmRNA", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = F, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel, itermax = 500)
  results_list[[indicator]] <- cbind(mod = "mod_3_dose_nonmRNA", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  ## Sex
  
  df_sel <- df %>% filter(sex == "Man") %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_men", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  df_sel <- df %>% filter(sex == "Woman") %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_women", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  ## Chronic disease
  
  df_sel <- df %>% filter(chronicdisease_datevac_cat == "No") %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_nocd", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  df_sel <- df %>% filter(chronicdisease_datevac_cat == "One") %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_onecd", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  df_sel <- df %>% filter(chronicdisease_datevac_cat == "Multiple") %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_morecd", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  ## Previous registered positive COVID-19 infection
  
  df_sel <- df %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_prevposinf, expo2 = expo2_prevposinf, expo3 = expo3_prevposinf) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_prevposinf", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  df_sel <- df %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_noprevposinf, expo2 = expo2_noprevposinf, expo3 = expo3_noprevposinf) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_noprevposinf", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  # Sensitivity analyses
  
  ## Model without adjustment background mortality rate
  
  df <- df %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = NULL, dataformat = "multi", verbose = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_3_noseason", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_3, sameexpopar = T, agegrp = NULL, dataformat = "multi", verbose = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_3_week_noseason", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  ## Model with exclusion of individuals with a positive registered SARS-CoV-2 tests in the 8 weeks before vaccination or during the exposure period
  
  df_sens <- df %>% filter(is.na(pos8weeks) & is.na(posexpperiod)) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sens)
  results_list[[indicator]] <- cbind(mod = "mod_3_sens1", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_3, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sens)
  results_list[[indicator]] <- cbind(mod = "mod_3_week_sens1", aedrug = "aedrug_3", expogrp = "expogrp_3", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  ## Model with calender time adjustment using a restricted cubic spline
  
  df <- df %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp_spline, dataformat = "multi", verbose = T, data = df, spline = T, X = 3)
  results_list[[indicator]] <- cbind(mod = "mod_3_spline", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = F, agegrp = agegrp_spline, dataformat = "multi", verbose = T, data = df, spline = T, X = 3)
  results_list[[indicator]] <- cbind(mod = "mod_3_spline_dose", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  df <- df %>% create_parameters()
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_3, sameexpopar = T, agegrp = agegrp_spline, dataformat = "multi", verbose = T, data = df, spline = T, X = 3)
  results_list[[indicator]] <- cbind(mod = "mod_3_spline_week", aedrug = "aedrug_3", expogrp = "expogrp_3", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_3, sameexpopar = F, agegrp = agegrp_spline, dataformat = "multi", verbose = T, data = df, spline = T, X = 3)
  results_list[[indicator]] <- cbind(mod = "mod_3_spline_dose_week", aedrug = "aedrug_3", expogrp = "expogrp_3", mod[["output"]])
  
  ## Model with a 12-week risk period
  
  df <- df %>% create_parameters()
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_12, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_12", aedrug = "aedrug_12", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  df <- df %>% create_parameters()
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_12, 
                            expogrp = expogrp_12, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df, itermax = 500)
  results_list[[indicator]] <- cbind(mod = "mod_12_week", aedrug = "aedrug_12", expogrp = "expogrp_12", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  # Stratified by administrator (GGD, GP or other)
  
  df_sel <- df %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_GGD, expo2 = expo2_GGD, expo3 = expo3_GGD) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_GGD", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1  
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  df_sel <- df %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_GP, expo2 = expo2_GP, expo3 = expo3_GP) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_GP", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1  
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  df_sel <- df %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_other, expo2 = expo2_other, expo3 = expo3_other) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_other", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  # Add variable to results_list with the outcome
  if (outcome == "alldeaths") {
    results_list <- lapply(results_list, function(list) transform(list, outcome = "alldeaths"))
  } else if (outcome == "noncovdeaths") {
    results_list <- lapply(results_list, function(list) transform(list, outcome = "noncovdeaths"))
  } 
  
  # Clean environment
  rm(df, adrug, aedrug_3, aedrug_12, envir = .GlobalEnv)
  
  return(results_list)
  
}

# Function to extract the information from the list with model results
export_outcome_primary <- function(output) {
  
  output <- do.call(rbind, results_list)
  row.names(output) <- NULL
  
  out <- c("outcome", "mod", "term", "ncases", "nexpo1", "nexpo2", "nexpo3", "IRR", "irr_low_ci", "irr_up_ci")
  output <- output %>% dplyr::select(all_of(out))
  
  output1 <- output %>% filter(outcome == "alldeaths")
  output2 <- output %>% filter(outcome == "noncovdeaths")
  
  rbind_results <- function(output) { 
    
    results <- rbind(
      output %>% filter(mod == "mod_3" & term == "adrug2"),
      output %>% filter(mod == "mod_3_week" & (term == "adrug2" | term == "adrug3" | term == "adrug4")),
      output %>% filter(mod == "mod_3_dose" & term == "adrug2"),
      output %>% filter(mod == "mod_3_dose_week" & (term == "adrug2" | term == "adrug3" | term == "adrug4")),
      output %>% filter(mod == "mod_3_dose" & term == "adrug4"),
      output %>% filter(mod == "mod_3_dose_week" & (term == "adrug6" | term == "adrug7" | term == "adrug8")),
      output %>% filter(mod == "mod_3_dose" & term == "adrug6"),
      output %>% filter(mod == "mod_3_dose_week" & (term == "adrug9" | term == "adrug10" | term == "adrug11")),
      
      output %>% filter(mod == "mod_3_mRNA" & term == "adrug2"),
      output %>% filter(mod == "mod_3_nonmRNA" & term == "adrug2"),
      output %>% filter(mod == "mod_3_dose_mRNA" & term == "adrug2"),
      output %>% filter(mod == "mod_3_dose_nonmRNA" & term == "adrug2"),
      output %>% filter(mod == "mod_3_dose_mRNA" & term == "adrug4"),
      output %>% filter(mod == "mod_3_dose_nonmRNA" & term == "adrug4"),
      output %>% filter(mod == "mod_3_dose_mRNA" & term == "adrug6"),
      output %>% filter(mod == "mod_3_dose_nonmRNA" & term == "adrug6"),
      
      output %>% filter(mod == "mod_3_men" & term == "adrug2"),
      output %>% filter(mod == "mod_3_women" & term == "adrug2"),
      
      output %>% filter(mod == "mod_3_nocd" & term == "adrug2"),
      output %>% filter(mod == "mod_3_onecd" & term == "adrug2"),
      output %>% filter(mod == "mod_3_morecd" & term == "adrug2"),
      
      output %>% filter(mod == "mod_3_prevposinf" & term == "adrug2"),
      output %>% filter(mod == "mod_3_noprevposinf" & term == "adrug2"),
      
      output %>% filter(mod == "mod_3_GGD" & term == "adrug2"),
      output %>% filter(mod == "mod_3_GP" & term == "adrug2"),
      output %>% filter(mod == "mod_3_other" & term == "adrug2"),
      
      output %>% filter(mod == "mod_3_noseason" & term == "adrug2"),
      output %>% filter(mod == "mod_3_week_noseason" & (term == "adrug2" | term == "adrug3" | term == "adrug4")),
      
      output %>% filter(mod == "mod_3_sens1" & term == "adrug2"),
      output %>% filter(mod == "mod_3_week_sens1" & (term == "adrug2" | term == "adrug3" | term == "adrug4")),
      
      output %>% filter(mod == "mod_3" & term == "adrug2"),
      
      output %>% filter(mod == "mod_12" & term == "adrug2"),
      output %>% filter(mod == "mod_12_week" & (term %in% c(paste0("adrug", 2:13)))),
      
      output %>% filter(mod == "mod_3_spline" & term == "adrug2"),
      output %>% filter(mod == "mod_3_spline_week" & (term == "adrug2" | term == "adrug3" | term == "adrug4")),
      output %>% filter(mod == "mod_3_spline_dose" & term == "adrug2"),
      output %>% filter(mod == "mod_3_spline_dose_week" & (term == "adrug2" | term == "adrug3" | term == "adrug4")),
      output %>% filter(mod == "mod_3_spline_dose" & term == "adrug4"),
      output %>% filter(mod == "mod_3_spline_dose_week" & (term == "adrug6" | term == "adrug7" | term == "adrug8")),
      output %>% filter(mod == "mod_3_spline_dose" & term == "adrug6"),
      output %>% filter(mod == "mod_3_spline_dose_week" & (term == "adrug9" | term == "adrug10" | term == "adrug11")))
    
    return(results)                   
  }
  
  results <- rbind(rbind_results(output1), rbind_results(output2))
  results$IRR <- round(results$IRR, 4)
  results$irr_low_ci <- round(results$irr_low_ci, 4)
  results$irr_up_ci <- round(results$irr_up_ci, 4)
  
  return(results)
  
}

# Function to run SCCS models for booster series vaccination
run_models_booster <- function(df, outcome) { 
  # Make list to store results
  results_list <- list()
  indicator = 1
  
  df <- df %>% create_parameters_booster()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, booster = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_3", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = F, agegrp = agegrp, dataformat = "multi", verbose = T, booster = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_3_dose", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1   
  save(results_list, file = "H:/Annemarijn/oversterfte_covidvac/Output/Modelresultaten/results_list")
  
  # Add variable to results_list with the outcome
  if (outcome == "alldeaths") {
    results_list <- lapply(results_list, function(list) transform(list, outcome = "alldeaths"))
  } else if (outcome == "noncovdeaths") {
    results_list <- lapply(results_list, function(list) transform(list, outcome = "noncovdeaths"))
  } 
  
  # Clean environment
  rm(df, adrug, aedrug_3, aedrug_12, envir = .GlobalEnv)
  
  return(results_list)
  
}

# Function to extract the information from the list with model results
export_outcome_booster <- function(output) {
  
  output <- do.call(rbind, results_list)
  row.names(output) <- NULL
  
  out <- c("outcome", "mod", "term", "ncases", "nexpo1", "nexpo2", "nexpo3", "nexpo4", "nexpo5", "nexpo6", "nexpo7", "IRR", "irr_low_ci", "irr_up_ci")
  output <- output %>% dplyr::select(all_of(out))
  
  output1 <- output %>% filter(outcome == "alldeaths")
  output2 <- output %>% filter(outcome == "noncovdeaths")
  
  rbind_results <- function(output) { 
    
    results <- rbind(
      output %>% filter(mod == "mod_3" & term == "adrug2"),
      output %>% filter(mod == "mod_3_dose" & (term == "adrug2")),
      output %>% filter(mod == "mod_3_dose" & (term == "adrug4")),
      output %>% filter(mod == "mod_3_dose" & (term == "adrug6")),
      output %>% filter(mod == "mod_3_dose" & (term == "adrug8")),
      output %>% filter(mod == "mod_3_dose" & (term == "adrug10")),
      output %>% filter(mod == "mod_3_dose" & (term == "adrug12")),
      output %>% filter(mod == "mod_3_dose" & (term == "adrug14")))
    
    return(results)                   
  }
  
  results <- rbind(rbind_results(output1), rbind_results(output2))
  results$IRR <- round(results$IRR, 4)
  results$irr_low_ci <- round(results$irr_low_ci, 4)
  results$irr_up_ci <- round(results$irr_up_ci, 4)
  
  return(results)
  
}

# Function to run SCCS models for positive registered SARS-CoV-2 tests
run_models_inf <- function(df, outcome) {  
  # Make list to store results
  results_list <- list()
  indicator = 1
  
  df <- df %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_inf, expo2 = expo2_inf, expo3 = expo3_inf) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_3", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_3, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_3_week", aedrug = "aedrug_3", expogrp = "expogrp_3", mod[["output"]])
  rm(mod)
  indicator = indicator + 1
  
  # Stratified analyses
  
  ## Registered vaccination within 6 months before infection
  ### Only ran if there are exposures
  
  df_sel <- df %>% filter(vaccin_status_vac_6m == 1 | (is.na(vaccin_status_vac_6m) & is.na(vaccin_status_unvac_6m))) %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_inf_vac_6m, expo2 = expo2_inf_vac_6m, expo3 = expo3_inf_vac_6m) %>% create_parameters()
  
  if(all(is.na(df_sel$expo1))) {
    output_empty <- data.frame(matrix(ncol = 17, nrow = 1))
    colnames(output_empty) <- c("ncases", "nevents", "nexpo1", "nexpo2", "nexpo3", "term", "coef", "IRR", "se", "irr_low_ci_mod", "irr_up_ci_mod", "p_value_mod", "ses", "irr_low_ci", "irr_up_ci", "z", "p_value")
    
    results_list[[indicator]] <- cbind(mod = "mod_3_vaccinated_6m", aedrug = "aedrug_3", expogrp = "expogrp_1", output_empty)
    indicator = indicator + 1
    results_list[[indicator]] <- cbind(mod = "mod_3_vaccinated_week_6m", aedrug = "aedrug_3", expogrp = "expogrp_1", output_empty)
    indicator = indicator + 1
    
  } else { 
    
    mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                              expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel, tolerance = 1e-05, itermax = 100)
    results_list[[indicator]] <- cbind(mod = "mod_3_vaccinated_6m", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
    rm(mod)
    indicator = indicator + 1
    
    mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                              expogrp = expogrp_3, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel, tolerance = 1e-05, itermax = 100)
    results_list[[indicator]] <- cbind(mod = "mod_3_vaccinated_week_6m", aedrug = "aedrug_3", expogrp = "expogrp_3", mod[["output"]])
    rm(mod)
    indicator = indicator + 1
    
  }
  
  ## No registered vaccination within 6 months before infection
  
  df_sel <- df %>% filter(vaccin_status_unvac_6m == 1 | (is.na(vaccin_status_vac_6m) & is.na(vaccin_status_unvac_6m))) %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_inf_unvac_6m, expo2 = expo2_inf_unvac_6m) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_notvaccinated_6m", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_3, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_notvaccinated_week_6m", aedrug = "aedrug_3", expogrp = "expogrp_3", mod[["output"]])
  rm(mod)
  indicator = indicator + 1
  
  ## Sex
  
  df_sel <- df %>% filter(sex == "Man") %>% create_parameters()
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_men", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1
  
  df_sel <- df %>% filter(sex == "Woman") %>% create_parameters()
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_women", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1
  
  ## Chronic disease
  
  df_sel <- df %>% filter(chronicdisease_dateinf_cat == "No") %>% create_parameters()
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_nocd", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1
  
  df_sel <- df %>% filter(chronicdisease_dateinf_cat == "One") %>% create_parameters()
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_onecd", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1
  
  df_sel <- df %>% filter(chronicdisease_dateinf_cat == "Multiple") %>% create_parameters()
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sel)
  results_list[[indicator]] <- cbind(mod = "mod_3_morecd", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sel)
  indicator = indicator + 1
  
  # Second mortality wave
  
  df_wave2 <- df %>% filter(wave == 2) %>% dplyr::select(!c(case, start, end, event, expo1, expo2, expo3)) %>%
    dplyr::rename(case = case_inf_wave, start = start_wave2, end = end_wave2, event = event_wave2, 
                  expo1 = expo1_inf_wave2, expo2 = expo2_inf_wave2, expo3 = expo3_inf_wave2) %>% create_parameters()
  
  agegrp_wave2 = seq(14,119,14)
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp_wave2, dataformat = "multi", verbose = T, data = df_wave2)
  results_list[[indicator]] <- cbind(mod = "mod_3_wave2", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_3, sameexpopar = T, agegrp = agegrp_wave2, dataformat = "multi", verbose = T, data = df_wave2, tolerance = 1e-06, itermax = 200)
  results_list[[indicator]] <- cbind(mod = "mod_3_week_wave2", aedrug = "aedrug_3", expogrp = "expogrp_3", mod[["output"]])
  rm(mod, df_wave2, agegrp_wave2)
  indicator = indicator + 1
  
  ## Third mortality wave
  
  df_wave3 <- df %>% filter(wave == 3) %>% dplyr::select(!c(case, start, end, event, expo1, expo2, expo3)) %>%
    dplyr::rename(case = case_inf_wave, start = start_wave3, end = end_wave3, event = event_wave3, 
                  expo1 = expo1_inf_wave3, expo2 = expo2_inf_wave3, expo3 = expo3_inf_wave3) %>% create_parameters()
  
  agegrp_wave3 = seq(14,133,14)
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp_wave3, dataformat = "multi", verbose = T, data = df_wave3)
  results_list[[indicator]] <- cbind(mod = "mod_3_wave3", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_3, sameexpopar = T, agegrp = agegrp_wave3, dataformat = "multi", verbose = T, data = df_wave3)
  results_list[[indicator]] <- cbind(mod = "mod_3_week_wave3", aedrug = "aedrug_3", expogrp = "expogrp_3", mod[["output"]])
  rm(mod, df_wave3, agegrp_wave3)
  indicator = indicator + 1
  
  # Sensitivity analyses
  
  ## Model without adjustment background mortality rate
  df <- df %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = NULL, dataformat = "multi", verbose = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_3_noseason", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_3, sameexpopar = T, agegrp = NULL, dataformat = "multi", verbose = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_3_week_noseason", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1
  
  ## Induction interval  2-, 4-, 7- days or 14 days before registration date
  
  df_sens2 <- df %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_inf_sens2, expo2 = expo2_inf_sens2, expo3 = expo3_inf_sens2) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sens2)
  results_list[[indicator]] <- cbind(mod = "mod_3_sens2", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sens2)
  indicator = indicator + 1
  
  df_sens4 <- df %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_inf_sens4, expo2 = expo2_inf_sens4, expo3 = expo3_inf_sens4) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sens4)
  results_list[[indicator]] <- cbind(mod = "mod_3_sens4", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sens4)
  indicator = indicator + 1
  
  df_sens7 <- df %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_inf_sens7, expo2 = expo2_inf_sens7, expo3 = expo3_inf_sens7) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sens7)
  results_list[[indicator]] <- cbind(mod = "mod_3_sens7", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sens7)
  indicator = indicator + 1
  
  df_sens14 <- df %>% dplyr::select(!c(expo1, expo2, expo3)) %>%
    dplyr::rename(expo1 = expo1_inf_sens14, expo2 = expo2_inf_sens14, expo3 = expo3_inf_sens14) %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df_sens14)
  
  results_list[[indicator]] <- cbind(mod = "mod_3_sens14", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod, df_sens14)
  indicator = indicator + 1
  
  ## Model with calender time adjustment using a restricted cubic spline
  
  df <- df %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_3, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp_spline, dataformat = "multi", verbose = T, data = df, spline = T, X = 3)
  results_list[[indicator]] <- cbind(mod = "mod_3_spline", aedrug = "aedrug_3", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1
  
  ## Model with a 12-week risk period
  
  df <- df %>% create_parameters()
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_12, 
                            expogrp = expogrp_1, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df)
  results_list[[indicator]] <- cbind(mod = "mod_12", aedrug = "aedrug_12", expogrp = "expogrp_1", mod[["output"]])
  rm(mod)
  indicator = indicator + 1
  
  mod <- eventdependexp_adj(indiv = case, astart = start, aend = end, aevent = event, adrug = adrug, aedrug = aedrug_12, 
                            expogrp = expogrp_12, sameexpopar = T, agegrp = agegrp, dataformat = "multi", verbose = T, data = df, itermax = 500)
  results_list[[indicator]] <- cbind(mod = "mod_12_week", aedrug = "aedrug_12", expogrp = "expogrp_12", mod[["output"]])
  rm(mod)
  indicator = indicator + 1
  
  # Add variable to results_list with the outcome
  if (outcome == "alldeaths") {
    results_list <- lapply(results_list, function(list) transform(list, outcome = "alldeaths"))
  } else if (outcome == "noncovdeaths") {
    results_list <- lapply(results_list, function(list) transform(list, outcome = "noncovdeaths"))
  } 
  
  # Clean environment
  rm(df, adrug, aedrug_3, aedrug_12, envir = .GlobalEnv)
  
  return(results_list)
  
}

export_outcome_inf <- function(output) {
  
  output <- do.call(rbind, results_list)
  row.names(output) <- NULL
  
  out <- c("outcome", "mod", "term", "ncases", "nexpo1", "nexpo2", "nexpo3", "IRR", "irr_low_ci", "irr_up_ci")
  output <- output %>% dplyr::select(all_of(out))
  
  output1 <- output %>% filter(outcome == "alldeaths")
  output2 <- output %>% filter(outcome == "noncovdeaths")
  
  rbind_results <- function(output) { 
    
    results <- rbind(
      output %>% filter(mod == "mod_3" & term == "adrug1"),
      output %>% filter(mod == "mod_3_week" & (term == "adrug1" | term == "adrug2" | term == "adrug3")),
      
      output %>% filter(mod == "mod_3_vaccinated_6m" & term == "adrug1"),
      output %>% filter(mod == "mod_3_vaccinated_week_6m" & (term == "adrug1" | term == "adrug2" | term == "adrug3")),
      
      output %>% filter(mod == "mod_3_notvaccinated_6m" & term == "adrug1"),
      output %>% filter(mod == "mod_3_notvaccinated_week_6m" & (term == "adrug1" | term == "adrug2" | term == "adrug3")),
      
      output %>% filter(mod == "mod_3_men" & term == "adrug1"),
      output %>% filter(mod == "mod_3_women" & term == "adrug1"),
      
      output %>% filter(mod == "mod_3_nocd" & term == "adrug1"),
      output %>% filter(mod == "mod_3_onecd" & term == "adrug1"),
      output %>% filter(mod == "mod_3_morecd" & term == "adrug1"),
      
      output %>% filter(mod == "mod_3_wave2" & term == "adrug1"),
      output %>% filter(mod == "mod_3_week_wave2" & (term == "adrug1" | term == "adrug2" | term == "adrug3")),
      
      output %>% filter(mod == "mod_3_noseason" & term == "adrug2"),
      output %>% filter(mod == "mod_3_week_noseason" & (term == "adrug2" | term == "adrug3" | term == "adrug4")),
      
      output %>% filter(mod == "mod_3_wave3" & term == "adrug1"),
      output %>% filter(mod == "mod_3_week_wave3" & (term == "adrug1" | term == "adrug2" | term == "adrug3")),
      
      output %>% filter(mod == "mod_3_sens2" & term == "adrug1"),
      output %>% filter(mod == "mod_3_sens4" & term == "adrug1"),
      output %>% filter(mod == "mod_3_sens7" & term == "adrug1"),
      output %>% filter(mod == "mod_3_sens14" & term == "adrug1"),
      
      output %>% filter(mod == "mod_3_spline" & term == "adrug1"),
      
      output %>% filter(mod == "mod_12" & term == "adrug1"),
      output %>% filter(mod == "mod_12_week" & (term %in% c(paste0("adrug", 1:12)))))
    
    return(results)           
    
  }
  
  results <- rbind(rbind_results(output1), rbind_results(output2))
  results$IRR <- round(results$IRR, 4)
  results$irr_low_ci <- round(results$irr_low_ci, 4)
  results$irr_up_ci <- round(results$irr_up_ci, 4)
  
  return(results)
  
}

# Function to run SCCS analysis with event-dependent exposures
## Function adapted from eventdepenexp from SCCS package
### Adaptions:
#### Option to include adjustment for spline (options spline = T) with specified number of knots (X = k)
#### Different outputs according to 1-3 exposures or 1-7 exposures (booster = T)
#### Instead of a matrix calculation to obtain SE, the SEs are obtained in a loop function. This takes much more time, but avoids memory issues. The progresses is printed every 100 iterations.

eventdependexp_adj <- function(indiv, astart, aend, aevent, adrug, aedrug, expogrp = 0, 
                               sameexpopar = T, agegrp = NULL, dataformat = "stack", verbose = F, 
                               tolerance = 1e-08, itermax = 100, data, spline = F, X = 3, booster = F) 
{
  yon <- deparse(substitute(adrug))
  yon1 <- as.formula(paste("z", "~", yon))
  adrugcolnames <- all.vars(yon1, functions = FALSE, unique = TRUE)[-1]
  adrug <- eval(substitute(adrug), data, parent.frame())
  if ((dataformat == "multi" & !is.null(ncol(adrug)))) {
    adrug <- data.frame(adrug)
    adrug <- list(adrug)
  }
  else if (dataformat == "stack" & !is.null(ncol(adrug))) {
    adrug <- data.frame(adrug)
    adrug1 <- list()
    for (i in 1:ncol(adrug)) {
      adrug1[[i]] <- adrug[, i]
    }
    adrug <- adrug1
  }
  else if (length(adrugcolnames) == 1 & length(adrug) != 1) {
    adrug <- list(adrug)
  }
  else {
    adrug <- adrug
  }
  for (i in 1:length(adrug)) {
    adrug[[i]] <- data.frame(adrug[[i]])
  }
  ncoladrug <- NULL
  for (i in 1:length(adrug)) {
    ncoladrug[i] <- ncol(adrug[[i]])
  }
  for (i in 1:length(adrug)) {
    colnames(adrug[[i]]) <- adrugcolnames[c(1, cumsum(ncoladrug) + 
                                              1)[-(length(ncoladrug) + 1)][i]:cumsum(ncoladrug)[i]]
  }
  colname <- adrugcolnames[1]
  indiv <- eval(substitute(indiv), data, parent.frame())
  astart <- eval(substitute(astart), data, parent.frame())
  aend <- eval(substitute(aend), data, parent.frame())
  aevent <- eval(substitute(aevent), data, parent.frame())
  aedrug <- eval(substitute(aedrug), data, parent.frame())
  if ((dataformat == "multi" & !is.null(ncol(aedrug)))) {
    aedrug <- data.frame(aedrug)
    aedrug <- list(aedrug)
  }
  else if (dataformat == "stack" & !is.null(ncol(aedrug))) {
    aedrug <- data.frame(aedrug)
    aedrug1 <- list()
    for (i in 1:ncol(aedrug)) {
      aedrug1[[i]] <- aedrug[, i]
    }
    aedrug <- aedrug1
  }
  else if (length(adrugcolnames) == 1 & length(aedrug) != 1) {
    aedrug <- list(aedrug)
  }
  else {
    aedrug <- aedrug
  }
  if (dataformat == "stack") {
    adrug_all <- list()
    for (i in 1:length(adrug)) {
      adrug_all[[i]] <- adrug_matrix(indiv, aevent, adrug[[i]])
    }
    aedrug_all <- list()
    for (i in 1:length(aedrug)) {
      aedrug_all[[i]] <- adrug_matrix(indiv, aevent, aedrug[[i]])
    }
  }
  else if (dataformat == "multi") {
    adrug_all <- adrug
    aedrug_all <- aedrug
    adrug_all <- list()
    for (i in 1:length(adrug)) {
      adrug_all[[i]] <- cbind(indiv, aevent, adrug[[i]])
      adrug_all[[i]] <- data.frame(adrug_all[[i]][order(adrug_all[[i]][, 
                                                                       1], adrug_all[[i]][, 3]), -c(1, 2)])
    }
    aedrug_all <- list()
    for (i in 1:length(aedrug)) {
      aedrug_all[[i]] <- cbind(indiv, aevent, aedrug[[i]])
      aedrug_all[[i]] <- data.frame(aedrug_all[[i]][order(aedrug_all[[i]][, 
                                                                          1], aedrug_all[[i]][, 3]), -c(1, 2)])
    }
  }
  else {
    stop("dataformat should be multi or stack")
  }
  adrug <- as.matrix(adrug_all[[1]])
  aedrug <- as.matrix(aedrug_all[[1]])
  data1 <- data.frame(unique(cbind(indiv, aevent, astart, aend)))
  data1 <- data1[order(data1$indiv), ]
  if (length(data1$indiv) != length(unique(data1$indiv))) {
    warning("Multiple events per case detected: analysis restricted to first events")
    ord <- order(data1$indiv)
    minev <- unlist(tapply(data1$aevent[ord], data1$indiv[ord], 
                           function(x) rep(min(x), length(x))), use.names = F)
    data1$aevent <- minev[order(ord)]
    first.event <- 1 - duplicated(data1$indiv)
  }
  else {
    first.event <- rep(1, length(data1$indiv))
  }
  if (expogrp[[1]][1] < 0) {
    stop("Pre-exposure risk periods not allowed")
  }
  aevent <- data1$aevent
  nrem <- 0
  if (is.matrix(adrug)) {
    for (i in 1:ncol(adrug)) {
      nrem <- nrem + sum(ifelse(first.event == 1, aevent < 
                                  adrug[, i], 0), na.rm = T)
      adrug[, i] <- ifelse(aevent < adrug[, i], NA, adrug[, 
                                                          i])
      aedrug[, i] <- ifelse(aevent < adrug[, i], NA, aedrug[, 
                                                            i])
    }
  }
  else {
    nrem <- nrem + sum(ifelse(first.event == 1, aevent < 
                                adrug, 0), na.rm = T)
    adrug <- ifelse(aevent < adrug, NA, adrug)
    aedrug <- ifelse(aevent < adrug, NA, aedrug)
  }
  if (verbose == T) {
    cat(paste0("No. exposures after first event (treated as missing): ", 
               nrem))
    cat("\n")
  }
  combinedoses <- ifelse(sameexpopar == T, 1, 0)
  riskstart <- ifelse(is.null(expogrp), 0, expogrp[[1]][1])
  if (riskstart > 0) {
    expogrP <- if (is.list(expogrp)) {
      c(0, as.vector(expogrp[[1]]))
    }
    else {
      c(0, expogrp)
    }
  }
  else {
    expogrP <- if (is.list(expogrp)) {
      as.vector(expogrp[[1]])
    }
    else {
      expogrp
    }
  }
  all.data <- data.frame(indiv = data1$indiv, astart = data1$astart, 
                         aend = data1$aend, aevent = aevent, adrug, aedrug, first.event = first.event)
  all.data_fe <- all.data[all.data$first.event == 1, ]
  adrug_fe <- adrug[first.event == 1, ]
  aedrug_fe <- aedrug[first.event == 1, ]
  base.dat <- subset(formatdata(indiv = indiv, astart = astart, 
                                aend = aend, aevent = aevent, adrug = list(adrug_fe), 
                                aedrug = list(aedrug_fe), expogrp = expogrP, sameexpopar = F, 
                                agegrp = agegrp, dataformat = "multi", data = all.data_fe), 
                     select = c(1, 2, 7, 8, 6))
  numrisk <- length(expogrP)
  base.dat$risk <- as.numeric(as.character(base.dat[, 4]))
  base.dat$dose <- ceiling(base.dat$risk/numrisk)
  if (riskstart > 0) {
    base.dat$risk <- ifelse(base.dat$risk%%numrisk == 1, 
                            0, base.dat$risk - base.dat$dose)
  }
  maxdose <- max(base.dat$dose)
  increment <- ceiling(log10(maxdose + 1))
  base.dat$indivL <- as.numeric(as.character(base.dat$indivL))
  base.dat$expo <- rep(0, length(base.dat$indivL))
  base.dat$weight <- ifelse(base.dat$event == 1, base.dat$risk, 
                            0)
  base.dat$maxw <- unlist(tapply(base.dat$weight, base.dat$indivL, 
                                 function(x) rep(max(x), length(x))), use.names = F)
  base.dat$previ <- rep(0, length(base.dat$indivL))
  temp <- unlist(tapply(base.dat[, 4], base.dat$indivL, function(x) rep(x[1], 
                                                                        length(x))), use.names = F)
  base.dat$previ <- ifelse(temp == 0, 0, -1)
  rm(temp)
  stack.dat <- subset(base.dat, base.dat$previ == 0)
  nlevel0 <- sum(unique(stack.dat$indivL) > 0)
  if (verbose == T) {
    cat(paste0("No. events included at stack level ", 0, 
               ": ", nlevel0))
    cat("\n")
  }
  for (i in 1:maxdose) {
    temp1 <- ifelse(base.dat$dose == i, i, 0)
    temp2 <- unlist(tapply(temp1, base.dat$indivL, cumsum), 
                    use.names = F)
    temp3 <- ifelse((temp1 == temp2) & temp1 > 0, i, 0)
    base.dat$previ <- unlist(tapply(temp3, base.dat$indivL, 
                                    cumsum), use.names = F)
    rm(temp1, temp2, temp3)
    stack.dat.i <- subset(base.dat, base.dat$previ == i)
    stack.dat.i$indivL <- stack.dat.i$indivL + i/(10^increment)
    stack.dat.i$expo <- ifelse(stack.dat.i$dose == i, stack.dat.i$risk, 
                               0)
    stack.dat.i$weight <- ifelse(stack.dat.i$event == 1 & 
                                   stack.dat.i$dose > i, stack.dat.i$risk, 0)
    stack.dat.i$maxw <- unlist(tapply(stack.dat.i$weight, 
                                      stack.dat.i$indivL, function(x) rep(max(x), length(x))), 
                               use.names = F)
    stack.dat <- rbind(stack.dat, stack.dat.i)
    nleveli <- sum(unique(stack.dat.i$indivL) > 0)
    if (verbose == T) {
      cat(paste0("No. events included at stack level ", 
                 i, ": ", nleveli))
      cat("\n")
    }
  }
  fitrisk <- ifelse(riskstart > 0, numrisk - 1, numrisk)
  if (combinedoses == 1) {
    stack.dat$expo <- as.factor(ifelse(stack.dat$expo > 0, 
                                       stack.dat$risk - fitrisk * (stack.dat$dose - 1), 
                                       0))
    stack.dat$weight <- ifelse(stack.dat$weight > 0, stack.dat$weight - 
                                 fitrisk * (stack.dat$dose - 1), 0)
    stack.dat$maxw <- ifelse(stack.dat$maxw > 0, stack.dat$maxw - 
                               (ceiling(stack.dat$maxw/fitrisk) - 1) * fitrisk, 
                             0)
  }
  else {
    stack.dat$expo <- as.factor(stack.dat$expo)
  }
  stack.dat$indivL <- as.factor(stack.dat$indivL)
  stack.dat$expo1 <- stack.dat$expo
  colnames(stack.dat)[which(names(stack.dat) == "expo1")] <- colname
  
  if (spline == T) {
    if (week_start != start_date) { 
      stack.dat <- stack.dat %>% mutate(week_length = ifelse(age == 1, as.numeric(week_start - start_date) + 7, 
                                                             ifelse(age == max(as.numeric(age)), as.numeric(end_date - max(week_start_days) + 1), 7)))
      stack.dat$age <- as.numeric(stack.dat$age)
    } 
    else {
      stack.dat <- stack.dat %>% mutate(week_length = 7)
      stack.dat$age <- as.numeric(stack.dat$age)
    }
  }
  
  lenbeta <- ifelse(combinedoses == 1, fitrisk, maxdose * fitrisk)
  
  if (spline == T) {
    lenalpha <- length(seq(1:X)-1)
  } else {
    lenalpha <- length(agegrp)
  }
  
  if (spline == T) {
    fmla1 = as.formula(paste("wevent ~", colname[1], "+ ns(age, df =", X, ", Boundary.knots = quantile(age, c(.10,.90)))"))
  }
  else if (lenalpha > 0) {
    fmla1 = as.formula(paste("wevent ~", colname[1], "+", 
                             "age"))
  }
  else {
    fmla1 <- as.formula(paste("wevent ~", colname[1]))
  }
  
  diff <- 1
  numiter <- 0
  beta <- rep(0, lenbeta)
  while (diff > tolerance * lenbeta & numiter < itermax) {
    stack.dat$wevent <- stack.dat$event
    for (i in 1:lenbeta) {
      stack.dat$wevent <- ifelse(stack.dat$weight == i, 
                                 exp(-beta[i]), stack.dat$wevent)
    }
    betaold <- beta
    indivL <- stack.dat$indivL
    interval <- stack.dat$interval
    
    if (spline == T) {
      mod <- suppressWarnings(gnm(fmla1, eliminate = indivL, 
                                  offset = log(week_length), family = poisson, data = stack.dat))
    } else {
      mod <- suppressWarnings(gnm(fmla1, eliminate = indivL, 
                                  offset = log(interval), family = poisson, data = stack.dat))
    }
    
    beta <- as.vector(mod$coefficients[1:lenbeta])
    diff <- as.numeric(sqrt(t(beta - betaold) %*% (beta - 
                                                     betaold)))
    numiter <- numiter + 1
    if (verbose == T) {
      cat(paste("iteration:", numiter, "diff:", diff))
      cat("\n")
      cat(paste("beta:", 1:lenbeta, beta))
      cat("\n")
    }
    if (numiter == itermax) 
      warning("Maximum number of iterations reached")
  }
  beta <- as.vector(mod$coefficients[1:lenbeta])
  alpha <- NULL
  if (lenalpha > 0) {
    a1 <- lenbeta + 1
    a2 <- lenbeta + lenalpha
    alpha <- as.vector(mod$coefficients[a1:a2])
    rm(a1, a2)
  }
  ncases <- sum(unique(floor(as.numeric(as.character(stack.dat$indivL)))) > 
                  0)
  
  residual <- (stack.dat$wevent - mod$fitted.values)
  
  Hmat <- (-1) * solve(vcov(mod))
  
  if (lenalpha > 0) {
    if (spline == T) {
      Salpha <- matrix(rep(1, length(stack.dat$indivL) * lenalpha), 
                       ncol = lenalpha)
      for (j in 1:lenalpha) {
        Salpha[, j] <- Salpha[, j] * (as.numeric(as.character(stack.dat$age)) == 
                                        j + 1)
      }
    }
    else {
      Salpha <- matrix(rep(1, as.numeric(length(stack.dat$indivL)) * as.numeric(lenalpha)), 
                       ncol = lenalpha)
      for (j in 1:lenalpha) {
        Salpha[, j] <- Salpha[, j] * (as.numeric(as.character(stack.dat$age)) == 
                                        j + 1)
      }
    }
  }
  
  Sbeta <- matrix(rep(1, length(stack.dat$indivL) * lenbeta), 
                  ncol = lenbeta)
  for (j in 1:lenbeta) {
    Sbeta[, j] <- Sbeta[, j] * (as.numeric(as.character(stack.dat$expo)) == 
                                  j) * (stack.dat$previ >= 1)
    
  }
  if (lenalpha > 0) {
    eeealpha <- matrix(0, nrow = ncases, ncol = lenalpha)
  }
  
  eeebeta <- matrix(0, nrow = ncases, ncol = lenbeta)
  
  indivL <- floor(as.numeric(as.character(stack.dat$indivL)))
  
  for (i in 1:ncases) {
    select = (indivL == i)
    tempresidual = residual[select]
    if (lenalpha > 0) {
      eeealpha[i, ] = tempresidual %*% Salpha[select, ]
    }
    eeebeta[i, ] = tempresidual %*% Sbeta[select, ]
    
    if (verbose == T & (i %% 100 == 0)) {
      print(i)
    }
  }
  
  if (lenalpha > 0) {
    Mmat <- cbind(eeebeta, eeealpha)
  } else {
    Mmat <- eeebeta
  }
  
  deebeta <- matrix(apply(expand.grid(1:lenbeta, 1:lenbeta), 
                          1, function(x) {
                            (-1) * sum(residual[(stack.dat$previ >= 1) & 
                                                  (as.numeric(as.character(stack.dat$expo)) == 
                                                     x[1]) & (stack.dat$maxw == x[2])])
                          }), ncol = lenbeta, byrow = F)
  
  if (lenalpha > 0) {
    deealpha <- matrix(apply(expand.grid(1:lenalpha, 1:lenbeta), 
                             1, function(x) {
                               (-1) * sum(residual[(as.numeric(as.character(stack.dat$age)) == 
                                                      x[1] + 1) & (stack.dat$maxw == x[2])])
                             }), ncol = lenbeta, byrow = F)
  }
  if (lenalpha > 0) {
    Amat <- cbind(rbind(deebeta, deealpha), matrix(rep(0, 
                                                       lenalpha * (lenbeta + lenalpha)), ncol = lenalpha, 
                                                   nrow = lenbeta + lenalpha))
  } else {
    Amat <- deebeta
  }
  
  bread <- Amat + Hmat
  
  temp <- Mmat %*% t(solve(bread))
  
  ses = rep(0, lenbeta + lenalpha) # SE based on sandwich model (ses)
  for (i in 1:(lenbeta + lenalpha)) {
    ses[i] = sqrt(sum(temp[, i]^2))
  }
  
  b <- coef(mod)
  se <- sqrt(diag(vcov(mod))) # SE based on model
  
  z_score <- b / ses # Calculation of p-values based on ses
  p_value <- 2 * (1 - pnorm(abs(z_score)))
  
  if (booster == T) {
    
    res <- data.frame(ncases = ncases,
                      nevents = nrow(data),
                      nexpo1 = unname(table(!is.na(data$expo1))[2]),
                      nexpo2 = unname(table(!is.na(data$expo2))[2]),
                      nexpo3 = unname(table(!is.na(data$expo3))[2]),
                      nexpo4 = unname(table(!is.na(data$expo4))[2]),
                      nexpo5 = unname(table(!is.na(data$expo5))[2]),
                      nexpo6 = unname(table(!is.na(data$expo6))[2]),
                      nexpo7 = unname(table(!is.na(data$expo7))[2]),
                      term = names(b),
                      coef = b,
                      IRR = exp(b),
                      se = se,
                      irr_low_ci_mod = exp(b - (1.96 * se)),
                      irr_up_ci_mod = exp(b + (1.96 * se)),
                      # p_value_mod = summary(mod)$coefficients[,4],
                      ses = ses,
                      irr_low_ci = exp(b - (1.96 * ses)),
                      irr_up_ci = exp(b + (1.96 * ses)),
                      z = z_score,
                      p_value = p_value)
    
  } else {
    
    res <- data.frame(ncases = ncases,
                      nevents = nrow(data1),
                      nexpo1 = data %>% filter(!is.na(expo1)) %>% nrow(),
                      nexpo2 = data %>% filter(!is.na(expo2)) %>% nrow(),
                      nexpo3 = data %>% filter(!is.na(expo3)) %>% nrow(),
                      term = names(b),
                      coef = b,
                      IRR = exp(b),
                      se = se,
                      irr_low_ci_mod = exp(b - (1.96 * se)),
                      irr_up_ci_mod = exp(b + (1.96 * se)),
                      p_value_mod = summary(mod)$coefficients[,4],
                      ses = ses,
                      irr_low_ci = exp(b - (1.96 * ses)),
                      irr_up_ci = exp(b + (1.96 * ses)),
                      z = z_score,
                      p_value = p_value)
    
  }
  
  result <- list(mod = mod, output = res, 
                 ncases = length(unique(data1$indiv)), nevents = nrow(data1))
  return(result)
  
}