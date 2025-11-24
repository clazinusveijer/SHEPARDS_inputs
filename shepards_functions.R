#### functions to be used in the SHEPARDS model ###
### Developed as part of the COVend project
## UMCG - Global Health unit © November 2025
## Correspondence: c.veijer@gmail.com

# Define Simulation Functions for the short-term model

## 1.1 Probability function

## The `Probs` function updates the transition probabilities of every cycle is shown below.

Probs <- function(M_t, df_X, t, Trt, n_i
                  , p_mild_moderate_soc, p_mild_moderate_trt
                  , p_moderate_severe_soc, p_moderate_severe_trt
                  , p_mild_severe_soc, p_mild_severe_trt
                  , p_MV_ICU_mild_soc, p_MV_ICU_mild_trt
                  , p_MV_ICU_moderate_soc, p_MV_ICU_moderate_trt
                  , p_MV_ICU_severe_soc, p_MV_ICU_severe_trt
                  , p_D_mild, p_ICU_GW_mild, p_GW_REH_mild
                  , p_D_moderate, p_ICU_GW_moderate, p_GW_REH_moderate
                  , p_D_severe, p_ICU_GW_severe, p_GW_REH_severe
                  , p_REH_REC
                  , df_hr_iculos_rtw, df_p_REH_D, df_p_REC_D) { 
  
  # Treatment specific transition probabilities
  if (Trt == "SOC") {
    p_mild_moderate <- p_mild_moderate_soc
  } else if (Trt == "TRT") {
    p_mild_moderate <- p_mild_moderate_trt
  }  else { 
    warning("Invalid treatment type (mild_moderate)") 
  }
  if (Trt == "SOC") {
    p_moderate_severe <- p_moderate_severe_soc
  } else if (Trt == "TRT") {
    p_moderate_severe <- p_moderate_severe_trt
  }  else {
    warning("Invalid treatment type (moderate_severe)") 
  }
  if (Trt == "SOC") {
    p_mild_severe <- p_mild_severe_soc
  } else if (Trt == "TRT") {
    p_mild_severe <- p_mild_severe_trt
  }  else {
    warning("Invalid treatment type (mild_severe)") 
  }
  if (Trt == "SOC") {
    p_MV_ICU_mild <- p_MV_ICU_mild_soc
  } else if (Trt == "TRT") {
    p_MV_ICU_mild <- p_MV_ICU_mild_trt
  } else {
    warning("Invalid treatment type (MV_ICU_mild)") 
  }
  if (Trt == "SOC") {
    p_MV_ICU_moderate <- p_MV_ICU_moderate_soc
  } else if (Trt == "TRT") {
    p_MV_ICU_moderate <- p_MV_ICU_moderate_trt
  } else { 
    warning("Invalid treatment type (MV_ICU_moderate") 
  }
  if (Trt == "SOC") {
    p_MV_ICU_severe <- p_MV_ICU_severe_soc
  } else if (Trt == "TRT") {
    p_MV_ICU_severe <- p_MV_ICU_severe_trt
  } else {
    warning("Invalid treatment type (MV_ICU_severe)") 
  }
  
  # create matrix of state transition probabilities
  m_p_t           <- matrix(0, nrow = n_states, ncol = n_i)  
  # give the state names to the rows
  rownames(m_p_t) <- v_names_states                               
  
  hr_RTW_all <- inner_join(df_X, df_hr_iculos_rtw, by = "n_cycles_ICU")
  hr_RTW     <- hr_RTW_all[M_t == "REH", "hr_RTW"]
  
  p_REH_D_all      <- inner_join(df_X, df_p_REH_D, by = "age_group_id")
  p_REH_D          <- p_REH_D_all[M_t == "REH", "p_REH_D"]
  
  p_REC_D_all      <- inner_join(df_X, df_p_REC_D, by = "Age")
  p_REC_D          <- p_REC_D_all[M_t == "REC", "p_mort_d"]
  
  
  #### mild #### 
  # transition probabilities when hospitalised in ICU on MV
  m_p_t["MV_mild",      M_t == "MV_mild"]  <- (1 - p_D_mild - max(p_MV_ICU_mild[df_X$n_cycles_MV]) - max(p_mild_moderate[df_X$mild_moderate_YN]) - max(p_mild_severe[df_X$mild_severe_YN]))
  m_p_t["ICU_mild",     M_t == "MV_mild"]  <- max(p_MV_ICU_mild[df_X$n_cycles_MV])
  m_p_t["D",            M_t == "MV_mild"]  <- p_D_mild    
  m_p_t["MV_moderate",  M_t == "MV_mild"] <- max(p_mild_moderate[df_X$mild_moderate_YN])
  m_p_t["MV_severe",    M_t == "MV_mild"] <- max(p_mild_severe[df_X$mild_severe_YN])
  
  # transition probabilities when hospitalised in ICU without MV
  m_p_t["ICU_mild",     M_t == "ICU_mild"] <- 1 - p_D_mild - p_ICU_GW_mild - max(p_mild_moderate[df_X$mild_moderate_YN]) - max(p_mild_severe[df_X$mild_severe_YN])
  m_p_t["GW_mild",      M_t == "ICU_mild"] <- p_ICU_GW_mild        
  m_p_t["D",            M_t == "ICU_mild"] <- p_D_mild    
  m_p_t["ICU_moderate", M_t == "ICU_mild"] <- max(p_mild_moderate[df_X$mild_moderate_YN])
  m_p_t["ICU_severe",   M_t == "ICU_mild"] <- max(p_mild_severe[df_X$mild_severe_YN])
  
  # transition probabilities when hospitalised outside ICU
  m_p_t["GW_mild",      M_t == "GW_mild"] <- (1 - p_D_mild - p_GW_REH_mild)
  m_p_t["REH",          M_t == "GW_mild"] <- p_GW_REH_mild
  m_p_t["D",            M_t == "GW_mild"] <- p_D_mild    
  
  
  #### moderate ####
  # transition probabilities when hospitalised in ICU on MV
  m_p_t["MV_moderate",  M_t == "MV_moderate"]  <- (1 - p_D_moderate - max(p_MV_ICU_moderate[df_X$n_cycles_MV]) - max(p_moderate_severe[df_X$moderate_severe_YN]))
  m_p_t["ICU_moderate", M_t == "MV_moderate"]  <- max(p_MV_ICU_moderate[df_X$n_cycles_MV])
  m_p_t["D",            M_t == "MV_moderate"]  <- p_D_moderate 
  m_p_t["MV_severe",    M_t == "MV_moderate"] <- max(p_moderate_severe[df_X$moderate_severe_YN])
  
  # transition probabilities when hospitalised in ICU without MV
  m_p_t["ICU_moderate", M_t == "ICU_moderate"] <- (1 - p_D_moderate - p_ICU_GW_moderate - max(p_moderate_severe[df_X$moderate_severe_YN]))
  m_p_t["GW_moderate",  M_t == "ICU_moderate"] <- p_ICU_GW_moderate 
  m_p_t["D",            M_t == "ICU_moderate"] <- p_D_moderate    
  m_p_t["ICU_severe",   M_t == "ICU_moderate"] <- max(p_moderate_severe[df_X$moderate_severe_YN])
  
  # transition probabilities when hospitalised outside ICU
  m_p_t["GW_moderate",  M_t == "GW_moderate"] <- (1 - p_D_moderate - p_GW_REH_moderate)
  m_p_t["REH",          M_t == "GW_moderate"] <- p_GW_REH_moderate
  m_p_t["D",            M_t == "GW_moderate"] <- p_D_moderate
  
  
  #### severe ####
  # transition probabilities when hospitalised in ICU on MV
  m_p_t["MV_severe",    M_t == "MV_severe"]  <- (1 - p_D_severe - max(p_MV_ICU_severe[df_X$n_cycles_MV]))
  m_p_t["ICU_severe",   M_t == "MV_severe"]  <- max(p_MV_ICU_severe[df_X$n_cycles_MV])
  m_p_t["D",            M_t == "MV_severe"]  <- p_D_severe    
  
  # transition probabilities when hospitalised in ICU without MV
  m_p_t["ICU_severe",   M_t == "ICU_severe"] <- (1 - p_D_severe - p_ICU_GW_severe) 
  m_p_t["GW_severe",    M_t == "ICU_severe"] <- p_ICU_GW_severe 
  m_p_t["D",            M_t == "ICU_severe"] <- p_D_severe    
  
  # transition probabilities when hospitalised outside ICU
  m_p_t["GW_severe",    M_t == "GW_severe"] <- (1 - p_D_severe - p_GW_REH_severe) 
  m_p_t["REH",   M_t == "GW_severe"] <- p_GW_REH_severe
  m_p_t["D",            M_t == "GW_severe"] <- p_D_severe    
  
  
  # transition probabilities when in rehabilitation
  m_p_t["REH",          M_t == "REH"] <- 1 - p_REH_D - (p_REH_REC * hr_RTW)
  m_p_t["REC",          M_t == "REH"] <- p_REH_REC * hr_RTW
  m_p_t["D",            M_t == "REH"] <- p_REH_D
  
  # transition probabilities when recovered
  m_p_t["REC",          M_t == "REC"] <- 1 - p_REC_D
  m_p_t["D",            M_t == "REC"] <- p_REC_D 
  
  # transition probabilities when dead
  m_p_t["D",            M_t == "D"]  <- 1  
  
  return(t(m_p_t))
}

## 1.2 Cost function

## The `Costs` function estimates the costs at every cycle.

Costs <- function (M_t, df_X, Trt, df_c) {
  # Arguments:
  # M_t: health state occupied at cycle t (character variable)
  # Returns: 
  # costs accrued in this cycle
  # Trt:  treatment
  
  # Treatment specific transition costs
  if (Trt == "SOC") {
    c_trt <- 0
  } else if (Trt == "TRT") {
    c_trt <- df_c$value[which(df_c$type == "c_trt")]
  } 
  
  c_t <- c()
  
  c_t[M_t %in% c("MV_mild", "MV_moderate", "MV_severe")] <- if_else(Trt == "TRT" & df_X$trt_YN == 1, df_c$value[which(df_c$type == "c_MV")] + df_c$value[which(df_c$type == "c_trt")], df_c$value[which(df_c$type == "c_MV")])
  c_t[M_t %in% c("ICU_mild", "ICU_moderate", "ICU_severe")] <- if_else(Trt == "TRT" & df_X$trt_YN == 1, df_c$value[which(df_c$type == "c_ICU")] + df_c$value[which(df_c$type == "c_trt")], df_c$value[which(df_c$type == "c_ICU")])
  c_t[M_t %in% c("GW_mild", "GW_moderate", "GW_severe")] <- if_else(Trt == "TRT" & df_X$trt_YN == 1, df_c$value[which(df_c$type == "c_GW")] + df_c$value[which(df_c$type == "c_trt")], df_c$value[which(df_c$type == "c_GW")])
  c_t[M_t == "REH"]  <- df_c$value[which(df_c$type == "c_REH")]
  c_t[M_t == "REC"]  <- df_c$value[which(df_c$type == "c_REC")]
  c_t[M_t == "D"]    <- 0  
  
  return(c_t)  # return costs accrued this cycle
}

## 1.3 Health outcome function

## The `Effs` function to update the utilities at every cycle.

Effs <- function (M_t, df_X, Trt, df_u) {
  # Arguments:
  # M_t: health state occupied at cycle t (character variable)
  # Returns: 
  # QALYs accrued this cycle
  
  u_t <- 0                          # by default the utility for everyone is zero
  
  # Treatment specific utilities
  
  q_t <- 0 
  q_t[M_t %in% c("MV_mild", "MV_moderate", "MV_severe")]    <- df_u$value[which(df_u$type == "u_MV")]
  q_t[M_t %in% c("ICU_mild", "ICU_moderate", "ICU_severe")]    <- df_u$value[which(df_u$type == "u_ICU")] 
  q_t[M_t %in% c("GW_mild", "GW_moderate", "GW_severe")]    <- df_u$value[which(df_u$type == "u_GW")]
  q_t[M_t == "REH"]   <- df_u$value[which(df_u$type == "u_REH")]
  q_t[M_t == "REC"]   <- case_when(df_X$Age < 25 ~ df_u$value[which(df_u$type == "u_REC_1")]
                                   , df_X$Age >= 25 & df_X$Age < 35 ~ df_u$value[which(df_u$type == "u_REC_2")]
                                   , df_X$Age >= 35 & df_X$Age < 45 ~ df_u$value[which(df_u$type == "u_REC_3")]
                                   , df_X$Age >= 45 & df_X$Age < 55 ~ df_u$value[which(df_u$type == "u_REC_4")]
                                   , df_X$Age >= 55 & df_X$Age < 65 ~ df_u$value[which(df_u$type == "u_REC_5")]
                                   , df_X$Age >= 65 & df_X$Age < 75 ~ df_u$value[which(df_u$type == "u_REC_6")]
                                   , df_X$Age >= 75 ~ df_u$value[which(df_u$type == "u_REC_7")])
  q_t[M_t == "D"]     <- 0
  
  return(q_t)  # return the QALYs accrued this cycle
}

## 1.4 The Microsimulation function for the short-term model

MicroSim <- function(n_i, df_X, seed = 234, Trt, df_c, df_u, 
                     p_mild_moderate_soc, p_mild_moderate_trt
                     , p_moderate_severe_soc, p_moderate_severe_trt
                     , p_mild_severe_soc, p_mild_severe_trt
                     , p_MV_ICU_mild_soc, p_MV_ICU_mild_trt
                     , p_MV_ICU_moderate_soc, p_MV_ICU_moderate_trt
                     , p_MV_ICU_severe_soc, p_MV_ICU_severe_trt
                     , p_D_mild, p_ICU_GW_mild, p_GW_REH_mild
                     , p_D_moderate, p_ICU_GW_moderate, p_GW_REH_moderate
                     , p_D_severe, p_ICU_GW_severe, p_GW_REH_severe
                     , p_REH_REC
                     , df_hr_iculos_rtw, df_p_REH_D, df_p_REC_D) {
  
  set.seed(seed) # set a seed to be able to reproduce the same results
  n_states <- length(v_names_states) # the number of health states
  # create three matrices called m_M, m_C and m_E
  # number of rows is equal to the n_i, the number of columns is equal to n_cycles 
  # (the initial state and all the n_cycles cycles)
  # m_M is used to store the health state information over time for every individual
  # m_C is used to store the costs information over time for every individual
  # m_E is used to store the effects information over time for every individual
  
  m_M <- m_C <- m_E <-  matrix(nrow = n_i, ncol = n_cycles + 1, 
                               dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                               paste("cycle", 0:n_cycles, sep = " ")))  
  
  m_M [, 1] <- as.character(df_X$M_init) # initial health state at cycle 0 for individual i
  # calculate costs per individual during cycle 0
  m_C[, 1]  <- Costs(m_M[, 1], df_X, Trt=Trt, df_c)     
  # calculate QALYs per individual during cycle 0
  m_E[, 1]  <- Effs (m_M[, 1], df_X, Trt=Trt, df_u)
  
  # open a loop for time running cycles 1 to n_cycles 
  for (t in 1:n_cycles) {
    # calculate the transition probabilities for the cycle based on health state t
    m_P <- Probs(m_M[, t], df_X, t, Trt = Trt, n_i
                 , p_mild_moderate_soc, p_mild_moderate_trt
                 , p_moderate_severe_soc, p_moderate_severe_trt
                 , p_mild_severe_soc, p_mild_severe_trt
                 , p_MV_ICU_mild_soc, p_MV_ICU_mild_trt
                 , p_MV_ICU_moderate_soc, p_MV_ICU_moderate_trt
                 , p_MV_ICU_severe_soc, p_MV_ICU_severe_trt
                 , p_D_mild, p_ICU_GW_mild, p_GW_REH_mild
                 , p_D_moderate, p_ICU_GW_moderate, p_GW_REH_moderate
                 , p_D_severe, p_ICU_GW_severe, p_GW_REH_severe
                 , p_REH_REC
                 , df_hr_iculos_rtw, df_p_REH_D, df_p_REC_D)
    # check if transition probabilities are between 0 and 1
    check_transition_probability(m_P, verbose = F)
    # check if each of the rows of the transition probabilities matrix sum to one
    check_sum_of_transition_array(m_P, n_rows = n_i, n_cycles = n_cycles, verbose = F)
    
    # sample the next health state and store that state in matrix m_M
    m_M[, t + 1]  <- samplev(m_P, 1)    
    # calculate costs per individual during cycle t + 1
    m_C[, t + 1]  <- Costs(m_M[, t + 1], df_X, Trt = Trt, df_c)  
    # calculate QALYs per individual during cycle t + 1
    m_E[, t + 1]  <- Effs (m_M[, t + 1], df_X, Trt = Trt, df_u)  
    
    # update time since illness onset for t + 1 
    df_X$mild_moderate_YN <- df_X$mild_severe_YN <- if_else(t+1 < 7 & df_X$Severity == "mild", 1, 2)
    df_X$moderate_severe_YN <- if_else(t+1 < 7 & df_X$Severity == "moderate", 1, 2)
    
    df_X$n_cycles_MV <- if_else(grepl("MV", m_M[, t + 1]), df_X$n_cycles_MV +1, df_X$n_cycles_MV)
    
    df_X$n_cycles_ICU <- if_else(m_M[, t + 1] %in% c("MV_mild", "MV_moderate", "MV_severe", "ICU_mild", "ICU_moderate", "ICU_severe"), 
                                 df_X$n_cycles_ICU +1,
                                 df_X$n_cycles_ICU)
    
    df_X$trt_YN <- if_else(m_M[, t + 1] %in% c("MV_mild", "MV_moderate", "MV_severe",
                                               "ICU_mild", "ICU_moderate", "ICU_severe",
                                               "GW_mild", "GW_moderate", "GW_severe")
                           & t+1 < 5
                           , 1
                           , 0)
    
    #df_X   <- inner_join(df_X, df_hr_iculos_rtw, by = "n_cycles_ICU")
    
    # Display simulation progress
    if(t/(n_cycles/10) == round(t/(n_cycles/10), 0)) { # display progress every 10%
      cat('\r', paste(t/n_cycles * 100, "% done", sep = " "))
    }
    
  } # close the loop for the time points 
  
  # store individual characteristics in dataframe for the long-term model and create a column for running age
  df_X   <- df_X %>% mutate(Age_end = Age)
  
  # add hazard ratios of the association between ICU LOS and recovery from rehabilitation
  df_X   <- inner_join(df_X, df_hr_iculos_rtw, by = "n_cycles_ICU")
  
  # calculate outcomes 
  tc      <- rowSums(m_C)   # total cost per individual
  te      <- rowSums(m_E) * cycle_length  # total QALYs per individual 
  tc_hat  <- mean(tc)       # average cost 
  te_hat  <- mean(te)       # average QALY  
  # store the results from the simulation in a list
  results <- list(df_X = df_X, m_M = m_M, m_C = m_C, m_E = m_E, tc = tc , te = te, 
                  tc_hat = tc_hat, te_hat = te_hat)   
  
  return(results)  # return the results
}

# 2 Define Simulation Functions for the long-term model

## 2.1 Probability function

Probs_lt <- function(M_t_lt, df_X_lt, t, Trt, p_LTS_REC_lt, df_LTS_D_lt, df_p_REC_D_lt, n_i, n_states_lt, v_names_states_lt) {
  
  # create matrix of state transition probabilities
  m_p_t_lt           <- matrix(0, nrow = n_states_lt, ncol = n_i)  
  # give the state names to the rows
  rownames(m_p_t_lt) <- v_names_states_lt   
  
  hr_RTW_lt      <- df_X_lt[M_t_lt == "LTS", "hr_RTW"]
  
  p_LTS_D_lt_all <- inner_join(df_X_lt, df_LTS_D_lt, by = "age_group_id")
  p_LTS_D_lt     <- p_LTS_D_lt_all[M_t_lt == "LTS", "p_LTS_D_lt"]
  
  p_REC_D_lt_all <- inner_join(df_X_lt, df_p_REC_D_lt, by = "Age_end")
  p_REC_D_lt     <- p_REC_D_lt_all[M_t_lt == "REC", "p_mort_m"]
  
  # transition probabilities
  m_p_t_lt["REC",  M_t_lt == "LTS"]   <- p_LTS_REC_lt * hr_RTW_lt
  m_p_t_lt["LTS",  M_t_lt == "LTS"]   <- 1 - p_LTS_D_lt - (p_LTS_REC_lt * hr_RTW_lt)
  m_p_t_lt["D",    M_t_lt == "LTS"]   <- p_LTS_D_lt
  
  m_p_t_lt["REC",  M_t_lt == "REC"]   <- (1 - p_REC_D_lt)
  m_p_t_lt["D",    M_t_lt == "REC"]   <- p_REC_D_lt
  #p_REH_D[df_X_lt$Age_end]
  
  m_p_t_lt["D",    M_t_lt == "D"]     <- 1
  
  return(t(m_p_t_lt))
}

## 2.2 Costs function
Costs_lt <- function (M_t_lt, df_X_lt, Trt, c_LTS, c_REC, c_D) {
  
  c_t_lt <- c()
  
  c_t_lt[M_t_lt == "LTS"] <- c_LTS
  c_t_lt[M_t_lt == "REC"] <- c_REC 
  c_t_lt[M_t_lt == "D"]   <- c_D
  
  return(c_t_lt)  # return costs accrued this cycle
}

## 2.3 Effects function
Effs_lt <- function (M_t_lt, df_X_lt, Trt, df_u) {
  # Arguments:
  # M_t: health state occupied at cycle t (character variable)
  
  u_t_lt <- 0                          # by default the utility for everyone is zero
  
  q_t_lt <- 0 
  q_t_lt[M_t_lt == "LTS"]   <- df_u$value[which(df_u$type == "u_LTS")]
  q_t_lt[M_t_lt == "REC"]   <- case_when(df_X_lt$Age_end < 25 ~ df_u$value[which(df_u$type == "u_REC_1")]
                                         , df_X_lt$Age_end >= 25 & df_X_lt$Age_end < 35 ~ df_u$value[which(df_u$type == "u_REC_2")]
                                         , df_X_lt$Age_end >= 35 & df_X_lt$Age_end < 45 ~ df_u$value[which(df_u$type == "u_REC_3")]
                                         , df_X_lt$Age_end >= 45 & df_X_lt$Age_end < 55 ~ df_u$value[which(df_u$type == "u_REC_4")]
                                         , df_X_lt$Age_end >= 55 & df_X_lt$Age_end < 65 ~ df_u$value[which(df_u$type == "u_REC_5")]
                                         , df_X_lt$Age_end >= 65 & df_X_lt$Age_end < 75 ~ df_u$value[which(df_u$type == "u_REC_6")]
                                         , df_X_lt$Age_end >= 75 ~ df_u$value[which(df_u$type == "u_REC_7")])
  q_t_lt[M_t_lt == "D"]     <- 0
  
  return(q_t_lt)  # return the QALYs accrued this cycle
}

## 2.4 Microsimulation function for the long-term model

MicroSim_lt <- function(n_i, v_names_states_lt, n_cycles_lt, cycle_length_lt, df_X_lt, df_m_M_lt, seed = 234, Trt, p_LTS_REC_lt, df_LTS_D_lt, df_p_REC_D_lt,
                        c_LTS, c_REC, c_D,
                        df_u,
                        d_c, d_e) {
  
  #if (Trt == "SOC") {
  #  df_X_lt <- df_X_soc
  #} else if (Trt == "TRT") {
  #  df_X_lt <- df_X_trt
  #}
  #if (Trt == "SOC") {
  #  df_m_M_lt <- df_m_M_soc
  #} else if (Trt == "TRT") {
  #  df_m_M_lt <- df_m_M_trt
  #}
  
  set.seed(seed) # set a seed to be able to reproduce the same results
  n_states_lt <- length(v_names_states_lt) # the number of health states
  
  m_M_lt <- m_C_lt <- m_E_lt <- matrix(nrow = n_i, ncol = n_cycles_lt + 1, 
                                       dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                                       paste("month", 0:n_cycles_lt, sep = " ")))
  
  m_M_lt [, 1] <- df_m_M_lt[, n_cycles + 1] # initial health state at cycle 0 for individual i
  m_M_lt [, 1] <- if_else(m_M_lt[, 1] %in% c("MV_mild", "MV_moderate", "MV_severe", "ICU_mild", "ICU_moderate", "ICU_severe", "GW_mild", "GW_moderate", "GW_severe", "REH"),
                          "LTS", m_M_lt[, 1])
  
  m_C_lt[, 1]  <- Costs_lt(m_M_lt[, 1], df_X_lt, Trt=Trt, c_LTS, c_REC, c_D)     
  
  m_E_lt[, 1]  <- Effs_lt(m_M_lt[, 1], df_X_lt, Trt=Trt, df_u)
  
  # open a loop for time running cycles 1 to n_cycles 
  for (t in 1:n_cycles_lt) {
    # calculate the transition probabilities for the cycle based on health state t
    m_P_lt <- Probs_lt(m_M_lt[, t], df_X_lt, t, Trt = Trt, p_LTS_REC_lt, df_LTS_D_lt, df_p_REC_D_lt, n_i, n_states_lt, v_names_states_lt)
    # check if transition probabilities are between 0 and 1
    check_transition_probability(m_P_lt, verbose = F)
    # check if each of the rows of the transition probabilities matrix sum to one
    check_sum_of_transition_array(m_P_lt, n_rows = n_i, n_cycles = n_cycles_lt, verbose = F)
    
    # sample the next health state and store that state in matrix m_M
    m_M_lt[, t + 1]  <- samplev(m_P_lt, 1) 
    # calculate costs per individual during cycle t + 1
    m_C_lt[, t + 1]  <- Costs_lt(m_M_lt[, t + 1], df_X_lt, Trt = Trt, c_LTS, c_REC, c_D)  
    # calculate QALYs per individual during cycle t + 1
    m_E_lt[, t + 1]  <- Effs_lt(m_M_lt[, t + 1], df_X_lt, Trt = Trt, df_u)  
    
    # update age column if alive
    df_X_lt$Age_end <- if_else(!m_M_lt[, t + 1] == "D", df_X_lt$Age + as.integer((t+1)*cycle_length_lt), df_X_lt$Age_end)
    df_X_lt$age_group_id <- if_else(df_X_lt$Age_end < 65, 1, 2)
    
    # Display simulation progress
    if(t/(n_cycles_lt/10) == round(t/(n_cycles_lt/10), 0)) { # display progress every 10%
      cat('\r', paste(t/n_cycles_lt * 100, "% done", sep = " "))
    }
  }  
  # calculate 
  ### Discount weight for costs and effects 
  v_dwc_lt   <- 1 / ((1 + (d_c * cycle_length_lt)) ^ (0:n_cycles_lt))
  v_dwe_lt   <- 1 / ((1 + (d_e * cycle_length_lt)) ^ (0:n_cycles_lt))
  
  tc_lt      <- m_C_lt %*% v_dwc_lt  # total (discounted) cost per individual
  te_lt      <- (m_E_lt * cycle_length_lt) %*% v_dwe_lt  # total (discounted) QALYs per individual 
  tc_lt_hat  <- mean(tc_lt)       # average (discounted) cost 
  te_lt_hat  <- mean(te_lt)       # average (discounted) QALY  
  
  results_lt <- list(df_X_lt = df_X_lt
                     , m_P_lt = m_P_lt
                     , m_M_lt = m_M_lt
                     , m_C_lt = m_C_lt
                     , m_E_lt = m_E_lt
                     , tc_lt = tc_lt 
                     , te_lt = te_lt
                     , tc_lt_hat = tc_lt_hat
                     , te_lt_hat = te_lt_hat)
  return(results_lt)
}


#------------------------------------------------------------------------------#
####             Adapted version of plot_trace_microsim                 ####
#------------------------------------------------------------------------------#
#' Generate parameter sets for the probabilistic sensitivity analysis (PSA)
#'
#' \code{plot_trace_microsim_lt} generates a trace plot by merging subgroup states 
#' @export
#' 

plot_trace_microsim_lt <- function(m_M) {
  # plot the distribution of the population across health states over time (trace)
  # count the number of individuals in each health state at each cycle
  m_TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_names_states_lt, ordered = TRUE))))
  # m_TR <- m_TR / n_i                                 # calculate the proportion of individuals
  m_TR <- m_TR / nrow(m_M)
  colnames(m_TR) <- v_names_states_lt                   # name the rows of the matrix
  rownames(m_TR) <- paste("Cycle", 0:(ncol(m_M)-1), sep = " ") # name the columns of the matrix
  # Plot trace of first health state
  plot(0:(ncol(m_M)-1), m_TR[, 1], type = "l", main = "Health state trace long-term model",
       ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "month")
  # add a line for each additional state
  for (n_states_lt in 2:length(v_names_states_lt)) {
    lines(0:(ncol(m_M)-1), m_TR[, n_states_lt], col = n_states_lt)   # adds a line to current plot
  }
  legend("topright", v_names_states_lt, col = 1:length(v_names_states_lt), # add a legend to current plot
         lty = rep(1, length(v_names_states_lt)), bty = "n", cex = 0.65)
  
}



#' Decision Model
#'
#' \code{decision_model} implements the decision model used.
#'
#' @param l_params_all List with all parameters of decision model
#' @param verbose Logical variable to indicate print out of messages
#' @return The transition probability array and the cohort trace matrix.
#' @export
decision_model <- function(l_params_all, verbose = FALSE) {
  with(as.list(l_params_all), {
    ########################### Process model inputs ###########################
    country_analysis <- 'Germany'
    ## Cycle length
    cycle_length   <- 1/365
    ## Number of cycles
    n_cycles       <- 60 # time horizon, number of cycles
    # Define health states of the hospitalisation model:
    v_names_states <- c("MV_mild",
                        "ICU_mild",       
                        "GW_mild",        
                        "MV_moderate",
                        "ICU_moderate",   
                        "GW_moderate",    
                        "MV_severe",
                        "ICU_severe",      
                        "GW_severe",
                        "REH",
                        "REC",
                        "D"
    )       
    n_states       <- length(v_names_states) # number of health states 
    ### Strategies 
    v_names_str   <- c("SOC",   # store the strategy names
                       "TRT") 
    n_str         <- length(v_names_str)     # number of strategies
    
    # total number of patients
    n_i            <- 1000
    # initial number of patients per disease severity
    n_i_mild       <- round(n_i * r_i_mild,0)
    n_i_moderate   <- round(n_i * r_i_moderate,0)
    n_i_severe     <- n_i - n_i_mild - n_i_moderate
    
    #n_i_mild       <- round(n_i_dsa * r_i_mild,0)
    #n_i_moderate   <- round(n_i_dsa * r_i_moderate,0)
    #n_i_severe     <- n_i_dsa - round(n_i_dsa * r_i_mild,0) - round(n_i_dsa * r_i_moderate,0)
    #n_i            <- round(n_i_dsa * r_i_mild,0) + round(n_i_dsa * r_i_moderate,0) + (n_i_dsa - round(n_i_dsa * r_i_mild,0) - round(n_i_dsa * r_i_moderate,0))
    
    r_i_ICU <- 1-r_i_MV # 0.144 
    n_i_mild_MV  <- round(n_i_mild*r_i_MV,0)
    n_i_mild_ICU <- n_i_mild - round(n_i_mild*r_i_MV,0)
    n_i_moderate_MV  <- round(n_i_moderate*r_i_MV,0)
    n_i_moderate_ICU <- n_i_moderate - round(n_i_moderate*r_i_MV,0)
    n_i_severe_MV  <- round(n_i_severe*r_i_MV,0)
    n_i_severe_ICU <- n_i_severe - round(n_i_severe*r_i_MV,0)
    
    #r_i_ICU          <- 1 - r_i_MV 
    #n_i_mild_MV      <- round(round(n_i_dsa * r_i_mild,0)*r_i_MV,0)
    #n_i_mild_ICU     <- round(n_i_dsa * r_i_mild,0) - round(round(n_i_dsa * r_i_mild,0)*r_i_MV,0) #round(n_i_mild*r_i_ICU,0)
    #n_i_moderate_MV  <- round(round(n_i_dsa * r_i_moderate,0)*r_i_MV,0)
    #n_i_moderate_ICU <- round(n_i_dsa * r_i_moderate,0) - round(round(n_i_dsa * r_i_moderate,0)*r_i_MV,0)
    #n_i_severe_MV    <- round((n_i_dsa - round(n_i_dsa * r_i_mild,0) - round(n_i_dsa * r_i_moderate,0))*r_i_MV,0)
    #n_i_severe_ICU   <- round(n_i_dsa - round(n_i_dsa * r_i_mild,0) - round(n_i_dsa * r_i_moderate,0)) - round((n_i_dsa - round(n_i_dsa * r_i_mild,0) - round(n_i_dsa * r_i_moderate,0))*r_i_MV,0)
    
    df_p_REC_D <- p_mort %>% filter(country == country_analysis)
    
    #df_hr_iculos_rtw <- rbind(df_hr_iculos_rtw, data.frame(n_cycles_ICU = rep(max(df_hr_iculos_rtw$n_cycles_ICU+1):(n_cycles+1)), hr_RTW = df_hr_iculos_rtw$hr_RTW[which(df_hr_iculos_rtw$n_cycles_ICU == 30)]))
    
    p_mild_moderate_soc     <- 1-((1-r_mild_moderate)^(1/d_progression_period))
    p_moderate_severe_soc   <- 1-((1-r_moderate_severe)^(1/d_progression_period))
    p_mild_severe_soc       <- 1-((1-r_mild_severe)^(1/d_progression_period))
    
    p_mild_moderate_soc     <- c(p_mild_moderate_soc, 0)
    p_moderate_severe_soc   <- c(p_moderate_severe_soc, 0)
    p_mild_severe_soc       <- c(p_mild_severe_soc, 0)
    
    r_mild_moderate_trt     <- r_mild_moderate * hr_mild_moderate
    r_moderate_severe_trt   <- r_moderate_severe *  hr_moderate_severe
    r_mild_severe_trt       <- r_mild_severe * hr_mild_severe
    
    p_mild_moderate_trt     <- 1-((1-r_mild_moderate_trt)^(1/d_progression_period))
    p_moderate_severe_trt   <- 1-((1-r_moderate_severe_trt)^(1/d_progression_period))
    p_mild_severe_trt       <- 1-((1-r_mild_severe_trt)^(1/d_progression_period))
    
    p_mild_moderate_trt     <- c(p_mild_moderate_trt, 0)
    p_moderate_severe_trt   <- c(p_moderate_severe_trt, 0)
    p_mild_severe_trt       <- c(p_mild_severe_trt, 0)
    
    p_MV_ICU_mild_soc         <- 1-exp(-df_p_mv_icu$probability[which(df_p_mv_icu$severity_group == 'mild')])
    p_MV_ICU_moderate_soc     <- 1-exp(-df_p_mv_icu$probability[which(df_p_mv_icu$severity_group == 'moderate')])
    p_MV_ICU_severe_soc       <- 1-exp(-df_p_mv_icu$probability[which(df_p_mv_icu$severity_group == 'severe')])
    
    r_MV_ICU_mild_trt         <- df_p_mv_icu$probability[which(df_p_mv_icu$severity_group == 'mild')]*hr_MV_ICU_mild
    r_MV_ICU_moderate_trt     <- df_p_mv_icu$probability[which(df_p_mv_icu$severity_group == 'moderate')]*hr_MV_ICU_moderate
    r_MV_ICU_severe_trt       <- df_p_mv_icu$probability[which(df_p_mv_icu$severity_group == 'severe')]*hr_MV_ICU_severe
    
    p_MV_ICU_mild_trt         <- 1-exp(-r_MV_ICU_mild_trt)
    p_MV_ICU_moderate_trt     <- 1-exp(-r_MV_ICU_moderate_trt)
    p_MV_ICU_severe_trt       <- 1-exp(-r_MV_ICU_severe_trt)
    
    ### MV LOS
    median_MV_LOS_mild      <- 6
    median_MV_LOS_moderate  <- 8
    median_MV_LOS_severe    <- 11
    
    # ICU LOS
    median_ICU_LOS_mild     <- 10
    median_ICU_LOS_moderate <- 12
    median_ICU_LOS_severe   <- 14
    median_ICU_LOS_mild_net     <- median_ICU_LOS_mild - median_MV_LOS_mild
    median_ICU_LOS_moderate_net <- median_ICU_LOS_moderate - median_MV_LOS_moderate
    median_ICU_LOS_severe_net   <- median_ICU_LOS_severe - median_MV_LOS_severe
    
    r_ICU_mild              <- -log(0.5)/median_ICU_LOS_mild_net
    r_ICU_moderate          <- -log(0.5)/median_ICU_LOS_moderate_net
    r_ICU_severe            <- -log(0.5)/median_ICU_LOS_severe_net
    
    p_ICU_GW_mild           <- 1 - exp(-r_ICU_mild)
    p_ICU_GW_moderate       <- 1 - exp(-r_ICU_moderate)
    p_ICU_GW_severe         <- 1 - exp(-r_ICU_severe)
    
    # GW LOS
    median_GW_LOS_mild      <- 23
    median_GW_LOS_moderate  <- 22
    median_GW_LOS_severe    <- 26
    median_GW_LOS_mild_net      <- median_GW_LOS_mild - median_ICU_LOS_mild
    median_GW_LOS_moderate_net  <- median_GW_LOS_moderate - median_ICU_LOS_moderate
    median_GW_LOS_severe_net    <- median_GW_LOS_severe - median_ICU_LOS_severe
    
    r_GW_mild               <- -log(0.5)/median_GW_LOS_mild_net
    r_GW_moderate           <- -log(0.5)/median_GW_LOS_moderate_net
    r_GW_severe             <- -log(0.5)/median_GW_LOS_severe_net
    
    p_GW_REH_mild           <- 1 - exp(-r_GW_mild)
    p_GW_REH_moderate       <- 1 - exp(-r_GW_moderate)
    p_GW_REH_severe         <- 1 - exp(-r_GW_severe)
    
    p_D_mild                 = 1 - (1 - r_mort_28_mild)^(1 / 28)
    p_D_moderate             = 1 - (1 - r_mort_28_moderate)^(1 / 28)
    p_D_severe               = 1 - (1 - r_mort_28_severe)^(1 / 28)
    
    LTS_RTW_median_d <- 203
    p_REH_REC        <- 1-(exp(-(-log(0.5)/LTS_RTW_median_d)))
    
    p_REH_D_st_1     <- 1-((1-0.042)^(1/(2*365))) #related to p_REH_REC and specific for age <65
    p_REH_D_st_2     <- 1-(1-0.356)^(1/(2*30)) #difference in mortality between day 30 and day 365, not related to p_REH_REC and specific for age >=65 
    df_p_REH_D       <- data.frame(age_group_id = c(1,2), p_REH_D = c(p_REH_D_st_1, p_REH_D_st_2))
    
    c_LTS     <- c_REH_LTS * n_rehab_visits_LTS
    df_c      <- data.frame(type = c("c_MV", "c_ICU", "c_GW", "c_REH", "c_LTS"), country = country_analysis, value = c(c_MV, c_ICU, c_GW, c_REH, c_LTS))
    df_c_add  <- data.frame(type = c("c_REC", "c_D", "c_trt"), country = country_analysis, value = c(c_REC, c_D, c_trt))
    df_c      <- rbind(df_c, df_c_add)
    
    df_u    <- data.frame(type = c("u_MV", "u_ICU", "u_GW", "u_REH", "u_LTS", "u_REC_1", "u_REC_2", "u_REC_3", "u_REC_4", "u_REC_5", "u_REC_6", "u_REC_7"), value = c(u_MV, u_ICU, u_GW, u_REH, u_LTS, u_REC_1, u_REC_2, u_REC_3, u_REC_4, u_REC_5, u_REC_6, u_REC_7))
    
    ## sample from age distribution an initial age for every individual
    set.seed(234)
    v_age_init <- rnorm(n_i, mean = mean_age, sd = sd_age)
    df_age <- data.frame(age = round(v_age_init, 0))
    df_age <- df_age %>% mutate(age = if_else(age<min_age, min_age, if_else(age>max_age, max_age, age)))
    df_age <- df_age %>% group_by(age) %>% summarize(prop = n()/n_i)
    v_age_init  <- sample(x = df_age$age, prob = df_age$prop, size = n_i, replace = TRUE) 
    
    v_M_init_mild <- rep(c("MV_mild", "ICU_mild"), times = c(n_i_mild_MV, n_i_mild_ICU))
    df_X_mild <- data.frame(ID = 1:n_i_mild, Severity = "mild", M_init = v_M_init_mild, n_cycles_MV = if_else(v_M_init_mild == "MV_mild", 1, 0), n_cycles_ICU = 1, trt_YN = 1, mild_moderate_YN = 1, mild_severe_YN = 1, moderate_severe_YN = 2)
    v_M_init_moderate <- rep(c("MV_moderate", "ICU_moderate"), times = c(n_i_moderate_MV, n_i_moderate_ICU))
    df_X_moderate <- data.frame(ID = (n_i_mild+1):(n_i_mild+n_i_moderate), Severity = "moderate", M_init = v_M_init_moderate, n_cycles_MV = if_else(v_M_init_moderate == "MV_moderate", 1, 0), n_cycles_ICU = 1, trt_YN = 1, mild_moderate_YN = 2, mild_severe_YN = 2, moderate_severe_YN = 1)
    v_M_init_severe <- rep(c("MV_severe", "ICU_severe"), times = c(n_i_severe_MV, n_i_severe_ICU))
    df_X_severe <- data.frame(ID = (n_i_mild+n_i_moderate+1):n_i, Severity = "severe", M_init = v_M_init_severe, n_cycles_MV = if_else(v_M_init_severe == "MV_severe", 1, 0), n_cycles_ICU = 1, trt_YN = 1, mild_moderate_YN = 2, mild_severe_YN = 2, moderate_severe_YN = 2)
    
    df_X <- rbind(df_X_mild, df_X_moderate, df_X_severe)
    df_X <- df_X %>% mutate(Age = v_age_init,
                            age_group_id = if_else(Age < 65, 1, 2))
    
    outcomes_SoC <- MicroSim(n_i = n_i, df_X = df_X, seed = 234, Trt="SOC", df_c = df_c, df_u = df_u, 
                             p_mild_moderate_soc, p_mild_moderate_trt
                             , p_moderate_severe_soc, p_moderate_severe_trt
                             , p_mild_severe_soc, p_mild_severe_trt
                             , p_MV_ICU_mild_soc, p_MV_ICU_mild_trt
                             , p_MV_ICU_moderate_soc, p_MV_ICU_moderate_trt
                             , p_MV_ICU_severe_soc, p_MV_ICU_severe_trt
                             , p_D_mild, p_ICU_GW_mild, p_GW_REH_mild
                             , p_D_moderate, p_ICU_GW_moderate, p_GW_REH_moderate
                             , p_D_severe, p_ICU_GW_severe, p_GW_REH_severe
                             , p_REH_REC
                             , df_hr_iculos_rtw = df_hr_iculos_rtw
                             , df_p_REH_D = df_p_REH_D
                             , df_p_REC_D = df_p_REC_D) 
    outcomes_trt <- MicroSim(n_i = n_i, df_X = df_X, seed = 234, Trt="TRT", df_c = df_c, df_u = df_u,
                             p_mild_moderate_soc, p_mild_moderate_trt
                             , p_moderate_severe_soc, p_moderate_severe_trt
                             , p_mild_severe_soc, p_mild_severe_trt
                             , p_MV_ICU_mild_soc, p_MV_ICU_mild_trt
                             , p_MV_ICU_moderate_soc, p_MV_ICU_moderate_trt
                             , p_MV_ICU_severe_soc, p_MV_ICU_severe_trt
                             , p_D_mild, p_ICU_GW_mild, p_GW_REH_mild
                             , p_D_moderate, p_ICU_GW_moderate, p_GW_REH_moderate
                             , p_D_severe, p_ICU_GW_severe, p_GW_REH_severe
                             , p_REH_REC
                             , df_hr_iculos_rtw = df_hr_iculos_rtw
                             , df_p_REH_D = df_p_REH_D
                             , df_p_REC_D = df_p_REC_D)
    
    #df_m_M_soc <- as.data.frame(outcomes_SoC$m_M)
    #df_m_M_trt <- as.data.frame(outcomes_trt$m_M)
    
    ### long term model
    #n_age_max         <- 120
    n_cycles_lt       <- round((n_age_max-mean(df_X$Age))*12, 0)
    cycle_length_lt   <- 1/12     # cycle length equal to one year (use 1/12 for monthly)
    v_names_states_lt <- c("LTS", "REC", "D")
    
    #c_D   <- 0
    
    #df_X_soc   <- outcomes_SoC$df_X
    #df_X_trt   <- outcomes_trt$df_X
    #df_X_soc   <- df_X_soc %>% mutate(Age_end = Age)
    #df_X_trt   <- df_X_trt %>% mutate(Age_end = Age)
    df_p_REC_D_lt <- df_p_REC_D %>% mutate(Age_end = Age)
    #df_X_soc   <- inner_join(df_X_soc, df_hr_iculos_rtw, by = "n_cycles_ICU")
    #df_X_trt   <- inner_join(df_X_trt, df_hr_iculos_rtw, by = "n_cycles_ICU")
    
    LTS_RTW_median_d <- 203
    p_LTS_REC_lt     <- 1-(exp(-(-log(0.5)/LTS_RTW_median_d)))^30
    
    p_LTS_D_lt_1     <- 1-((1-0.042)^(1/24)) #related to p_LTS_REC and specific for age <65
    p_LTS_D_lt_2     <- 1-((1-0.077)^(1/10)) #not related to p_LTS_REC and specific for age >=65 
    df_LTS_D_lt <- data.frame(age_group_id = c(1,2), p_LTS_D_lt = c(p_LTS_D_lt_1, p_LTS_D_lt_2))
    
    outcomes_SoC_lt   <- MicroSim_lt(n_i = n_i, v_names_states_lt, n_cycles_lt, cycle_length_lt, df_X_lt = outcomes_SoC$df_X, df_m_M_lt = as.data.frame(outcomes_SoC$m_M), seed = 234, Trt="SOC", 
                                     p_LTS_REC_lt, df_LTS_D_lt, df_p_REC_D_lt,
                                     c_LTS, c_REC, c_D,
                                     df_u,
                                     d_c, d_e)
    outcomes_trt_lt   <- MicroSim_lt(n_i = n_i, v_names_states_lt, n_cycles_lt, cycle_length_lt, df_X_lt = outcomes_trt$df_X, df_m_M_lt = as.data.frame(outcomes_trt$m_M), seed = 234, Trt="TRT", 
                                     p_LTS_REC_lt, df_LTS_D_lt, df_p_REC_D_lt,
                                     c_LTS, c_REC, c_D,
                                     df_u,
                                     d_c, d_e)
    
    df_SoC_tc_st <- as.data.frame(outcomes_SoC$tc)
    df_SoC_tc_lt <- as.data.frame(outcomes_SoC_lt$tc_lt)
    df_SoC_tc_st_lt <- cbind(df_SoC_tc_st, df_SoC_tc_lt) 
    df_SoC_tc_st_lt <- df_SoC_tc_st_lt %>% mutate(total_cost = rowSums(df_SoC_tc_st_lt)) %>% select(total_cost)
    
    df_SoC_te_st <- as.data.frame(outcomes_SoC$te)
    df_SoC_te_lt <- as.data.frame(outcomes_SoC_lt$te_lt)
    df_SoC_te_st_lt <- cbind(df_SoC_te_st, df_SoC_te_lt) 
    df_SoC_te_st_lt <- df_SoC_te_st_lt %>% mutate(total_qaly = rowSums(df_SoC_te_st_lt)) %>% select(total_qaly)
    
    df_trt_tc_st <- as.data.frame(outcomes_trt$tc)
    df_trt_tc_lt <- as.data.frame(outcomes_trt_lt$tc_lt)
    df_trt_tc_st_lt <- cbind(df_trt_tc_st, df_trt_tc_lt) 
    df_trt_tc_st_lt <- df_trt_tc_st_lt %>% mutate(total_cost = rowSums(df_trt_tc_st_lt)) %>% select(total_cost)
    
    df_trt_te_st <- as.data.frame(outcomes_trt$te)
    df_trt_te_lt <- as.data.frame(outcomes_trt_lt$te_lt)
    df_trt_te_st_lt <- cbind(df_trt_te_st, df_trt_te_lt) 
    df_trt_te_st_lt <- df_trt_te_st_lt %>% mutate(total_qaly = rowSums(df_trt_te_st_lt)) %>% select(total_qaly)
    
    # store the mean costs of each strategy in a new variable C (vector of costs)
    v_C <- c(mean(df_SoC_tc_st_lt$total_cost), mean(df_trt_tc_st_lt$total_cost))
    # store the mean QALYs of each strategy in a new variable E (vector of effects)
    v_E <- c(mean(df_SoC_te_st_lt$total_qaly), mean(df_trt_te_st_lt$total_qaly))
    
    #store the ICER in a new variable icer
    v_icer <- (mean(df_trt_tc_st_lt$total_cost) - mean(df_SoC_tc_st_lt$total_cost))/(mean(df_trt_te_st_lt$total_qaly) - mean(df_SoC_te_st_lt$total_qaly))
    
    ## Vector with discounted net monetary benefits (NMB)
    #n_wtp <- 80000
    v_nmb <- v_E * n_wtp - v_C
    
    # use dampack to calculate the ICER
    #cea <- calculate_icers(cost       = v_C,
    #                       effect     = v_E,
    #                       strategies = v_names_str)
    #df_cea <- as.data.frame(cea)
    
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_C,
                        Effect   = v_E,
                        ICER     = v_icer,
                        NMB      = v_nmb
    )
    
    return(df_ce)
    
  }
  )
}