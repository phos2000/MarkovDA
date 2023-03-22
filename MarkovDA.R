## ----setup, include=FALSE----------------------------------------
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
library(tidyverse)
library(knitr)
library(kableExtra)
library(DiagrammeR)
library(here)
library(glue)
library(Matrix)
library(expm)
library(ggrepel)    # For plotting
library(ellipse)    # For plotting
library(scales)     # For dollar signs and commas
library(randtoolbox)
library(igraph)
library(tidygraph)
library(ggraph)
library(progress)
library(hrbrthemes)
library(ggsci)
library(directlabels)
library(dampack)
library(tinytex)
source(here('exercise/functions_markov.r'))


## ----------------------------------------------------------------
knitr::kable(
  dplyr::tibble(
    type = c("Probability","Rate","Matrix","Cost","Utility","Hazard Ratio"),
    prefix = c("p_","r_","m_","c_","u_","hr_")
  )
)


## ----input-excel-------------------------------------------------
input_file = normalizePath(here("exercise/markov_exercise.xlsx"))
input_raw = readxl::read_xlsx(input_file,sheet="Parameters")[,1:7]
params_sc = input_raw %>% 
  select(`Variable name`, Value) %>%
  na.omit() %>%
  deframe() %>%
  as.list()
# list2env(params_sc, envir=.GlobalEnv)

v_names_states = str_sub(names(as_tibble(params_sc) %>% select(starts_with("u"))), 2, -1)
v_n_states <- length(v_names_states)
v_names_str = c("quo", str_c("t",str_sub(names(as_tibble(params_sc) %>% select(starts_with("cTrt"))), 3, -1)))
n_strategies <- length(v_names_str)

options("scipen"=1000, "digits"=2)

t(as_tibble(params_sc)) %>%
  kable(caption = "input params from excel") %>%
  kable_styling() %>%
  scroll_box(height = "400px")


## ----build-prob-matrix-------------------------------------------
build_matrices <- function(params, mtype) {
  m_X = list()
  for (m in v_names_str){
    m_X_ = matrix(0, nrow = v_n_states, ncol = v_n_states, dimnames = list(v_names_states,v_names_states))
    for (i in rownames(m_X_)) {
      for (j in colnames(m_X_)) {
        if (paste0(mtype,i,j) %in% names(params)) {
          m_X_[i,j] = as.numeric(params[paste0(mtype,i,j)])
        }
      }
    }
    
    if (mtype == "p") {
      ## if drawing the probabilities from the given parameters and building a probability matrix
      
      if (paste0("pRemissionStage1_",m) %in% names(params)) {
        ## the relapse probs change over the strategies. A special rule here. 
        m_X_["Remission","Stage1"] = as.numeric(params[paste0("pRemissionStage1_",m)])
        m_X_["Remission","Stage2"] = as.numeric(params[paste0("pRemissionStage2_",m)])
      }
        
      for (n in 1:v_n_states) {
        m_X_[n,n] = 1-sum(m_X_[n,])
      }
    } else {
      ## if drawing the rates from the given paramters and building a rate matrix
      
      if (paste0("hrRemissionStage1_",m) %in% names(params)) {
        ## use the hazard ratios to adjust the relapse rates over the strategies. 
        m_X_["Remission","Stage1"] = m_X_["Remission","Stage1"] * as.numeric(params[paste0("hrRemissionStage1_",m)])
        m_X_["Remission","Stage2"] = m_X_["Remission","Stage2"] * as.numeric(params[paste0("hrRemissionStage2_",m)])
      }
        
      for (n in 1:v_n_states) {
        m_X_[n,n] = -sum(m_X_[n,])
      }
    }
    m_X[[m]] = m_X_
  }
  m_X
}
    
m_P = build_matrices(params_sc, "p")

m_P


## ----formula-----------------------------------------------------
m_P_form = list()
m_P_ = matrix(0, nrow = v_n_states, ncol = v_n_states, dimnames = list(v_names_states,v_names_states))

for (i in v_names_states) {
  tempsum = sum(as_tibble(params_sc) %>% select(starts_with(paste0("r",i))))
  for (j in v_names_states) {
    if (paste0("r",i,j) %in% names(params_sc)) {
      ptemp = params_sc[[paste0("r",i,j)]] * (1 - exp(-tempsum)) / tempsum
    } else if (j == i) {
      ptemp = exp(-tempsum)
    } else {
      ptemp = 0
    }
    m_P_[i,j] = ptemp
  }
}
m_P_form[[v_names_str[1]]] = m_P_

for (str in v_names_str[2:4]) {
  params_sc[[paste0("rRemissionStage1_",str)]] = params_sc$rRemissionStage1 * params_sc[[paste0("hrRemissionStage1_",str)]]
  params_sc[[paste0("rRemissionStage2_",str)]] = params_sc$rRemissionStage2 * params_sc[[paste0("hrRemissionStage2_",str)]]
  m_P_form[[str]] = m_P_
  tempsum = params_sc[[paste0("rRemissionStage1_",str)]] + params_sc[[paste0("rRemissionStage2_",str)]]
  m_P_form[[str]]["Remission", "Stage1"] = params_sc[[paste0("rRemissionStage1_",str)]] * (1 - exp(-tempsum)) / tempsum
  m_P_form[[str]]["Remission", "Stage2"] = params_sc[[paste0("rRemissionStage2_",str)]] * (1 - exp(-tempsum)) / tempsum
  m_P_form[[str]]["Remission", "Remission"] = exp(-tempsum)
}

m_P_form


## ----expm--------------------------------------------------------
m_R = build_matrices(params_sc, "r")
m_P_expm = list()
for (m in v_names_str){
  m_P_expm[[m]] = expm(m_R[[m]])
}
m_P_expm


## ----payoffs-----------------------------------------------------
build_payoffs = function(params){
  payoffs = list()
  cost_basic = c()
  cost_basic["Healthy"] = 0
  for (i in v_names_states[2:(v_n_states-1)]) {
    cost_basic[i] = as.numeric(params[paste0("cTreatment",i)]) + as.numeric(params[paste0("nVisits",i)]) * params$cVisit
  }
  cost_basic["Death"] = 0

  payoffs[["quo"]] = matrix(c(cost_basic, as_vector(as_tibble(params) %>% select(starts_with("u")))),
                      nrow = 2, ncol = 6, byrow = TRUE, dimnames = list(c("costs", "qalys"), v_names_states))
  for (m in v_names_str[2:n_strategies]) {
    payoffs[[m]] = payoffs[["quo"]]
    payoffs[[m]]["costs","Remission"] = payoffs[[m]]["costs","Remission"] + as.numeric(params[paste0("cT",str_sub(m,2,-1))]) - as.numeric(params[paste0("cTreatment","Remission")])
  }
  payoffs
}

payoffs = build_payoffs(params_sc)
payoffs


## ----trace-------------------------------------------------------
build_traces = function(params){
  m_P = build_matrices(params, "p")
  
  l_m_M <-
    m_P %>%
    map(~({
      tmp <- .x[1,] %>% as.matrix() %>% t()
      tmp <- matrix(0,ncol=(ncol(tmp)),nrow=(params$n_cycles+1),dimnames = list(paste0(0:params$n_cycles),colnames(.x)))
      tmp[1,1] <- params$nPop
      tmp
    }))

  for (t in 1:params$n_cycles) {
    res <-
      map2(l_m_M,m_P,~({
        .x[paste0(t-1),] %*% .y
      }))
    l_m_M <-
      map2(l_m_M,res,~({
        .x[paste0(t),] <- .y
        .x
      }))
  }
  l_m_M
}

l_m_M = build_traces(params_sc)
l_m_M


## ----discount----------------------------------------------------
### Discount weight for costs and effects ----

get_ce = function(params){
  v_dwc  <- 1 / ((1 + (params$discount_rate * params$cycle_length)) ^ (0:params$n_cycles))
  v_dwe  <- 1 / ((1 + (params$discount_rate * params$cycle_length)) ^ (0:params$n_cycles))
  ## Within-cycle correction (WCC) using Simpson's 1/3 rule ----
  # Function included in "Functions_markov.R"
  v_wcc <- gen_wcc(n_cycles = params$n_cycles,
                 method = "Simpson1/3") # vector of wcc
  
  payoffs = build_payoffs(params)
  l_m_M = build_traces(params)
  
  v_tot_qaly <-
    map2(l_m_M, payoffs,~({
      v_u_str <- .y["qalys",v_names_states] %>% as.vector()
      t(.x[,v_names_states] %*% v_u_str) %*% (v_dwe * v_wcc)
    })) %>%
    unlist()

  v_tot_cost <-
    map2(l_m_M, payoffs,~({
      v_c_str <- .y["costs",v_names_states] %>% as.vector()
      t(.x[,v_names_states] %*% v_c_str) %*% (v_dwc * v_wcc)
    })) %>%
    unlist()
  
  ce=list()
  ce[["discounted_weight_of_cost"]] = v_dwc
  ce[["discounted_weight_of_qaly"]] = v_dwe
  ce[["Within_cycle_correction"]] = v_dwc
  ce[["cost"]] = v_tot_cost
  ce[["qaly"]] = v_tot_qaly
  
  ce
}

tot_ce = get_ce(params_sc)
tot_ce


## ----cea---------------------------------------------------------
## Compare different strategies with icers
# Function also included in "Functions_markov.R"
calculate_icers <- function(cost, effect, strategies) {
  # checks on input
  n_cost <- length(cost)
  n_eff <- length(effect)
  n_strat <- length(strategies)
  if (n_cost != n_eff | n_eff != n_strat) {
    stop("cost, effect, and strategies must all be vectors of the same length", call. = FALSE)
  }

  # coerce to character, in case they are provided as numeric
  char_strat <- as.character(strategies)

  # create data frame to hold data
  df <- data.frame("Strategy" = char_strat,
                   "Cost" = cost,
                   "Effect" = effect,
                   stringsAsFactors = FALSE)
  nstrat <- nrow(df)

  # if only one strategy was provided, return df with NAs for incremental
  if (nstrat == 1) {
    df[, c("ICER", "Inc_Cost", "Inc_Effect")] <- NA
    return(df)
  }

  # three statuses: dominated, extended dominated, and non-dominated
  d <- NULL

  # detect dominated strategies
  # dominated strategies have a higher cost and lower effect
  df <- df %>%
    arrange(.data$Cost, desc(.data$Effect))

  # iterate over strategies and detect (strongly) dominated strategies
  # those with higher cost and equal or lower effect
  for (i in 1:(nstrat - 1)) {
    ith_effect <- df[i, "Effect"]
    for (j in (i + 1):nstrat) {
      jth_effect <- df[j, "Effect"]
      if (jth_effect <= ith_effect) {
        # append dominated strategies to vector
        d <- c(d, df[j, "Strategy"])
      }
    }
  }

  # detect weakly dominated strategies (extended dominance)
  # this needs to be repeated until there are no more ED strategies
  ed <- vector()
  continue <- TRUE  # ensure that the loop is run at least once
  while (continue) {
    # vector of all dominated strategies (strong or weak)
    dom <- union(d, ed)

    # strategies declared to be non-dominated at this point
    nd <- setdiff(strategies, dom)

    # compute icers for nd strategies
    nd_df <- df[df$Strategy %in% nd, ] %>%
      compute_icers()

    # number non-d
    n_non_d <- nrow(nd_df)

    # if only two strategies left, we're done
    if (n_non_d <= 2) {
      break
    }

    # strategy identifiers for non-d
    nd_strat <- nd_df$Strategy

    # now, go through non-d strategies and detect any
    # with higher ICER than following strategy
    ## keep track of whether any ED strategies are picked up
    # if not, we're done - exit the loop
    new_ed <- 0
    for (i in 2:(n_non_d - 1)) {
      if (nd_df[i, "ICER"] > nd_df[i + 1, "ICER"]) {
        ed <- c(ed, nd_strat[i])
        new_ed <- new_ed + 1
      }
    }
    if (new_ed == 0) {
      continue <- FALSE
    }
  }

  # recompute icers without weakly dominated strategies
  nd_df_icers <- nd_df[!(nd_df$Strategy %in% dom), ] %>%
    mutate(Status = "ND") %>%
    compute_icers()

  # dominated and weakly dominated
  d_df <- df[df$Strategy %in% d, ] %>%
    mutate(ICER = NA, Status = "D")

  ed_df <- df[df$Strategy %in% ed, ] %>%
    mutate(ICER = NA, Status = "ED")

  # when combining, sort so we have ref,ND,ED,D
  results <- bind_rows(d_df, ed_df, nd_df_icers) %>%
    arrange(desc(.data$Status), .data$Cost, desc(.data$Effect))

  # re-arrange columns
  results <- results %>%
    select(.data$Strategy, .data$Cost, .data$Effect,
           .data$Inc_Cost, .data$Inc_Effect, .data$ICER, .data$Status)

  # declare class of results
  class(results) <- c("icers", "data.frame")
  return(results)
}

df_cea <- calculate_icers(cost       = tot_ce$cost,
                          effect     = tot_ce$qaly,
                          strategies = v_names_str)
df_cea


## ----comments----------------------------------------------------
knitr::kable(
  dplyr::tibble(
    Status = c("D","ND","ED"),
    Descriptiion = c("Dominated","Not Dominated","Extended Dominated")
  )
)


## ----cea-output--------------------------------------------------
## display a table for CEA outputs
# Function included in "Functions_markov.R"
table_cea <- format_table_cea(df_cea) 
table_cea

## CEA frontier -----
# depends on the `ggplot2`  and `ggrepel` packages.
# Function included in "Functions_markov.R"
plot(df_cea, label = "all", txtsize = 16) +
  expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.8, 0.2))


## ----owdsa-------------------------------------------------------
# set lower & upper case for cTrtA
# params_lower = params_upper  = params_sc
# params_lower$cTrtA = 8000
# params_upper$cTrtA = 10500
# 
# result8000 = calculate_icers(cost       = get_ce(params_lower)$cost,
#                 effect     = get_ce(params_lower)$qaly,
#                 strategies = v_names_str)
# 
# result10500 = calculate_icers(cost       = get_ce(params_upper)$cost,
#                 effect     = get_ce(params_upper)$qaly,
#                 strategies = v_names_str)
# 
# # set lower & upper case for pStage3Death
# params_lower = params_upper  = params_sc
# params_lower$pStage3Death = 0.44
# params_upper$pStage3Death = 0.5
# 
# calculate_icers(cost       = get_ce(params_lower)$cost,
#                 effect     = get_ce(params_lower)$qaly,
#                 strategies = v_names_str)
# 
# calculate_icers(cost       = get_ce(params_upper)$cost,
#                 effect     = get_ce(params_upper)$qaly,
#                 strategies = v_names_str)


#dampack_dsa

simulate_strategies = function(params, wtp = params_sc$wtp) {
  # Create cost-effectiveness results data frame
  # LY is the total dicounted life-years, here will not use. 
  df_ce <- data.frame(Strategy = v_names_str,
                      Cost = numeric(n_strategies),
                      QALY = numeric(n_strategies),
                      stringsAsFactors = FALSE)

  result = calculate_icers(cost = get_ce(params)$cost,
                           effect = get_ce(params)$qaly,
                           strategies = v_names_str)
  
  df_ce[c("Cost", "QALY","ICER")] <- c(result$Cost,
                                result$Effect,
                                result$ICER)
  df_ce["NMB"] <- result$Effect * wtp - result$Cost
  
  return(df_ce)
}

my_owsa_params_range <- data.frame(pars = c("cTrtA", "pStage3Death"),
                                   min = c(8000, 0.3),
                                   max = c(10500, 0.9))

l_owsa_det = run_owsa_det(params_range = my_owsa_params_range,
                           params_basecase = params_sc,
                           nsamp = 2,
                           FUN = simulate_strategies,
                           outcomes = c("Cost", "QALY", "ICER", "NMB"),
                           strategies = v_names_str,
                           progress = FALSE)

owsaNMBtrtA = l_owsa_det$owsa_NMB[l_owsa_det$owsa_NMB[,"strategy"] == "trtA",]
owsaICERtrtA = l_owsa_det$owsa_ICER[l_owsa_det$owsa_ICER[,"strategy"] == "trtA",]
owsa_tornado(owsaNMBtrtA)
owsa_tornado(owsaICERtrtA)

# # Select the net monetary benefit (NMB) owsa object
# my_owsa_NMB <- l_owsa_det$owsa_NMB
# 
# # Plot outcome of each strategy over each parameter range
# plot(my_owsa_NMB,
#      n_x_ticks = 4)
# 
# owsa_tornado(l_owsa_det$owsa_NMB)
# 
# owsa_opt_strat(l_owsa_det$owsa_NMB)

