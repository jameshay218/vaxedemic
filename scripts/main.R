library(reshape2)
library(ggplot2)
library(Matrix)
library(data.table)

wd <- "~/Documents/vaxedemic" 
devtools::load_all(wd)

## LIFE HISTORY PARAMETER INPUTS
## R_0, recovery time and latent period
life_history_params <- list(R0=1.8, TR=2.6, LP = 1.5, case_fatality_ratio = c(1e-3,1e-2))

## travel parameters: scaling of off-diagonals
travel_params <- list(epsilon = 1e-3)

## vaccine efficacy and initial vaccinated proportion
# this example roughly brings effective R to 1.2
vax_params <- list(efficacy = 1 - 1.2/1.8, propn_vax0 = 0)

## example parameters for vaccine production (see cum_vax_pool_func_closure)
vax_production_params <- list(detection_delay = 0, production_delay = 0, 
                              production_rate = 0, max_vax = 5e9)

## example parameters for vaccine allocation
vax_allocation_params <- list(priorities = NULL)
