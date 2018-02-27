######################################################################
# Functions for manipulating simulation output
######################################################################

# this is tend... I presume the final time?
time_end <- function(results){
  ncol(results$S)
}

# vector of deaths... What does each entry corresponds to? I presume split by
#  coutry, risk group, age group?
deaths <- function(results){
  tend <- time_end(results)
  X - results$S[,tend] - results$SV[,tend] - results$E[,tend] - results$EV[,tend] -
    results$I[,tend] - results$IV[,tend] - results$R[,tend] - results$RV[,tend]
}

# total number of deaths
worldwide_deaths <- function(results){
  sum(deaths(results))
}

# global attack rate
global_attack <- function(results){
  pop_total <- sum(X)
  tend <- time_end(results)
  sum(results$R[,tend] + results$RV[,tend] + deaths(results))/pop_total
}
