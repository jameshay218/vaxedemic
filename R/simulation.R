#' function to set up simulation from inputs, then run
#' 
#' @param simulation_flags named vector (or list) with the logical elements "normaliseTravel"
#' and "seasonal"
#' @param life_history_params named vector (or list) with the numeric elements
#' "R0", "TR" (time to recovery) and "LP" (latent period)
#' @param vax_params named vector (or list) with the numeric element "efficacy"
#' @param seasonality_params list of seasonality parameters.
#' contains the elements tdelay (0 <= tdelay <= 364): shifts the seasonality function - changing this effectively changes the seed time.
#' tdelay = 0 is seed at t = 0 in sinusoidal curve, roughly start of autumn in Northern hemisphere
#' days_per_block: average seasonality over blocks of this many days
#' amp: amplitude of seasonality
#' @param time_params list of parameters to do with time steps in simulation.
# contains the elements tmax (Maximum time of simulation), tdiv (Number of time steps per day)
#' @param vax_alloc_period allocate vaccines once every this many timesteps
#' @param processed_inputs list. output of setup_inputs function
#' @param n_runs number of repeats to run
#' @param calculate_summaries_func string specifying what function to use to calculate summaries.
#' @param other_info list which provides information needed to calculate the summaries
#' @return a list of lists containing results of main simulation loop (see main_simulation). The length of this list is equal to n_runs
#' @import foreach
#' @importFrom Matrix Matrix
#' @export
run_simulation <- function(simulation_flags, life_history_params, vax_params, seasonality_params,
                           time_params, vax_alloc_period, processed_inputs,
                           n_runs=1,
                           calculate_summaries_func="return_all_res", other_info){
  
  #travelMatrix <- diag(n_countries) ##DH debug - decouples countries, keeping seed
  normaliseTravel <- simulation_flags[["normaliseTravel"]]
  seasonal <- simulation_flags[["seasonal"]]
  
  R0 <- life_history_params[["R0"]]
  TR <- life_history_params[["TR"]] # time to recovery
  LP <- life_history_params[["LP"]] # latent period
  efficacy <- vax_params[["efficacy"]]
  gamma <- 1/TR # recovery rate
  sigma <- 1/LP # rate of exposed -> infectious
  tdelay <- seasonality_params[["tdelay"]] #Delay from peak summer in northern hemisphere######## ##DH
  amp <- seasonality_params[["amp"]] #Amplitude of seasonality########
  tmax <- time_params[["tmax"]]
  tdiv <- time_params[["tdiv"]]
  
  n_countries <- processed_inputs[["n_countries"]]
  n_ages <- processed_inputs[["n_ages"]]
  n_riskgroups <- processed_inputs[["n_riskgroups"]]
  if(seasonal) {
    n_divisions_in_block <- seasonality_params[["days_per_block"]] * tdiv
    n_blocks <- tmax / seasonality_params[["days_per_block"]]
    if(round(n_blocks) != n_blocks) {
      stop("time_params$tmax must be a multiple of seasonality_params$days_per_block")
    }
    if(round(n_divisions_in_block) != n_divisions_in_block){
      stop("time_params$tdiv must be a multiple of seasonality_params$days_per_block")
    }
  }
  travelMatrix <- processed_inputs[["travelMatrix"]]
  contactMatrix <- processed_inputs[["contactMatrix"]]
  X <- processed_inputs[["popns"]]
  labels <- processed_inputs[["labels"]]
  cum_vax_pool_func <- processed_inputs[["cum_vax_pool_func"]]
  vax_allocation_func <- processed_inputs[["vax_allocation_func"]]
  seed_vec <- processed_inputs[["seed_vec"]]
  
  maxIndex <- n_countries*n_ages*n_riskgroups
  groupsPerLoc <- n_ages*n_riskgroups
  
  #### calculation of force of infection matrix (LD) starts here
  
  ##Seasonality:
  Beta1 <- function(lats,n_countries){
    trop <- 23.5 #Tropics (in degrees)
    Y=matrix(1,n_countries,1) #Flat outside of tropics - could modify (DH/SR)
    
    in_tropics <- abs(lats) < trop # find countries in tropics
    Y[in_tropics] <- lats[in_tropics]/trop # linear with abs(latitude) in tropics
    Y[lats < (-trop)] <- -1 # opposite sign in southern hemisphere
    return(Y)
  }
  
  ## calculate Phi, the amplifying factor due to seasonality at each simulation time,
  ## for each location, age, risk group
  if(seasonal){
    beta1 <- Beta1(processed_inputs[["latitudes"]],n_countries)
    B <- kronecker(matrix(1,groupsPerLoc,n_blocks),beta1)
    timevec <- seq(1/tdiv,tmax,by=1/tdiv)
    wave <- sin((timevec-tdelay)*2*pi/365)
    wave <- t(colMeans(matrix(wave, nrow=n_divisions_in_block)))
    Phi <- 1+amp*B*kronecker(matrix(1,maxIndex,1),wave)#One column of B per time step
  } else {
    Phi <- 1
  }
  ## row normalise travel matrix
  if(normaliseTravel){
    Krow <- rowSums(travelMatrix)
    Knorm <- kronecker(matrix(1,1,n_countries),matrix(Krow,n_countries,1))
    travelMatrix <- travelMatrix/Knorm
  }
  Kdelta <- kronecker(diag(groupsPerLoc),travelMatrix)
  K1 <- kronecker(matrix(1,groupsPerLoc,groupsPerLoc),t(travelMatrix))
  
  ## combine contact and travel matrices
  if(is.list(contactMatrix)) { # if age-risk mixing is country-specific
    # row concatenation of age/risk matrices by country
    Cvec <- do.call(rbind, contactMatrix)
    # next three lines: getting the row order right
    idx_vec <- matrix(seq_len(nrow(Cvec)), n_ages * n_riskgroups, n_countries)
    idx_vec <- matrix(t(idx_vec), nrow(Cvec), 1)
    Cvec <- Cvec[idx_vec,]
    
    CT <- t(Cvec)
    Cnew <- matrix(CT,groupsPerLoc,groupsPerLoc*n_countries)
    Chome <- kronecker(t(Cnew),matrix(1,1,n_countries))#Mix as though were home
    KC <- kronecker(matrix(1,n_ages * n_riskgroups,n_ages * n_riskgroups),travelMatrix)*Chome
  } else {
    KC <- kronecker(contactMatrix,t(travelMatrix))
  }
  
  ## Denominator of force of infection term
  M <- K1%*%X
  M1 <- 1/M
  M1[M==0] <- 0
  
  M1 <- kronecker(matrix(1,maxIndex,1),t(M1))
  
  Xover=1/X
  Xover[X==0] <- 0
  
  #Country-specific betas: ##DH added
  
  n_classes <- n_ages * n_riskgroups
  c_names <- levels(labels[,2])
  c_pops <- labels[order(labels$Location),]
  c_pops <- c_pops[,1]
  c_pops <- matrix(c_pops, n_classes, n_countries, byrow=FALSE)
  betaV <- matrix(0, n_countries, 1)
  for (i in 1:n_countries){
    NiClasses <- c_pops[,i]#Pop in each class
    Ni <- sum(NiClasses)#Country population
    if(is.list(contactMatrix)) {
      Ci <- contactMatrix[[i]]#Country age-mix - =C for all?
    } else {
      Ci <- contactMatrix
    }
    Li <- 1/(Ni*gamma)*rep(NiClasses,1,n_classes)*Ci##DH travelMatrix[i,i]/...
    ev <- eigen(Li)
    ev <- ev$values
    Rstar <- max(abs(ev))
    betaV[i] <- R0/Rstar
  }
  betaV <- kronecker(matrix(1,n_classes,1), betaV)
  
  beta <- betaV
  
  ## Convert to smaller time step rates
  beta <- beta/tdiv
  gamma <- gamma/tdiv
  
  LD <- (Kdelta*M1)%*%KC ##DH: no beta
  
  #### calculation of force of infection matrix (LD) ends here
  
  ## set initial conditions
  
  # intial exposed are distributed among vaccinated and unvaccinated proportionally
  E <- seed_vec
  S <- X - seed_vec
  
  SV <- EV <- I <- R <- IV <- RV <- matrix(0, maxIndex)
  
  ## Pre-compute seaonal contribution to FOI
  if(seasonal){
    LD <- foreach(i=1:ncol(Phi)) %do% {
      M1Phi <- M1*kronecker(matrix(1,n_ages * n_riskgroups *n_countries,1),t(Phi[,i]))
      as.matrix(Matrix(Kdelta*M1Phi)%*%KC)
    }
  }
  
  ## gather model parameters
  modelParameters <- list("gamma"=gamma, "sigma" = sigma, "efficacy" = efficacy,
                          "beta" = beta, "M1" = M1, "Kdelta" = Kdelta, "KC"= KC,
                          "Phi" = Phi, "seasonal" = seasonal, 
                          "case_fatality_ratio" = processed_inputs[["case_fatality_ratio_vec"]])
  
  # put arguments of other_info into global environment, so that they can be passed
  # to calculate_summaries_func.  Warns if we're overwriting anything.
  
  list2here(other_info)
  calculate_summaries_arg_names <- formalArgs(calculate_summaries_func)
  calculate_summaries_arg_names <- calculate_summaries_arg_names[calculate_summaries_arg_names != "res"]
  calculate_summaries_args <- list_vars_from_environment(calculate_summaries_arg_names)
  # run simulation
  run_parallel <- TRUE
  if(run_parallel) {
    result <- foreach(i = 1:n_runs) %dopar% {
      res <- main_simulation(tmax,tdiv, vax_alloc_period, LD, S, E, I, R,
                             SV, EV, IV, RV, modelParameters, cum_vax_pool_func,
                             vax_allocation_func)
      res <- c(list(res = res), calculate_summaries_args)
      res <- do.call(calculate_summaries_func, res)
      res
    }
  } else {
    # series version for debugging
    result <- list(n_runs)
    for(i in 1:n_runs) {
      res <- main_simulation(tmax,tdiv, vax_alloc_period, LD, S, E, I, R,
                             SV, EV, IV, RV, modelParameters, cum_vax_pool_func,
                             vax_allocation_func)
      # put arguments of other_info into global environment, so that they can be passed
      # to calculate_summaries_func.  Warns if we're overwriting anything.
      list2here(other_info)
      calculate_summaries_args <- list_vars_from_environment(formalArgs(calculate_summaries_func))
      res <- do.call(calculate_summaries_func, calculate_summaries_args)
      result[[i]] <- res
    }
  }
  result
}


#' main simulation loop
#' 
#' @param tmax number of days for which to run simulation
#' @param tdiv number of timesteps per day
#' @param vax_alloc_period allocate vaccines every this number of timesteps
#' @param LD force of infection matrix: multiply by the I vector to get force of infection.
#' a square matrix of side length maxIndex  = n_countries * n_ages * n_riskgroups
#' @param S0 numeric vector of length max_Index: initial number of unvaccinated
#' susceptibles in each location, age, risk group
#' @param E0 numeric vector of length max_Index: initial number of unvaccinated
#' exposed
#' @param I0 numeric vector of length max_Index: initial number of unvaccinated
#' infectious
#' @param R0 numeric vector of length max_Index: initial number of unvaccinated
#' removed
#' @param SV0 numeric vector of length max_Index: initial number of vaccinated
#' susceptibles
#' @param EV0 numeric vector of length max_Index: initial number of vaccinated
#' exposed
#' @param IV0 numeric vector of length max_Index: initial number of vaccinated
#' infectious
#' @param RV0 numeric vector of length max_Index: initial number of vaccinated
#' removed
#' @param params list of model parameters
#' @param cum_vax_pool_func a function of time which gives the numebr of vaccines
#' ever produced at that time
#' @param vax_allocation_func a function of the current state of the simulation which
#' gives the number of vaccines to allocate to each location, age, risk group
#' @return list containing the elements:
#' beta: infectivity parameter
#' incidence: matrix containing the incidence in each
#' location, age, risk group at each timestep.
#' nrow(S) = n_groups = n_countres * n_ages * n_riskgroups
#' ncol(S) = tend = tmax * tdiv
#' deaths: same as incidence but for deaths
main_simulation <- function(tmax, tdiv, vax_alloc_period, LD, S0, E0, I0, R0, 
                            SV0, EV0, IV0, RV0, params,
                            cum_vax_pool_func, vax_allocation_func){
  ## extract model parameters  
  gamma <- params[["gamma"]]
  sigma <- params[["sigma"]]
  case_fatality_ratio <- params[["case_fatality_ratio"]]
  efficacy <- params[["efficacy"]]
  beta <- params[["beta"]]
  M1 <- params[["M1"]]
  Kdelta <- params[["Kdelta"]]
  KC <- Matrix(params[["KC"]])
  Phi <- params[["Phi"]]
  seasonal <- params[["seasonal"]]
  
  ## initialisevectors for current state of simulation
  S <- S0
  E <- E0
  I <- I0
  R <- R0
  SP <- SV0
  EP <- EV0
  IP <- IV0
  RP <- RV0
  incidence <- SV2 <- EV2 <- IV2 <- RV2 <- SV1 <- EV1 <- IV1 <- RV1 <- 0 * S0
  
  n_groups <- length(I0)
  
  ## calculate number of time steps
  tend <- tmax*tdiv
  ## make vector of simulation times
  times <- seq(0,tmax,by=1/tdiv)
  
  ## If we have seasonality, find the resolution
  ## of the seasonality vector. We re-calculate seasonal impact on FOI
  ## every n_blocks iterations
  switch_freq <- tend+1
  index <- 1
  LD1 <- LD
  if(is.list(LD)){
    n_blocks <- length(LD)
    switch_freq <- tend/n_blocks
    LD1 <- LD[[index]]
  }
  
  ## initialise matrices to store simulation outputs
  incidence_mat <-  deaths_mat <- matrix(0, n_groups, length(times))
  
  ## initialise vector to store number of vaccines over time
  vax_pool_vec <- double(length(times))
  ## calculate the number of vaccines ever produced at each timestep
  cum_vax_pool <- vapply(times, cum_vax_pool_func, double(1))
  
  if(any(diff(cum_vax_pool) < 0)) {
    stop("cumulative number of vaccines not a monotonically non-decreasing function")
  }
  
  for(i in 2:(tend+1)){      
    ## check that current state is sensible
    stopifnot(all(S >= 0),all(E >= 0), all(I >= 0), all(R >= 0), 
              all(SV1 >= 0), all(EV1 >= 0), all(IV1 >= 0), all(RV1 >= 0),
              all(SV2 >= 0),all(EV2 >= 0), all(IV2 >= 0), all(RV2 >= 0), 
              all(SP >= 0), all(EP >= 0), all(IP >= 0), all(RP >= 0))
    
    if(i %% switch_freq == 0 && i < tend){
      index <- index + 1
      LD1 <- LD[[index]]
    }
    
    #################
    ## VAX PRODUCTION
    #################
    ## current number of vaccines
    if(i == 2) {
      vax_pool <- cum_vax_pool[i] # if allocating vaccines for the first time
    } else {
      vax_pool <- vax_pool + cum_vax_pool[i] - cum_vax_pool[i - 1]
    }
    
    #################
    ## VAX ALLOCATION
    #################       
    if(i %% vax_alloc_period == 2) { # allocate at start of simulation, and
      # every vax_alloc_period timesteps thereafter
    
    # move individuals vaccinated two weeks ago to protected class
    SP <- SP + SV2
    EP <- EP + EV2
    IP <- IP + IV2
    RP <- RP + RV2

    # move individuals vaccinated last week to intermediate vaccination class
    SV2 <- SV1
    EV2 <- EV1
    IV2 <- IV1
    RV2 <- RV1

    SV1 <- EV1 <- IV1 <- RV1 <- 0

      # initialisation
      vax_alloc <- 1
      actual_alloc <- 0
      sum_vax_alloc <- S * 0
      alloc_this_round <- Inf
      
      # notes on while loop to follow
      while(alloc_this_round > 0 && vax_pool >= 1) {
        
        # first pass at allocating vaccines.
        # may allocate more vaccines to a location / age / risk group/ infection status
        # combination than there are individuals in that combination
        vax_alloc <- vax_allocation_func(S, E, I, R, 
                                         SV1, EV1, IV1, RV1, incidence,
                                         vax_pool)
        # the actual number of vaccines allocated is the smaller of the
        # number of vaccines according to the algorithm and the actual
        # number of individuals in the combination
        # if the acutal number of vaccines allocated is less than that
        # allocated by the algorithm, we will try to allocate the 
        # remaining vaccines in the next iteration of the while loop
        actual_alloc <- Map(pmin, vax_alloc, list(S, E, I, R))
        sum_vax_alloc <- sum_vax_alloc + actual_alloc$S + actual_alloc$E +
          actual_alloc$I + actual_alloc$R
        
        ## update current vax pool
        alloc_this_round <- sum(unlist(actual_alloc))
        # if (alloc_this_round == 0 && vax_pool >= 1) {
        #   break
        # }
        
        vax_pool <- vax_pool - alloc_this_round
        ## update the vaccination status of individuals
        S <- S - actual_alloc$S
        E <- E - actual_alloc$E
        I <- I - actual_alloc$I
        R <- R - actual_alloc$R
        
        SV1 <- SV1 + actual_alloc$S
        EV1 <- EV1 + actual_alloc$E
        IV1 <- IV1 + actual_alloc$I
        RV1 <- RV1 + actual_alloc$R
      }
    }
    
    # numeric errors can be introduced, so make sure numbers of individuals
    # are still integers before the transmission process
    
    S <- round(S)
    E <- round(E)
    I <- round(I)
    R <- round(R)
    SV1 <- round(SV1)
    EV1 <- round(EV1)
    IV1 <- round(IV1)
    RV1 <- round(RV1)
    
    #################
    ## INFECTIONS
    #################
    ## Generate force of infection on each group/location
    lambda <- beta*LD1%*%(I + IV1 + IV2 + IP)
    
    ## Generate probability of infection from this
    P_infection <- 1 - exp(-lambda)
    P_infection_vax <- 1 - exp(-lambda*(1 - efficacy))
    
    ## Simulate new infections for each location
    newInfections <- rbinom(n_groups, S, P_infection)
    newInfectionsV1 <- rbinom(n_groups, SV1, P_infection)
    newInfectionsV2 <- rbinom(n_groups, SV2, P_infection)
    newInfectionsP <- rbinom(n_groups, SP, P_infection_vax)
    
    ## Update populations
    S <- S - newInfections
    E <- E + newInfections
    SV1 <- SV1 - newInfectionsV1
    EV1 <- EV1 + newInfectionsV1
    SV2 <- SV2 - newInfectionsV2
    EV2 <- EV2 + newInfectionsV2
    SP <- SP - newInfectionsP
    EP <- EP + newInfectionsP
    
    #################
    ## EXPOSED BECOMING INFECTIOUS
    #################
    ## Generate probability of exposed becoming infectious
    P_infectious <- 1 - exp(-sigma)
    
    ## Simulate new becoming infectious
    newInfectious <- rbinom(n_groups, E, P_infectious)
    newInfectiousV1 <- rbinom(n_groups, EV1, P_infectious)
    newInfectiousV2 <- rbinom(n_groups, EV2, P_infectious)
    newInfectiousP <- rbinom(n_groups, EP, P_infectious)
    
    ## Update populations
    E <- E - newInfectious
    I <- I + newInfectious
    EV1 <- EV1 - newInfectiousV1
    IV1 <- IV1 + newInfectiousV1
    EV2 <- EV2 - newInfectiousV2
    IV2 <- IV2 + newInfectiousV2
    EP <- EP - newInfectiousP
    IP <- IP + newInfectiousP
    incidence <- newInfectious + newInfectiousV1 +
        newInfectiousV2 + newInfectiousP
    
    #################
    ## RECOVERIES AND DEATHS
    #################
    ## Generate probability of removal
    P_removal <- 1 - exp(-gamma)
    
    ## Simulate new recoveries
    newRemovals <- rbinom(n_groups, I, P_removal)
    newRemovalsV1 <- rbinom(n_groups, IV1, P_removal)
    newRemovalsV2 <- rbinom(n_groups, IV2, P_removal)
    newRemovalsP <- rbinom(n_groups, IP, P_removal)
    
    newRecovered <- rbinom(n_groups, newRemovals, 1 - case_fatality_ratio)
    newRecoveredV1 <- rbinom(n_groups, newRemovalsV1, 1 - case_fatality_ratio)
    newRecoveredV2 <- rbinom(n_groups, newRemovalsV2, 1 - case_fatality_ratio)
    newRecoveredP <- rbinom(n_groups, newRemovalsP, 1 - case_fatality_ratio)
    ## Update populations
    I <- I - newRemovals
    R <- R + newRecovered
    IV1 <- IV1 - newRemovalsV1
    RV1 <- RV1 + newRecoveredV1
    IV2 <- IV2 - newRemovalsV2
    RV2 <- RV2 + newRecoveredV2
    IP <- IP - newRemovalsP
    RP <- RP + newRecoveredP
    
    #################
    ## SAVE RESULTS
    #################
    incidence_mat[,i] <- incidence
    deaths_mat[,i] <- newRemovals + newRemovalsV1 + 
        newRemovalsV2 + newRemovalsP -
        newRecovered - newRecoveredV1 - 
        newRecoveredV2 - newRecoveredP
    
    ## stop simulation if there are no more exposed/infectious individuals
    
    if(i < (tend + 1) && sum(E + I + EV1 + IV1 + EV2 + IV2 + EP + IP) == 0) {
      break
    }
  }
  
  ## put times as column names for readability of output
  colnames(incidence_mat) <- colnames(deaths_mat) <- times
  
  return(list(incidence = incidence_mat, deaths = deaths_mat))
}
