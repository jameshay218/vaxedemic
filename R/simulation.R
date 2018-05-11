#' function to set up simulation from inputs, then run
#' 
#' @param simulation_flags named vector (or list) with the logical elements "normaliseTravel"
#' and "seasonal"
#' @param life_history_params named vector (or list) with the numeric elements
#' "R0", "TR" (time to recovery) and "LP" (latent period)
#' @param vax_params named vector (or list) with the numeric elements "efficacy"
#' and "propn_vax0" (initial proportion of vaccinated individuals; assumed constant
#' across location, age and risk groups)
#' @param seasonality_params list of seasonality parameters.
#' contains the elements tdelay (0 <= tdelay <= 364): shifts the seasonality function - changing this effectively changes the seed time.
#' tdelay = 0 is seed at t = 0 in sinusoidal curve, roughly start of autumn in Northern hemisphere
#' division: Average seasonality into this many blocks of time
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
    propn_vax0 <- vax_params[["propn_vax0"]] # initial proportion of vaccinated individuals
    gamma <- 1/TR # recovery rate
    sigma <- 1/LP # rate of exposed -> infectious
    tdelay <- seasonality_params[["tdelay"]] #Delay from peak summer in northern hemisphere######## ##DH
    amp <- seasonality_params[["amp"]] #Amplitude of seasonality########
    tmax <- time_params[["tmax"]]
    tdiv <- time_params[["tdiv"]]

    n_countries <- processed_inputs[["n_countries"]]
    n_ages <- processed_inputs[["n_ages"]]
    n_riskgroups <- processed_inputs[["n_riskgroups"]]
    season_res <- time_params[["tmax"]]*time_params[["tdiv"]]/seasonality_params[["division"]]
    
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
      B <- kronecker(matrix(1,groupsPerLoc,tmax*tdiv/season_res),beta1)
      timevec <- seq(1/tdiv,tmax,by=1/tdiv)
      wave <- sin((timevec-tdelay)*2*pi/365)
      wave <- t(colMeans(matrix(wave, nrow=season_res)))
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
    EV <- round(seed_vec * propn_vax0)
    E <- seed_vec - EV
    SV <- round((X - seed_vec) * propn_vax0) ##DH
    S <- X - seed_vec - SV

    I <- R <- IV <- RV <- matrix(0, maxIndex)

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
    result <- foreach(i = 1:n_runs) %dopar% {
        res <- main_simulation(tmax,tdiv, vax_alloc_period, LD, S, E, I, R,
                               SV, EV, IV, RV, modelParameters, cum_vax_pool_func,
                               vax_allocation_func)
        res <- c(list(res = res), calculate_summaries_args)
        res <- do.call(calculate_summaries_func, res)
        res
    }
    # series version for debugging
    # result <- list(n_runs)
    # for(i in 1:n_runs) {
    #       res <- main_simulation(tmax,tdiv, vax_alloc_period, LD, S, E, I, R,
    #                              SV, EV, IV, RV, modelParameters, cum_vax_pool_func,
    #                              vax_allocation_func)
    #       # put arguments of other_info into global environment, so that they can be passed
    #       # to calculate_summaries_func.  Warns if we're overwriting anything.
    #       list2here(other_info)
    #       calculate_summaries_args <- list_vars_from_environment(formalArgs(calculate_summaries_func))
    #       res <- do.call(calculate_summaries_func, calculate_summaries_args)
    #       result[[i]] <- res
    # }
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
#' S: matrix containing the number of unvaccinated susceptibles in each
#' location, age, risk group at each timestep.
#' nrow(S) = n_groups = n_countres * n_ages * n_riskgroups
#' ncol(S) = tend = tmax * tdiv
#' E: same as Smat but for unvaccinated exposed
#' I: same as Smat but for unvaccinated infectious
#' R: same as Smat but for unvaccinated removed
#' SV: same as Smat but for vaccinated susceptible
#' EV: same as Smat but for vaccinated exposed
#' IV: same as Smat but for vaccinated infectious
#' RV: same as Smat but for vaccinated removed
#' vax_pool: numeric vector containing number of vaccines available at each timestep.
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
    SV <- SV0
    EV <- EV0
    IV <- IV0
    RV <- RV0
    
    n_groups <- length(I0)
    
    ## calculate number of time steps
    tend <- tmax*tdiv
    ## make vector of simulation times
    times <- seq(0,tmax,by=1/tdiv)

    ## If we have seasonality, find the resolution
    ## of the seasonality vector. We re-calculate seasonal impact on FOI
    ## every season_res iterations
    switch_freq <- tend+1
    index <- 1
    LD1 <- LD
    if(is.list(LD)){
        season_res <- length(LD)
        switch_freq <- tend/season_res
        LD1 <- LD[[index]]
    }
    
    ## initialise matrices to store simulation outputs
    Smat <- Emat <- Imat <- Rmat <- matrix(0, n_groups, length(times))
    SVmat <- EVmat <- IVmat <- RVmat <- vax_alloc_mat <- Smat
    
    Smat[,1] <- S0
    Emat[,1] <- E0
    Imat[,1] <- I0
    Rmat[,1] <- R0
    SVmat[,1] <- SV0
    EVmat[,1] <- EV0
    IVmat[,1] <- IV0
    RVmat[,1] <- RV0
    
    ## current number of vaccines
    vax_pool <- 0
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
                  all(SV >= 0), all(EV >= 0), all(IV >= 0), all(RV >= 0))

         if(i %% switch_freq == 0 && i < tend){
             index <- index + 1
             LD1 <- LD[[index]]
        }
        
        #################
        ## VAX PRODUCTION
        #################            
        vax_pool <- vax_pool + cum_vax_pool[i] - cum_vax_pool[i - 1]
        
        #################
        ## VAX ALLOCATION
        #################       
        if(i %% vax_alloc_period == 0) {
            
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
                vax_alloc <- vax_allocation_func(S, E, I, R, vax_pool)
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
                
                SV <- SV + actual_alloc$S
                EV <- EV + actual_alloc$E
                IV <- IV + actual_alloc$I
                RV <- RV + actual_alloc$R
            }
            
            vax_alloc_mat[,i] <- sum_vax_alloc
            
        }
        
        # numeric errors can be introduced, so make sure numbers of individuals
        # are still integers before the transmission process
        
        S <- round(S)
        E <- round(E)
        I <- round(I)
        R <- round(R)
        SV <- round(SV)
        EV <- round(EV)
        IV <- round(IV)
        RV <- round(RV)
        
        #################
        ## INFECTIONS
        #################
        ## Generate force of infection on each group/location
        lambda <- beta*LD1%*%(I + IV)

        ## Generate probability of infection from this
        P_infection <- 1 - exp(-lambda)
        P_infection_vax <- 1 - exp(-lambda*(1 - efficacy))
        
        ## Simulate new infections for each location
        newInfections <- rbinom(n_groups, S, P_infection)
        newInfectionsVax <- rbinom(n_groups, SV, P_infection_vax)
        
        ## Update populations
        S <- S - newInfections
        E <- E + newInfections
        SV <- SV - newInfectionsVax
        EV <- EV + newInfectionsVax
        
        #################
        ## EXPOSED BECOMING INFECTIOUS
        #################
        ## Generate probability of exposed becoming infectious
        P_infectious <- 1 - exp(-sigma)
        
        ## Simulate new recoveries
        newInfectious <- rbinom(n_groups, E, P_infectious)
        newInfectiousVax <- rbinom(n_groups, EV, P_infectious)
        
        ## Update populations
        E <- E - newInfectious
        I <- I + newInfectious
        EV <- EV - newInfectiousVax
        IV <- IV + newInfectiousVax
        
        #################
        ## RECOVERIES AND DEATHS
        #################
        ## Generate probability of removal
        P_removal <- 1 - exp(-gamma)
        
        ## Simulate new recoveries
        newRemovals <- rbinom(n_groups, I, P_removal)
        newRemovalsVax <- rbinom(n_groups, IV, P_removal)

        newRecovered <- rbinom(n_groups, newRemovals, 1 - case_fatality_ratio)
        newRecoveredVax <- rbinom(n_groups, newRemovalsVax, 1 - case_fatality_ratio)
        ## Update populations
        I <- I - newRemovals
        R <- R + newRecovered
        IV <- IV - newRemovalsVax
        RV <- RV + newRecoveredVax
        
        #################
        ## SAVE RESULTS
        #################
        Smat[,i] <- S
        Emat[,i] <- E
        Imat[,i] <- I
        Rmat[,i] <- R
        SVmat[,i] <- SV
        EVmat[,i] <- EV
        IVmat[,i] <- IV
        RVmat[,i] <- RV
        vax_pool_vec[i] <- vax_pool
        
        ## stop simulation if there are no more exposed/infectious individuals

        if(i < (tend + 1) && sum(E + I + EV + IV) == 0) {
            remaining_idx <- seq((i+1), ncol(Smat))
                                        # make a function to fill in the remaining parts of the state matrix
                                        # once we stop the simulation
            
            fill_remaining_closure <- function(n_groups, remaining_idx) {
                f <- function(vec) {
                    matrix(vec, n_groups, length(remaining_idx))
              }
              f
            }
            
            fill_remaining <- fill_remaining_closure(n_groups, remaining_idx)
            
            Smat[,remaining_idx] <- fill_remaining(S)
            Rmat[,remaining_idx] <- fill_remaining(R)
            SVmat[,remaining_idx] <- fill_remaining(SV)
            RVmat[,remaining_idx] <- fill_remaining(RV)
            vax_pool_vec[remaining_idx] <- vax_pool
            break
        }
    }
    
    ## put times as column names for readability of output
    colnames(Smat) <- colnames(Emat) <- colnames(Imat) <- colnames(Rmat) <- times
    colnames(SVmat) <- colnames(EVmat) <- colnames(IVmat) <- colnames(RVmat) <- times
    colnames(vax_alloc_mat) <- times
    
    return(list(beta=beta,S=Smat,E = Emat, I=Imat,R=Rmat,
                SV = SVmat, EV = EVmat, IV = IVmat, RV = RVmat, vax_pool = vax_pool_vec,
                vax_alloc = vax_alloc_mat))
}
