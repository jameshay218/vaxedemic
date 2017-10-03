#' function to set up simulation from inputs, then run
#' 
#' @param simulation_flags named vector (or list) with the logical elements "normaliseTravel"
#' and "seasonal"
#' @param life_history_params named vector (or list) with the numeric elements
#' "R0", "TR" (time to recovery) and "LP" (latent period)
#' @param vax_params named vector (or list) with the numeric elements "efficacy"
#' and "propn_vax0" (initial proportion of vaccinated individuals; assumed constant
#' across location, age and risk groups)
#' @param sim_params list with the numeric elements n_countries, n_ages,
#' n_riskgroups and seed_vec.  The last of these is a vector specifying the 
#' initial number of exposed individuals in each location, age and risk group
#' @param X vector of population size in each location, age, risk group:
#' length n_countries * n_ages * n_riskgroups
#' @param contactMatrix either a square matrix with side length n_ages * n_riskgroups,
#' or a list of such matrices, of length n_countries.  
#' If a single matrix, the contact matrix is assumed to be constant across countries.
#' contactMatrix[i,j] denotes the amount of influence an individual of type i
#' has on an individual of type j, where type includes age and risk groups.
#' @param travelMatrix a square matrix with side length n_countries.
#' travelMatrix[i,j] denotes the proportion of time an individual in country i
#' spends in country j.  If simulation_flags[["normaliseTravel"]] == TRUE,
#' the code will row normalise travelMatrix for you.
#' @param latitudes a matrix with 1 column and n_countries rows, giving the latitude
#' of each country in degrees.
#' @param cum_vax_pool_func a function of time which gives the numebr of vaccines
#' ever produced at that time
#' @param vax_allocation_func a function of the current state of the simulation which
#' gives the number of vaccines to allocate to each location, age, risk group
#' @param tmax the number of days for which to run the simulation
#' @param tdiv the number of timesteps per day
#' @param vax_alloc_period allocate vaccines once every this many timesteps
#' @return list containing results of main simulation loop (see main_simulation)
#' @export
run_simulation <- function(simulation_flags, life_history_params,
                           vax_params, sim_params,
                           X,
                           contactMatrix,
                           travelMatrix,
                           latitudes,
                           cum_vax_pool_func,
                           vax_allocation_func,
                           tmax=100,tdiv=24, vax_alloc_period = 24 * 7){
    
    normaliseTravel <- simulation_flags[["normaliseTravel"]]
    seasonal <- simulation_flags[["seasonal"]]
    
    R0 <- life_history_params[["R0"]]
    TR <- life_history_params[["TR"]] # time to recovery
    LP <- life_history_params[["LP"]] # latent period
    efficacy <- vax_params[["efficacy"]]
    propn_vax0 <- vax_params[["propn_vax0"]] # initial proportion of vaccinated individuals
    gamma <- 1/TR # recovery rate
    sigma <- 1/LP # rate of exposed -> infectious
    tdelay <- 90 #Delay from peak summer in northern hemisphere########
    amp <- .5 #Amplitude of seasonality########

    n_countries <- sim_params[["n_countries"]]
    n_ages <- sim_params[["n_ages"]]
    n_riskgroups <- sim_params[["n_riskgroups"]]

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
      beta1 <- Beta1(latitudes,n_countries)
      B <- kronecker(matrix(1,groupsPerLoc,tmax*tdiv),beta1)
      timevec <- seq(1/tdiv,tmax,by=1/tdiv)
      timevec <- t(timevec)
      wave <- sin((timevec-tdelay)*2*pi/365)
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

    L1 <- (kronecker(matrix(1,1,maxIndex),X)*Kdelta*M1)%*%KC

    ## Beta calcualtion (global):
    XX <- 1/gamma*L1
    ev <- eigen(XX)
    ev <- ev$values
    Rstar <- max(abs(ev))
    beta <- R0/Rstar

    ## Convert to smaller time step rates
    beta <- beta/tdiv
    gamma <- gamma/tdiv
    
    LD <- beta*(Kdelta*M1)%*%KC

    #### calculation of force of infection matrix (LD) ends here
    
    ## set initial conditions

    seed_vec <- sim_params[["seed_vec"]]

    # intial exposed are distributed among vaccinated and unvaccinated proportionally
    EV <- round(seed_vec * propn_vax0)
    E <- seed_vec - EV
    SV <- round((X - seed_vec) * propn_vax0)
    S <- X - seed_vec - SV

    I <- R <- IV <- RV <- matrix(0, maxIndex)

    ## gather model parameters
    modelParameters <- list("gamma"=gamma, "sigma" = sigma, "efficacy" = efficacy,
                         "beta" = beta, "M1" = M1, "Kdelta" = Kdelta, "KC"= KC,
                         "Phi" = Phi, "seasonal" = seasonal)
    ## run simulation
    result <- main_simulation(tmax,tdiv, vax_alloc_period, LD, S, E, I, R, 
                              SV, EV, IV, RV, modelParameters, cum_vax_pool_func,
                              vax_allocation_func)
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
    efficacy <- params[["efficacy"]]
    beta <- params[["beta"]]
    M1 <- params[["M1"]]
    Kdelta <- params[["Kdelta"]]
    KC <- params[["KC"]]
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
    
    ## initialise matrices to store simulation outputs
    Smat <- Emat <- Imat <- Rmat <- matrix(0, n_groups, length(times))
    SVmat <- EVmat <- IVmat <- RVmat <- Smat

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
    
    ## function that allocates vaccines proportionally to S, E, I, or R and rounds
    alloc_minifunc_closure <- function(sum_SEIR, vax_alloc){
      function(comp) {
        alloc <- pmin(comp, round(comp / sum_SEIR * vax_alloc))
        alloc[is.na(alloc)] <- 0
        return(alloc)
      }
    }

    for(i in 2:(tend+1)){
      
      ## check that current state is sensible
      stopifnot(all(S >= 0),all(E >= 0), all(I >= 0), all(R >= 0), 
                all(SV >= 0), all(EV >= 0), all(IV >= 0), all(RV >= 0))

#################
## VAX PRODUCTION
#################            
        vax_pool <- vax_pool + cum_vax_pool[i] - cum_vax_pool[i - 1]
        
#################
## VAX ALLOCATION
#################       
        if(i %% vax_alloc_period == 0) {
          # allocate vaccines
          vax_alloc <- round(vax_allocation_func(S, E, I, R, SV, EV, IV, RV, vax_pool))
          if(any(vax_alloc > 0)) {
            # distribute vaccines proportionally among S, E, I, R
            sum_SEIR <- S + E + I + R
            
            alloc_minifunc <- alloc_minifunc_closure(sum_SEIR, vax_alloc)
            E_alloc <- alloc_minifunc(E)
            I_alloc <- alloc_minifunc(I)
            R_alloc <- alloc_minifunc(R)
            S_alloc <- pmin(S, vax_alloc - E_alloc - I_alloc - R_alloc)
            S_alloc[is.na(S_alloc)] <- 0
            
            # actual number allocated after all the rounding
            vax_alloc <- S_alloc + E_alloc + I_alloc + R_alloc
            
            ## update current vax pool
            vax_pool <- vax_pool - sum(vax_alloc)
            
            S <- S - S_alloc
            E <- E - E_alloc
            I <- I - I_alloc
            R <- R - R_alloc
            
            SV <- SV + S_alloc
            EV <- EV + E_alloc
            IV <- IV + I_alloc
            RV <- RV + R_alloc
          }

        }
#################
        ## INFECTIONS
#################

        ## recalculate force of infection matrix if necessary
        if(seasonal){
          M1Phi <- M1*kronecker(matrix(1,n_ages * n_riskgroups *n_countries,1),t(Phi[,i-1]))
          LD <- beta*(Kdelta*M1Phi)%*%KC
        }
        
        ## Generate force of infection on each group/location
        lambda <- LD%*%(I + IV)

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
        ## RECOVERIES
#################
        ## Generate probability of recoveries
        P_recover <- 1 - exp(-gamma)

## Simulate new recoveries
        newRecoveries <- rbinom(n_groups, I, P_recover)
        newRecoveriesVax <- rbinom(n_groups, IV, P_recover)

        ## Update populations
        I <- I - newRecoveries
        R <- R + newRecoveries
        IV <- IV - newRecoveriesVax
        RV <- RV + newRecoveriesVax

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
        if(sum(E + I + EV + IV) == 0) {
          remaining_idx <- seq((i+1), ncol(Smat))
          Smat[,remaining_idx] <- matrix(S, n_groups, length(remaining_idx))
          Rmat[,remaining_idx] <- matrix(R, n_groups, length(remaining_idx))
          SVmat[,remaining_idx] <- matrix(SV, n_groups, length(remaining_idx))
          RVmat[,remaining_idx] <- matrix(RV, n_groups, length(remaining_idx))
          vax_pool_vec[remaining_idx] <- vax_pool
          break
        }
    }
    
    ## put times as column names for readability of output
    colnames(Smat) <- colnames(Emat) <- colnames(Imat) <- colnames(Rmat) <- times
    colnames(SVmat) <- colnames(EVmat) <- colnames(IVmat) <- colnames(RVmat) <- times

    return(list(beta=beta,S=Smat,E = Emat, I=Imat,R=Rmat,
                SV = SVmat, EV = EVmat, IV = IVmat, RV = RVmat, vax_pool = vax_pool_vec))
}

#' function to setup synthetic population sizes for each location, age, risk group, 
#' and labels associated with these
#' 
#' @param popn_size numeric vector of length 1: population size of a country.
#' assumed to be same across countries.
#' @param n_countries numeric vector of length 1: number of countries
#' @param age_propns numeric vector: proportion of individuals in each age group.
#' assumed to be the same across countries.
#' @param risk_propns numeric vector: proportion of individuals in each age group.
#' assumed to be the same across countries and age groups.
#' @param risk_factors numeric vector: degree of increased susceptibility 
#' in each risk group.
#' assumed to be the same across countries and age groups.
#' @return list with the following elements:
#' X: numeric vector of length n_groups = n_countries * n_ages * n_riskgroups
#' containing the population sizes in each country, age, risk group
#' labels: data frame containing three columns: the location, age, risk group
#' corresponding to each element of X.
#' nrow(labels) = length(X)
#' pop_size: numeric vector of length n_countries containing the total population
#' size in each country.
#' @export
setup_populations <- function(popn_size, n_countries,
                              age_propns, risk_propns, risk_factors){

    stopifnot(sum(age_propns) == 1, all(rowSums(risk_propns) == 1))
  
    n_ages <- length(age_propns)
    n_riskgroups <- ncol(risk_propns)
    
    ## Assuming the same population size for each country, the same age
    ## distribution and the same risk proportions. Once we input real data,
    ## we would input these matrices directly
    country_popns <- rep(popn_size, n_countries)
    country_popns <- t(matrix(rep(country_popns,n_ages),nrow=n_ages,byrow=TRUE)) ## Copy population size for each age group

    ## Get proportion of population in each age group
    age_propns_country <- matrix(rep(age_propns, n_countries),
                                 nrow=n_countries,ncol=n_ages,byrow=TRUE)

    ## Use proportion to get actual population size
    age_groups <- country_popns * age_propns_country

    ## Generate risk propns for each age group
    ## Enumerate out risk groups to same dimension as countries
    ## Risk proportions are already enumerated out for ages (ie. n_riskgroups*n_ages)
    risk_matrix <- kronecker(c(risk_propns),matrix(1,n_countries,1))

    ## Enumerate out country/age propns to same dimensions as risk groups (ie. n_riskgroups*n_ages*n_countries x 1)
    age_group_risk <- c(kronecker(matrix(1,n_riskgroups,1), age_groups))

    ## Multiply together to get proportion of population in each location/age/risk group combo
    X <- trunc(risk_matrix*age_group_risk)
    
    labels <- cbind(X,expand.grid("Location"=1:n_countries, "RiskGroup"=1:n_riskgroups,"Age"=1:n_ages))
    return(list(X=X,labels=labels, "pop_size" = rep(popn_size, n_countries)))
}

#' function to setup population sizes for each location, age, risk group, 
#' and labels associated with these, from data
#' 
#' @param demography_filename character vector of length 1: file where demographic
#' data is located
#' @param risk_propns numeric vector: proportion of individuals in each age group.
#' assumed to be the same across countries and age groups.
#' @param risk_factors numeric vector: degree of increased susceptibility 
#' in each risk group.
#' assumed to be the same across countries and age groups.
#' @return list with the following elements:
#' X: numeric vector of length n_groups = n_countries * n_ages * n_riskgroups
#' containing the population sizes in each country, age, risk group
#' labels: data frame containing three columns: the location, age, risk group
#' corresponding to each element of X.
#' nrow(labels) = length(X)
#' pop_size: numeric vector of length n_countries containing the total population
#' size in each country.
#' @export
setup_populations_real_data <- function(demography_filename,
                              risk_propns, risk_factors,
                              n_riskgroups){
  
  n_riskgroups <- ncol(risk_propns)
  
  ## read in demographic data
  demographic_data <- read.csv(demography_filename,sep = ",")
  pop_size <- demographic_data[,"N"]
  age_groups <- as.matrix(pop_size * demographic_data[,seq(3,ncol(demographic_data))])

  ## Generate risk propns for each age group
  ## Enumerate out risk groups to same dimension as countries
  ## Risk proportions are already enumerated out for ages (ie. n_riskgroups*n_ages)
  risk_matrix <- kronecker(c(risk_propns),matrix(1,n_countries,1))
  
  ## Enumerate out country/age propns to same dimensions as risk groups (ie. n_riskgroups*n_ages*n_countries x 1)
  age_group_risk <- c(kronecker(matrix(1,n_riskgroups,1), age_groups))
  
  ## Multiply together to get proportion of population in each location/age/risk group combo
  X <- round(risk_matrix*age_group_risk)
  
  location_names <- demographic_data[,"countryID"]
  labels <- cbind(X,expand.grid("Location"=location_names, "RiskGroup"=1:n_riskgroups,"Age"=1:n_ages))
  return(list("X"=X,"labels"=labels, "pop_size" = pop_size))
}

#' function to construct contact matrix from data
#' 
#' @param contact_filename character vector of length 1: file where contact
#' data is located
#' @return list of length n_countries, of square contact matrices of side length
#' n_ages * n_riskgroups.
#' for each matrix, contactMatrix[i,j] denotes the amount of influence 
#' an individual of type i has on an individual of type j, 
#' where type includes age and risk groups.
#' @export
read_contact_data <- function(contact_filename){
  
  contact_data <- read.table(contact_filename,sep = ",",stringsAsFactors = FALSE,row.names = 1,header = TRUE)
  contact_data <- as.data.frame(t(contact_data))
  contact_data <- lapply(contact_data, function(x) matrix(x,nrow = sqrt(nrow(contact_data))))
  return(contact_data)
  
}

#' function to construct travel matrix from data
#' 
#' @param travel_filename character vector of length 1: file where travel
#' data is located
#' @param pop_size numeric vector of length n_countries: population size in each
#' country
#' @param travel_params named vector/list containing the element "epsilon":
#' a numeric vector of length 1 which scales the off-diagonal terms of the 
#' travel matrix.  Roughly, the ratio of off-diagonal to on-diagonal term size.
#' @return square travel matrix of side length n_countries (unnormalised)
#' @export
setup_travel_real_data <- function(travel_filename, pop_size, travel_params) {

  travel_data <- read.table(travel_filename,sep = ",",stringsAsFactors = FALSE,header = TRUE)
  travel_data <- as.matrix(travel_data)
  colnames(travel_data) <- NULL
  # normalise travel data to mean so that epsilon is more meaningful
  travel_data <- travel_data / mean(travel_data[travel_data != 0]) * mean(pop_size)
  travel_matrix <- diag(pop_size) + travel_params[["epsilon"]] * travel_data
  return(travel_matrix)
  
}

#' function to construct latitude matrix from data
#' 
#' @param latitude_filename character vector of length 1: file where latitude
#' data is located
#' @return matrix with 1 column and nrow = n_countries.
#' The latitude of each country in degrees.
#' @export
read_latitude_data <- function(latitude_filename){
  
  latitude_data <- read.table(latitude_filename, sep = ",", header = TRUE)
  latitudes <- latitude_data$latitude
  return(matrix(latitudes, ncol = 1))
  
}
