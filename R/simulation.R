run_simulation <- function(simulation_flags, life_history_params,
                           vax_params, sim_params,
                           X,
                           contactMatrix,
                           travelMatrix,
                           cum_vax_pool_func,
                           vax_allocation_func,
                           tmax=100,tdiv=24, vax_alloc_period = 24 * 7){
    
    ageMixing <- simulation_flags[["ageMixing"]]
    riskGroups <- simulation_flags[["riskGroups"]]
    spatialCoupling <- simulation_flags[["spatialCoupling"]]
    normaliseTravel <- simulation_flags[["normaliseTravel"]]
    
    R0 <- life_history_params[["R0"]]
    TR <- life_history_params[["TR"]]
    LP <- life_history_params[["LP"]]
    efficacy <- vax_params[["efficacy"]]
    propn_vax0 <- vax_params[["propn_vax0"]]
    gamma <- 1/TR
    alpha <- 1/LP

    n_countries <- sim_params[["n_countries"]]
    n_ages <- sim_params[["n_ages"]]
    n_riskgroups <- sim_params[["n_riskgroups"]]
    seed_countries <- sim_params[["seedCs"]]
    seed_ns <- sim_params[["seedNs"]]
    seed_ages <- sim_params[["seedAges"]]

    maxIndex <- n_countries*n_ages*n_riskgroups
    groupsPerLoc <- n_ages*n_riskgroups

    if(normaliseTravel){
        Krow <- rowSums(travelMatrix)
        Knorm <- kronecker(matrix(1,1,n_countries),matrix(Krow,n_countries,1))
        travelMatrix <- travelMatrix/Knorm
    }
    Kdelta <- kronecker(diag(groupsPerLoc),travelMatrix)
    K1 <- kronecker(matrix(1,groupsPerLoc,groupsPerLoc),t(travelMatrix))
    
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

    sigma <- sim_params[["seed_vec"]]

    # intial exposed are distributed among vaccinated and unvaccinated proportionally
    EV <- round(sigma * propn_vax0)
    E <- sigma - EV
    SV <- round((X - sigma) * propn_vax0)
    S <- X - sigma - SV

    I <- R <- IV <- RV <- matrix(0, maxIndex)

    modelParameters <- c("gamma"=gamma, "alpha" = alpha, "efficacy" = efficacy,
                         "beta" = beta)
    
    result <- main_simulation(tmax,tdiv, vax_alloc_period, LD, S, E, I, R, 
                              SV, EV, IV, RV, modelParameters, cum_vax_pool_func,
                              vax_allocation_func)
    result
}



## This could be moved to C

main_simulation <- function(tmax, tdiv, vax_alloc_period, LD, S0, E0, I0, R0, 
                            SV0, EV0, IV0, RV0, params,
                            cum_vax_pool_func, vax_allocation_func){
  
    gamma <- params["gamma"]
    alpha <- params["alpha"]
    efficacy <- params["efficacy"]
    beta <- params["beta"]

    S <- S0    
    E <- E0
    I <- I0
    R <- R0
    SV <- SV0    
    EV <- EV0
    IV <- IV0
    RV <- RV0

    n_groups <- length(I0)
    
    tend <- tmax*tdiv
    times <- seq(0,tmax,by=1/tdiv)
    
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
    
    vax_pool_vec <- double(length(times))
    vax_pool <- 0
    cum_vax_pool <- vapply(times, cum_vax_pool_func, double(1))
    
    alloc_minifunc_closure <- function(sum_SEIR, vax_alloc){
      function(comp) {
        alloc <- pmin(comp, round(comp / sum_SEIR * vax_alloc))
        alloc[is.na(alloc)] <- 0
        return(alloc)
      }
    }

    for(i in 2:(tend+1)){
      
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
            
            ## update vax pool
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
        P_infectious <- 1 - exp(-alpha)
        
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
        if(sum(E + I + EV + IV) == 0) {
          # sorry for terrible coding -- will fix later (Ada)
          Smat[,(i+1):ncol(Smat)] <- matrix(S, length(S), length((i+1):ncol(Smat)))
          Emat[,(i+1):ncol(Smat)] <- matrix(E, length(S), length((i+1):ncol(Smat)))
          Imat[,(i+1):ncol(Smat)] <- matrix(I, length(S), length((i+1):ncol(Smat)))
          Rmat[,(i+1):ncol(Smat)] <- matrix(R, length(S), length((i+1):ncol(Smat)))
          SVmat[,(i+1):ncol(Smat)] <- matrix(SV, length(S), length((i+1):ncol(Smat)))
          EVmat[,(i+1):ncol(Smat)] <- matrix(EV, length(S), length((i+1):ncol(Smat)))
          IVmat[,(i+1):ncol(Smat)] <- matrix(IV, length(S), length((i+1):ncol(Smat)))
          RVmat[,(i+1):ncol(Smat)] <- matrix(RV, length(S), length((i+1):ncol(Smat)))
          vax_pool_vec[(i+1):ncol(Smat)] <- vax_pool
          break
        }
    }
    colnames(Smat) <- colnames(Emat) <- colnames(Imat) <- colnames(Rmat) <- times
    colnames(SVmat) <- colnames(EVmat) <- colnames(IVmat) <- colnames(RVmat) <- times

    return(list(beta=beta,S=Smat,E = Emat, I=Imat,R=Rmat,
                SV = SVmat, EV = EVmat, IV = IVmat, RV = RVmat, vax_pool = vax_pool_vec))
}


setup_populations <- function(popn_size, n_countries,
                              age_propns, n_ages,
                              risk_propns, risk_factors,
                              n_riskgroups){
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
    age_group_risk <- rep(c(age_groups), n_riskgroups)

    ## Multiply together to get proportion of population in each location/age/risk group combo
    X <- trunc(risk_matrix*age_group_risk)
    
    labels <- cbind(X,expand.grid("Location"=1:n_countries, "RiskGroup"=1:n_riskgroups,"Age"=1:n_ages))
    return(list(X=X,labels=labels))
}

setup_populations_real_data <- function(demography_filename,
                              risk_propns, risk_factors,
                              n_riskgroups){
  
  demographic_data <- read.csv(demography_filename,sep = ",")
  pop_size <- demographic_data[,"N"]
  age_groups <- as.matrix(pop_size * demographic_data[,seq(3,ncol(demographic_data))])

  ## Generate risk propns for each age group
  ## Enumerate out risk groups to same dimension as countries
  ## Risk proportions are already enumerated out for ages (ie. n_riskgroups*n_ages)
  risk_matrix <- kronecker(c(risk_propns),matrix(1,n_countries,1))
  
  ## Enumerate out country/age propns to same dimensions as risk groups (ie. n_riskgroups*n_ages*n_countries x 1)
  age_group_risk <- rep(c(age_groups), n_riskgroups)
  
  ## Multiply together to get proportion of population in each location/age/risk group combo
  X <- round(risk_matrix*age_group_risk)
  
  location_names <- demographic_data[,"countryID"]
  labels <- cbind(X,expand.grid("Location"=location_names, "RiskGroup"=1:n_riskgroups,"Age"=1:n_ages))
  return(list("X"=X,"labels"=labels, "pop_size" = pop_size))
}

read_contact_data <- function(contact_filename){
  
  contact_data <- read.table(contact_filename,sep = ",",stringsAsFactors = FALSE,row.names = 1,header = TRUE)
  contact_data <- as.data.frame(t(contact_data))
  contact_data <- lapply(contact_data, function(x) matrix(x,nrow = sqrt(nrow(contact_data))))
  return(contact_data)
  
}

setup_travel_real_data <- function(travel_filename, pop_size, travel_params) {

  travel_data <- read.table(travel_filename,sep = ",",stringsAsFactors = FALSE,header = TRUE)
  travel_data <- as.matrix(travel_data)
  colnames(travel_data) <- NULL
  # normalise travel data to mean so that epsilon is more meaningful
  travel_data <- travel_data / mean(travel_data[travel_data != 0]) * mean(pop_size)
  travel_matrix <- diag(pop_size) + travel_params[["epsilon"]] * travel_data
  return(travel_matrix)
  
}
