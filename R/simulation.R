run_simulation <- function(simulation_flags, life_history_params,
                           sim_params,
                           X,
                           contactMatrix,
                           travelMatrix,
                           tmax=100,tdiv=24){
    
    ageMixing <- simulation_flags[["ageMixing"]]
    riskGroups <- simulation_flags[["riskGroups"]]
    spatialCoupling <- simulation_flags[["spatialCoupling"]]
    normaliseTravel <- simulation_flags[["normaliseTravel"]]
    
    R0 <- life_history_params[["R0"]]
    TR <- life_history_params[["TR"]]
    gamma <- 1/TR

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

    sigma <- matrix(0,maxIndex,1)
    sigma[(seed_countries-1)*groupsPerLoc + seed_ages,1] <- seed_ns

    I <- sigma
    S <- X - I
    R <- matrix(0, maxIndex)

    modelParameters <- c("gamma"=gamma)
    
    result <- main_simulation(tmax,tdiv,LD, I, S, R, modelParameters)
    result
}



## This could be moved to C

main_simulation <- function(tmax, tdiv, LD, I0, S0, R0, params){
    gamma <- params["gamma"]

    I <- I0
    S <- S0
    R <- R0

    n_groups <- length(I0)
    
    tend <- tmax*tdiv
    times <- seq(0,tmax,by=1/tdiv)
    
    Imat <- matrix(0, length(I0), length(times))
    Smat <- matrix(0, length(S0), length(times))
    Rmat <- matrix(0, length(R0), length(times))

    Imat[,1] <- I0
    Smat[,1] <- S0
    Rmat[,1] <- R0

    for(i in 2:(tend+1)){
#################
        ## INFECTIONS
#################
        ## Generate force of infection on each group/location
        lambda <- LD%*%I

        ## Generate probability of infection from this
        P_infection <- 1 - exp(-lambda)

        ## Simulate new infections for each location
        newInfections <- rbinom(n_groups, S, P_infection)

        ## Update populations
        S <- S - newInfections
        I <- I + newInfections

#################
        ## RECOVERIES
#################
        ## Generate probability of recoveries
        P_recover <- 1 - exp(-gamma)

## Simulate new recoveries
        newRecoveries <- rbinom(n_groups, I, P_recover)

        ## Update populations
        I <- I - newRecoveries
        R <- R + newRecoveries

#################
        ## SAVE RESULTS
#################
        Smat[,i] <- S
        Imat[,i] <- I
        Rmat[,i] <- R
    }
    colnames(Smat) <- colnames(Imat) <- colnames(Rmat) <- times

    return(list(beta=beta,S=Smat,I=Imat,R=Rmat))
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
  age_groups <- as.matrix(demographic_data[,2] * demographic_data[,seq(3,ncol(demographic_data))])

  ## Generate risk propns for each age group
  ## Enumerate out risk groups to same dimension as countries
  ## Risk proportions are already enumerated out for ages (ie. n_riskgroups*n_ages)
  risk_matrix <- kronecker(c(risk_propns),matrix(1,n_countries,1))
  
  ## Enumerate out country/age propns to same dimensions as risk groups (ie. n_riskgroups*n_ages*n_countries x 1)
  age_group_risk <- rep(c(age_groups), n_riskgroups)
  
  ## Multiply together to get proportion of population in each location/age/risk group combo
  X <- round(risk_matrix*age_group_risk)
  
  labels <- cbind(X,expand.grid("Location"=1:n_countries, "RiskGroup"=1:n_riskgroups,"Age"=1:n_ages))
  return(list(X=X,labels=labels))
}
