#' @export
setup_inputs <- function(wd, simulation_flags, life_params, travel_params){
    if(simulation_flags[["real_data"]]){
        setup <- setup_real_data(wd=wd, simulation_flags=simulation_flags, travel_params=travel_params)
    } else {
        setup <- setup_sim_data(simulation_flags=simulation_flags, travel_params=travel_params)
    }
    return(setup)
}

#' @export
setup_real_data <- function(simulation_flags,
                            life_params=NULL,
                            travel_params=NULL,
                            wd="~/Documents/vaxedemic/",
                            demography_filename="data/demographic_data_intersect.csv",
                            contact_filename="data/contact_data_intersect.csv",
                            travel_filename="data/flight_data_intersect.csv",
                            latitude_filename="data/latitudes_intersect.csv",
                            risk_filename="data/risk_group_data.csv"){   
    ## Read in demography data
    dem_tmp <- read.csv(paste0(wd,demography_filename), sep=",")
    n_countries <- nrow(dem_tmp)
    n_ages <- ncol(dem_tmp)-2

    ## Risk group data
    risk_propns <- as.matrix(read.csv(paste0(wd,risk_filename),sep=",",header=FALSE))
    n_riskgroups <- ncol(risk_propns)
    
    if(nrow(risk_propns) != n_ages) {
        stop("number of age groups inconsistent between data sets")
    }
    if(ncol(risk_propns) !=n_riskgroups) {
        stop("number of risk groups inconsistent between data sets")
    }

    ############
    ## THIS IS WHERE WE'D CHANGE TO READ IN RISK FACTORS
    ############
    risk_factors <- rep(1, n_riskgroups) ## Assume that each risk group has same modifier

    ## Enumerate out risk factors for each age group
    age_specific_riskgroup_factors <- matrix(rep(risk_factors,each=n_ages),
                                             ncol=n_riskgroups)   
    
    ## construct demography matrix
    demography_matrix <- setup_populations_real_data(paste0(wd,demography_filename),
                                                     risk_propns, risk_factors)

    ## construct travel matrix
    travelMatrix <- setup_travel_real_data(paste0(wd,travel_filename), demography_matrix$pop_size, travel_params)
    
    ## construct latitude vector
    latitudes <- read_latitude_data(paste0(wd, latitude_filename))
    
    ## Generate risk factor modifier. ie. modifier for each age/risk group pair, same dimensions as C2
    ## The risk factor modifier modifies the susceptibility of age/risk groups
    risk <- c(t(age_specific_riskgroup_factors))
    risk_matrix <- t(kronecker(risk,matrix(1,1,n_riskgroups*n_ages)))
    
    ## Generate a contact matrix with dimensions (n_ages*n_riskgroups) * (n_ages*n_riskgroups).
    ## ie. get age specific, then enumerate out by risk group.
    ## If we have country specific contact rates, we get a list
    ## of these matrices of length n_countries
    C1 <- read_contact_data(paste0(wd,contact_filename))
    if(is.list(C1)) { # for country specific contact rates
        C2 <- lapply(C1, function(x) kronecker(x, matrix(1,n_riskgroups,n_riskgroups)))
        C3 <- lapply(C2, function(x) x*risk_matrix)
    } else {
        C2 <- kronecker(C1, matrix(1,n_riskgroups,n_riskgroups))
        C3 <- C2*risk_matrix
    }
    return(list("popns"=demography_matrix$X,"labels"=demography_matrix$labels,
                "contactMatrix"=C3,"travelMatrix"=travelMatrix,"latitudes"=latitudes,
                "n_countries"=n_countries,"n_ages"=n_ages,"n_riskgroups"=n_riskgroups))
}

#' @export
setup_sim_data <- function(simulation_flags,
                           life_params=NULL,
                           travel_params=NULL,
                           popn_size=100000,
                           n_riskgroups=2,
                           n_ages=4,
                           age_propns=c(5,14,45,16)/80,
                           n_countries=10){
    ## Setup risk group data
    risk_propns <- rep(1/n_riskgroups,n_riskgroups) ## Assume risk groups are uniformly distributed
    risk_propns <- matrix(rep(risk_propns,each=n_ages),ncol=n_riskgroups) ## Assume that proportion of ages in each risk group are the same for all ages
    risk_factors <- rep(1, n_riskgroups) ## Assume that each risk group has same modifier
    ## Enumerate out risk factors for each age group
    age_specific_riskgroup_factors <- matrix(rep(risk_factors,each=n_ages),ncol=n_riskgroups)
    
    ## construct demography matrix
    tmp <- setup_populations(popn_size,n_countries,age_propns, 
                             risk_propns, risk_factors)
    ## construct travel matrix
    travelMatrix <- matrix(1,n_countries,n_countries)+999*diag(n_countries) #Travel coupling - assumed independent of age (but can be changed)
    ## construct latitude vector
    latitudes <- matrix(seq_len(n_countries), n_countries, 1)

    ## Generate a contact matrix with dimensions (n_ages*n_riskgroups) * (n_ages*n_riskgroups).
    ## ie. get age specific, then enumerate out by risk group.
    ## If we have country specific contact rates, we get a list
    ## of these matrices of length n_countries
    ## Contact rates
    contactRates <- c(6.92,.25,.77,.45,.19,3.51,.57,.2,.42,.38,1.4,.17,.36,.44,1.03,1.83)
    contactDur <- c(3.88,.28,1.04,.49,.53,2.51,.75,.5,1.31,.8,1.14,.47,1,.85,.88,1.73)

    ## Generate risk factor modifier. ie. modifier for each age/risk group pair, same dimensions as C2
    ## The risk factor modifier modifies the susceptibility of age/risk groups
    risk <- c(t(age_specific_riskgroup_factors))
    risk_matrix <- t(kronecker(risk,matrix(1,1,n_riskgroups*n_ages)))
    
    ## for now, make contact matrices same for all countries
    C1 <- generate_contact_matrix(contactRates, contactDur, n_ages, simulation_flags[["ageMixing"]])
    if(simulation_flags[["country_specific_contact"]]) {
        C1 <- rep(list(C1), n_countries)
        C2 <- lapply(C1, function(x) kronecker(x, matrix(1,n_riskgroups,n_riskgroups)))
        C3 <- lapply(C2, function(x) x*risk_matrix)
    } else {
        C2 <- kronecker(C1, matrix(1,n_riskgroups,n_riskgroups))
        C3 <- C2*risk_matrix
    }
    return(list("popns"=demography_matrix$X,"labels"=demography_matrix$labels,
                "contactMatrix"=C3,"travelMatrix"=travelMatrix,"latitudes"=latitudes,
                "n_countries"=n_countries,"n_ages"=n_ages,"n_riskgroups"=n_riskgroups))
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
    X <- round(risk_matrix*age_group_risk)
    
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
                                        risk_propns, risk_factors){
    
    n_riskgroups <- ncol(risk_propns)
    
    ## read in demographic data
    demographic_data <- read.csv(demography_filename,sep = ",")
    n_countries <- nrow(demographic_data)
    n_ages <- ncol(demographic_data)-2
    
    ## ensure proportions sum to 1 -- eliminate rounding errors
    demographic_data[,ncol(demographic_data)] <- 1 - rowSums(demographic_data[,seq(3,ncol(demographic_data) - 1)])
    ## alphabetise countries for consistency across data sets
    demographic_data <- demographic_data[order(demographic_data$countryID),]
    pop_size <- demographic_data[,"N"]
    age_groups <- as.matrix(pop_size * demographic_data[,seq(3,ncol(demographic_data))])
    
    ## Generate risk propns for each age group
    ## Enumerate out risk groups to same dimension as countries
    ## Risk proportions are already enumerated out for ages (ie. n_riskgroups*n_ages)
    risk_matrix <- kronecker(c(t(risk_propns)),matrix(1,n_countries,1))
    
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
    # alphabetise
    contact_data <- contact_data[order(rownames(contact_data)),]
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
    # alphabetise for consistency between data sets
    travel_data <- travel_data[order(colnames(travel_data)),order(colnames(travel_data))]
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
    # alphabetise
    latitudes <- latitudes[order(latitude_data$Location)]
    return(matrix(latitudes, ncol = 1))
    
}

#' Uniform location age matrix
#'
#' Given a vector of age boundaries and a maximum age, returns a matrix with a uniform age distribution
#' @param ages the vector of age boundaries
#' @param maxAge the maximum age
#' @return the matrix of proportions in each age goup
#' @export
generate_age_matrix_uniform <- function(ages,maxAge){
    ageProp <- matrix(ages/maxAge, length(ages),1)
    ageProp
}

#' Age matrix proportions
#'
#' Given a vector of proportions in each age group, converts this to a matrix
#' @param ages the vector of ages
#' @return the same vector but as an nx1 matrix
#' @export
generate_age_matrix <- function(ages){
    ageProp <- matrix(ages, nrow=length(ages),ncol=1)
    ageProp    
}


#' Age mixing manipulatoin
#'
#' Converts a vector of age-mixing rates or contact durations into a matrix. If a matrix is passed, just returns the original matrix. Can be switched off to return a matrix of 1s
#' @param contactVector vector of length n_ages*n_ages giving age specific contact rates/durations
#' @param n_ages the number of age groups
#' @param ON bool, if FALSE, turns off age mixing
#' @param TRANSPOSE bool, if TRUE, returns the transpose of the created matrix
#' @return the matrix of age-specific contact rates/durations
#' @export
generate_age_mixing <- function(contactVector, n_ages, ON=TRUE, TRANSPOSE=TRUE){
    if(!ON) return(matrix(1,n_ages,n_ages))
    if(is.vector(contactVector)){
        Cnum <- matrix(contactVector, n_ages, n_ages)
        if(TRANSPOSE) Cnum <- t(Cnum)
    } else {
        Cnum <- contactVector
    }
    return(Cnum)
}

#' Age mixing matrix generation
#' 
#' Generates the age-mixing contact matrix taking into account contact rates and durations
#' @param contactRates the vector or matrix of age-specific contact rates
#' @param contactDur the vector or matrix of age-specific contact durations
#' @param n_ages the number of age groups
#' @param ON bool, if FALSE, turns off age mixing
#' @param TRANSPOSE bool, if TRUE, returns the transpose of the created matrix
#' @return the matrix of age-specific contact rates
#' @export
generate_contact_matrix <- function(contactRates, contactDur, n_ages, ON=TRUE, TRANSPOSE=TRUE){
    Cnum <- generate_age_mixing(contactRates, n_ages, ON, TRANSPOSE)
    Cdur <- generate_age_mixing(contactDur, n_ages, ON, TRANSPOSE)
    C1 <- Cnum * Cdur
    return(C1)   
}

kronecker_by_group <- function(n_groups, y){
    x <- matrix(1, n_groups, n_groups)
    y <- kronecker(x, y)
    y    
}

generate_risk_matrix <- function(propRisk, n_riskgroups, transpose=TRUE){
    non_risk <- 1 - propRisk
    risk_mat <- matrix(c(non_risk, propRisk), length(propRisk), n_riskgroups)
    if(TRANSPOSE) risk_mat <- t(risk_mat)
}

#' Read and process seasonal vaccine coverage data
#' 
#' Calculates the proportion of seasonal vaccine doses distributed to each country in 2013
#' @param coverage_filename the vector or matrix of age-specific contact rates
#' @param labels a data frame containing the number of individuals in each location, age, risk group
#' @return the proportion of doses distributed to each country
#' @export
read_coverage_data <- function(coverage_filename, labels) {
  sum_age_risk_func <- sum_age_risk_closure(labels)
  pop_size <- sum_age_risk_func(labels$X)
  coverage_df <- read.table(coverage_filename, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  browser()
  stopifnot(all(coverage_df$country == levels(labels$Location)))
  coverage <- coverage_df$dose_per_1000 * pop_size
  normalise(coverage)
}
