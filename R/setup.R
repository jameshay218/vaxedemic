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
    
    ## ensure proportions sum to 1 -- eliminate rounding errors
    demographic_data[,ncol(demographic_data)] <- 1 - rowSums(demographic_data[,seq(3,ncol(demographic_data) - 1)])
    ## alphabetise countries for consistency across data sets
    demographic_data <- demographic_data[order(demographic_data$countryID),]
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
    # alphebetise
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
