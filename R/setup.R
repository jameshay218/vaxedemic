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
