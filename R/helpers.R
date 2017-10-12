#' returns a function to sum a state vector (e.g. the current number of susceptibles in each
#' country, age and risk group) across age and risk groups
#' 
#' @param labels data frame with n rows, and columns "Location", "RiskGroup" and "Age"
#' @return a function to sum a state vector (e.g. the current number of susceptibles in each
#' country, age and risk group) across age and risk groups
#' @export
sum_age_risk_closure <- function(labels) {
    function(state) {
        aggregate(state, by = list(labels$Location), FUN = sum)[,2]
    }
}

#' returns a function which, 
#' given the number of vaccines allocated to each country, distributes vaccines
#'  between age and risk groups
#'  
#'  @param priorities data frame with rows n_groups = n_ages * n_riskgroups.
#'  Contains the columns "RiskGroup","Age" and "Priority".
#'  If x vaccines are allocated to a country, and the unvaccinated population in
#'  each age/risk group in that country is given by a vector p of length n_groups,
#'  then the number of vaccines allocated to each age/risk group is
#'  x * p * y. where y is a normalised version of priorities$Priority.
#'  For example, to prioritise risk/age group 1, y = [1,0,0,...].
#'  For non-discriminatory distribution, y = [1,1,1,...] / n_groups.
#'  For consistency with labels, 
#'  @param labels data frame with n rows, and columns "Location", "RiskGroup" and "Age".
#'  n = n_countries * n_ages * n_riskgroups.
#'  @return a function f.
#'  f takes the arguments
#'  n_vax_allocated: numeric vector of length n_countries. 
#'  The number of vaccines allocated to each country.
#'  n_unvaccinated: numeric vector of length n.
#'  The number of unvaccinated individuals in each location, age, risk group.
#'  f returns a numeric vector of length n: the
#'  number of vaccines allocated to each country, age, risk group
distribute_vax_among_age_risk_closure <- function(priorities, labels) {
    ## if priorities not given, default behaviour is non-discriminatory distribution
    if(is.null(priorities)) {
        n_groups <- length(unique(labels$Age)) * length(unique(labels$RiskGroup))
        y_long_vec <- rep(1 / n_groups, nrow(labels))
    } else {
        ## sort priorities
        priorities <- priorities[order(priorities$RiskGroup, priorities$Age),]
        n_groups <- nrow(priorities)
        n_countries <- round(nrow(labels) / n_groups)
        y <- priorities$Priority / sum(priorities$Priority)
        y_long_vec <- rep(y, each = n_countries)
    }
    
    f <- function(n_vax_allocated, n_unvaccinated) {
        vax_alloc_long_vec <- rep(n_vax_allocated, times = n_groups)
        vax_alloc_loc_age_risk <- vax_alloc_long_vec * n_unvaccinated * y_long_vec
        vax_alloc_loc_age_risk
    }
    f
}

#' rounds a numeric vector probabilistically such that its sum (an integer) is preserved
#' 
#' @param vec numeric vector which sums to an integer
#' @return numeric vector of the same length as vec: rounded vector
round_preserve_sum <- function(vec) {
    sum_vec <- sum(vec)
    rounded_sum <- round(sum_vec)
    if(!all.equal(rounded_sum, sum_vec)) {
        stop("vector passed to round_preserve_sum does not sum to integer")
    }
    
    n_round_up <- rounded_sum - sum(floor(vec))
    round_up_idx <- sample.int(length(vec), size = n_round_up, replace = FALSE,
                               prob = vec - floor(vec))
    vec_out <- floor(vec)
    vec_out[round_up_idx] <- ceiling(vec[round_up_idx])
    vec_out
}