#' returns a function to sum a state vector across age and risk groups
#' 
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

#' returns a function which distributes vaccines between age and risk groups
#' 
#' returns a function which, 
#' given the number of vaccines allocated to each country, distributes vaccines
#'  between age and risk groups
#'  
#'  @param priorities data frame with rows n_groups = n_ages * n_riskgroups.
#'  Contains the columns "RiskGroup","Age" and "Priority".
#'  If x vaccines are allocated to a country, and the unvaccinated population in
#'  each age/risk group in that country is given by a vector p of length n_groups,
#'  then the number of vaccines allocated to each age/risk group is
#'  x * p * y / sum(p * y) where y is a normalised version of priorities$Priority.
#'  For example, to give vaccines only to risk/age group 1, y = [1,0,0,...].
#'  For non-discriminatory distribution, y = [1,1,1,...] / n_groups.
#'  Note that if y_i = 0 for any age/risk group i, vaccines will never be allocated
#'  to that age/risk group even if all other age/risk groups have been fully vaccinated --
#'  best to use y_i = epsilon << 1 instead.
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
    n_countries <- round(nrow(labels) / n_groups)
    y_long_vec <- rep(1 / n_groups, nrow(labels))
  } else {
    stop("code with specified priorities doesn't work yet")
    ## sort priorities
    # priorities <- priorities[order(priorities$Age, priorities$RiskGroup),]
    # n_groups <- nrow(priorities)
    # n_countries <- round(nrow(labels) / n_groups)
    # y <- normalise(priorities$Priority)
    # y_long_vec <- rep(y, each = n_countries)
  }
  
  state_names <- c("S", "E", "I", "R")
  n_states <- length(state_names)
  ## non-discriminatory distribution among S, E, I, R
  y_long_vec <- rep(y_long_vec, times = n_states)

  sum_age_risk <- sum_age_risk_closure(labels)
  
  # index of each location in y_long_vec
  loc_idx <- lapply(seq_len(n_countries), function(x) seq(x, length(y_long_vec), by = n_countries))
  
  apply_by_loc_closure <- function(loc_idx, n_countries) {
    f <- function(f, vec) {
      list_by_loc <- lapply(loc_idx, function(x) f(vec[x]))
      return(as.numeric(matrix(unlist(list_by_loc), nrow = n_countries, byrow = TRUE)))
    }
    f
  }
  
  apply_by_loc <- apply_by_loc_closure(loc_idx, n_countries)
  
  f <- function(n_vax_allocated, S, E, I, R) {
    
    # ensure that an integer number of vaccines is allocated to each country
    n_vax_allocated <- round_preserve_sum(n_vax_allocated)

    # distribute vaccines between age and risk groups
    vax_alloc <- rep(n_vax_allocated, times = n_groups * n_states)
    distribution_factor <- c(S, E, I, R) * y_long_vec
    distribution_factor <- apply_by_loc(normalise, distribution_factor)
    vax_alloc <- vax_alloc * distribution_factor
    # round ensuring that the number of vaccines allocated to each country stays the same  
    vax_alloc <- apply_by_loc(round_preserve_sum, vax_alloc)

    # split the vaccine allocation vector into S, E, I, R
    vax_alloc <- split(vax_alloc, ceiling(seq_along(vax_alloc) / n_groups / n_countries))
    names(vax_alloc) <- state_names
    return(vax_alloc)
  }
  f
}

#' rounds a numeric vector probabilistically such that its sum (an integer) is preserved
#' 
#' @param vec numeric vector which sums to an integer
#' @return numeric vector of the same length as vec: rounded vector
round_preserve_sum <- function(vec) {
  
  # if vector already integers, do nothing
  # old implementation -- all.equal(vec - round(vec)) -- fails for large values in vec
  if(isTRUE(all.equal(vec - round(vec), numeric(length(vec))))) {
    return(vec)
  }
  
  # if vector doesn't sum to integer, 
  # round down so it does sum to an integer (ensuring positivity)
  sum_vec <- sum(vec)
  rounded_sum <- round(sum_vec)
  fractional_part <- rounded_sum - sum_vec
  if(!isTRUE(all.equal(rounded_sum - sum_vec, 0))) {
    larger_idx <- which(vec > fractional_part)
    # subtract fractional_part from one of those groups, chosen randomly
    subtract_idx <- sample(larger_idx, 1)
    vec[subtract_idx] <- vec[subtract_idx] - fractional_part
  }
  
  n_round_up <- rounded_sum - sum(floor(vec))
  round_up_idx <- sample.int(length(vec), size = n_round_up, replace = FALSE,
                             prob = vec - floor(vec))
  vec_out <- floor(vec)
  vec_out[round_up_idx] <- ceiling(vec[round_up_idx])
  vec_out
}

#' normalises a vector so that it sums to 1
#' 
#' @param vec numeric vector to be normalised
#' @return normalised vector of same length as vec
normalise <- function(vec) {
  
  if(all(vec == 0)) {
    return(vec)
  } else if(any(vec < 0)) {
    stop("vector to be normalised must have all non-negative elements")
  }
  
  vec / sum(vec)
  
}

#' generate n_int random integers that sum to sum_int
#' 
#' @param n_int number of random integers to generate
#' @param sum_int integer to which those integers should sum to
#' @return numeric vector of length n_int
gen_int_sum_int <- function(n_int, sum_int) {
  round_preserve_sum(normalise(runif(n_int)) * sum_int)
}

#' closure to make vaccine allocation function
#' 
#' a closure to make a function which returns a vector of length
#' n = n_countries * n_ages * n_risk_groups:
#' the number of vaccines to be allocated to each location, age, risk group
#'
#' @param user_specified_vax_alloc_func function to allocate vaccines
#' given the arguments
#' #' @param sum_age_risk_func: a function which sums a vector across age and risk groups
#' #' (created in the closure)
#' #' @param travel_matrix: see below
#' #' @param vax_allocation_params: see below
#' #' @param S: numeric vector of length n. number of unvaccinated susceptibles
#' #' @param E: numeric vector of length n. number of unvaccinated exposed
#' #' @param I: numeric vector of length n. number of unvaccinated infectious
#' #' @param R: numeric vector of length n. number of unvaccinated recovered
#' #' @param vax_pool: numeric vector of length 1. number of vaccines available
#' @param travel_matrix square matrix with side length n_countries. 
#' travel_matrix[x,y] is the proportion of time an individual in location x spends
#' in location y
#' @param vax_allocation_params list of parameters for vaccine allocation
#' @param labels data frame containing location, age and risk group corresponding to a vector element
#' @return a function which returns a vector of length
# n_countries * n_ages * n_risk_groups:
# the number of vaccines to be allocated to each location, age, risk group
#' @export
vaccine_allocation_closure <- function(user_specified_vax_alloc_func, 
                                       travel_matrix, vax_allocation_params, labels) {
  
  ## create a function which sums a state vector across age and risk groups
  sum_age_risk_func <- sum_age_risk_closure(labels)
  distribute_vax_among_age_risk <- 
    distribute_vax_among_age_risk_closure(vax_allocation_params$priorities, labels)
  
  function(S, E, I, R, SV, EV, IV, RV, incidence, vax_pool) {
    
    # vax_pool can be not an integer -- round down
    vax_pool <- floor(vax_pool)
    
    # allocate vaccines to countries
    n_vax_allocated <- user_specified_vax_alloc_func(sum_age_risk_func, travel_matrix,
                                                     vax_allocation_params,
                                                     S, E, I, R, 
                                                     SV, EV, IV, RV, 
                                                     incidence, vax_pool)
    n_vax_allocated <- round_preserve_sum(n_vax_allocated)
    
    # ensure that number of vaccines allocated to each country is fewer than or
    # equal to the number of people eligible for vaccination
    total_individuals <- sum_age_risk_func(S) +
      sum_age_risk_func(E) +
      sum_age_risk_func(I) +
      sum_age_risk_func(R)
    n_vax_allocated <- pmin(n_vax_allocated, total_individuals)
    
    n_vax_allocated <- round_preserve_sum(n_vax_allocated)
    
    # distribute vaccines in each country between individuals of different
    # ages, risk groups, and infection statuses
    split_allocation <- distribute_vax_among_age_risk(n_vax_allocated, S, E, I, R)
    return(split_allocation)
  }
}

#' closure to make vaccine production function
#' 
#' a closure to make a function which returns a scalar:
# the number of vaccines which have ever existed at time t
# (may involve integrating the production rate curve if you want,
# otherwise provide directly)
#' @param user_specified_cum_vax_pool_func: function of vax_production_params and t
#' which returns the number of vaccines which have ever existed at time t
#' @param vax_production_params: list of parameters for vaccine production
#' @return function of t which returns
#' the number of vaccines which have ever existed at time t
#' @export
cum_vax_pool_func_closure <- function(user_specified_cum_vax_pool_func, vax_production_params) {
  f <- function(t) {
    user_specified_cum_vax_pool_func(vax_production_params, t)
  }
  f
}

#' run the simulation and get a vector of the number of deaths and recovereds
#' in each country/age/risk group
#' 
#' for code testing purposes.
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
#' @param case_fatality_ratio_vec vector of case fatality ratio in each location, age, risk group:
#' length n_countries * n_ages * n_riskgroups
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
#' @return list containing the number of deaths and recovereds
#' in each country/age/risk group
run_and_get_deaths_recovered <- function(simulation_flags, life_history_params, vax_params, sim_params,
                                         case_fatality_ratio_vec, X, labels, C3, K, latitudes, 
                                         cum_vax_pool_func, vax_allocation_func, tmax, tdiv, vax_alloc_period) {
  
  res <- run_simulation(simulation_flags, life_history_params, vax_params, sim_params,
                        case_fatality_ratio_vec, X, labels, C3, K, latitudes, 
                        cum_vax_pool_func, vax_allocation_func, tmax, tdiv, vax_alloc_period)
  
  tend <- ncol(res$S)
  deaths <- as.numeric(X - res$S[,tend] - res$SV[,tend] - res$E[,tend] - res$EV[,tend] -
                         res$I[,tend] - res$IV[,tend] - res$R[,tend] - res$RV[,tend])
  recovered <- unname(res$R[,tend] + res$RV[,tend])
  list("deaths" = deaths, "recovered" = recovered)
  
}

#' gather variables with given names from an environment into a list
#' 
#' @param var_names character vector of variable names to gather
#' @param envir environment in which to find parameters
#' @return list of variables with names var_names
#' @export
list_vars_from_environment <- function(var_names, envir = parent.frame()) {
  env_vars <- as.list(envir)
  env_vars <- get_vars_from_list_with_check(env_vars, var_names)
  env_vars
}

#' extract variables from list, throwing an error if they are not found
#' 
#' @param x list of variables
#' @param var_names character vector containing names of variables to extract
#' @return list of selected variables
#' @export
get_vars_from_list_with_check <- function(x, var_names) {
  missing_vars <- var_names[!(var_names %in% names (x))]
  if(length(missing_vars) > 0) {
    stop(paste0("variables missing from list: ", paste0(missing_vars, collapse = " ")))
  }
  x <- x[var_names]
  x
}

#' dump variables from a list to the parent frame
#' 
#' @param x list containing variables
#' @param var_names character vector containing names of variables to dump into 
#' parent frame
#' @param overwrite optional parameter controlling the behaviour if any variables
#' with the same name(s) alraedy exist in the parent frame. If set to "warn",
#' throws a warning; if set to "error", throws an error
#' @return NULL
list2here <- function(x, var_names, overwrite) {
  if(!missing(var_names)) {
    x <- get_vars_from_list_with_check(x, var_names)
  }
  
  
  if(!missing(overwrite)) {
    parent_frame_vars <- ls(parent.frame)
    overwriting_vars <- intersect(names(x), parent_frame_vars)
    if(length(overwriting_vars) > 0) {
      if(overwrite == "warn") {
        lapply(overwriting_vars, function(x) warning(cat("overwriting", x)))
      } else if(overwrite == "error") {
        stop(cat("Attempting to overwrite", x[1]))
      }
    }
  }
  
  list2env(x, envir = parent.frame())
  invisible(NULL)
}

#' return names of all functions in a .R file
#' 
#' @param filename location of .R file
#' @return character vector where each entry is the name of a function in the file
get_func_names <- function(filename) {
  suppress_final_line_warning <- function(w) {
    if(any(grepl("incomplete final line found on", w, fixed = TRUE))) {
      invokeRestart("muffleWarning")
    }
  }
  # read the .R file, suppressing warnings about incomplete final line
  code <- withCallingHandlers(readLines(filename), warning = suppress_final_line_warning)
  # remove spaces
  code <- gsub(" ", "", code)
  # identify lines which define functions
  func_str <- "<-function("
  potential_func_lines <- grep(func_str, code, fixed = TRUE)
  
  # determine whether each function is an outermost function
  is_outside_func <- function(code, potential_func_line) {
    # assume functions on first line are outside functions
    if(potential_func_line == 1) {
      return(TRUE)
    }
    # paste all code up to name of potential function
    pasted_code <- paste0(code[seq_len(potential_func_line - 1)], collapse = "")
    potential_func_subline <- strsplit(code, split = func_str, fixed = TRUE)
    potential_func_subline <- potential_func_subline[1]
    pasted_code <- paste0(pasted_code, potential_func_subline, collapse = "")
    # count number of open and close curly brackets up to potential function name
    count_characters <- function(pasted_code, char_in) {
      n <- gregexpr(char_in, pasted_code, fixed = TRUE)
      length(n[[1]])
    }
    n_brackets <- vapply(c("{", "}"), function(x) count_characters(pasted_code, x), numeric(1))
    # the functino is an outermost function if the number of open and close brackets is the same
    n_brackets[1] == n_brackets[2]
  }
  
  func_lines <- potential_func_lines[vapply(potential_func_lines, 
                                            function(x) is_outside_func(code, x),
                                            logical(1))]
  # split off the function names
  func_lines <- strsplit(code[func_lines], split = func_str, fixed = TRUE)
  func_lines <- vapply(func_lines, function(x) x[1], character(1))
  func_lines
}

#' return function options for vaxedemic
#' 
#' user_specified_cum_vax_pool_func: options for vaccine production functions
#' user_specified_vax_alloc_func: options for vaccine allocation functions
#' calculate_summaries_func: options for functions to calculate summary statistics
#' postprocessing_func: options for postprocessing functions
#' run_func: options for functions to run simulation
#' this function relies on the above functions being in the right .R files 
#' @param package_dir vaxedemic repository directory
#' @return a list with elements user_specified_cum_vax_pool_func, 
#' user_specified_vax_alloc_func etc. where each element is a character vector
#' and each element in that vector is an option for that function
#' @export
get_vaxedemic_func_options <- function(package_dir = getwd()) {
  # where the functions live
  filenames <- c("user_specified_cum_vax_pool_func" = "vax_production_funcs.R",
                 "user_specified_vax_alloc_func" = "vax_alloc_funcs.R",
                 "calculate_summaries_func" = "calc_summaries_funcs.R",
                 "postprocessing_func" = "postprocessing_funcs.R",
                 "run_func" = "run_funcs.R")
  
  func_options_names <- names(filenames)
  
  # append directory
  filenames <- paste0(package_dir, "/R/", filenames)
  filenames <- gsub("//", "/", filenames, fixed = TRUE)
  # read the function names in the .R files
  func_options <- lapply(filenames, get_func_names)
  names(func_options) <- func_options_names
  func_options
}

#' make argument list with shorter number of runs
#' 
#' @param args_list an argument list to run simulations on or off the cluster
#' @param n_runs_test shorter number of runs
#' @return an argument list to run simulations on or off the cluster
shorten_runs <- function(args_list, n_runs_test) {
  args_list$n_runs <- n_runs_test
  args_list$output_prefix <- paste0("test_", args_list$output_prefix)
  short_base_dir <- dirname(args_list$output_prefix)
  if(!file.exists(short_base_dir)) dir.create(short_base_dir)
  args_list
}

#' thin a time series matrix
#' 
#' @param time_series_matrix matrix where each column corresponds to a given
#' time, and the column names are the times
#' @param thin_every an integer.  If 0, ignore.  Otherwise, thin the matrix
#' every this many columns
#' @param thin_integer a logical. If TRUE, thin the matrix to integer values of 
#' times.  If FALSE, ignore
#' @param thin_by_sum logical.  IF TRUE, thin by adding together columns (e.g.
#' if thin_every = 2, add togeher columns 1-2, columns 3-4 etc.)  This is useful
#' for thinning quantities such as incidence.  If FALSE, ignore.
#' @return a thinned matrix with the same number of rows as time_series_matrix,
#' and fewer columns
thin_time_series <- function(time_series_matrix, 
                             thin_every = 0, 
                             thin_integer = FALSE,
                             thin_by_sum = FALSE) {
  if(thin_every == 0 && thin_integer == FALSE) {
    stop("thinning interval not specified")
  }
  
  if(thin_every > 0 && thin_integer) {
    stop("thinning interval specified in two different ways")
  }
  
  stopifnot(is_integer_like(thin_every) && thin_every >= 0)
  
  t_vec <- as.numeric(colnames(time_series_matrix))
  
  if(!missing(thin_every)) {
    idx <- seq(1, ncol(time_series_matrix), by = thin_every)
  } else {
    idx <- vlapply(t_vec, is_integer_like)
    idx <- which(idx)
  }
  
  if(thin_by_sum) {
    idx_end <- c(idx[-1] - 1, ncol(time_series_matrix))
    sum_submatrix <- function(idx, idx_end) {
      submatrix <- time_series_matrix[,seq(idx, idx_end)]
      if(idx != idx_end) {
        submatrix <- apply(submatrix, 1, sum)
      }
      submatrix
    }
    time_series_matrix <- Map(sum_submatrix, idx, idx_end)
    time_series_matrix <- t(do.call(rbind, time_series_matrix))
    # add a t = 0 column
    time_series_matrix <- cbind(matrix(0, ncol = 1, nrow = nrow(time_series_matrix)), 
                                time_series_matrix)
    colnames(time_series_matrix) <- t_vec[c(1, idx_end)]
  } else {
    time_series_matrix <- time_series_matrix[,idx]
  }

  time_series_matrix
}

#' vapply for logicals
#' 
#' @param X X in vapply
#' @param FUN FUN in vapply
#' @return vapply(X, FUN, logical(1), ...)
vlapply <- function(X, FUN, ...) {
  vapply(X, FUN, logical(1), ...)
}

#' vapply for integers
#' 
#' @param X X in vapply
#' @param FUN FUN in vapply
#' @return vapply(X, FUN, integer(1), ...)
viapply <- function(X, FUN, ...) {
  vapply(X, FUN, integer(1), ...)
}

#' vapply for numerics
#' 
#' @param X X in vapply
#' @param FUN FUN in vapply
#' @return vapply(X, FUN, numeric(1), ...)
vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}

#' vapply for characters
#' 
#' @param X X in vapply
#' @param FUN FUN in vapply
#' @return vapply(X, FUN, character(1), ...)
vcapply <- function(X, FUN, ...) {
  vapply(X, FUN, character(1), ...)
}

#' check if a number is integer-type or close enough to an integer
#' 
#' @param x numeric: the number
#' @param tol numeric: tolerance
#' @return TRUE if x is integer-type or within tol of an integer, FALSE otherwise
is_integer_like <- function(x, tol = sqrt(.Machine$double.eps)) {
  is.integer(x) || (is.numeric(x) && abs(x - round(x)) < tol)
}

#' add two things
#' 
#' @param x something that can be added to y
#' @param y something that can be added to x
#' @return x + y
add <- function(x, y) x + y

#' sum over a list
#' 
#' example use case: element-by-element sum over matrices in a list
#' @param xs list of objects that can be summed together
#' @return the sum over the list
sum_list <- function(xs) {
  Reduce(function(x, y) add(x, y), xs)
}

#' dump variables from a list to the parent frame
#' 
#' @param x list containing variables
#' @param var_names character vector containing names of variables to dump into 
#' parent frame
#' @param overwrite optional parameter controlling the behaviour if any variables
#' with the same name(s) alraedy exist in the parent frame. If set to "warn",
#' throws a warning; if set to "error", throws an error
#' @return NULL
#' @export
list2here <- function(x, var_names, overwrite) {
  if(!missing(var_names)) {
    x <- get_vars_from_list_with_check(x, var_names)
  }
  
  
  if(!missing(overwrite)) {
    parent_frame_vars <- ls(parent.frame)
    overwriting_vars <- intersect(names(x), parent_frame_vars)
    if(length(overwriting_vars) > 0) {
      if(overwrite == "warn") {
        lapply(overwriting_vars, function(x) warning(cat("overwriting", x)))
      } else if(overwrite == "error") {
        stop(cat("Attempting to overwrite", x[1]))
      }
    }
  }
  
  list2env(x, envir = parent.frame())
  invisible(NULL)
}