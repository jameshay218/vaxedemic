% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{get_vaxedemic_func_options}
\alias{get_vaxedemic_func_options}
\title{return function options for vaxedemic}
\usage{
get_vaxedemic_func_options(package_dir = getwd())
}
\arguments{
\item{package_dir}{vaxedemic repository directory}
}
\value{
a list with elements user_specified_cum_vax_pool_func, 
user_specified_vax_alloc_func etc. where each element is a character vector
and each element in that vector is an option for that function
}
\description{
user_specified_cum_vax_pool_func: options for vaccine production functions
user_specified_vax_alloc_func: options for vaccine allocation functions
calculate_summaries_func: options for functions to calculate summary statistics
postprocessing_func: options for postprocessing functions
run_func: options for functions to run simulation
this function relies on the above functions being in the right .R files
}
