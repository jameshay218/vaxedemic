% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{list_vars_from_environment}
\alias{list_vars_from_environment}
\title{gather variables with given names from an environment into a list}
\usage{
list_vars_from_environment(var_names, envir = parent.frame())
}
\arguments{
\item{var_names}{character vector of variable names to gather}

\item{envir}{environment in which to find parameters}
}
\value{
list of variables with names var_names
}
\description{
gather variables with given names from an environment into a list
}
