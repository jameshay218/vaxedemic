% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postprocessing_funcs.R
\name{postprocessing_deaths}
\alias{postprocessing_deaths}
\title{postprocess simulations to save the cumulative number of dead individuals for each country}
\usage{
postprocessing_deaths(res_list, runName, other_info, output_prefix)
}
\arguments{
\item{res_list}{summary statistics for each run and processed inputs produced by run_fixed_params}

\item{runName}{character string to make filename out of}

\item{other_info}{list which provides any other information needed, such as to calculate the summaries}

\item{output_prefix}{character vector of length 1.  Prefix for output filenames}
}
\description{
postprocess simulations to save the cumulative number of dead individuals for each country
}
