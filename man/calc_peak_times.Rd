% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_summaries_funcs.R
\name{calc_peak_times}
\alias{calc_peak_times}
\title{calculate the peak time by country}
\usage{
calc_peak_times(res, labels)
}
\arguments{
\item{res}{output of main_simulation}

\item{labels}{data frame outputted by setup_inputs containing the location/
age/risk group for each row in res}
}
\value{
a list containing a numeric vector of length n_countries: 
the peak time for each country
}
\description{
calculate the peak time by country
}
