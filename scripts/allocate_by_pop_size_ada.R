library(magrittr)

# allocate vaccines to top n countries, in proportion to those countries' population
# sizes.

allocate_top_n_countries <- function(n_countries_allocated) {
  new_coverage <- current_coverage <- 
    read.csv("data/coverage_data_intersect.csv",stringsAsFactors = FALSE)
  pop_size <- read.csv("data/demographic_data_intersect.csv", stringsAsFactors = FALSE)
  pop_size <- pop_size[order(pop_size$N, decreasing = TRUE),]
  countries_allocated <- pop_size[seq_len(n_countries_allocated), "countryID"]
  if(n_countries_allocated == 1) {
    new_coverage$dose_per_1000 <- 1e-6
  } else {
    new_coverage$dose_per_1000 <- 0
  }

  new_coverage[new_coverage$country %in% countries_allocated, "dose_per_1000"] <- 1000
  write.table(new_coverage, paste0("data/coverage_tables_by_popn/coverage_data_",n_countries_allocated,".csv"),sep=",",row.names=FALSE)
  invisible(new_coverage)
}

if(FALSE) {
  n_countries_allocated <- c(seq_len(5), seq(10, 120, by = 10), 127)
  invisible(lapply(n_countries_allocated, allocate_top_n_countries))
}

# allocate vaccines per capita in a country as a functino of the population size
# in that country. alpha is the power coefficient, where 0 means equal
# allocation per capita, and 1 means applcation per capita as a linear function
# of population size in the country
allocate_fn_pop_size <- function(alpha) {
  new_coverage <- current_coverage <- 
    read.csv("data/coverage_data_intersect.csv",stringsAsFactors = FALSE)
  pop_size <- read.csv("data/demographic_data_intersect.csv", stringsAsFactors = FALSE)
  pop_size <- pop_size[order(pop_size$countryID),]
  new_coverage$dose_per_1000 <- 0
  pop_size_normalised_to_max <- pop_size$N / max(pop_size$N)
  new_coverage$dose_per_1000 <- pop_size_normalised_to_max^alpha
  write.table(new_coverage, paste0("data/coverage_tables_fn_pop_size/coverage_data_",num2str(alpha),".csv"),sep=",",row.names=FALSE)
  invisible(new_coverage)
}

# alpha <- seq(0, 1, by = .1)
alpha <- seq(1.5, 5, by = .5)
invisible(lapply(alpha, allocate_fn_pop_size))
