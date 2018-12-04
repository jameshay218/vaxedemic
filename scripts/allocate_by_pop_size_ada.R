# allocate vaccines to top n countries, in proportion to those countries' population
# sizes.

allocate_top_n_countries <- function(n_countries_allocated) {
  new_coverage <- current_coverage <- 
    read.csv("data/coverage_data_intersect.csv",stringsAsFactors = FALSE)
  pop_size <- read.csv("data/demographic_data_intersect.csv", stringsAsFactors = FALSE) %>%
    .[order(pop_size$N, decreasing = TRUE),]
  countries_allocated <- pop_size[seq_len(n_countries_allocated), "countryID"]
  new_coverage$dose_per_1000 <- 0
  new_coverage[new_coverage$country %in% countries_allocated, "dose_per_1000"] <- 1000
  write.table(new_coverage, paste0("data/coverage_tables_by_popn/coverage_data_",n_countries_allocated,".csv"),sep=",",row.names=FALSE)
  invisible(new_coverage)
}

n_countries_allocated <- c(seq_len(5), seq(10, 120, by = 10), 127)
invisible(lapply(n_countries_allocated, allocate_top_n_countries))