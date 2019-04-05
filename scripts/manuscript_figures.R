make_Fig1 <- function() {
  
}

# calc_median_deaths_all <- function() {
#   strategy <- c("incidence", 
#                 "curr_alloc",
#                 "top_n_countries",
#                 "fn_pop_size")
#   production_delay <- c(7, 90, 180)
#   pars <- expand.grid(strategy = strategy, production_delay = production_delay)
#   pars <- rbind(pars, data.frame(strategy = "no_vaccination",
#                                  production_delay = 0))
#   pars$dir_name <- paste0("outputs/pd", pars$production_delay, pars$strategy)
#   
#   filenames <- sapply(pars$dir_name, list.files, pattern = "deaths")
#   filenames <- unlist(Map(function(x, y) paste(x, y, sep = "/"), names(filenames), filenames))
#   median_global_deaths <- vnapply(filenames, calc_median_deaths)
#   median_global_deaths_filenames <- 
# }

read_median_global_deaths <- function(strategy, production_delay, other_info) {
  dir_name <- paste0("outputs/pd", pars$production_delay, pars$strategy, "/")
  if(strategy %in% c("incidence", "curr_alloc", "no_vaccination")) {
    filename <- paste0(dir_names, pars$strategy, "_fixed_median_deaths.rds")
  } else {
    
  }
  if(file.exists(filename)) {
    readRDS(filename)
  } else {
    source_filename <- sub("_median", "", filename, fixed = TRUE)
    median_global_deaths <- calc_median_global_deaths(source_filename)
    saveRDS(median_global_deaths, filename)
    median_global_deaths
  }
}

make_Fig2 <- function() {
    strategy <- c("incidence",
                  "curr_alloc")
    production_delay <- c(7, 90, 180)
    pars <- expand.grid(strategy = strategy, production_delay = production_delay)
    pars <- rbind(pars, data.frame(strategy = "no_vaccination",
                                   production_delay = 0))
    pars$median_global_deaths <- Map(pars$strategy, pars$production_delay, read_median_global_deaths)
    pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_global_deaths"] - pars$median_global_deaths
    
    pars <- pars[pars$strategy != "no_vaccination",]
    pars$production_delay <- factor(pars$production_delay)
    plot1 <- ggplot(pars, aes(x = production_delay, y = deaths_averted, fill = strategy)) +
      geom_bar(stat = "identity") +
      facet_wrap(~strategy, nrow = 1) +
      expand_limits(y = 0)
    plot2 <- ggplot(pars, aes(x = production_delay, y = deaths_averted, fill = strategy)) +
      geom_bar(stat = "identity", position = "dodge") +
      expand_limits(y = 0)
    list(table = pars, plot1 = plot1, plot2 = plot2)
}

make_Fig3 <- function() {
  n_countries <- c(1, 2, 5, 10, 127)
  alpha <- c(.5, 1, 2, 3)

  production_delay <- c(7, 90, 180)
  pars <- expand.grid(strategy = strategy, production_delay = production_delay)
  pars <- rbind(pars, data.frame(strategy = "no_vaccination",
                                 production_delay = 0))
  pars$median_global_deaths <- Map(pars$strategy, pars$production_delay, read_median_global_deaths)
  pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_global_deaths"] - pars$median_global_deaths
  
  pars <- pars[pars$strategy != "no_vaccination",]
  pars$production_delay <- factor(pars$production_delay)
  plot1 <- ggplot(pars, aes(x = production_delay, y = deaths_averted, fill = strategy)) +
    geom_bar(stat = "identity") +
    facet_wrap(~strategy, nrow = 1) +
    expand_limits(y = 0)
  plot2 <- ggplot(pars, aes(x = production_delay, y = deaths_averted, fill = strategy)) +
    geom_bar(stat = "identity", position = "dodge") +
    expand_limits(y = 0)
  list(table = pars, plot1 = plot1, plot2 = plot2)
}