# How to make figures in manuscript
# 1. Download files from https://drive.google.com/open?id=150nv86_WLibdEAZDurlOe1X-TWV0sc8S
# (this is Data/vaxedemic_manuscript_data.zip in the vaxedemic google drive folder)
# 2. make a directory in your vaxedemic code main folder (so one down from the folder containing this script) called "outputs"
# 3. unzip the files there, maintaining data structure.  there should now be two directories outputs/deaths_only/ 
# and outputs/pd0no_vaccination/ containing
# e.g. outputs/deaths_only/pd0no_vaccination/no_vaccination_fixed_global_deaths.rds
# 4. set your working directory to the vaxedemic code main folder and run devtools::load_all()
# 5. change save_dir on line 20 of this script to a folder on your drive
# 6. source this script
# 7. run 
# make_Fig1(TRUE, TRUE)
# make_Fig2(TRUE, TRUE)
# make_Fig3(TRUE, TRUE)
# make_Fig3(FALSE, FALSE, TRUE, TRUE)
# make_FigS1(TRUE, TRUE)
# make_FigS2(TRUE, TRUE)
# the raw data for Fig 1 is in https://drive.google.com/open?id=1540gR4un0Dz-m3P3TB9JhdnFIR3iAqJi
# in case you want to revisualise
# (this is Data/no_vax_manuscript_data.zip in the vaxedemic google drive folder)
# contact Ada with any questions

# directory in which to save plots and tables
save_dir <- "~/overleaf/vaxedemic_manuscript/figs_delayed_protection/"

make_dirname <- function(strategy, production_delay, stockpile_size, free_param, seedCountries) {
  # no vaccination results are independent of production delay and stockpile size, so we only ran them once
  if(strategy == "no_vaccination") {
    production_delay <- stockpile_size <- 0 
  }
  dir_name <- paste0("no_seasonality/pd", production_delay, strategy,
                     "_stockpile", num2str(stockpile_size), "_", seedCountries, "/")
  # dir_name <- paste0("outputs_delayed_protection/pd", production_delay, strategy,
                     # "_stockpile", num2str(stockpile_size), "_", seedCountries, "/")
  # if(!dir.exists(dir_name)) {
  #   dir_name <- paste0("outputs/deaths_only/pd", production_delay, strategy, "/")
  # }
  dir_name
}

read_median_global_deaths <- function(strategy, production_delay, stockpile_size, free_param, seedCountries) {
  
  dir_name <- make_dirname(strategy, production_delay, stockpile_size, free_param, seedCountries)
  print(dir_name)
  if(strategy %in% c("incidence", "curr_alloc", "no_vaccination")) {
    filename <- paste0(dir_name, strategy, "_fixed_median_deaths.rds")
  } else {
    filename <- paste0(dir_name, strategy, "_coverage_data_", 
                       num2str(free_param), "_median_deaths.rds")
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

read_global_deaths <- function(strategy, production_delay, stockpile_size, free_param, seedCountries,
                               alloc_split = "global") {
  dir_name <- make_dirname(strategy, production_delay, stockpile_size, free_param, seedCountries)
  
  if(strategy %in% c("incidence", "curr_alloc", "no_vaccination")) {
    filename <- paste0(dir_name, strategy, "_fixed_",
                       alloc_split, "_deaths.rds")
  } else {
    filename <- paste0(dir_name, strategy, "_coverage_data_", 
                       num2str(free_param), "_", alloc_split, "_deaths.rds")
  }
  if(file.exists(filename)) {
    readRDS(filename)
  } else {
    source_filename <- sub(paste0("_", alloc_split), "", filename, fixed = TRUE)
    global_deaths <- calc_global_deaths(source_filename, alloc_split)
    saveRDS(global_deaths, filename)
    global_deaths
  }
}

bootstrap_deaths_averted <- function(strategies, production_delays, 
                                     stockpile_sizes, free_params, 
                                     seedCountries, alloc_split = "global") {
  # set seed for reproducible bootstrapping
  set.seed(1)
  label_with_strategy <- function(global_deaths, strategy_label) {
    data.frame(deaths = global_deaths, 
               strategy = strategy_label)
  }
  
  global_deaths <- Map(read_global_deaths, strategies, production_delays, 
                       stockpile_sizes, free_params, seedCountries,
                       alloc_split) %>%
    Map(label_with_strategy, ., c(1,2)) %>%
    do.call(rbind, .)
  
  difference_medians <- function(data, indices) {
    d <- data[indices,1]
    d1 <- d[data$strategy == 1]
    d2 <- d[data$strategy == 2]
    median(d1) - median(d2)
  }
  
  boot::boot(global_deaths, difference_medians, 500, strata = global_deaths$strategy)
}

bootstrap_median_deaths <- function(strategy, production_delay, stockpile_size, 
                                    free_param, seedCountries, alloc_split = "global") {
  # set seed for reproducible bootstrapping
  set.seed(1)  
  global_deaths <- read_global_deaths(strategy, production_delay, stockpile_size,
                                      free_param, seedCountries, alloc_split)
  my_median <- function(data, indices) {
    median(data[indices])
  }
  
  boot::boot(global_deaths, my_median, 500)
}

bootstrap_median_deaths_wrapper <- function(strategy, production_delay, 
                                            stockpile_size = 0, free_param = 0, 
                                            seedCountries, alloc_split = "global") {
  production_delay <- as.numeric(production_delay)
  stockpile_size <- as.numeric(stockpile_size)
  free_param <- as.numeric(free_param)
  bootstrap_sol <- bootstrap_median_deaths(strategy, production_delay, 
                                           stockpile_size, free_param, 
                                           seedCountries, alloc_split)
  quantile(bootstrap_sol$t, probs = c(.025, .5, .975))
}

bootstrap_deaths_averted_wrapper <- function(strategy, production_delay, 
                                             stockpile_size = 0, free_param = 0, 
                                             seedCountries, alloc_split = "global") {
  production_delay <- as.numeric(production_delay)
  stockpile_size <- as.numeric(stockpile_size)
  free_param <- as.numeric(free_param)
  strategies <- c("no_vaccination", strategy)
  production_delays <- c(0, production_delay)
  stockpile_sizes <- rep(stockpile_size, 2)
  free_params <- c(0, free_param)
  seedCountries <- rep(seedCountries, 2)
  bootstrap_sol <- bootstrap_deaths_averted(strategies, production_delays, 
                                            stockpile_sizes, free_params, 
                                            seedCountries, alloc_split)
  quantile(bootstrap_sol$t, probs = c(.025, .5, .975))
}

print_table <- function(pars, filename, strategy) {
  bootstrap <- any(grepl("lower", colnames(pars)))
  
  if(bootstrap) {
    print(xtable::xtable(pars),
          sanitize.text.function = identity,
          file = filename)
    pars <- pars[,-grep("lower", colnames(pars))]
    pars <- pars[,-grep("upper", colnames(pars))]
  }
  
  if(missing(strategy)) {
    pars$strategy <- factor(pars$strategy, 
                            levels = c("No vaccination", "Current allocation", "Incidence"))
    pars <- pars[order(pars$strategy),]
  } else {
    pars <- pars[pars$strategy == strategy, -1]
  }  
  
  
  levels(pars$production_delay)[1] <- 0
  levels(pars$production_delay) <- c(levels(pars$production_delay), "NA")
  
  pars$stockpile_size <- pars$stockpile_size / 1e6
  
  round_colnames <- c("stockpile_size", "median_deaths", "deaths_averted")
  pars[,round_colnames] <- 
    t(apply(pars[, round_colnames], 1, function(x) as.character(round(x))))
  
  pars[pars$strategy == "No vaccination", c("production_delay", "stockpile_size", "deaths_averted")] <- "NA"
  
  alignment <- rep("r", ncol(pars) + 1) 
  if(missing(strategy)) {
    colnames(pars) <- c("Strategy", 
                        "Production delay (days)", 
                        "Stockpile size (million)",
                        "Median deaths",
                        "Deaths averted")
    alignment[2] <- "l"
  } else {
    if(strategy == "Top n countries") {
      free_param_name <- "$n$"
    } else {
      free_param_name <- "$\\alpha$"
    }
    colnames(pars) <- c("Production delay (days)", 
                        free_param_name,
                        "Stockpile size (million)",
                        "Median deaths",
                        "Deaths averted")
  }
  
  rownames(pars) <- NULL
  
  print(xtable::xtable(pars, align = alignment),
        sanitize.text.function = identity,
        file = filename, include.rownames=FALSE, append = bootstrap)
}

#' make Fig. 1 in the manuscript: Simulation results without vaccination: attack rates and peak times by country.
#' @param save_output logical: whether to save the output as a pdf as well as outputting the ggplot object
#' @output a list of two ggplot objects: attack_rates and peak_times
make_Fig1 <- function(save_output) {
  
  read_and_sort_data <- function(filename) {
    dat <- read.csv(filename)
    # sort countries by latitude
    dat$Location <- as.character(dat$Location)
    dat$Location <- factor(dat$Location, 
                           levels = as.character(dat$Location[order(dat$latitude)]))
    dat
  }
  
  # read in median and 95% CI for attack rates for each country
  attack_rates_data <- read_and_sort_data("outputs/pd0no_vaccination/no_vaccination_country_attack_rates_fixed_median_ci.csv")
  # plot attack rates
  attack_rates <- plot_country_attack_rates(attack_rates_data)
  
  # read in median and 95% CI for peak times for each country
  peak_times_data <- read_and_sort_data("outputs/pd0no_vaccination/no_vaccination_peakTimes_fixed_median_ci.csv")
  # plot peak times
  peak_times <- plot_peak_times(peak_times_data)
  
  # save output
  if(save_output) {
    png(paste0(save_dir, "no_vaccination_country_attack_rates_fixed_plot.png"),width=800,height=1200)
    plot(attack_rates)
    dev.off()
    png(paste0(save_dir, "no_vaccination_peakTimes_fixed_plot.png"),width=800,height=1200)
    plot(peak_times)
    dev.off()
  }
  
  list(attack_rates = attack_rates, peak_times = peak_times)
}

#' make Fig. 2 in the manuscript: comparing deaths averted for allocation by incidence vs current allocation,
#' for different production delays
#' 
#' deaths averted = median deaths over n simulations of no vaccination,
#' minus median deaths over n simulations of vaccination
#' 
#' @param save_output logical: whether to save the output as a pdf as well as outputting the ggplot object
#' @param bootstrap logical: if TRUE, plot median and 95% CI for deaths averted, calculated using
#' bootstrap replicates of the stochastic simulations.  if FALSE, plot point
#' estimate of deaths averted only.
#' @output a list of three objects:
#' pars is a data frame containing the median deaths and deaths averted under each
#' vaccination strategy and production delay.  if bootsrap = TRUE, also contains 95% CI.
#' plot1 and plot2 are ggplot objects: two different visualisations of the deaths averted
make_Fig2 <- function(save_output = FALSE, bootstrap = TRUE) {
  # define vaccination strategies to plot
  strategy <- c("incidence",
                "curr_alloc")
  # define production delays to plot
  production_delay <- c(0, 7, 90, 180)
  # make grid of all combinations of vaccination strategies and production delays
  pars <- expand.grid(strategy = strategy, production_delay = production_delay,
                      stringsAsFactors = FALSE)
  pars <- rbind(pars, data.frame(strategy = "no_vaccination", production_delay = 0,
                                 stringsAsFactors = FALSE))
  # define other parameters (used to construct filename from which to read the data)
  pars$stockpile_size <- 550e6
  pars[pars$production_delay > 0,"stockpile_size"] <- 0
  pars$seedCountries <- "China"
  pars$free_param <- 0
  if(bootstrap) {
    # calculate median and 95% CI for the median deaths under each vaccination strategy
    # and production delay
    median_deaths <- t(apply_named_args(pars, 1, bootstrap_median_deaths_wrapper)) %>%
      as.data.frame
    # calculate median and 95% CI for the deaths averted under each vaccination strategy
    # and production delay
    deaths_averted <- t(apply_named_args(pars, 1, bootstrap_deaths_averted_wrapper)) %>%
      as.data.frame
    # do.call(rbind, .)
    
    # format stuff
    colnames(median_deaths) <- c("lower_deaths", "median_deaths", "upper_deaths")
    colnames(deaths_averted) <- c("lower_deaths_averted", "deaths_averted", "upper_deaths_averted")
    pars <- cbind(pars, as.data.frame(median_deaths), as.data.frame(deaths_averted))
    pars$strategy <- factor(pars$strategy, 
                            levels = c("no_vaccination", "curr_alloc", "incidence"))
    levels(pars$strategy) <- c("No vaccination", "Current allocation", "Incidence")
  } else {
    pars <- rbind(pars, data.frame(strategy = "no_vaccination",
                                   production_delay = 0,
                                   stockpile_size = 0,
                                   seedCountries = pars$seedCountries[1],
                                   stringsAsFactors = FALSE))
    # read in median deaths under each vaccination strategy
    # and production delay, as well as median deaths under no vaccination
    pars$median_deaths <- Map_vapply(read_median_global_deaths, numeric(1), 
                                     pars$strategy, 
                                     pars$production_delay,
                                     pars$stockpile_size,
                                     seedCountries = pars$seedCountries)
    # calculate deaths averted = median deaths with no vaccination - median deaths with vaccination
    pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_deaths"] - pars$median_deaths
    # format stuff
    levels(pars$strategy) <- c("Incidence", "Current allocation", "No vaccination")
    pars$strategy <- factor(pars$strategy, 
                            levels = c("Current allocation", "Incidence", "No vaccination"))
    
  }
  pars$production_delay <- factor(pars$production_delay)
  levels(pars$production_delay) <- c("stockpile", production_delay[production_delay > 0])
  
  # make plots
  plot_wrapper <- function(production_delay_subset) {
    pars <- pars[pars$production_delay %in% production_delay_subset,]
    y_breaks <- seq(0, max(pars$deaths_averted) * 1.1, by = 1e5)
    plot1 <- ggplot(pars[pars$strategy != "No vaccination",], aes(x = production_delay, y = deaths_averted, fill = strategy)) +
      geom_bar(stat = "identity") +
      facet_wrap(~strategy, nrow = 1) +
      scale_y_continuous(breaks = y_breaks, 
                         labels = formatC(y_breaks / 1e5, digits = 2),
                         expand = expand_scale(mult = c(0, .1))) +
      xlab("Production delay (days)") +
      ylab("Deaths averted (x100,000)") +
      theme_bw() +
      theme(legend.position = "none",
            text = element_text(size = 12)) +
      scale_fill_manual(values = gg_color_hue(4)[c(1,4)])
    if(bootstrap) {
      plot1 <- plot1 + geom_errorbar(aes(ymin = lower_deaths_averted, ymax = upper_deaths_averted))
    }
    if(save_output) {
      width <- 20
      height <- 10
      filename_temp <- paste0(save_dir, "Fig2_", production_delay_subset[1])
      ggsave(paste0(filename_temp, ".pdf"),
             plot1, width = width, height = height, units = "cm")
      png(paste0(filename_temp, ".png"),width=width,height=height, units = "cm", res = 300)
      plot(plot1)
      dev.off()
      saveRDS(plot1, file = paste0(filename_temp, ".rds"))
    }
    plot1
  }
  production_delay_subset <- list(c("stockpile", 7),
                                  c(90, 180))
  plot1 <- lapply(production_delay_subset, plot_wrapper)
  # save output as files
  if(save_output) {
    print_table(pars, filename = paste0(save_dir, "Fig2.tex"))
  }
  
  list(table = pars, plot1 = plot1)
}

#' make Fig. 3 in the manuscript: comparing deaths averted for allocation strategies by population size,
#' for different production delays
#' 
#' deaths averted = median deaths over m simulations of no vaccination,
#' minus median deaths over m simulations of vaccination
#' strategy 1 = vaccinate top n countries by population size,
#' with equal allocation to each of those countries by capita
#' strategy 2 = allocate vaccines per capita to each country
#' proportional to x^alpha where x in the population size and alpha is a parameter
#' 
#' @inheritParams make_Fig2
#' 
#' @output a list of three objects:
#' pars is a data frame containing the median deaths and deaths averted under each
#' vaccination strategy and production delay.  if bootsrap = TRUE, also contains 95% CI.
#' plot1 and plot2 are ggplot objects: two different visualisations of the deaths averted
make_Fig3 <- function(plot_reference_curr_alloc = FALSE, 
                      plot_reference_curr_alloc_separate_panel = FALSE,
                      save_output = FALSE, bootstrap = TRUE) {
  # values of n for strategy 1
  n_countries <- c(1, 2, 5, 10, 127)
  # values of alpha for strategy 2
  alpha <- c(.5, 1, 2, 3)
  
  # define production delays to plot
  production_delay <- c(0, 7, 90, 180)
  # make grid of all combinations of vaccination strategies and production delays
  pars <- expand.grid(strategy = "top_n_countries", 
                      production_delay = production_delay,
                      free_param = n_countries,
                      stringsAsFactors = FALSE)
  pars <- rbind(pars, expand.grid(strategy = "fn_pop_size", 
                                  production_delay = production_delay,
                                  free_param = alpha,
                                  stringsAsFactors = FALSE))
  # define other parameters (used to construct filename from which to read the data)
  pars$stockpile_size <- 550e6
  pars[pars$production_delay > 0,"stockpile_size"] <- 0
  pars$seedCountries <- "China"
  
  pars_reference <- pars[seq_len(4),]
  pars_reference$strategy <- "curr_alloc"
  pars_reference$free_param <- 0
  
  if(bootstrap) {
    # calculate median and 95% CI for the median deaths under each vaccination strategy
    # and production delay
    median_deaths <- t(apply_named_args(pars, 1, bootstrap_median_deaths_wrapper)) %>%
      as.data.frame
    # calculate median and 95% CI for the deaths averted under each vaccination strategy
    # and production delay
    deaths_averted <- t(apply_named_args(pars, 1, bootstrap_deaths_averted_wrapper)) %>%
      as.data.frame
    # format stuff
    colnames(median_deaths) <- c("lower_deaths", "median_deaths", "upper_deaths")
    colnames(deaths_averted) <- c("lower_deaths_averted", "deaths_averted", "upper_deaths_averted")
    pars <- cbind(pars, median_deaths, deaths_averted)
    pars$strategy <- factor(pars$strategy, levels = c("top_n_countries", "fn_pop_size"))
    levels(pars$strategy) <- c("Top n countries", "Function of population size")
  } else {
    pars <- rbind(pars, data.frame(strategy = "no_vaccination",
                                   production_delay = 0,
                                   stockpile_size = 0,
                                   free_param = NA,
                                   seedCountries = "China"))
    # read in median deaths under each vaccination strategy
    # and production delay, as well as median deaths under no vaccination
    pars$median_deaths <- Map_vapply(read_median_global_deaths, 
                                     numeric(1), 
                                     pars$strategy, 
                                     pars$production_delay,
                                     pars$stockpile_size,
                                     pars$free_param,
                                     pars$seedCountries)
    # calculate deaths averted = median deaths with no vaccination - median deaths with vaccination
    pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_deaths"] - pars$median_deaths
    
    # format stuff
    levels(pars$strategy) <- c("Top n countries", "Function of population size", "No vaccination")
  }
  
  pars$production_delay <- factor(pars$production_delay)
  levels(pars$production_delay) <- c("stockpile", production_delay[production_delay > 0])
  pars$free_param <- factor(pars$free_param)
  
  # get reference values
  if(plot_reference_curr_alloc) {
    if(bootstrap) {
      # calculate median and 95% CI for the deaths averted under each vaccination strategy
      # and production delay
      deaths_averted_curr_alloc <- t(apply_named_args(pars_reference, 1, bootstrap_deaths_averted_wrapper)) %>%
        as.data.frame
      # format stuff
      colnames(deaths_averted_curr_alloc) <- c("lower_deaths_averted", "deaths_averted", "upper_deaths_averted")
      deaths_averted_curr_alloc$production_delay <- production_delay
    } else {
      stop("Not yet implemented")
      # pars <- rbind(pars, data.frame(strategy = "no_vaccination",
      #                                production_delay = 0,
      #                                stockpile_size = 0,
      #                                seedCountries = pars$seedCountries[1]))
      # # read in median deaths under each vaccination strategy
      # # and production delay, as well as median deaths under no vaccination
      # pars$median_deaths <- Map_vapply(read_median_global_deaths, numeric(1), 
      #                                  pars$strategy, 
      #                                  pars$production_delay,
      #                                  pars$stockpile_size,
      #                                  seedCountries = pars$seedCountries)
      # # calculate deaths averted = median deaths with no vaccination - median deaths with vaccination
      # pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_deaths"] - pars$median_deaths
      # 
      # # format stuff
      # levels(pars$strategy) <- c("Incidence", "Current allocation", "No vaccination")
      # pars$strategy <- factor(pars$strategy, 
      #                         levels = c("Current allocation", "Incidence", "No vaccination"))
      
    }
    
    deaths_averted_curr_alloc$production_delay <- factor(deaths_averted_curr_alloc$production_delay)
    levels(deaths_averted_curr_alloc$production_delay) <- c("stockpile", production_delay[production_delay > 0])
    if(plot_reference_curr_alloc_separate_panel) {
      deaths_averted_curr_alloc <- data.frame(strategy = factor("Current Allocation",
                                                                levels = c(levels(pars$strategy), "Current Allocation")),
                                              production_delay = deaths_averted_curr_alloc$production_delay,
                                              free_param = factor(1, levels = levels(pars$free_param)),
                                              stockpile_size = 0,
                                              seedCountries = "China",
                                              lower_deaths = 0,
                                              median_deaths = 0,
                                              upper_deaths = 0,
                                              lower_deaths_averted = deaths_averted_curr_alloc$lower_deaths_averted,
                                              deaths_averted = deaths_averted_curr_alloc$deaths_averted,
                                              upper_deaths_averted = deaths_averted_curr_alloc$upper_deaths_averted)
      levels(pars$strategy) <- levels(deaths_averted_curr_alloc$strategy)
      
      pars <- rbind(pars, deaths_averted_curr_alloc)
      
    }
  }
  
  # make plots
  plot_wrapper <- function(production_delay_subset) {
    pars <- pars[pars$production_delay %in% production_delay_subset,]
    pars$production_delay <- droplevels(pars$production_delay)
    deaths_averted_curr_alloc <- deaths_averted_curr_alloc[deaths_averted_curr_alloc$production_delay %in% production_delay_subset,]
    deaths_averted_curr_alloc$production_delay <- droplevels(deaths_averted_curr_alloc$production_delay)
    if("stockpile" %in% production_delay_subset) {
      interval <- 5e5
    } else {
      interval <- 1e5
    }
    y_breaks <- seq(0, max(pars$deaths_averted) * 1.1, by = interval)

    plot1 <- ggplot(pars[pars$strategy != "No vaccination",], aes(x = free_param, y = deaths_averted, fill = strategy)) +
      geom_bar(stat = "identity")
    if(plot_reference_curr_alloc && !plot_reference_curr_alloc_separate_panel) {
      plot1 <- plot1 + geom_hline(data = deaths_averted_curr_alloc, 
                                  aes(yintercept = deaths_averted),
                                  colour = gg_color_hue(4)[4],
                                  size = 1)
    }
    
    if(plot_reference_curr_alloc && plot_reference_curr_alloc_separate_panel) {
      xlab_string <- "n                                       alpha                            "
    } else {
      xlab_string <- "n                                             alpha"
    }
    plot1 <- plot1 + 
      facet_grid(production_delay~strategy, scales = "free_x", drop = TRUE) +
      scale_y_continuous(breaks = y_breaks, 
                         labels = formatC(y_breaks / 1e5, digits = 2),
                         expand = expand_scale(mult = c(0, .1))) +
      xlab(xlab_string) +
      ylab("Deaths averted (x100,000)") +
      theme_bw() +
      theme(legend.position = "none",
            text = element_text(size = 12)) +
      scale_fill_manual(values = gg_color_hue(4)[2:3])
    
    if(bootstrap) {
      plot1 <- plot1 + geom_errorbar(aes(ymin = lower_deaths_averted, ymax = upper_deaths_averted))
    }
    
    if(save_output) {
      width <- 15
      height <- 15
      filename_temp <- paste0(save_dir, "Fig3_", production_delay_subset[1])
      ggsave(paste0(filename_temp, ".pdf"),
             plot1, width = width, height = height, units = "cm")
      png(paste0(filename_temp, ".png"),width=width,height=height, units = "cm", res = 300)
      plot(plot1)
      dev.off()
      saveRDS(plot1, file = paste0(filename_temp, ".rds"))
    }
    plot1
  }
  
  production_delay_subset <- list(c("stockpile", 7),
                                  c(90, 180))
  plot1 <- lapply(production_delay_subset, plot_wrapper)
  
  # save output as files
  if(save_output) {
    print_table(pars, filename = paste0(save_dir, "deaths_averted_top_n_countries.tex"), "Top n countries")
    print_table(pars, filename = paste0(save_dir, "deaths_averted_fn_pop_size.tex"), 
                "Function of population size")
  }
  
  list(table = pars, plot1 = plot1)
}

#' make Fig S1 in manuscript: death averted for vaccination allocation
#' accoring to incidence or corrent allocation, for different
#' seeding countries of the epidemic, and no production delay, with a stockpile
#' 
#' @inheritParams make_Fig2
#' 
#' @output a list of two objects:
#' pars is a data frame containing the median deaths and deaths averted under each
#' vaccination strategy and seed country.  if bootsrap = TRUE, also contains 95% CI.
#' plot1 is a ggplot object: visualisation of the deaths averted
make_FigS1 <- function(save_output = FALSE, bootstrap = TRUE) {
  # define vaccination strategies to plot
  strategy <- c("incidence",
                "curr_alloc")
  # define seed countries to plot
  seedCountries <- c("Belgium", "Sao_Tome_and_Principe", "Singapore", "Uganda")
  # make grid of all combinations of vaccination strategies and seed countries
  pars <- expand.grid(strategy = strategy, seedCountries = seedCountries)
  # define other parameters (used to construct filename from which to read the data)
  pars$stockpile_size <- 550e6
  pars$production_delay <- 0
  pars <- pars[,c("strategy", "production_delay", "stockpile_size", "seedCountries")]
  
  if(bootstrap) {
    # calculate median and 95% CI for the median deaths under each vaccination strategy
    # and seed country
    median_deaths <- apply_named_args(pars, 1, bootstrap_median_deaths_wrapper)
    # calculate median and 95% CI for the deaths averted under each vaccination strategy
    # and seed country
    deaths_averted <- apply_named_args(pars, 1, bootstrap_deaths_averted_wrapper)
    # format stuff
    rownames(median_deaths) <- c("lower_deaths", "median_deaths", "upper_deaths")
    rownames(deaths_averted) <- c("lower_deaths_averted", "deaths_averted", "upper_deaths_averted")
    pars <- cbind(pars, as.data.frame(t(median_deaths)), as.data.frame(t(deaths_averted)))
    levels(pars$strategy) <- c("Incidence", "Current allocation")
    pars$strategy <- factor(pars$strategy,
                            levels = c("Current allocation", "Incidence"))
  } else {
    pars <- rbind(pars, data.frame(strategy = "no_vaccination",
                                   production_delay = 0,
                                   stockpile_size = pars$stockpile_size[1],
                                   seedCountries = seedCountries))
    # read in median deaths under each vaccination strategy
    # and seed country, as well as median deaths under no vaccination
    pars$median_deaths <- Map_vapply(read_median_global_deaths, numeric(1),
                                     pars$strategy,
                                     pars$production_delay,
                                     pars$stockpile_size,
                                     seedCountries = pars$seedCountries)
    # calculate deaths averted = median deaths with no vaccination - median deaths with vaccination
    pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_deaths"] - pars$median_deaths
    
    # format stuff
    levels(pars$strategy) <- c("Incidence", "Current allocation", "No vaccination")
    pars$strategy <- factor(pars$strategy,
                            levels = c("Current allocation", "Incidence", "No vaccination"))
    
  }
  
  y_breaks <- seq(0, max(pars$deaths_averted) * 1.1, by = 2e5)
  
  # make plot
  
  plot1 <- ggplot(pars[pars$strategy != "No vaccination",], aes(x = strategy, y = deaths_averted)) +
    geom_bar(stat = "identity") +
    facet_wrap(~seedCountries, nrow = 1, drop = TRUE) +
    scale_y_continuous(breaks = y_breaks,
                       labels = formatC(y_breaks / 1e6, digits = 2),
                       expand = expand_scale(mult = c(0, .1))) +
    xlab("Vaccination strategy") +
    ylab("Deaths averted (million)") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  if(bootstrap) {
    plot1 <- plot1 + geom_errorbar(aes(ymin = lower_deaths_averted, ymax = upper_deaths_averted))
  }
  
  # save output as files
  if(save_output) {
    ggsave(paste0(save_dir, "deaths_averted_seedCountries_simple_allocation.pdf"),
           plot1, width = 20, height = 10, units = "cm")
    print_table(pars, filename = paste0(save_dir, "deaths_averted_seedCountries_simple_allocation.tex"))
  }
  
  list(table = pars, plot1 = plot1)
}

#' make Fig. S2 in the manuscript: comparing deaths averted for allocation strategies by population size,
#' for different seed countries and no production delay, with a stockpile
#' 
#' deaths averted = median deaths over m simulations of no vaccination,
#' minus median deaths over m simulations of vaccination
#' strategy 1 = vaccinate top n countries by population size,
#' with equal allocation to each of those countries by capita
#' strategy 2 = allocate vaccines per capita to each country
#' proportional to x^alpha where x in the population size and alpha is a parameter
#' 
#' @inheritParams make_Fig2
#' 
#' @output a list of three objects:
#' pars is a data frame containing the median deaths and deaths averted under each
#' vaccination strategy and seed country.  if bootsrap = TRUE, also contains 95% CI.
#' plot1 is a ggplot object: visualisation of the deaths averted
make_FigS2 <- function(save_output = FALSE, bootstrap = TRUE) {
  # values of n for strategy 1
  n_countries <- c(1, 2, 5, 10, 127)
  # values of alpha for strategy 2
  alpha <- c(.5, 1, 2, 3)
  
  stockpile_size <- 550e6
  # define seed countries to plot
  seedCountries <- c("Belgium", "Sao_Tome_and_Principe", "Singapore", "Uganda")
  # make grid of all combinations of vaccination strategies andseed countries
  pars <- expand.grid(strategy = "top_n_countries", 
                      production_delay = 0,
                      stockpile_size = stockpile_size,
                      free_param = n_countries,
                      seedCountries = seedCountries)
  pars <- rbind(pars, expand.grid(strategy = "fn_pop_size", 
                                  production_delay = 0,
                                  stockpile_size = stockpile_size,
                                  free_param = alpha,
                                  seedCountries = seedCountries))
  
  if(bootstrap) {
    # calculate median and 95% CI for the median deaths under each vaccination strategy
    # and seed country
    median_deaths <- apply_named_args(pars, 1, bootstrap_median_deaths_wrapper)
    # calculate median and 95% CI for the deaths averted under each vaccination strategy
    # and seed country
    deaths_averted <- apply_named_args(pars, 1, bootstrap_deaths_averted_wrapper)
    # format stuff
    rownames(median_deaths) <- c("lower_deaths", "median_deaths", "upper_deaths")
    rownames(deaths_averted) <- c("lower_deaths_averted", "deaths_averted", "upper_deaths_averted")
    pars <- cbind(pars, as.data.frame(t(median_deaths)), as.data.frame(t(deaths_averted)))
    levels(pars$strategy) <- c("Top n countries", "Function of population size")
  } else {
    # haven't written the non-bootstrap version for this yet
    stop("Not yet implemented")
    # pars <- rbind(pars, data.frame(strategy = "no_vaccination",
    #                                production_delay = 0,
    #                                stockpile_size = 0,
    #                                free_param = NA,
    #                                seedCountries = "China"))
    # pars$median_deaths <- Map_vapply(read_median_global_deaths, 
    #                                  numeric(1), 
    #                                  pars$strategy, 
    #                                  pars$production_delay,
    #                                  pars$stockpile_size,
    #                                  pars$free_param,
    #                                  pars$seedCountries)
    # pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_deaths"] - pars$median_deaths
    # 
    # levels(pars$strategy) <- c("Top n countries", "Function of population size", "No vaccination")
  }
  
  # make plot
  pars$free_param <- factor(pars$free_param)
  y_breaks <- seq(0, max(pars$deaths_averted) * 1.1, by = 2e5)
  plot1 <- ggplot(pars[pars$strategy != "No vaccination",], aes(x = free_param, y = deaths_averted, fill = strategy)) +
    geom_bar(stat = "identity") +
    facet_grid(seedCountries~strategy, scales = "free_x") +
    scale_y_continuous(breaks = y_breaks, 
                       labels = formatC(y_breaks / 1e6, digits = 2),
                       expand = expand_scale(mult = c(0, .1))) +
    xlab("n                                             alpha") +
    ylab("Deaths averted (million)") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 12))
  
  if(bootstrap) {
    plot1 <- plot1 + geom_errorbar(aes(ymin = lower_deaths_averted, ymax = upper_deaths_averted))
  }
  # save output as files
  if(save_output) {
    ggsave(paste0(save_dir, "deaths_averted_seedCountries_complex_allocation.pdf"),
           plot1, width = 15, height = 15, units = "cm")
    print_table(pars, filename = paste0(save_dir, "deaths_averted_seedCountries_top_n_countries.tex"), "Top n countries")
    print_table(pars, filename = paste0(save_dir, "deaths_averted_seedCountries_fn_pop_size.tex"), 
                "Function of population size")
  }
  
  list(table = pars, plot1 = plot1)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' make Fig. 2 in the manuscript: comparing deaths averted for allocation by incidence vs current allocation,
#' for different production delays
#' 
#' deaths averted = median deaths over n simulations of no vaccination,
#' minus median deaths over n simulations of vaccination
#' 
#' @param save_output logical: whether to save the output as a pdf as well as outputting the ggplot object
#' @param bootstrap logical: if TRUE, plot median and 95% CI for deaths averted, calculated using
#' bootstrap replicates of the stochastic simulations.  if FALSE, plot point
#' estimate of deaths averted only.
#' @output a list of three objects:
#' pars is a data frame containing the median deaths and deaths averted under each
#' vaccination strategy and production delay.  if bootsrap = TRUE, also contains 95% CI.
#' plot1 and plot2 are ggplot objects: two different visualisations of the deaths averted
make_Fig2_split_alloc <- function(save_output = FALSE, bootstrap = TRUE) {
  # define production delays to plot
  production_delay <- c(0, 7, 90, 180)
  alloc_split <- c("with", "without")
  # make grid of all combinations of vaccination strategies and production delays
  pars <- expand.grid(strategy = "curr_alloc", production_delay = production_delay,
                      alloc_split = alloc_split,
                      stringsAsFactors = FALSE)
  pars <- rbind(pars, data.frame(strategy = "no_vaccination", production_delay = 0,
                                 alloc_split = alloc_split,
                                 stringsAsFactors = FALSE))
  # define other parameters (used to construct filename from which to read the data)
  pars$stockpile_size <- 550e6
  pars[pars$production_delay > 0,"stockpile_size"] <- 0
  pars$seedCountries <- "China"
  pars$free_param <- 0
  if(bootstrap) {
    # calculate median and 95% CI for the median deaths under each vaccination strategy
    # and production delay
    median_deaths <- t(apply_named_args(pars, 1, bootstrap_median_deaths_wrapper)) %>%
      as.data.frame
    # calculate median and 95% CI for the deaths averted under each vaccination strategy
    # and production delay
    deaths_averted <- t(apply_named_args(pars, 1, bootstrap_deaths_averted_wrapper)) %>%
      as.data.frame
    # do.call(rbind, .)
    
    # format stuff
    colnames(median_deaths) <- c("lower_deaths", "median_deaths", "upper_deaths")
    colnames(deaths_averted) <- c("lower_deaths_averted", "deaths_averted", "upper_deaths_averted")
    pars <- cbind(pars, as.data.frame(median_deaths), as.data.frame(deaths_averted))
    pars$strategy <- factor(pars$strategy, 
                            levels = c("no_vaccination", "curr_alloc", "incidence"))
    levels(pars$strategy) <- c("No vaccination", "Current allocation", "Incidence")
  } else {
    pars <- rbind(pars, data.frame(strategy = "no_vaccination",
                                   production_delay = 0,
                                   stockpile_size = 0,
                                   seedCountries = pars$seedCountries[1],
                                   stringsAsFactors = FALSE))
    # read in median deaths under each vaccination strategy
    # and production delay, as well as median deaths under no vaccination
    pars$median_deaths <- Map_vapply(read_median_global_deaths, numeric(1), 
                                     pars$strategy, 
                                     pars$production_delay,
                                     pars$stockpile_size,
                                     seedCountries = pars$seedCountries)
    # calculate deaths averted = median deaths with no vaccination - median deaths with vaccination
    pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_deaths"] - pars$median_deaths
    # format stuff
    levels(pars$strategy) <- c("Incidence", "Current allocation", "No vaccination")
    pars$strategy <- factor(pars$strategy, 
                            levels = c("Current allocation", "Incidence", "No vaccination"))
    
  }
  pars$production_delay <- factor(pars$production_delay)
  levels(pars$production_delay) <- c("stockpile", production_delay[production_delay > 0])
  
  # make plots
  plot_wrapper <- function(production_delay_subset) {
    pars <- pars[pars$production_delay %in% production_delay_subset,]
    y_breaks <- seq(0, max(pars$deaths_averted) * 1.1, by = 1e5)
    plot1 <- ggplot(pars[pars$strategy != "No vaccination",], aes(x = production_delay, y = deaths_averted, fill = strategy)) +
      geom_bar(stat = "identity") +
      facet_wrap(~strategy, nrow = 1) +
      scale_y_continuous(breaks = y_breaks, 
                         labels = formatC(y_breaks / 1e5, digits = 2),
                         expand = expand_scale(mult = c(0, .1))) +
      xlab("Production delay (days)") +
      ylab("Deaths averted (x100,000)") +
      theme_bw() +
      theme(legend.position = "none",
            text = element_text(size = 12)) +
      scale_fill_manual(values = gg_color_hue(4)[c(1,4)])
    if(bootstrap) {
      plot1 <- plot1 + geom_errorbar(aes(ymin = lower_deaths_averted, ymax = upper_deaths_averted))
    }
    if(save_output) {
      width <- 20
      height <- 10
      filename_temp <- paste0(save_dir, "Fig2_", production_delay_subset[1])
      ggsave(paste0(filename_temp, ".pdf"),
             plot1, width = width, height = height, units = "cm")
      png(paste0(filename_temp, ".png"),width=width,height=height, units = "cm", res = 300)
      plot(plot1)
      dev.off()
      saveRDS(plot1, file = paste0(filename_temp, ".rds"))
    }
    plot1
  }
  production_delay_subset <- list(c("stockpile", 7),
                                  c(90, 180))
  plot1 <- lapply(production_delay_subset, plot_wrapper)
  # save output as files
  if(save_output) {
    print_table(pars, filename = paste0(save_dir, "Fig2.tex"))
  }
  
  list(table = pars, plot1 = plot1)
}

get_deaths_countries_without_alloc <- function(strategy, production_delay, stockpile_size, free_param, seedCountries) {
  dir_name <- make_dirname(strategy, production_delay, stockpile_size, free_param, seedCountries)
  filename <- paste0(dir_name, strategy, "_fixed_deaths.rds")
  deaths <- readRDS(filename)
  coverage <- read.csv("data/coverage_data_intersect.csv")
  deaths <- lapply(deaths,
                   function(x) x[coverage$dose_per_1000 == 0,] )
  deaths <- vapply(deaths, colSums, numeric(ncol(deaths[[1]])))
  deaths
}

plot_deaths_countries_without_alloc_stockpile <- function() {
  strategy <- c("no_vaccination", "curr_alloc")
  get_quantile_deaths_countries_without_alloc <- function(strategy) {
    deaths <- get_deaths_countries_without_alloc(strategy, 0, 550e6, 0, "China")%>%
      apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
      t %>%
      as.data.frame
    deaths$t <- as.numeric(rownames(deaths))
    deaths$strategy <- strategy
    deaths
  }
  deaths <- lapply(strategy, get_quantile_deaths_countries_without_alloc) %>%
    do.call(rbind, .)
  g <- ggplot(deaths, aes(x = t, group = strategy)) +
    geom_line(aes(y = `50%`, color = strategy)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = strategy), alpha = 0.3) +
    theme_bw() +
    xlab("Time (days)") +
    ylab("Cumulative deaths")
}