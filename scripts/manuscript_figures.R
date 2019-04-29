make_dirname <- function(strategy, production_delay, stockpile_size, free_param, seedCountries) {
  if(seedCountries == "China") {
    dir_name <- paste0("outputs/deaths_only/pd", production_delay, strategy, 
                       "_stockpile", num2str(stockpile_size), "/")
  } else {
    dir_name <- paste0("outputs/deaths_only/pd", production_delay, strategy, 
                       "_stockpile", num2str(stockpile_size), "_", seedCountries, "/")
  }
  if(!dir.exists(dir_name)) {
    dir_name <- paste0("outputs/deaths_only/pd", production_delay, strategy, "/")
  }
  dir_name
}

read_median_global_deaths <- function(strategy, production_delay, stockpile_size, free_param, seedCountries) {

  dir_name <- make_dirname(strategy, production_delay, stockpile_size, free_param, seedCountries)
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

read_global_deaths <- function(strategy, production_delay, stockpile_size, free_param, seedCountries) {
  dir_name <- make_dirname(strategy, production_delay, stockpile_size, free_param, seedCountries)

  if(strategy %in% c("incidence", "curr_alloc", "no_vaccination")) {
    filename <- paste0(dir_name, strategy, "_fixed_global_deaths.rds")
  } else {
    filename <- paste0(dir_name, strategy, "_coverage_data_", 
                       num2str(free_param), "_global_deaths.rds")
  }
  if(file.exists(filename)) {
    readRDS(filename)
  } else {
    source_filename <- sub("_global", "", filename, fixed = TRUE)
    global_deaths <- calc_global_deaths(source_filename)
    saveRDS(global_deaths, filename)
    global_deaths
  }
}

bootstrap_deaths_averted <- function(strategies, production_delays, stockpile_sizes, free_params, seedCountries) {

  label_with_strategy <- function(global_deaths, strategy_label) {
    data.frame(deaths = global_deaths, 
                                strategy = strategy_label)
  }
  
  global_deaths <- Map(read_global_deaths, strategies, production_delays, stockpile_sizes, free_params, seedCountries) %>%
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

bootstrap_median_deaths <- function(strategy, production_delay, stockpile_size, free_param, seedCountries) {
  
  global_deaths <- read_global_deaths(strategy, production_delay, stockpile_size, free_param, seedCountries)
  
  my_median <- function(data, indices) {
    median(data[indices])
  }
  
  boot::boot(global_deaths, my_median, 500)
}

bootstrap_median_deaths_wrapper <- function(strategy, production_delay, stockpile_size = 0, free_param = 0, seedCountries) {
  production_delay <- as.numeric(production_delay)
  stockpile_size <- as.numeric(stockpile_size)
  free_param <- as.numeric(free_param)
  bootstrap_sol <- bootstrap_median_deaths(strategy, production_delay, stockpile_size, free_param, seedCountries)
  quantile(bootstrap_sol$t, probs = c(.025, .5, .975))
}

bootstrap_deaths_averted_wrapper <- function(strategy, production_delay, stockpile_size = 0, free_param = 0, seedCountries) {
  production_delay <- as.numeric(production_delay)
  stockpile_size <- as.numeric(stockpile_size)
  free_param <- as.numeric(free_param)
  strategies <- c("no_vaccination", strategy)
  production_delays <- c(0, production_delay)
  stockpile_sizes <- c(0, stockpile_size)
  free_params <- c(0, free_param)
  seedCountries <- rep(seedCountries, 2)
  bootstrap_sol <- bootstrap_deaths_averted(strategies, production_delays, stockpile_sizes, free_params, seedCountries)
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

make_Fig2 <- function(save_output = FALSE, bootstrap = FALSE) {
    strategy <- c("incidence",
                  "curr_alloc")
    production_delay <- c(0, 7, 90, 180)
    pars <- expand.grid(strategy = strategy, production_delay = production_delay)
    pars$stockpile_size <- 550e6
    pars[pars$production_delay > 0,"stockpile_size"] <- 0
    pars$seedCountries <- "China"
    
    if(bootstrap) {
      median_deaths <- apply_named_args(pars, 1, bootstrap_median_deaths_wrapper)
      deaths_averted <- apply_named_args(pars, 1, bootstrap_deaths_averted_wrapper)
      rownames(median_deaths) <- c("lower_deaths", "median_deaths", "upper_deaths")
      rownames(deaths_averted) <- c("lower_deaths_averted", "deaths_averted", "upper_deaths_averted")
      pars <- cbind(pars, as.data.frame(t(median_deaths)), as.data.frame(t(deaths_averted)))
      levels(pars$strategy) <- c("Incidence", "Current allocation")
      pars$strategy <- factor(pars$strategy, 
                              levels = c("Current allocation", "Incidence"))
    } else {
      pars <- rbind(pars, data.frame(strategy = "no_vaccination",
                                     production_delay = 0,
                                     stockpile_size = 0,
                                     seedCountries = "China"))
      pars$median_deaths <- Map_vapply(read_median_global_deaths, numeric(1), 
                                       pars$strategy, 
                                       pars$production_delay,
                                       pars$stockpile_size,
                                       seedCountries = pars$seedCountries)
      pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_deaths"] - pars$median_deaths
      
      levels(pars$strategy) <- c("Incidence", "Current allocation", "No vaccination")
      pars$strategy <- factor(pars$strategy, 
                              levels = c("Current allocation", "Incidence", "No vaccination"))

    }
    pars$production_delay <- factor(pars$production_delay)
    levels(pars$production_delay) <- c("stockpile", production_delay[production_delay > 0])
    y_breaks <- seq(0, max(pars$deaths_averted) * 1.1, by = 2e5)

    plot1 <- ggplot(pars[pars$strategy != "No vaccination",], aes(x = production_delay, y = deaths_averted, fill = strategy)) +
      geom_bar(stat = "identity") +
      facet_wrap(~strategy, nrow = 1) +
      scale_y_continuous(breaks = y_breaks, 
                         labels = formatC(y_breaks / 1e6, digits = 2),
                         expand = expand_scale(mult = c(0, .1))) +
      xlab("Production delay (days)") +
      ylab("Deaths averted (million)") +
      theme_bw() +
      theme(legend.position = "none",
            text = element_text(size = 12))

    plot2 <- ggplot(pars[pars$strategy != "No vaccination",], aes(x = production_delay, y = deaths_averted, fill = strategy)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_y_continuous(breaks = y_breaks, 
                         labels = formatC(y_breaks / 1e6, digits = 2),
                         expand = expand_scale(mult = c(0, .1))) +
      xlab("Production delay (days)") +
      ylab("Deaths averted (million)") +
      theme_bw() +
      theme(legend.position = c(.8,.8),
            text = element_text(size = 12))
    
    if(bootstrap) {
      plot1 <- plot1 + geom_errorbar(aes(ymin = lower_deaths_averted, ymax = upper_deaths_averted))
      plot2 <- plot2 + geom_errorbar(aes(ymin = lower_deaths_averted, ymax = upper_deaths_averted), position = "dodge")
    }
    
    if(save_output) {
      save_dir <- "~/overleaf/vaxedemic_manuscript/figs/"
      ggsave(paste0(save_dir, "deaths_averted_simple_allocation.pdf"),
             plot1, width = 20, height = 10, units = "cm")
      ggsave(paste0(save_dir, "deaths_averted_simple_allocation_alt.pdf"),
             plot2, width = 10, height = 10, units = "cm")
      print_table(pars, filename = paste0(save_dir, "deaths_averted_simple_allocation.tex"))
    }
    
    list(table = pars, plot1 = plot1, plot2 = plot2)
}

make_Fig3 <- function(save_output = FALSE, bootstrap = FALSE) {
  n_countries <- c(1, 2, 5, 10, 127)
  alpha <- c(.5, 1, 2, 3)

  production_delay <- c(0, 7, 90, 180)
  pars <- expand.grid(strategy = "top_n_countries", 
                      production_delay = production_delay,
                      free_param = n_countries)
  pars <- rbind(pars, expand.grid(strategy = "fn_pop_size", 
                                  production_delay = production_delay,
                                  free_param = alpha))
  pars$stockpile_size <- 550e6
  pars[pars$production_delay > 0,"stockpile_size"] <- 0
  pars$seedCountries <- "China"

  if(bootstrap) {
    median_deaths <- apply_named_args(pars, 1, bootstrap_median_deaths_wrapper)
    deaths_averted <- apply_named_args(pars, 1, bootstrap_deaths_averted_wrapper)
    rownames(median_deaths) <- c("lower_deaths", "median_deaths", "upper_deaths")
    rownames(deaths_averted) <- c("lower_deaths_averted", "deaths_averted", "upper_deaths_averted")
    pars <- cbind(pars, as.data.frame(t(median_deaths)), as.data.frame(t(deaths_averted)))
    levels(pars$strategy) <- c("Top n countries", "Function of population size")
  } else {
    pars <- rbind(pars, data.frame(strategy = "no_vaccination",
                                   production_delay = 0,
                                   stockpile_size = 0,
                                   free_param = NA,
                                   seedCountries = "China"))
    pars$median_deaths <- Map_vapply(read_median_global_deaths, 
                                            numeric(1), 
                                            pars$strategy, 
                                            pars$production_delay,
                                            pars$stockpile_size,
                                            pars$free_param,
                                     pars$seedCountries)
    pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_deaths"] - pars$median_deaths
    
    levels(pars$strategy) <- c("Top n countries", "Function of population size", "No vaccination")
  }

  pars$production_delay <- factor(pars$production_delay)
  levels(pars$production_delay) <- c("stockpile", production_delay[production_delay > 0])
  pars$free_param <- factor(pars$free_param)
  y_breaks <- seq(0, max(pars$deaths_averted) * 1.1, by = 2e5)
  plot1 <- ggplot(pars[pars$strategy != "No vaccination",], aes(x = free_param, y = deaths_averted, fill = strategy)) +
    geom_bar(stat = "identity") +
    facet_grid(production_delay~strategy, scales = "free_x") +
    scale_y_continuous(breaks = y_breaks, 
                       labels = formatC(y_breaks / 1e6, digits = 2),
                       expand = expand_scale(mult = c(0, .1))) +
    xlab("n                                             alpha") +
    ylab("Deaths averted (million)") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 12))
  plot2 <- ggplot(pars[pars$strategy != "No vaccination",], aes(x = free_param, y = deaths_averted, fill = production_delay, group = production_delay)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(~strategy, scales = "free_x") +
    scale_y_continuous(breaks = y_breaks, 
                       labels = formatC(y_breaks / 1e6, digits = 2),
                       expand = expand_scale(mult = c(0, .1))) +
    xlab("n                                            alpha") +
    ylab("Deaths averted (million)") +
    theme_bw() +
    theme(text = element_text(size = 12)) +
    scale_fill_manual(name = "Production delay (days)",
                      values = rev(RColorBrewer::brewer.pal(4,"Blues")))
  
  if(bootstrap) {
    plot1 <- plot1 + geom_errorbar(aes(ymin = lower_deaths_averted, ymax = upper_deaths_averted))
    plot2 <- plot2 + geom_errorbar(aes(ymin = lower_deaths_averted, ymax = upper_deaths_averted), position = "dodge")
  }
  
  if(save_output) {
    save_dir <- "~/overleaf/vaxedemic_manuscript/figs/"
    ggsave(paste0(save_dir, "deaths_averted_complex_allocation.pdf"),
           plot1, width = 15, height = 15, units = "cm")
    ggsave(paste0(save_dir, "deaths_averted_complex_allocation_alt.pdf"),
           plot2, width = 20, height = 10, units = "cm")
    print_table(pars, filename = paste0(save_dir, "deaths_averted_top_n_countries.tex"), "Top n countries")
    print_table(pars, filename = paste0(save_dir, "deaths_averted_fn_pop_size.tex"), 
                "Function of population size")
  }
  
  list(table = pars, plot1 = plot1, plot2 = plot2)
}

make_Fig2_seedCountries <- function(save_output = FALSE, bootstrap = FALSE) {
  strategy <- c("incidence",
                "curr_alloc")
  seedCountries <- c("Belgium", "Sao_Tome_and_Principe", "Singapore", "Uganda")
  pars <- expand.grid(strategy = strategy, seedCountries = seedCountries)
  pars$stockpile_size <- 550e6
  pars$production_delay <- 0
  pars <- pars[,c("strategy", "production_delay", "stockpile_size", "seedCountries")]

  if(bootstrap) {
    median_deaths <- apply_named_args(pars, 1, bootstrap_median_deaths_wrapper)
    deaths_averted <- apply_named_args(pars, 1, bootstrap_deaths_averted_wrapper)
    rownames(median_deaths) <- c("lower_deaths", "median_deaths", "upper_deaths")
    rownames(deaths_averted) <- c("lower_deaths_averted", "deaths_averted", "upper_deaths_averted")
    pars <- cbind(pars, as.data.frame(t(median_deaths)), as.data.frame(t(deaths_averted)))
    levels(pars$strategy) <- c("Incidence", "Current allocation")
    pars$strategy <- factor(pars$strategy,
                            levels = c("Current allocation", "Incidence"))
  } else {
    pars <- rbind(pars, data.frame(strategy = "no_vaccination",
                                   production_delay = 0,
                                   stockpile_size = 0,
                                   seedCountries = seedCountries))
    pars$median_deaths <- Map_vapply(read_median_global_deaths, numeric(1),
                                     pars$strategy,
                                     pars$production_delay,
                                     pars$stockpile_size,
                                     seedCountries = pars$seedCountries)
    pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_deaths"] - pars$median_deaths

    levels(pars$strategy) <- c("Incidence", "Current allocation", "No vaccination")
    pars$strategy <- factor(pars$strategy,
                            levels = c("Current allocation", "Incidence", "No vaccination"))

  }

  y_breaks <- seq(0, max(pars$deaths_averted) * 1.1, by = 2e5)

  plot1 <- ggplot(pars[pars$strategy != "No vaccination",], aes(x = strategy, y = deaths_averted)) +
    geom_bar(stat = "identity") +
    facet_wrap(~seedCountries, nrow = 1) +
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

  if(save_output) {
    save_dir <- "~/overleaf/vaxedemic_manuscript/figs/"
    ggsave(paste0(save_dir, "deaths_averted_seedCountries_simple_allocation.pdf"),
           plot1, width = 20, height = 10, units = "cm")
    print_table(pars, filename = paste0(save_dir, "deaths_averted_seedCountries_simple_allocation.tex"))
  }

  list(table = pars, plot1 = plot1)
}

make_Fig3_seedCountries <- function(save_output = FALSE, bootstrap = FALSE) {
  n_countries <- c(1, 2, 5, 10, 127)
  alpha <- c(.5, 1, 2, 3)
  
  stockpile_size <- 550e6
  
  seedCountries <- c("Belgium", "Sao_Tome_and_Principe", "Singapore", "Uganda")
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
    median_deaths <- apply_named_args(pars, 1, bootstrap_median_deaths_wrapper)
    deaths_averted <- apply_named_args(pars, 1, bootstrap_deaths_averted_wrapper)
    rownames(median_deaths) <- c("lower_deaths", "median_deaths", "upper_deaths")
    rownames(deaths_averted) <- c("lower_deaths_averted", "deaths_averted", "upper_deaths_averted")
    pars <- cbind(pars, as.data.frame(t(median_deaths)), as.data.frame(t(deaths_averted)))
    levels(pars$strategy) <- c("Top n countries", "Function of population size")
  } else {
    stop("Not yet implemented")
    pars <- rbind(pars, data.frame(strategy = "no_vaccination",
                                   production_delay = 0,
                                   stockpile_size = 0,
                                   free_param = NA,
                                   seedCountries = "China"))
    pars$median_deaths <- Map_vapply(read_median_global_deaths, 
                                     numeric(1), 
                                     pars$strategy, 
                                     pars$production_delay,
                                     pars$stockpile_size,
                                     pars$free_param,
                                     pars$seedCountries)
    pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_deaths"] - pars$median_deaths
    
    levels(pars$strategy) <- c("Top n countries", "Function of population size", "No vaccination")
  }
  
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
  
  if(save_output) {
    save_dir <- "~/overleaf/vaxedemic_manuscript/figs/"
    ggsave(paste0(save_dir, "deaths_averted_seedCountries_complex_allocation.pdf"),
           plot1, width = 15, height = 15, units = "cm")
    print_table(pars, filename = paste0(save_dir, "deaths_averted_seedCountries_top_n_countries.tex"), "Top n countries")
    print_table(pars, filename = paste0(save_dir, "deaths_averted_seedCountries_fn_pop_size.tex"), 
                "Function of population size")
  }
  
  list(table = pars, plot1 = plot1)
}