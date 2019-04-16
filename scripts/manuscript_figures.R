read_median_global_deaths <- function(strategy, production_delay, free_param) {
  dir_name <- paste0("outputs/pd", production_delay, strategy, "/")
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

read_global_deaths <- function(strategy, production_delay, free_param) {
  dir_name <- paste0("outputs/deaths_only/pd", production_delay, strategy, "/")
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

bootstrap_deaths_averted <- function(strategies, production_delays, free_params) {

  label_with_strategy <- function(global_deaths, strategy_label) {
    data.frame(deaths = global_deaths, 
                                strategy = strategy_label)
  }
  
  global_deaths <- Map(read_global_deaths, strategies, production_delays, free_params) %>%
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

make_Fig2 <- function(save_output = FALSE, bootstrap = FALSE) {
    strategy <- c("incidence",
                  "curr_alloc")
    production_delay <- c(7, 90, 180)
    pars <- expand.grid(strategy = strategy, production_delay = production_delay)
    
    if(bootstrap) {
      bootstrap_deaths_averted_wrapper <- function(strategy, production_delay) {
        production_delay <- as.numeric(production_delay)
        strategies <- c("no_vaccination", strategy)
        production_delays <- c(0, production_delay)
        free_params <- c(0, 0)
        bootstrap_sol <- bootstrap_deaths_averted(strategies, production_delays, free_params)
        quantile(bootstrap_sol$t, probs = c(.025, .5, .975))
      }
      deaths_averted <- apply_named_args(pars, 1, bootstrap_deaths_averted_wrapper)
      rownames(deaths_averted) <- c("lower", "deaths_averted", "upper")
      pars <- cbind(pars, as.data.frame(t(deaths_averted)))
      levels(pars$strategy) <- c("Incidence", "Current allocation")
      pars$strategy <- factor(pars$strategy, 
                              levels = c("Current allocation", "Incidence"))
    } else {
      pars <- rbind(pars, data.frame(strategy = "no_vaccination",
                                     production_delay = 0))
      pars$median_global_deaths <- Map_vapply(read_median_global_deaths, numeric(1), pars$strategy, pars$production_delay)
      pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_global_deaths"] - pars$median_global_deaths
      
      pars <- pars[pars$strategy != "no_vaccination",]
      levels(pars$strategy) <- c("Incidence", "Current allocation", "No vaccination")
      pars$strategy <- factor(pars$strategy, 
                              levels = c("Current allocation", "Incidence", "No vaccination"))

    }
    pars$production_delay <- factor(pars$production_delay)
    y_breaks <- seq(0, max(pars$deaths_averted) * 1.1, by = 2e5)

    plot1 <- ggplot(pars, aes(x = production_delay, y = deaths_averted, fill = strategy)) +
      geom_bar(stat = "identity") +
      facet_wrap(~strategy, nrow = 1) +
      scale_y_continuous(breaks = y_breaks, 
                         labels = scientific_10x(y_breaks, digits = 0),
                         expand = expand_scale(mult = c(0, .1))) +
      xlab("Production delay (days)") +
      ylab("Deaths averted") +
      theme_bw() +
      theme(legend.position = "none",
            text = element_text(size = 12))

    plot2 <- ggplot(pars, aes(x = production_delay, y = deaths_averted, fill = strategy)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_y_continuous(breaks = y_breaks, 
                         labels = scientific_10x(y_breaks, digits = 0),
                         expand = expand_scale(mult = c(0, .1))) +
      xlab("Production delay (days)") +
      ylab("Deaths averted") +
      theme_bw() +
      theme(legend.position = c(.8,.8),
            text = element_text(size = 12))
    
    if(bootstrap) {
      plot1 <- plot1 + geom_errorbar(aes(ymin = lower, ymax = upper))
      plot2 <- plot2 + geom_errorbar(aes(ymin = lower, ymax = upper), position = "dodge")
    }
    
    if(save_output) {
      save_dir <- "~/overleaf/vaxedemic_manuscript/figs/"
      ggsave(paste0(save_dir, "deaths_averted_simple_allocation.pdf"), 
             plot1, width = 20, height = 10, units = "cm")
      ggsave(paste0(save_dir, "deaths_averted_simple_allocation_alt.pdf"), 
             plot2, width = 10, height = 10, units = "cm")
      print(xtable::xtable(pars),
            sanitize.text.function = identity,
            file = paste0(save_dir, "deaths_averted_simple_allocation.tex"))
    }
    
    list(table = pars, plot1 = plot1, plot2 = plot2)
}

make_Fig3 <- function(save_output, bootstrap = FALSE) {
  n_countries <- c(1, 2, 5, 10, 127)
  alpha <- c(.5, 1, 2, 3)

  production_delay <- c(7, 90, 180)
  pars <- expand.grid(strategy = "top_n_countries", 
                                      production_delay = production_delay,
                                      free_param = n_countries)
  pars <- rbind(pars, expand.grid(strategy = "fn_pop_size", 
                      production_delay = production_delay,
                      free_param = alpha))
  
  if(bootstrap) {
    bootstrap_deaths_averted_wrapper <- function(strategy, production_delay, free_param) {
      production_delay <- as.numeric(production_delay)
      free_param <- as.numeric(free_param)
      strategies <- c("no_vaccination", strategy)
      production_delays <- c(0, production_delay)
      free_params <- c(0, free_param)
      bootstrap_sol <- bootstrap_deaths_averted(strategies, production_delays, free_params)
      quantile(bootstrap_sol$t, probs = c(.025, .5, .975))
    }
    deaths_averted <- apply_named_args(pars, 1, bootstrap_deaths_averted_wrapper)
    rownames(deaths_averted) <- c("lower", "deaths_averted", "upper")
    pars <- cbind(pars, as.data.frame(t(deaths_averted)))
    levels(pars$strategy) <- c("Top n countries", "Function of population size")
  } else {
    pars <- rbind(pars, data.frame(strategy = "no_vaccination",
                                   production_delay = 0,
                                   free_param = NA))
    pars$median_global_deaths <- Map_vapply(read_median_global_deaths, 
                                            numeric(1), 
                                            pars$strategy, 
                                            pars$production_delay,
                                            pars$free_param)
    pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_global_deaths"] - pars$median_global_deaths
    
    pars <- pars[pars$strategy != "no_vaccination",]
    levels(pars$strategy) <- c("Top n countries", "Function of population size", "No vaccination")
  }

  pars$production_delay <- factor(pars$production_delay)
  pars$free_param <- factor(pars$free_param)
  y_breaks <- seq(0, max(pars$deaths_averted) * 1.1, by = 2e5)
  plot1 <- ggplot(pars, aes(x = free_param, y = deaths_averted, fill = strategy)) +
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
  plot2 <- ggplot(pars, aes(x = free_param, y = deaths_averted, fill = production_delay, group = production_delay)) +
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
                      values = rev(RColorBrewer::brewer.pal(3,"Blues")))
  
  if(bootstrap) {
    plot1 <- plot1 + geom_errorbar(aes(ymin = lower, ymax = upper))
    plot2 <- plot2 + geom_errorbar(aes(ymin = lower, ymax = upper), position = "dodge")
  }
  
  if(save_output) {
    save_dir <- "~/overleaf/vaxedemic_manuscript/figs/"
    ggsave(paste0(save_dir, "deaths_averted_complex_allocation.pdf"), 
           plot1, width = 15, height = 15, units = "cm")
    ggsave(paste0(save_dir, "deaths_averted_complex_allocation_alt.pdf"), 
           plot2, width = 20, height = 10, units = "cm")
    print(xtable::xtable(pars),
          sanitize.text.function = identity,
          file = paste0(save_dir, "deaths_averted_complex_allocation.tex"))
  }
  
  list(table = pars, plot1 = plot1, plot2 = plot2)
}