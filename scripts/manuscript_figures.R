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

make_Fig2 <- function(save_output = FALSE) {
    strategy <- c("incidence",
                  "curr_alloc")
    production_delay <- c(7, 90, 180)
    pars <- expand.grid(strategy = strategy, production_delay = production_delay)
    pars <- rbind(pars, data.frame(strategy = "no_vaccination",
                                   production_delay = 0))
    pars$median_global_deaths <- Map_vapply(read_median_global_deaths, numeric(1), pars$strategy, pars$production_delay)
    pars$deaths_averted <- pars[pars$strategy == "no_vaccination", "median_global_deaths"] - pars$median_global_deaths
    
    pars <- pars[pars$strategy != "no_vaccination",]
    levels(pars$strategy) <- c("Incidence", "Current allocation", "No vaccination")
    pars$strategy <- factor(pars$strategy, 
                            levels = c("Current allocation", "Incidence", "No vaccination"))
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

make_Fig3 <- function(save_output) {
  n_countries <- c(1, 2, 5, 10, 127)
  alpha <- c(.5, 1, 2, 3)

  production_delay <- c(7, 90, 180)
  pars <- expand.grid(strategy = "top_n_countries", 
                                      production_delay = production_delay,
                                      free_param = n_countries)
  pars <- rbind(pars, expand.grid(strategy = "fn_pop_size", 
                      production_delay = production_delay,
                      free_param = alpha))
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