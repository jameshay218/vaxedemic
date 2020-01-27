calc_peak_global_incidence <- function(filename) {
  incidence <- readRDS(filename)
  find_peak_time <- function(incidence) {
    colSums(incidence) %>%
      which.max
  }
  vnapply(incidence, find_peak_time) + 1
}

plot_global_incidence <- function(filename) {
  grey <- "#D8D8D8"
  incidence <- readRDS(filename)
  global_incidence <- lapply(incidence, colSums) %>%
    do.call(rbind, .)
  quantiles <- apply(global_incidence, 2, quantile, probs = c(0.025, 0.975))
  plot_df <- cbind(data.frame(t = seq_len(ncol(quantiles)) - 1),
                   as.data.frame(t(quantiles)))
  g <- ggplot(plot_df, aes(x = t, ymin = `2.5%`, ymax = `97.5%`)) +
    geom_ribbon(fill = grey) +
    theme_bw() +
    coord_cartesian(expand = FALSE) +
    xlab("Time (days)") +
    ylab("Incidence")
  g
}

plot_global_incidence2 <- function(filenames) {
  wrapper <- function(filename, intervention) {
    incidence <- readRDS(filename)
    global_incidence <- lapply(incidence, colSums) %>%
      do.call(rbind, .)
    quantiles <- apply(global_incidence, 2, quantile, probs = c(0.025, 0.975))
    
    plot_df <- cbind(data.frame(t = seq_len(ncol(quantiles)) - 1),
                     as.data.frame(t(quantiles)))
    plot_df$intervention <- intervention
    plot_df
  }
  plot_df <- Map(wrapper, filenames, c("no intervention", "intervention")) %>%
    do.call(rbind, .)

  g <- ggplot(plot_df, aes(x = t, ymin = `2.5%`, ymax = `97.5%`)) +
    geom_ribbon(aes(fill = intervention, group = intervention)) +
    theme_bw() +
    coord_cartesian(expand = FALSE) +
    xlab("Time (days)") +
    ylab("Incidence")
  g
}

calc_global_attack_rate <- function(filename) {
  incidence <- readRDS(filename)
  sum_global_incidence <- vnapply(incidence, sum)
  world_popn <- sum(setup_demography_real_data("")$popns)
  return(list(sum_global_incidence = sum_global_incidence,
         global_attack_rate = sum_global_incidence / world_popn))
}

if(FALSE) {
  source("scripts/run_decrease_between_country_travel.R")
  source("scripts/run_decrease_between_country_travel_intervention.R")
  source("scripts/run_decrease_between_country_travel_seasonality.R")
  source("scripts/run_decrease_between_country_travel_intervention_seasonality.R")
}
if(FALSE) {
  library(ggplot2)
  library(magrittr)
  postprocess_all_interventions <- function(filename) {
    global_peak_incidence <- calc_peak_global_incidence(filename) %>%
      quantile(probs = c(0.025, 0.5, 0.975))
    global_attack_rate <- calc_global_attack_rate(filename)$global_attack_rate %>%
      quantile(probs = c(0.025, 0.5, 0.975))
    print(global_peak_incidence)
    print(global_attack_rate)
    g <- plot_global_incidence(filename)
    ggsave(sub("rds", "pdf", filename), g)
  }
  # filenames <- c("decrease_between_country_travel/no_intervention_fixed_incidence.rds",
  #                "decrease_between_country_travel/intervention_fixed_incidence.rds",
  #                "decrease_between_country_travel/no_intervention_seasonality_fixed_incidence.rds",
  #                "decrease_between_country_travel/intervention_seasonality_fixed_incidence.rds")
  
  filenames <- c("decrease_between_country_travel/no_intervention_fixed_incidence.rds",
                 "decrease_between_country_travel/intervention_fixed_incidence.rds",
                 "decrease_between_country_travel/intervention_0point1_fixed_incidence.rds",
                 "decrease_between_country_travel/intervention_0point05_fixed_incidence.rds",
                 "decrease_between_country_travel/intervention_0point01_fixed_incidence.rds")
  # lapply(filenames, postprocess_all_interventions)
  concat_peak_time <- function(filename) {
    calc_peak_global_incidence(filename) %>%
      quantile(probs = c(0.025, 0.5, 0.975))
  }
  peak_time_table <- lapply(filenames, concat_peak_time) %>%
    do.call(rbind, .)
  rownames(peak_time_table) <- c(0, 50, 90, 95, 99)
}

if(FALSE) {
  source("scripts/run_intervention.R")
  lapply(c(.1, .05, .01), run_intervention)
}
