library(ggplot2)
library(magrittr)

postprocess <- function(n_countries) {
  paste0("outputs_vax_by_pop_size/by_pop_size_coverage_data_",
         n_countries, "_incidence.rds") %>%
    readRDS %>%
    vnapply(., sum) %>%
    data.frame(n_countries = n_countries, attack_rate = .)
  
}

n_countries <- c(seq_len(5), seq(10, 120, by = 10), 127)
if(FALSE) {
  attack_rate <- lapply(n_countries, postprocess) %>%
    do.call(rbind, .)
  
  # g <- ggplot(attack_rate, aes(x = n_countries, y = attack_rate, group = n_countries)) +
  # geom_boxplot()
  
  attack_rate_summary <- tapply(attack_rate$attack_rate, 
                                attack_rate$n_countries, 
                                quantile, 
                                probs = c(.025, .5, .975)) %>%
    do.call(rbind, .) %>%
    as.data.frame %>%
    cbind(data.frame(n_countries = n_countries))
  
  g <- ggplot(attack_rate_summary, aes(x = n_countries)) +
    geom_point(aes(y = `50%`)) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
    theme_bw() +
    ylab("Attack rate") +
    theme(text = element_text(size=20))
}


vaccinated_by_country <- function(filename) {
  sum_stat <- readRDS(filename)
  
  sum_stat <- lapply(sum_stat, function(x) apply(x, 1, function(x) x[length(x)])) %>%
    do.call(rbind, .)
  
  pop_size <- read.csv("data/demographic_data_intersect.csv", stringsAsFactors = FALSE)
  colnames(pop_size)[1] <- "Location"
  
  sum_stat <- apply(sum_stat, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
  sum_stat <- as.data.frame(t(sum_stat))
  colnames(sum_stat) <- c("lower", "median", "upper")
  sum_stat$Location <- rownames(sum_stat)
  sum_stat <- merge(sum_stat, pop_size[,c("Location", "N")])
  sum_stat$lower <- sum_stat$lower / sum_stat$N
  sum_stat$median <- sum_stat$median / sum_stat$N
  sum_stat$upper <- sum_stat$upper / sum_stat$N
  sum_stat <- merge(sum_stat, other_info$other_info$regionDat[,c("Location", "region")])%>%
    merge(., other_info$other_info$latitudeDat[,c("Location", "latitude")])
  
  sum_stat$Location <- factor(sum_stat$Location, levels = sum_stat$Location[order(sum_stat$latitude)])
  sum_stat
}
other_info <- readRDS("outputs/sim1000/pd7_vacbycurralloc_args_list.rds")
plot_vaccinated_by_country <- function(n_countries) {
  filename <- paste0("outputs_vax_by_pop_size/by_pop_size_coverage_data_",
                     n_countries, "_vaccinated.rds")
  my_df <- vaccinated_by_country(filename)
  g <- ggplot(my_df) +
    geom_errorbarh(aes(y=Location,x=median,xmax=upper,xmin=lower)) +
    geom_point(aes(y=Location,x=median),size=0.5) +
    coord_cartesian(xlim=c(0,1)) +
    xlab("Vaccination rate, median and 95% quantiles") +
    facet_grid(region~.,scales="free_y", space="free",switch="both") +
    theme(axis.text.y=element_text(size=6)) + theme_bw()
  png(paste0("outputs_vax_by_pop_size/n_vaccinated_by_country_", n_countries, ".png"),width=800,height=1200)
  plot(g)
  dev.off()
  g
}
if(FALSE) {
  g <- lapply(n_countries, plot_vaccinated_by_country)
}


postprocess <- function(alpha) {
  paste0("outputs_vax_fn_pop_size/fn_pop_size_coverage_data_",
         num2str(alpha), "_incidence.rds") %>%
    readRDS %>%
    vnapply(., sum) %>%
    data.frame(alpha = alpha, attack_rate = .)
  
}

alpha <- seq(0, 1, by = .1)
if(FALSE) {
  attack_rate <- lapply(alpha, postprocess) %>%
    do.call(rbind, .)
  
  # g <- ggplot(attack_rate, aes(x = n_countries, y = attack_rate, group = n_countries)) +
  # geom_boxplot()
  
  attack_rate_summary <- tapply(attack_rate$attack_rate, 
                                attack_rate$alpha, 
                                quantile, 
                                probs = c(.025, .5, .975)) %>%
    do.call(rbind, .) %>%
    as.data.frame %>%
    cbind(data.frame(alpha = alpha))
  
  g <- ggplot(attack_rate_summary, aes(x = alpha)) +
    geom_point(aes(y = `50%`)) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
    theme_bw() +
    ylab("Attack rate") +
    theme(text = element_text(size=20))
}

postprocess <- function(n_countries, country) {
  filename <- paste0("outputs_vax_by_pop_size_seed/fn_pop_size_country_attack_rates_",
                     country, 
                     "_coverage_data_",
                     n_countries,
                     "_all_runs.csv")

  attack_rate <- read.csv(filename)
  attack_rate$Location
  pop_size <- read.csv("data/demographic_data_intersect.csv", stringsAsFactors = FALSE)
  colnames(pop_size)[1] <- "Location"
  
  attack_rate <- merge(attack_rate, pop_size[,c("Location", "N")])

  run_cols <-  1 + seq_len(ncol(attack_rate) - 9)
  attack_rate[, run_cols] <- attack_rate[, run_cols] * attack_rate$N
  attack_rate <- attack_rate[, run_cols] %>%
    colSums %>%
    quantile(probs = c(0.025, 0.5, 0.975))
  c(alpha = alpha, attack_rate)
}

plot_attack_rate <- function(country) {
  n_countries <- c(seq_len(5), seq(10, 120, by = 10), 127)
  attack_rate <- lapply(n_countries, postprocess, country = country) %>%
    do.call(rbind, .) %>%
    as.data.frame
  
  g <- ggplot(attack_rate, aes(x = n_countries)) +
    geom_point(aes(y = `50%`)) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
    theme_bw() +
    ylab("Attack rate") +
    theme(text = element_text(size=20))
  ggsave(paste0("outputs_vax_by_pop_size_seed/", country, ".png"), g)
  g
}

postprocess <- function(alpha, country) {
  filename <- paste0("outputs_vax_fn_pop_size_seed/fn_pop_size_country_attack_rates_",
                     country, 
                     "_coverage_data_",
                     num2str(alpha),
                     "_all_runs.csv")
  
  attack_rate <- read.csv(filename)
  attack_rate$Location
  pop_size <- read.csv("data/demographic_data_intersect.csv", stringsAsFactors = FALSE)
  colnames(pop_size)[1] <- "Location"
  
  attack_rate <- merge(attack_rate, pop_size[,c("Location", "N")])
  
  run_cols <-  1 + seq_len(ncol(attack_rate) - 9)
  attack_rate[, run_cols] <- attack_rate[, run_cols] * attack_rate$N
  attack_rate <- attack_rate[, run_cols] %>%
    colSums %>%
    quantile(probs = c(0.025, 0.5, 0.975))
  c(n_countries = n_countries, attack_rate)
}

if(FALSE) {
  seedCountries <- c("Sao_Tome_and_Principe", "Belgium",
                     "Singapore", "Uganda")
  lapply(seedCountries, plot_attack_rate)
}

plot_attack_rate <- function(country) {
  alpha <- seq(0, 1, by = .1)
  attack_rate <- lapply(alpha, postprocess, country = country) %>%
    do.call(rbind, .) %>%
    as.data.frame

  g <- ggplot(attack_rate, aes(x = alpha)) +
    geom_point(aes(y = `50%`)) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
    theme_bw() +
    ylab("Attack rate") +
    theme(text = element_text(size=20))
  ggsave(paste0("outputs_vax_fn_pop_size_seed/", country, ".png"), g)
  g
}

if(FALSE) {
  seedCountries <- c("Sao_Tome_and_Principe", "Belgium",
                     "Singapore", "Uganda")
  lapply(seedCountries, plot_attack_rate)
}