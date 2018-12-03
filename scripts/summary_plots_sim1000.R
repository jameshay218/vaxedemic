library(ggplot2)
library(magrittr)

postprocess <- function(strategy, vax_day, compartment = "incidence", proportion_produced = FALSE) {
  
  make_filename <- function(compartment) {
    paste0("outputs/sim1000/pd", vax_day, "_vacby",
           strategy, "_", compartment, ".rds")
  }
  
  filename <- make_filename(compartment)
  summary_func <- ifelse(compartment == "incidence", sum, function(x) x[length(x)])
  sum_stat <- readRDS(filename)
  
  if(proportion_produced) {
    find_vax_end_time <- function(vax_by_country) {
      global_cumul_vax <- apply(vax_by_country, 2, sum) %>%
        as.numeric
      
      last_consecutive_zeros <- rle(global_cumul_vax)
      last_consecutive_zeros <- last_consecutive_zeros$length
      last_consecutive_zeros <- last_consecutive_zeros[length(last_consecutive_zeros)]
      length(global_cumul_vax) - last_consecutive_zeros
    }
    
    calc_n_vax_produced <- function(vax_end_time, vax_day) {
      pmax((vax_end_time - vax_day) * other_info$vax_production_params$production_rate, 0)
    }
    
    n_vax_produced <- vapply(sum_stat, find_vax_end_time, numeric(1)) %>%
      calc_n_vax_produced(vax_day = vax_day)
  }
  sum_stat <- lapply(sum_stat, function(x) apply(x, 1, summary_func)) %>%
    do.call(rbind, .)
  sum_stat <- cbind(sum_stat, matrix(apply(sum_stat, 1, sum), ncol = 1, dimnames = list(NULL, "all")))

  if(proportion_produced) {
    proportion <- apply(sum_stat, 2, function(x) x / n_vax_produced)
  } else {
    pop_size <- read.csv("data/demographic_data_intersect.csv")$N
    pop_size <- c(pop_size, sum(as.numeric(pop_size)))
    proportion <- t(apply(sum_stat, 1, function(x) x / pop_size))
  }
  proportion <-  as.data.frame(proportion) %>%
    cbind(., strategy, vax_day)

  colnames(proportion)[ncol(proportion) + c(-1, 0)] <- c("strategy", "vax_day")
  proportion
}

vax_day <- c(7, 84, 272)
strategy <- c("inc", "curralloc")

other_info <- readRDS("outputs/sim1000/pd7_vacbycurralloc_args_list.rds")
pars <- expand.grid(strategy = strategy, vax_day = vax_day)

# attack_rate <- Map(postprocess, pars$strategy, pars$vax_day) %>%
#   do.call(rbind, .) %>%
#   as.data.frame
# 
# attack_rate$vax_day <- factor(attack_rate$vax_day, levels = vax_day)
# g <- ggplot(attack_rate, aes(x = vax_day, y = all)) + geom_boxplot() +
#   facet_wrap(~strategy) +
#   coord_cartesian(ylim = c(0,1), expand = FALSE) +
#   xlab("Vax production starting day") +
#   ylab("Global attack rate")
# ggsave("outputs/sim1000/summary.pdf", g)
# 
# vaccinated <- Map(postprocess, pars$strategy, pars$vax_day, compartment = "vaccinated") %>%
#   do.call(rbind, .) %>%
#   as.data.frame
# 
# vaccinated$vax_day <- factor(vaccinated$vax_day, levels = vax_day)
# g <- ggplot(vaccinated, aes(x = vax_day, y = all)) + geom_boxplot() +
#   facet_wrap(~strategy) +
#   coord_cartesian(ylim = c(0,1), expand = FALSE) +
#   xlab("Vax production starting day") +
#   ylab("Global vaccinated")
# ggsave("outputs/sim1000/vaccinated.pdf", g)

# vaccinated <- Map(postprocess, pars$strategy, pars$vax_day, 
#                   compartment = "vaccinated", proportion_produced = TRUE) %>%
#   do.call(rbind, .) %>%
#   as.data.frame

# vaccinated$vax_day <- factor(vaccinated$vax_day, levels = vax_day)
# g <- ggplot(vaccinated, aes(x = vax_day, y = all)) + geom_boxplot() +
#   facet_wrap(~strategy) +
#   coord_cartesian(ylim = c(0,1), expand = FALSE) +
#   xlab("Vax production starting day") +
#   ylab("Proportion of produced vaccines which are allocated")
# ggsave("outputs/sim1000/proportion_allocated.pdf", g)
# 
# incidence <- Map(postprocess, pars$strategy, pars$vax_day, 
#                  compartment = "incidence") %>%
#   do.call(rbind, .) %>%
#   as.data.frame
# 
# prop_vax_vs_incidence <- cbind(incidence[,"all", drop = FALSE], vaccinated[,c("all", "vax_day", "strategy")])
# colnames(prop_vax_vs_incidence)[1:2] <- c("attack_rate", "prop_vax")
# prop_vax_vs_incidence$vax_day_strategy <- paste0(prop_vax_vs_incidence$strategy, prop_vax_vs_incidence$vax_day)
# prop_vax_vs_incidence$vax_day_strategy <- 
#   factor(prop_vax_vs_incidence$vax_day_strategy, 
#          levels = c(paste0("curralloc", vax_day), paste0("inc", vax_day)))
# g <- ggplot(prop_vax_vs_incidence, aes(x = prop_vax, y = attack_rate)) + 
#   geom_point() + 
#   facet_wrap(~vax_day_strategy) + 
#   coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

vaccinated_by_country <- function(strategy, vax_day) {
  
  make_filename <- function(compartment) {
    paste0("outputs/sim1000/pd", vax_day, "_vacby",
           strategy, "_", "vaccinated", ".rds")
  }
  
  filename <- make_filename(compartment)
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

plot_vaccinated_by_country <- function(strategy, vax_day) {
  my_df <- vaccinated_by_country(strategy, vax_day)
  g <- ggplot(my_df) +
    geom_errorbarh(aes(y=Location,x=median,xmax=upper,xmin=lower)) +
    geom_point(aes(y=Location,x=median),size=0.5) +
    coord_cartesian(xlim=c(0,1)) +
    xlab("Attack rate, median and 95% quantiles") +
    facet_grid(region~.,scales="free_y", space="free",switch="both") +
    theme(axis.text.y=element_text(size=6)) + theme_bw()
  png(paste0("outputs/sim1000/", strategy, vax_day, "_vax_by_country.png"),width=800,height=1200)
  plot(g)
  dev.off()
  my_df
}
# Map(plot_vaccinated_by_country, pars$strategy, pars$vax_day)

attack_rate_by_country <- function(strategy, vax_day) {
  
  make_filename <- function(compartment) {
    paste0("outputs/sim1000/pd", vax_day, "_vacby",
           strategy, "_", "incidence", ".rds")
  }
  
  filename <- make_filename(compartment)
  sum_stat <- readRDS(filename)
  
  sum_stat <- lapply(sum_stat, function(x) apply(x, 1, sum)) %>%
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
  
  sum_stat$Location <- factor(sum_stat$Location, levels = sum_stat$Location[order(sum_stat$N)])
  sum_stat
}

plot_attack_rate_by_country <- function(strategy, vax_day) {
  my_df <- attack_rate_by_country(strategy, vax_day)
  g <- ggplot(my_df) +
    geom_errorbarh(aes(y=Location,x=median,xmax=upper,xmin=lower)) +
    geom_point(aes(y=Location,x=median),size=0.5) +
    coord_cartesian(xlim=c(0,1)) +
    xlab("Attack rate, median and 95% quantiles") +
    # facet_grid(region~.,scales="free_y", space="free",switch="both") +
    theme(axis.text.y=element_text(size=6)) + theme_bw()
  png(paste0("outputs/sim1000/", strategy, vax_day, "_attack_by_country_pop_size2.png"),width=800,height=1200)
  plot(g)
  dev.off()
  my_df
}

Map(plot_attack_rate_by_country, pars$strategy, pars$vax_day)