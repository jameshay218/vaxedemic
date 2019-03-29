library(ggplot2)
library(magrittr)
filenames <- c("outputs_no_vax_seed/no_vax_country_attack_rates_China_all_runs.csv",
               "outputs/sim1000/pd7_vacbycurralloc_country_attack_rates_fixed_all_runs.csv")

strategy <- c("no_vax", "curr_alloc")
demog_filename <- "data/demographic_data_intersect.csv"
pop_size <- read.csv(demog_filename)
pop_size <- pop_size[order(pop_size$countryID), "N"]

postprocess <- function(filename, strategy) {
  attack_rate <- read.csv(filename)
  countries <- attack_rate$Location
  attack_rate <- attack_rate[,seq(2, 501)]
  total_infected <- apply(attack_rate, 2, function(x) sum(x * as.numeric(pop_size)))
  attack_rate <- total_infected / sum(as.numeric(pop_size))
  data.frame(attack_rate = attack_rate, strategy = strategy)
}

attack_rate <- Map(postprocess, filenames, strategy) %>%
  do.call(rbind, .)

ggplot(attack_rate, aes(strategy, attack_rate)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Attack rate") +
  theme(text = element_text(size=20))
