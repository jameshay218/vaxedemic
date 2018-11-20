library(ggplot2)
library(magrittr)

postprocess_incidence <- function(strategy, vax_day) {
  filename <- paste0("outputs/sim1000/pd", vax_day, "_vacby",
                      strategy, "_incidence.rds")
  inc <- readRDS(filename) %>%
    lapply(., function(x) apply(x, 1, sum)) %>%
    do.call(rbind, .)
  inc <- cbind(inc, matrix(apply(inc, 1, sum), ncol = 1, dimnames = list(NULL, "all")))
  attack_rate <- t(apply(inc, 1, function(x) x / pop_size)) %>%
    as.data.frame %>%
    cbind(., strategy, vax_day)
  colnames(attack_rate)[ncol(attack_rate) + c(-1, 0)] <- c("strategy", "vax_day")
  attack_rate
}

vax_day <- c(7, 84, 272)
strategy <- c("inc", "curralloc")
pars <- expand.grid(strategy = strategy, vax_day = vax_day)

pop_size <- read.csv("data/demographic_data_intersect.csv")$N
pop_size <- c(pop_size, sum(as.numeric(pop_size)))

attack_rate <- Map(postprocess_incidence, pars$strategy, pars$vax_day) %>%
  do.call(rbind, .) %>%
  as.data.frame

attack_rate$vax_day <- factor(attack_rate$vax_day, levels = sort(unique(attack_rate$vax_day)))
g <- ggplot(attack_rate, aes(x = vax_day, y = all)) + geom_boxplot() +
  facet_wrap(~strategy) +
  coord_cartesian(ylim = c(0,1), expand = FALSE) +
  xlab("Vax production starting day") +
  ylab("Global attack rate")
ggsave("outputs/sim1000/summary.pdf", g)
