get_attack_rates_by_region <- function(filename) {
  demography_filename="data/demographic_data_intersect.csv"
  demographic_data <- read.csv(demography_filename,sep = ",")
  demographic_data[,ncol(demographic_data)] <- 1 - rowSums(demographic_data[,seq(3,ncol(demographic_data) - 1)])
  demographic_data <- demographic_data[order(demographic_data$countryID),]
  attack_rates_data <- read.csv(filename,
                                stringsAsFactors = FALSE)
  attack_rates_data$N <- demographic_data$N
  # retrieve numbers of infections from attack rates (attack rates are a proportion)
  attack_rates_data[,c("lower95", "median", "upper95")] <- 
    attack_rates_data[,c("lower95", "median", "upper95")] * attack_rates_data$N
  # sum infections over regions
  attack_rates_region <- apply(attack_rates_data[,c("lower95", "median", "upper95", "N")],
                               2,
                               function(x) tapply(x, attack_rates_data$region, sum)) %>%
    as.data.frame
  attack_rates_region[,c("lower95", "median", "upper95")] <- 
    attack_rates_region[,c("lower95", "median", "upper95")] / attack_rates_region$N
  attack_rates_region$region <- rownames(attack_rates_region)
  attack_rates_region
}

plot_attack_rates_region <- function(attack_rates_region) {
  ggplot(attack_rates_region) +
    geom_errorbarh(aes(y=region,xmax=upper95,xmin=lower95)) +
    geom_point(aes(y=region,x=median),size=0.5) +
    scale_x_continuous(limits=c(0,1)) +
    xlab("Attack rate, median and 95% quantiles") +
    theme(axis.text.y=element_text(size=6)) + theme_bw()
}

plot_attack_rates_region_strategy <- function(attack_rates_region, 
                                              facet_strategy_flag = FALSE,
                                               facet_flip_flag = FALSE, 
                                               coord_flip_flag = FALSE,
                                              save_plot = FALSE,
                                              delay) {
  g <- ggplot(attack_rates_region)
  if(facet_strategy_flag) {
    g <- g + geom_errorbar(aes(x=region,ymax=upper95,ymin=lower95, color = region)) +
      geom_point(aes(x=region,y=median, color = region)) #+
      # scale_color_manual(values = c("black", gg_color_hue(4)[c(1,4)]))
    if(facet_flip_flag) {
      g <- g + facet_grid(strategy~., switch = "both")
    } else {
      g <- g + facet_grid(.~strategy)
    }
  } else {
    g <- g + geom_errorbar(aes(x=strategy,ymax=upper95,ymin=lower95)) +
      geom_point(aes(x=strategy,y=median),size=0.5)
    if(facet_flip_flag) {
      g <- g + facet_grid(region~., switch = "both")
    } else {
      g <- g + facet_grid(.~region)
    }
  }

  if(coord_flip_flag) {
    g <- g + coord_flip()
  }
   g <- g + theme_bw()
   if(!coord_flip_flag) {
     g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
   }
    g <- g + scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
     ylab("Attack rate, median and 95% quantiles") +
      xlab("") +
      theme(legend.position = "none")
    if(save_plot) {
      filename <- paste0("plots/attack_rates_region_strategy",
                         delay,
                         facet_strategy_flag,
                         facet_flip_flag,
                         coord_flip_flag,
                         ".pdf")
      ggsave(filename, g, width = 20, height = 13, units = "cm")
    }
  g
}

read_and_plot_attack_rates_region_all <- function(delay) {
  filenames <- c("No vaccination" = "outputs_delayed_protection/pd0no_vaccination_stockpile0_China/no_vaccination_country_attack_rates_fixed_median_ci.csv",
                 "Current allocation" = paste0("outputs_delayed_protection/pd", delay, "curr_alloc_stockpile0_China/curr_alloc_country_attack_rates_fixed_median_ci.csv"),
                 "Incidence" = paste0("outputs_delayed_protection/pd", delay, "incidence_stockpile0_China/incidence_country_attack_rates_fixed_median_ci.csv"))
  attack_rates_region_all <- lapply(filenames, get_attack_rates_by_region) %>%
    do.call(rbind, .)
  strategies <- names(filenames)
  attack_rates_region_all$strategy <- factor(rep(strategies, each = nrow(attack_rates_region_all) / length(strategies)),
                                             levels = strategies)
  # plotting_options <- expand.grid(facet_strategy_flag = c(FALSE, TRUE),
  #                                 facet_flip_flag = c(FALSE, TRUE),
  #                                 coord_flip_flag = c(FALSE, TRUE))
  # attack_rates_region_all_plot <- Map(function(x, y, z)plot_attack_rates_region_strategy(attack_rates_region_all,
  #                                                                                     x, y, z, TRUE, delay),
  #                                     plotting_options$facet_strategy_flag,
  #                                     plotting_options$facet_flip_flag,
  #                                     plotting_options$coord_flip_flag)
  
  attack_rates_region_all_plot <- plot_attack_rates_region_strategy(attack_rates_region_all,
                                                                                         TRUE, FALSE, FALSE, FALSE, delay)
  save_dir <- "~/overleaf/vaxedemic_manuscript/figs_delayed_protection/"
  filename <- paste0(save_dir, "Fig1")
  saveRDS(attack_rates_region_all_plot,
          file = paste0(filename, ".rds"))
  width <- 20
  height <- 13
  png(paste0(filename, ".png"),width=width,height=height, units = "cm", res = 300)
  plot(attack_rates_region_all_plot)
  dev.off()
  ggsave(paste0(filename, ".pdf"), attack_rates_region_all_plot,
         width = width, height = height, units = "cm")
  
  attack_rates_region_all_plot
}