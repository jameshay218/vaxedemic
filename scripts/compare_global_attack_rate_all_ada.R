# library(ggplot2)
# library(magrittr)
# library(data.table)
# library(plyr)
# devtools::load_all()

make_prefix <- function(production_delay) {
  if(production_delay != 7) {
    paste0("outputs_pd", production_delay, "_")
  } else {
    "outputs_"
  }
}

load_GAR <- function(run_name, production_delay = 0) {
  output_prefix <- make_prefix(production_delay)
  if(run_name == "no_vax") {
    dir_name <- "outputs_no_vax/"
  } else {
    dir_name <- paste0(output_prefix, run_name, "/")
  }
  filename <- paste0(dir_name, "GAR.rds")
  if(!file.exists(filename)) {
    compile_GAR(run_name, production_delay = production_delay)
  } else {
    readRDS(filename)
  }
}

compile_GAR <- function(run_name, production_delay = 0) {
  pop_size_table <- read.csv("~/Documents/vaxedemic/data/demographic_data_intersect.csv")
  world_pop_size <- sum(as.numeric(pop_size_table$N))
  output_prefix <- make_prefix(production_delay)
  
  if(run_name == "top_n_countries") {
    dir_name <- paste0(output_prefix, run_name, "/")
    postprocess <- function(n_countries) {
      GAR <- paste0(dir_name, "by_pop_size_coverage_data_",
                    n_countries, "_incidence.rds") %>%
        readRDS %>%
        vnapply(., sum)
      data.frame(all = GAR / world_pop_size,
                 table_no = paste0("top_", n_countries),
                 run_name = run_name)
    }
    
    if(production_delay == 7) {
      n_countries <- c(seq_len(5), seq(10, 120, by = 10), 127)
    } else {
      n_countries <- c(1, 2, 5, 10, 127)
    }
    GAR <- lapply(n_countries, postprocess) %>%
      do.call(rbind, .)
  } else if (run_name == "fn_pop_size") {
    dir_name <- paste0(output_prefix, run_name, "/")
    postprocess <- function(alpha) {
      GAR <- paste0(dir_name, "fn_pop_size_coverage_data_",
                    num2str(alpha), "_incidence.rds") %>%
        readRDS %>%
        vnapply(., sum)
      data.frame(all = GAR / world_pop_size,
                 table_no = paste0("alpha_", alpha),
                 run_name = run_name)
    }
    
    if(production_delay == 7) {
      alpha <- c(seq(0, 1, by = .1), seq(1.5, 5, by = .5))
    } else {
      alpha <- c(.5, 1, 2, 3)
    }
    
    GAR <- lapply(alpha, postprocess) %>%
      do.call(rbind, .)
  } else if (run_name == "random") {
    pop_size_table <- pop_size_table[order(pop_size_table$countryID),]
    pop_size <- pop_size_table$N
    pop_size <- c(pop_size, sum(as.numeric(pop_size)))
    
    dir_name <- paste0(output_prefix, run_name, "/")
    completed_runs <- c("maximum_vacc_fixed_incidence.rds",
                        paste0("random_vaccinations_coverage_table_", seq_len(103), "_incidence.rds"))
    
    postprocess_incidence <- function(filename, table_no) {
      print(table_no)
      inc <- readRDS(paste0(dir_name, "/", filename)) %>%
        lapply(., function(x) apply(x, 1, sum)) %>%
        do.call(rbind, .)
      inc <- cbind(inc, matrix(apply(inc, 1, sum), ncol = 1, dimnames = list(NULL, "all")))
      attack_rate <- t(apply(inc, 1, function(x) x / pop_size)) %>%
        as.data.frame %>%
        cbind(., table_no, filename)
      attack_rate
    }
    
    file_labels <- data.frame(filename=completed_runs)
    
    labels <- seq(-1,length(completed_runs)-2,by=1)
    labels[1] <- "maxed_vaccination"
    labels[which(completed_runs == "random_vaccinations_coverage_table_103_incidence.rds")] <- "proportional_to_popn"
    labels[which(completed_runs == "random_vaccinations_coverage_table_102_incidence.rds")] <- "current_allocation"
    labels[which(completed_runs == "random_vaccinations_coverage_table_1_incidence.rds")] <- "no_vaccination"
    
    attack_rate <- Map(postprocess_incidence, completed_runs, labels) %>%
      do.call(rbind, .) %>%
      as.data.frame
    attack_rate$samp <- 1:nrow(attack_rate)
    GAR <- attack_rate[,c("all", "table_no")]
    GAR$run_name <- run_name
    rownames(GAR) <- NULL
    GAR$table_no <- as.factor(GAR$table_no)
    median_ar <- plyr::ddply(GAR, ~table_no, function(x) median(x[,"all"]))

    GAR$table_no <- factor(GAR$table_no, levels=median_ar$table_no[order(median_ar$V1)])
    
    # remove case with ridiculously high vaccination and proportional_to_popn (as repeated)
    GAR <- GAR[!(GAR$table_no %in% c("maxed_vaccination", "proportional_to_popn")),]
    
    # separate off special cases
    special_cases <- "current_allocation"
    GAR[GAR$table_no %in% special_cases, "run_name"] <- "special"
    GAR[GAR$table_no == "no_vaccination", "run_name"] <- "no_vaccination"
  } else if(run_name == "incidence"){
    dir_name <- paste0(output_prefix, run_name, "/")
    if(production_delay == 7) {
      filename <- paste0(dir_name, "pd7_vacbyinc_incidence.rds")
    } else {
      filename <- paste0(dir_name, "incidence_fixed_incidence.rds")
    }
    GAR <- readRDS(filename) %>%
      vnapply(., sum)
    GAR <- data.frame(all = GAR / world_pop_size,
                      table_no = "incidence",
                      run_name = "special")
  } else if(run_name == "curr_alloc"){
    dir_name <- paste0(output_prefix, run_name, "/")
    if(production_delay == 7) {
      filename <- paste0(dir_name, "pd7_vacbycurr_alloc_incidence.rds")
    } else {
      filename <- paste0(dir_name, "curr_alloc_fixed_incidence.rds")
    }
    GAR <- readRDS(filename) %>%
      vnapply(., sum)
    GAR <- data.frame(all = GAR / world_pop_size,
                      table_no = "curr_alloc",
                      run_name = "special")
  } else if(run_name == "no_vax_no_seasonality") {
    dir_name <- "outputs_no_vax_no_seasonality/"
    GAR <- readRDS(paste0(dir_name, "no_vax_fixed_incidence.rds")) %>%
      vnapply(., sum)
    GAR <- data.frame(all = GAR / world_pop_size,
                      table_no = "no_vaccination_no_seasonality",
                      run_name = "no_vaccination")
  } else if(run_name == "no_vax") {
    dir_name <- "outputs_no_vax/"
    pop_size_table <- pop_size_table[order(pop_size_table$countryID),]
    pop_size <- pop_size_table$N
    pop_size <- c(pop_size, sum(as.numeric(pop_size)))
    
    filename <- "outputs_random/random_vaccinations_coverage_table_1_incidence.rds"
    
    GAR <- readRDS(filename)  %>%
      vnapply(., sum)
    GAR <- data.frame(all = GAR / world_pop_size,
                      table_no = "no_vaccination",
                      run_name = "no_vaccination")
  }
  saveRDS(GAR, paste0(dir_name, "GAR.rds"))
  GAR
}



plot_GAR <- function(production_delay) {
  if(production_delay == 7) {
    run_name <- c("top_n_countries", "fn_pop_size", "random", "incidence",
                  "no_vax_no_seasonality")
  } else {
    run_name <- c("top_n_countries", "fn_pop_size", "incidence", "curr_alloc", "no_vax")
  }
  
  plot_dir <- paste0("outputs_pd_", production_delay, "/")
  dir.create(plot_dir, showWarnings = FALSE)
  GAR <- lapply(run_name, load_GAR, production_delay = production_delay) %>%
    do.call(rbind, .)

  if(production_delay == 7) {
    GAR$run_name <- factor(GAR$run_name, levels = c("top_n_countries",
                                                    "fn_pop_size",
                                                    "random",
                                                    "special",
                                                    "no_vaccination"))
  } else {
    GAR$run_name <- factor(GAR$run_name, levels = c("top_n_countries",
                                                    "fn_pop_size",
                                                    "special",
                                                    "no_vaccination"))
  }

  g <- ggplot(GAR, aes(x=table_no,y = all)) +
    geom_boxplot(outlier.size=0.1,lwd=0.3) +
    coord_flip(ylim=c(0.3,0.6))+
    facet_wrap(~run_name, ncol = 1, scales = "free_y") +
    ylab("Global attack rate") +
    xlab("Run ID") +
    theme_bw() +
    theme(axis.text.x=element_text(size=8))
  ggsave(paste0(plot_dir, "all_GAR.pdf"), g, width = 20, height = 35, units = "cm")
  ggsave(paste0(plot_dir, "all_GAR.png"), g, width = 20, height = 35, units = "cm")
  g_landscape <- ggplot(GAR, aes(x=table_no,y = all)) +
    geom_boxplot(outlier.size=0.1,lwd=0.3) +
    coord_cartesian(ylim=c(0.3,0.6))+
    facet_wrap(~run_name, nrow = 1, scales = "free_x") +
    ylab("Global attack rate") +
    xlab("Run ID") +
    theme_bw() +
    theme(axis.text.x=element_text(size=8, angle = 90, hjust = 1))
  ggsave(paste0(plot_dir, "all_GAR_landscape.pdf"), g_landscape, width = 35, height = 20, units = "cm")
  ggsave(paste0(plot_dir, "all_GAR_landscape.png"), g_landscape, width = 35, height = 20, units = "cm")
  list(GAR = GAR, g = g, g_landscape = g_landscape)
  
  # g <- ggplot(GAR[GAR$run_name == "no_vaccination",], aes(x=table_no,y = all)) +
  # geom_boxplot(outlier.size=0.1,lwd=0.3) +
  # coord_cartesian(ylim = c(.4, .6)) +
  # ylab("Global attack rate") +
  # xlab("Run ID") +
  # theme_bw() +
  # theme(axis.text.x=element_text(size=8))
  # ggsave("GAR_no_vax.pdf", g, width = 10, height = 10, units = "cm")
}