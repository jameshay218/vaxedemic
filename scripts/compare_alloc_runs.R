library(ggplot2)
library(magrittr)
setwd("~/net/home/vaxedemic")
postprocess_incidence <- function(filename, table_no) {
  print(table_no)
  inc <- readRDS(filename) %>%
    lapply(., function(x) apply(x, 1, sum)) %>%
    do.call(rbind, .)
  inc <- cbind(inc, matrix(apply(inc, 1, sum), ncol = 1, dimnames = list(NULL, "all")))
  attack_rate <- t(apply(inc, 1, function(x) x / pop_size)) %>%
    as.data.frame %>%
    cbind(., table_no, filename)
  attack_rate
}

postprocess_vaccinated <- function(filename, table_no) {
  print(table_no)
  allocd <- readRDS(filename) %>%
    lapply(., function(x) apply(x, 1, max)) %>%
    do.call(rbind, .)
  allocd <- cbind(allocd, matrix(apply(allocd, 1, sum), ncol = 1, dimnames = list(NULL, "all")))
  vaccinated_overall <- t(apply(allocd, 1, function(x) x / 1)) %>%
    as.data.frame %>%
    cbind(., table_no, filename)
  vaccinated_overall
}

pop_size_table <- read.csv("~/Documents/vaxedemic/data/demographic_data_intersect.csv")
pop_size_table <- pop_size_table[order(pop_size_table$countryID),]
pop_size <- pop_size_table$N
pop_size <- c(pop_size, sum(as.numeric(pop_size)))

setwd("outputs_random_vaccinations/")
completed_runs <- list.files(pattern="_incidence.rds")
vaccinated_runs <- list.files(pattern="_vaccinated.rds")
file_labels <- data.frame(filename=completed_runs, vacc_file=vaccinated_runs)

labels <- seq(0,length(completed_runs)-1,by=1)
labels[1] <- "maxed_vaccination"
labels[which(completed_runs == "random_vaccinations_coverage_table_103_incidence.rds")] <- "proportional_to_popn"
labels[which(completed_runs == "random_vaccinations_coverage_table_102_incidence.rds")] <- "current_allocation"
labels[which(completed_runs == "random_vaccinations_coverage_table_1_incidence.rds")] <- "no_vaccination"

attack_rate <- Map(postprocess_incidence, completed_runs, labels) %>%
  do.call(rbind, .) %>%
  as.data.frame

vaccinated_overall <- Map(postprocess_vaccinated, vaccinated_runs, labels) %>%
  do.call(rbind, .) %>%
  as.data.frame

attack_rate$samp <- 1:nrow(attack_rate)
vaccinated_overall$samp <- 1:nrow(vaccinated_overall)

attack_rate1 <- attack_rate
attack_rate1 <- data.table(attack_rate1)


## Branch to cor plot now
attack_rate1 <- merge(attack_rate1, data.table(file_labels),by="filename")
attack_rate1$table_no <- as.factor(attack_rate1$table_no)
median_ar <- ddply(attack_rate1, ~table_no, function(x) median(x[,"all"]))
attack_rate1$table_no <- factor(attack_rate1$table_no, levels=median_ar$table_no[order(median_ar$V1)])
tmp <- vaccinated_overall[,c("China","India","table_no","samp")]
colnames(tmp) <- c("China_vacc","India_vacc","table_no","samp")
tmp <- ddply(tmp, ~table_no, function(x) c(median(x$China_vacc), median(x$India_vacc)))
colnames(tmp) <- c("table_no","China_vacc","India_vacc")
tmp$table_no <- as.character(tmp$table_no)
attack_rate1$table_no <- as.character(attack_rate1$table_no)
attack_rate3 <- merge(attack_rate1, tmp)
attack_rate3$table_no <- factor(attack_rate3$table_no, levels=tmp[order(tmp$China_vacc),"table_no"])


g <- ggplot(attack_rate1, aes(x=table_no,y = all)) + 
  coord_cartesian(ylim=c(0.3,0.6))+
  geom_boxplot(outlier.size=0.1,lwd=0.3) +
  coord_flip(ylim=c(0.3,0.6))+
  ylab("Global attack rate") +
  xlab("Run ID") + 
  theme_bw() +
  theme(axis.text.x=element_text(size=8))
g
ggsave("~/net/home/vaxedemic/plots/different_attack_rates.pdf", g_china,height = 28,width=20,units="cm")

attack_rate3 <- merge(attack_rate1, tmp)
attack_rate3$sum_China_India_vacc <- attack_rate3$China_vacc + attack_rate3$India_vacc
tmp$sum_China_India_vacc <- tmp$China_vacc + tmp$India_vacc
attack_rate3$table_no <- factor(attack_rate3$table_no, levels=tmp[order(tmp$sum_China_India_vacc),"table_no"])

g_china <- ggplot(attack_rate3, aes(x=table_no,y = all, col=sum_China_India_vacc)) + 
  coord_cartesian(ylim=c(0.3,0.6))+
  geom_boxplot(outlier.size=0.1,lwd=0.3) +
  coord_flip(ylim=c(0.3,0.6))+
  scale_color_gradient2(low="red",high="green",mid="blue",midpoint = 3*max(attack_rate3$sum_China_India_vacc)/4)+
  ylab("Global attack rate") +
  xlab("Run ID") + 
  theme_bw() +
  theme(axis.text.x=element_text(size=8))
g_china
ggsave("~/net/home/vaxedemic/plots/different_attack_rates_by_china_plus_india.pdf", g_china,height = 28,width=20,units="cm")

setwd("~/net/home/vaxedemic/data/random_coverage_tables/")
alloc_files <- list.files()

alloc_tots <- NULL
index <- 1
for(file in alloc_files){
  tmp <- read.csv(file,stringsAsFactors=FALSE)
  tmp <- cbind(tmp, index, data_file=paste0("data/random_coverage_tables/",file))
  index <- index + 1
  alloc_tots <- rbind(alloc_tots, tmp)
}
alloc_tots$data_file 
unsplit_names <- strsplit(unique(as.character(alloc_tots$file)),split='[._]')
recomb <- lapply(unsplit_names, function(x) paste0("random_vaccinations_coverage_table_"))


runs1 <- read.csv("~/net/home/vaxedemic/data/run_key.csv",stringsAsFactors=FALSE)
runs1 <- rbind(runs1, c("maximum_vacc_fixed_incidence.rds","data/random_coverage_tables/random_coverage_data_MAX.csv"))
colnames(runs1) <- c("filename","data_file")

attack_rate1$filename <- as.character(attack_rate1$filename)
attack_rate2 <- merge(runs1, attack_rate1, by=c("filename"))
attack_rate2 <- melt(attack_rate2, id.vars=c("filename","data_file","table_no","vacc_file"))
attack_rate2$table_no <- as.character(attack_rate2$table_no)
attack_rate2$data_file <- as.character(attack_rate2$data_file)
colnames(attack_rate2)[5] <- "country"
attack_rate2$country <- as.character(attack_rate2$country)
alloc_tots$country <- as.character(alloc_tots$country)
alloc_tots$data_file <- as.character(alloc_tots$data_file)

attack_rate2 <- merge(data.table(attack_rate2),data.table(alloc_tots),by=c("data_file","country"))
#alloc_tots_USA <- alloc_tots[alloc_tots$country=="USA",]

attack_rate_USA <- attack_rate2[attack_rate2$country=="USA",]
median_ar_USA <- ddply(attack_rate_USA,~table_no,function(x) median(x$dose_per_1000))
levels <- order(median_ar_USA$V1)
attack_rate_USA$table_no <- factor(attack_rate_USA$table_no, levels=median_ar_USA[levels,"table_no"])

vaccinated_overall_melt <- melt(vaccinated_overall, id.vars=c("table_no","filename"))
vaccinated_overall_melt_USA <- vaccinated_overall_melt[vaccinated_overall_melt$variable == "USA",]
median_vacc_USA <- ddply(vaccinated_overall_melt_USA,~table_no,function(x) median(x$value))
levels <- order(median_vacc_USA$V1)
colnames(median_vacc_USA) <- c("table_no","number_vacc")
median_vacc_USA$table_no <- as.character(median_vacc_USA$table_no)
attack_rate_USA$table_no <- as.character(attack_rate_USA$table_no)
attack_rate_USA <- merge(attack_rate_USA, data.table(median_vacc_USA), by="table_no")
attack_rate_USA$table_no <- factor(attack_rate_USA$table_no, levels=median_vacc_USA[levels,"table_no"])

#median_ar_usa <- ddply(attack_rate1, ~table_no, function(x) median(x[,"USA"]))
#attack_rate_USA$table_no <- factor(attack_rate_USA$table_no, levels=median_ar_usa$table_no[order(median_ar_usa$V1)])
g_USA <- ggplot(attack_rate_USA, aes(x=table_no,y = as.numeric(value),col=number_vacc)) + 
  geom_boxplot(outlier.size=0.1,lwd=0.3) +
  scale_color_gradient2(low="red",high="green",mid="blue",midpoint = 1.5e08)+
  # facet_wrap(~x) +
  coord_flip()+
  #coord_cartesian(ylim = c(0,1), expand = FALSE) + 
  ylab("Attack rate in USA") +
  xlab("Run ID") + 
  theme_bw() +
  theme(axis.text.x=element_text(size=8))
g_USA
ggsave("~/net/home/vaxedemic/plots/attack_rate_USA.pdf", g_USA, height = 28,width=20,units="cm")
#attack_rate1$filename <- sapply(attack_rate1$filename, 
#         function(x) paste0(unlist(strsplit(x, split="[.]"))[1],".",unlist(strsplit(x, split="[.]"))[2]))
