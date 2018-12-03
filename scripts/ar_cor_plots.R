library(dplyr)

ar_melt <- melt(attack_rate1,id.vars=c("table_no","filename","samp","all"))
colnames(ar_melt) <- c("table_no","filename","samp","all","country","AR")
vacc_melt <- melt(data.table(vaccinated_overall),id.vars=c("table_no","filename","samp","all"))
colnames(vacc_melt) <- c("table_no","filename","samp","all_vacc","country","vacc_no")

ar_melt_median <- ddply(ar_melt, .(table_no,country), function(x) c(median(x$AR),median(x$all)))
colnames(ar_melt_median) <- c("table_no","country","median_AR","median_overall_AR")
vacc_melt_median <- ddply(vacc_melt, .(table_no,country), function(x) c(median(x$vacc_no),median(x$all_vacc)))
colnames(vacc_melt_median) <- c("table_no","country","median_vacc","median_overall_vacc")

ar_melt_median$table_no <- as.character(ar_melt_median$table_no)
ar_melt_median$country <- as.character(ar_melt_median$country)

vacc_melt_median$table_no <- as.character(vacc_melt_median$table_no)
vacc_melt_median$country <- as.character(vacc_melt_median$country)

all_melt <- merge(vacc_melt_median[,c("table_no","country","median_vacc","median_overall_vacc")],
                  ar_melt_median[,c("table_no","country","median_AR","median_overall_AR")],
                  id.vars=c("table_no","country"))

all_cors <- ddply(all_melt, ~country, function(x) cor(x$median_AR, x$median_overall_AR))
pop_size_table$countryID <- as.character(pop_size_table$countryID)
colnames(pop_size_table)[1] <- "country"
all_cors <- merge(all_cors, pop_size_table[,c("country","N")])
all_cors$country <- factor(all_cors$country, levels=all_cors[order(all_cors$N),"country"])
g_corr <- ggplot(all_cors) + 
  geom_point(aes(x=country,y=V1), size=0.75) + 
  coord_flip() +
  theme_bw() +
  xlab("Correlation between country AR and global AR") +
  ylab("Country")
ggsave("~/net/home/vaxedemic/plots/ar_corrs.pdf", g_corr, height = 28,width=20,units="cm")

vacc_corrs <- ddply(all_melt, ~country, function(x) cor(x$median_vacc, x$median_overall_AR))
vacc_corrs <- merge(vacc_corrs, pop_size_table[,c("country","N")])
vacc_corrs$country <- factor(vacc_corrs$country, levels=all_cors[order(vacc_corrs$N),"country"])
g_vacc <- ggplot(vacc_corrs) + 
  geom_point(aes(x=country,y=V1), size=0.75) + 
  coord_flip() +
  theme_bw() +
  xlab("Correlation between number vaccinated in\n country and global AR") +
  ylab("Country")
ggsave("~/net/home/vaxedemic/plots/vacc_corrs.pdf", g_corr, height = 28,width=20,units="cm")
