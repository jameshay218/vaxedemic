devtools::load_all("~/Documents/vaxedemic/")
attack_rates_data <- read.csv("outputs/pd0no_vaccination/no_vaccination_country_attack_rates_fixed_median_ci.csv",
                              stringsAsFactors = FALSE)
library(maps)
map_data_dir <- Sys.getenv("R_MAP_DATA_DIR_WORLD")
map_data_filename <- paste0(map_data_dir, "world.N")
map_data <- read.csv(map_data_filename, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
map_countries <- map_data[,1] %>%
  vcapply(function(x) strsplit(x, ":", fixed = TRUE)[[1]][1]) %>%
  gsub(" ", "_", .)
our_countries <- attack_rates_data$Location
our_countries[!(our_countries %in% map_countries)]
map_countries[map_countries == "UK"] <- "United_Kingdom"
map_countries [map_countries == "North_Korea"] <- "Republic_of_Korea"
map_countries [map_countries == "Syria"] <- "Syrian_Arab_Republic"
median_attack_rates <- attack_rates_data$median
attack_rates_map <- median_attack_rates[match(map_countries, our_countries)]
attack_rates_buckets <- as.numeric(cut(attack_rates_map, seq(0, 1, by = .01)))
colors <- colorRampPalette(brewer.pal(9,"Reds")[-(1:2)])(100)


png("plots/attack_rate_no_vaccination.png",width=1200,height=800)
map("world", col = colors[attack_rates_buckets], fill = TRUE)
# legend("topright", leg.txt, horiz = TRUE, fill = colors)
dev.off()

attack_rates_data <- read.csv("outputs/pd180incidence/incidence_country_attack_rates_fixed_median_ci.csv",
                              stringsAsFactors = FALSE)
median_attack_rates <- attack_rates_data$median
attack_rates_map <- median_attack_rates[match(map_countries, our_countries)]
attack_rates_buckets <- as.numeric(cut(attack_rates_map, seq(0, 1, by = .01)))
png("plots/attack_rate_incidence180.png",width=1200,height=800)
map("world", col = colors[attack_rates_buckets], fill = TRUE)
# legend("topright", leg.txt, horiz = TRUE, fill = colors)
dev.off()

attack_rates_data <- read.csv("outputs/pd180curr_alloc/curr_alloc_country_attack_rates_fixed_median_ci.csv",
                              stringsAsFactors = FALSE)
median_attack_rates <- attack_rates_data$median
attack_rates_map <- median_attack_rates[match(map_countries, our_countries)]
attack_rates_buckets <- as.numeric(cut(attack_rates_map, seq(0, 1, by = .01)))
png("plots/attack_rate_curr_alloc180.png",width=1200,height=800)
map("world", col = colors[attack_rates_buckets], fill = TRUE)
# legend("topright", leg.txt, horiz = TRUE, fill = colors)
dev.off()

# manually plot colour bar
my_df <- data.frame(x = c(0,1), y = c(0,1))
g <- ggplot(my_df, aes(x = x, y = y)) + 
  geom_point(aes(color = y)) + 
  scale_color_distiller(palette = "Reds", direction = 1, limits = c(-2/7,1),
                        name = "Attack rate")
ggsave("plots/colorbar.pdf", g, width = 10, height = 10, units = "cm")
