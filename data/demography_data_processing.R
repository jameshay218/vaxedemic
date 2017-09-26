# script to process demography csv file and turn into matrix

my_data <- read.csv("demography_data_world_population_prospects.csv",
                    stringsAsFactors = FALSE)
my_data <- my_data[my_data$Reference.date..as.of.1.July. == 2015 &
                     my_data$Country.code < 900,
                   c("Region..subregion..country.or.area..","X0.4","X5.14","X15.64","X65.")]
colnames(my_data) <- c("Country","0-4","15-64","65+")
numeric_cols <- seq(2,5)
my_data[,numeric_cols] <- as.data.frame(lapply(my_data[,numeric_cols], function(x) as.numeric(gsub(" ","",x))),
                               stringsAsFactors = FALSE)
my_data[,numeric_cols] <- my_data[,numeric_cols]*1000
saveRDS(my_data,"demography_data_world_population_prospects.rds")
write.table(my_data,file = "demography_data_world_population_prospects_clean.csv",sep = ",")