####################################################
## Combines all distinct countries and their IDs from the original
## dataset, and saves this to a file of original IDs and original
## corresponding location names
####################################################

topwd <- "~/Documents/vaxedemic/raw_data/JH_cleaning/"
setwd(topwd)

## All the datasets we consider
data_names <- c("WTO", "OAG", "huang2013", "demography", "contact", "coverage")

read_data <- function(data_name) {
  filename <- paste0(data_name, "_countries.csv")
  my_df <- read.csv(filename,stringsAsFactors = FALSE)
  colnames(my_df) <- c("originalID","Location","Data")
  my_df
}

all_countries <- lapply(data_names, read_data)
all_countries <- do.call(rbind, all_countries)
all_countries <- all_countries[order(all_countries$originalID),]
colnames(all_countries) <- c("originalID","Location","Data")

write.table(all_countries,paste0(topwd,"original_ids.csv"),sep=",",row.names=FALSE)
