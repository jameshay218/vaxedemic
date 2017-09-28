##################
## From a some specified data files, extracts the unique country names and their original IDs and 
## save these to their own file. These will then be used by "combine_all_countries.R" to combine into one
## data set
##################
topwd <- "~/Documents/vaxedemic/data/JH_cleaning/"
setwd(topwd)
#######
## Get all country names in different datasets
contacts <- read.csv("contact_data_clean.csv",stringsAsFactors=FALSE)

countries_contacts <- unique(contacts$contact_country_names)
countries_contacts <- countries_contacts[order(countries_contacts)]
countries_contacts <- data.frame("originalID"=1:length(countries_contacts),"Location"=countries_contacts,"Data"="contact")
write.table(countries_contacts,paste0(topwd,"contact_countries.csv"),sep=",",row.names=FALSE)

demography <- read.csv("demography_data_clean.csv",stringsAsFactors=FALSE)

countries_demography <- unique(demography$countryID)
countries_demography <- countries_demography[order(countries_demography)]
countries_demography <- data.frame("originalID"=1:length(countries_demography),"Location"=countries_demography,"Data"="demography")
write.table(countries_demography,paste0(topwd,"demography_countries.csv"),sep=",",row.names=FALSE)
