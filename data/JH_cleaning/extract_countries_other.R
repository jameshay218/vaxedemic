topwd <- "~/Documents/vaxedemic_tmp/"
setwd(topwd)
#######
## Get all country names in different datasets
contacts <- read.csv("contact_data_clean.csv",stringsAsFactors=FALSE)

countries_contacts <- unique(contacts$contact_country_names)
countries_contacts <- countries_contacts[order(countries_contacts)]
countries_contacts <- data.frame("ID"=1:length(countries_contacts),"Location"=countries_contacts)
write.table(countries_contacts,"~/Documents/vaxedemic_tmp/contact_countries.csv",sep=",",row.names=FALSE)

demography <- read.csv("demography_data_clean.csv",stringsAsFactors=FALSE)

countries_demography <- unique(demography$countryID)
countries_demography <- countries_demography[order(countries_demography)]
countries_demography <- data.frame("ID"=1:length(countries_demography),"Location"=countries_demography)
write.table(countries_demography,"~/Documents/vaxedemic_tmp/demography_countries.csv",sep=",",row.names=FALSE)

