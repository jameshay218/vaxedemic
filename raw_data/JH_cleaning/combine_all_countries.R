####################################################
## Combines all distinct countries and their IDs from the original
## dataset, and saves this to a file of original IDs and original
## corresponding location names
####################################################

topwd <- "~/Documents/vaxedemic/data/JH_cleaning/"
setwd(topwd)

## All the datasets we consider
wto <- read.csv("WTO_countries.csv",stringsAsFactors = FALSE)
oag <- read.csv("OAG_countries.csv",stringsAsFactors = FALSE)
huang2013 <- read.csv("huang2013_countries.csv",stringsAsFactors = FALSE)
demography <- read.csv("demography_countries.csv",stringsAsFactors = FALSE)
contacts <- read.csv("contact_countries.csv",stringsAsFactors = FALSE)

colnames(wto) <- colnames(oag) <- colnames(huang2013) <-
  colnames(demography) <- colnames(contacts) <-c("originalID","Location","Data")

all_countries <- rbind(wto,oag,huang2013,demography,contacts)
all_countries <- all_countries[order(all_countries$originalID),]
colnames(all_countries) <- c("originalID","Location","Data")

write.table(all_countries,paste0(topwd,"original_ids.csv"),sep=",",row.names=FALSE)
