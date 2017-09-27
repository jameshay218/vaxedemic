topwd <- "~/Documents/vaxedemic_tmp/"
setwd(topwd)

wto <- read.csv("WTO_countries.csv",stringsAsFactors = FALSE)[,2]
oag <- read.csv("OAG_countries.csv",stringsAsFactors = FALSE)[,2]
huang2013 <- read.csv("huang2013_countries.csv",stringsAsFactors = FALSE)[,2]
demography <- read.csv("demography_countries.csv",stringsAsFactors = FALSE)[,2]
contacts <- read.csv("contact_countries.csv",stringsAsFactors = FALSE)[,2]

all_countries <- unique(c(wto,oag,huang2013,demography,contacts))
all_countries <- all_countries[order(all_countries)]
all_countries <- data.frame("Location"=all_countries,"ID"=1:length(all_countries))

