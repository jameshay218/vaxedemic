## This script takes the WTO nonresident travel data and converts it to a form useable
## as a travel contact matrix between any two countries
codeFile <- "airports_extended_huang2013.csv"
flightFile <- "prediction_huang2013.csv"
SAVE <- TRUE
ADD_RETURNS <- FALSE
topwd <- "~/Documents/vaxedemic/data/JH_cleaning/Huang2013"
setwd(topwd)
source("../travel_plotting_funcs.R")

predict_data <- read.csv(flightFile,stringsAsFactors=FALSE)

## https://openflights.org/data.html
airport_data <- read.csv(codeFile,header = FALSE,stringsAsFactors=FALSE)
airport_data <- airport_data[,c(4,3)]

## else
#airport_data <- read.csv("airports_alternative.csv",stringsAsFactors=FALSE)
#airport_data <- airport_data[,c("iata_code","iso_country")]

colnames(airport_data) <- c("code","country")

#################################
## EXTRA CODE TO GET AIRPORTS THAT WEREN'T IN THE LIST
city_list <- read.csv("AirportInfo.csv",stringsAsFactors = FALSE)

## 12 airport codes in Huang 2013 not in the database
## These should be the same for origin and destination
origins_missing <- unique(predict_data$Origin)[which(!(unique(predict_data$Origin) %in% unique(airport_data$code)))]
destinations_missing <- unique(predict_data$Destination)[which(!(unique(predict_data$Destination) %in% unique(airport_data$code)))]

missing_cities <- city_list[city_list$NodeName %in% origins_missing,]

## These are ABC, CFK, KHY, MKX, MLH, NBW, NZE, OSM, PHG, RUD, SXV, SYJ
## Manual extraction of countries
additional_countries <- c("ABC"="Spain","CFK"="Algeria","KHY"="Iran","MKX"="Yemen",
                             "MLH"="France","NBW"="Cuba","NZE"="Guinea","OSM"="Iraq",
                             "PHG"="Nigeria","RUD"="Iran","SXV"="India","SYJ"="Iran")
## Only need airports that are in the data
airport_data <- airport_data[airport_data$code %in% predict_data$Origin | airport_data$code %in% predict_data$Destination,]

airports <- airport_data$country
names(airports) <- airport_data$code
airports <- c(airports,additional_countries)
##################################

###########################################
## Should be the same rows before and after
print(nrow(predict_data))
predict_data <- predict_data[predict_data$Origin %in% names(airports) & 
                               predict_data$Destination %in% names(airports),]
print(nrow(predict_data))
##########################################

## convert codes to counrties
predict_data$origin_country <- airports[predict_data$Origin]
predict_data$destination_country <- airports[predict_data$Destination]
#View(predict_data[predict_data$Origin=="LHR",])

## Aggregate by country
summed <- aggregate(predict_data[,"PredMu"], 
                    predict_data[,c("origin_country","destination_country")], FUN=sum)

## Origin as rows, destination as columns
final <- data.table::dcast(summed, origin_country~destination_country, value.var="x")
final[is.na(final)] <- 0


## Expand grid to get NA rather than actually missing in data set for missing pairs
filled_countries <- expand.grid(unique(summed$destination_country),unique(summed$origin_country),stringsAsFactors=FALSE)
colnames(filled_countries) <- c("destination_country","origin_country")

summed <- merge(filled_countries,summed,by=c("destination_country","origin_country"),all=TRUE)
summed$logVolume <- log(summed$x)
summed$logVolume[summed$x == 0] <- 0
summed[is.na(summed)] <- 0
colnames(summed) <- c("Destination","Origin","Volume","logVolume")

if(SAVE) write.table(summed,paste0(topwd,"/huang2013_melted.csv"),row.names=FALSE,sep=",")
if(SAVE) write.table(final,paste0(topwd,"/huang2013_travel_mat.csv"),row.names=FALSE,sep=",")


###########
## PLOTS

p1 <- plot_heatmap_travel(summed)
p2 <- plot_total_journies(summed,FALSE)


if(SAVE){
  png(paste0(topwd,"/huang2013_heatmap.png"),width=15,height=15,res=600,units="in")
  print(p1)
  dev.off()
  
  png(paste0(topwd,"/huang2013_numbers"),width=16,height=8,res=600,units="in")
  print(p2)
  dev.off()
  }



##########
## For country labels
countries <- unique(summed$Origin)
countries <- countries[order(countries)]
countries <- data.frame("originalID"=1:length(countries),"Location"=countries,"Data"="Huang2013")
if(SAVE) write.table(countries,"~/Documents/vaxedemic_tmp/huang2013_countries.csv",sep=",",row.names=FALSE)