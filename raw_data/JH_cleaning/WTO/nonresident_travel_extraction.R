## This script takes the WTO nonresident travel data and converts it to a form useable
## as a travel contact matrix between any two countries
codeFile <- "flight_ids_WTO.csv"
flightFile <- "raw_flights_WTO.csv"
SAVE <- TRUE
ADD_RETURNS <- FALSE
topwd <- "~/Documents/vaxedemic/data/JH_cleaning/WTO"
setwd(topwd)

## Read in country IDs and names
## Optionally filter if we only want to look at some
#tmp_ids <- c("China","USA","United_Kingdom","France","Germany","Brazil")
ids <- read.csv(codeFile,stringsAsFactors=FALSE)
tmp_ids <- unique(ids$Location)
some_ids <- ids[ids$Location %in% tmp_ids,]

## Read in arrivals data as full matrix
data <- read.csv(flightFile,stringsAsFactors=FALSE)
data <- data[,colnames(data)!="Location"]

## Everyone who arrives at country X from country Y returns to country Y from country X that year
if(ADD_RETURNS){
  for(i in 1:nrow(data)){
    for(j in 1:ncol(data)){
      data[i,j] <- data[i,j] + data[j,i]
    }
  }
}

## Extract locations for which we have on arrivals data
no_arrivals_data <- ids[(which(colSums(data) == 0)),]

## Melt data for clearer analyses
melted_data <- reshape2::melt(as.matrix(data))

## Label as 
colnames(melted_data) <- c("CountrySourceRow","CountryDestinationCol","Volume")
melted_data$CountrySourceRow <- as.numeric(melted_data$CountrySourceRow)

melted_data$CountrySourceRow <- ids$Location[melted_data$CountrySourceRow]
melted_data$CountryDestinationCol <- ids$Location[melted_data$CountryDestinationCol]

melted_data <- melted_data[melted_data$CountrySourceRow %in% some_ids$Location &
                             melted_data$CountryDestinationCol %in% some_ids$Location,]
melted_data$logVolume <- log(melted_data$Volume)
melted_data$logVolume[melted_data$Volume == 0] <- 0

colnames(melted_data) <- c("Origin","Destination","Volume")
melted_data$logVolume <- log(melted_data$Volume)
melted_data$logVolume[melted_data$Volume == 0] <- 0

if(SAVE) write.table(no_arrivals_data,paste0(topwd,"/no_arrivals_locations.csv"),sep=",",row.names=FALSE)
if(SAVE) write.table(melted_data,paste0(topwd,"/wto_melted.csv"),sep=",",row.names=FALSE)
  
p1 <- plot_heatmap_travel(melted_data)
p2 <- plot_total_journies(melted_data,FALSE)

if(SAVE){
  png(paste0(topwd,"/wto_heatmap.png"),width=15,height=15,res=600,units="in")
  print(p1)
  dev.off()
  
  png(paste0(topwd,"/wto_numbers.png"),width=16,height=8,res=600,units="in")
  print(p2)
  dev.off()
  
}

##########
## For country labels
countries <- unique(melted_data$Origin)
countries <- countries[order(countries)]
countries <- data.frame("originalID"=1:length(countries),"Location"=countries,"Data"="WTO")
if(SAVE) write.table(countries,"~/Documents/vaxedemic_tmp/WTO_countries.csv",sep=",",row.names=FALSE)
