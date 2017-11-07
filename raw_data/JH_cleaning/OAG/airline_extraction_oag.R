## Cleaning flight data
SAVE <- TRUE
topwd <- "~/Documents/vaxedemic/data/JH_cleaning/OAG"
setwd(topwd)
codeFile <- "countries.txt"
flight_ver <- "f2008_07"
flightFile <- paste0(flight_ver,".txt")

## Read in country IDs and names
## Optionally filter if we only want to look at some
#tmp_ids <- c("China","USA","United_Kingdom","France","Germany","Brazil")
ids <- read.table(codeFile,stringsAsFactors=FALSE)
colnames(ids) <- c("ID","Location")
tmp_ids <- unique(ids$Location)
some_ids <- ids[ids$Location %in% tmp_ids,]

## Read in arrivals data as full matrix
data <- read.table(flightFile,stringsAsFactors=FALSE)
data <- data[,colnames(data)!="Location"]

## Extract locations for which we have on arrivals data
no_arrivals_data <- ids[(which(rowSums(data) == 0)),]
if(SAVE_OUTPUTS) write.table(no_arrivals_data,paste0(topwd,"/no_arrivals_locations.csv"),sep=",",row.names=FALSE)

#message("Internal flight numbers:")
#for(i in 1:nrow(data)) message(paste0(ids$Location[i],": ",data[i,i]))

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

colnames(melted_data) <- c("Origin","Destination","Volume","logVolume")

if(SAVE) write.table(melted_data,paste0(topwd,"/",flight_ver,"_melted.csv"),sep=",",row.names=FALSE)

p1 <- plot_heatmap_travel(melted_data)
p2 <- plot_total_journies(melted_data,FALSE)


if(SAVE){
  png(paste0(topwd,"/",flight_ver,"_heatmap.png"),width=15,height=15,res=600,units="in")
  print(p1)
  dev.off()
  
  png(paste0(topwd,"/",flight_ver,"_numbers.png"),width=16,height=8,res=600,units="in")
  print(p2)
  dev.off()
  
}

##########
## For country labels
countries <- unique(melted_data$Origin)
countries <- countries[order(countries)]
countries <- data.frame("originalID"=1:length(countries),"Location"=countries,"Data"="OAG")
if(SAVE) write.table(countries,paste0(topwd,"_countries.csv"),sep=",",row.names=FALSE)

