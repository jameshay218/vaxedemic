#' Check for new countries
#'
#' Given a matrix of country names and IDs, checks which of these are not in the full country list and returns the 
#' appended list
#' @param newLocations a data frame of new country names, with columns "originalID","Location","Data"
#' @param countryFileLoc the file location for the correct matrix of country names
#' @param return_diff if TRUE, returns only the new countries (ie. removes countries that are already in the correct table)
#' @param save if TRUE, appends the new countries to the correct table and saves it
#' @return a data frame of country names and their correct IDs along with the original information
#' @export
check_new_countries <- function(newLocations, countryFileLoc="all_countries_FINAL.csv",return_diff=TRUE,save=TRUE){
  countries <- read.csv(countryFileLoc,stringsAsFactors = FALSE)
  
  ############################################
  missing_countries <- dplyr::anti_join(newLocations,countries,by=c("originalID","Location"))
  
  if(!is.null(missing_countries) && nrow(missing_countries) > 0){ 
    new_countries <- merge(missing_countries,unique(countries[,c("Location","Correct","correctID")]),all.x=TRUE)
    new_countries <- new_countries[,colnames(countries)]
    countries <- rbind(countries,new_countries)
  }
  ############################################
  if(save) write.table(countries, countryFileLoc, sep=",",row.names=FALSE)
  if(return_diff){
    if(!is.null(missing_countries) && nrow(missing_countries) > 0) return(missing_countries)
    else return(NULL)
  }
  return(countries)
}


## Make sure country list is self consistent. This country list should have a column for
## country name as in its original database, and a column for the correct name
clean_country_database <- function(fileLoc, topwd, SAVE){
  countries <- read.csv(fileLoc,stringsAsFactors=FALSE)
  if(length(unique(countries$Correct[!(countries$Correct %in% countries$Location)])) > 0){
    missing <- data.frame("Correct"=unique(countries$Correct[!(countries$Correct %in% countries$Location)]),
                        "Location"=unique(countries$Correct[!(countries$Correct %in% countries$Location)]),
                        "Data"="database","originalID"=NA)
  }
  original_ids <- read.csv(paste0(topwd,"original_ids.csv"),stringsAsFactors=FALSE)
  correct_ids <- read.csv(paste0(topwd,"correct_countries.csv"),stringsAsFactors=FALSE)

  if(!is.null(missing) && nrow(missing) > 0){ 
    missing <- merge(missing, correct_ids,by="Correct")
    countries <- rbind(countries,missing[,colnames(countries)])
  
  }
  if(SAVE) write.csv(countries,fileLoc,row.names=FALSE)
  return(countries)
}

## Generate correct country IDS
generate_correct_country_ids <- function(countries, SAVE=FALSE,topwd="~/Documents/vaxedemic/data/JH_cleaning/"){
  ## Get unique correct countries and generate IDs for these
  correct_countries <- unique(countries$Correct)
  correct_countries <- correct_countries[order(correct_countries)]
  correct <- data.frame("Correct"=correct_countries,"correctID"=1:length(correct_countries))
  
  if(SAVE) write.table(correct,paste0(topwd,"correct_countries.csv"),sep=",",row.names=FALSE)
  return(correct)
}

## Given a data frame with an originalID and Location, returns the same data frame but 
## with a correctID and Correct column
correct_original_id <- function(dataset, countries){
  corrected <- merge(dataset[,colnames(dataset) != "Data"],unique(countries[,colnames(countries) != "Data"]),
                     by=c("originalID","Location"),all.x=TRUE)
  return(corrected)
}

correct_original_id_col <- function(dataset, countries, col){
  colnames(dataset)[colnames(dataset) == col] <- "Location"
  correctNames <- unique(merge(dataset,countries[,c("Location","Correct","correctID")],by="Location"))
  colnames(correctNames)[colnames(correctNames) %in% c("Correct","correctID")] <- c(col,paste0(col,"ID"))
  colnames(correctNames)[colnames(correctNames) == "Location"] <- paste0("Old",col)
  return(correctNames)
}

#SAVE <- FALSE
#countries <- clean_country_database("~/Documents/vaxedemic/data/JH_cleaning/all_countries_FINAL.csv",
#                                    topwd="~/Documents/vaxedemic/data/JH_cleaning/",TRUE)
#correct <- generate_correct_country_ids(countries,TRUE)

# locations <- read.csv("LocationData/latitudes.csv",stringsAsFactors=FALSE)
# locations <- locations[,c("Location","originalID","Data")]
# wow <- check_new_countries(locations,return_diff=TRUE,save=FALSE)
