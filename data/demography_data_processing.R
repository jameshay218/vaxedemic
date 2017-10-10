# script to process demography and contact matrices csv files
# the age groups are 0-4, 5-14, 15-64, 65+

setwd("~/Documents/vaxedemic/data") ## change working directory as necessary
source("JH_cleaning/generate_correct_country_list.R")

## read in list of original and new country IDs, original and new country names
## to enforce consistency between data sets
countries <- read.table("JH_cleaning/all_countries_FINAL.csv", sep = ",", header = TRUE)

## demographic data
## read in raw data
demographic_data <- read.csv("demography_data_world_population_prospects.csv",
                    stringsAsFactors = FALSE)
## find latest data sets (2015) and relevant age columns
demographic_data <- demographic_data[demographic_data$Reference.date..as.of.1.July. == 2015 &
                     demographic_data$Country.code < 900,
                   c("Region..subregion..country.or.area..","Total","X0.4","X5.14","X15.64","X65.")]
## make total population size and population size in age group columns numeric
numeric_cols <- seq(2,ncol(demographic_data))
demographic_data[,numeric_cols] <- as.data.frame(lapply(demographic_data[,numeric_cols], function(x) as.numeric(gsub(" ","",x))),
                               stringsAsFactors = FALSE)
## population size reported as x 10^3 -- scale accordingly
demographic_data[,numeric_cols] <- demographic_data[,numeric_cols]*1000
## find proportion of population in each age group
demographic_data[,3:6] <- lapply(demographic_data[,3:6], function(x) x/demographic_data$Total)

demographic_colnames <- c("countryID", "N", "propn_age1", "propn_age2", "propn_age3", "propn_age4")
colnames(demographic_data) <- demographic_colnames
## enforce consistent country names between data sets
demographic_data <- demographic_data[order(demographic_data$countryID),]
demographic_data <- correct_original_id_col(demographic_data,countries,"countryID")
demographic_data <- demographic_data[,demographic_colnames]
## find list of countries in this data set
demographic_country_names <- unlist(demographic_data$countryID)

## get contact data filenames
n_files <- 2 ## contact data spread across two files
contact_filenames <- paste0("MUestimates_all_locations_", seq_len(n_files), ".xlsx")

## get country names from names of spreadsheets
workbooks <- lapply(contact_filenames, xlsx::loadWorkbook)
sheets <- lapply(workbooks, xlsx::getSheets)
contact_country_names <- lapply(sheets, names)

## get first country name in sheet 2 to figure out which countries are in which
## spreadsheet (they are in alphabetical order)
first_name_in_sheet_2 <- contact_country_names[[2]][1]
contact_country_names <- unlist(contact_country_names)

## enforce consistent country names between data sets
contact_country_names_df <- data.frame(country_name = contact_country_names)
contact_country_names_df <- correct_original_id_col(contact_country_names_df, countries, "country_name")
correct_contact_country_names <- as.character(contact_country_names_df$country_name)

read_contact_data_closure <- function(first_name_in_sheet_2, contact_filenames){

  
  f <- function(country_name){
    name_in_2 <- (country_name == first_name_in_sheet_2) ||
      (!(sort(c(country_name, first_name_in_sheet_2))[1] == country_name))
    filename <- contact_filenames[as.numeric(name_in_2) + 1]
    contact_data <- xlsx::read.xlsx(filename, sheetName = country_name, header = !name_in_2)
  }
  f
}

condense_age_groups <- function(contact_data, propn_ages){
  age_group_spacing <- 5
  nrow_contact_data <- 16
  age_groups <- age_group_spacing * (seq_len(nrow_contact_data) - 1)
  age_groups_wanted <- c(0, 5, 15, 65)
  age_groups_wanted_idx <- c(which(age_groups %in% age_groups_wanted),
                             length(age_groups) + 1)

  condense_age_contacts <- outer(age_groups, age_groups_wanted) * 0
  for(i in seq_along(age_groups)) {
    for(j in seq_along(age_groups_wanted)) {
      condense_age_contacts[i,j] <- sum(contact_data[i,
                                                     seq(age_groups_wanted_idx[j],age_groups_wanted_idx[j + 1] - 1)])
    }
  }
  condense_age_individual <- outer(age_groups_wanted, age_groups_wanted) * 0
  for(i in seq_along(age_groups_wanted)) {
    for(j in seq_along(age_groups_wanted)) {
      condense_age_individual[i,j] <- mean(condense_age_contacts[seq(age_groups_wanted_idx[i],age_groups_wanted_idx[i + 1] - 1),
                                                                 j]) / propn_ages[j]
    }
  }
  condense_age_individual
}

read_contact_data <- read_contact_data_closure(first_name_in_sheet_2, 
                            contact_filenames)

contact_data <- lapply(contact_country_names, read_contact_data)

flight_data <- read.table("JH_cleaning/Huang2013/huang2013_melted.csv", sep = ",", header = TRUE)
flight_data <- correct_original_id_col(flight_data, countries, "Destination")
flight_data <- correct_original_id_col(flight_data,countries,"Origin")
flight_data <- flight_data[,c("Destination", "Origin", "Volume")]
flight_data_names <- unique(flight_data$Destination)

countries_in_all_datasets <- intersect(demographic_country_names, 
                                       intersect(correct_contact_country_names,
                                                 flight_data_names))

demographic_data <- demographic_data[demographic_data$countryID %in% countries_in_all_datasets,]
contact_data <- lapply(which(correct_contact_country_names %in% countries_in_all_datasets), function(x) contact_data[[x]])
flight_data <- flight_data[flight_data$Destination %in% countries_in_all_datasets &
                             flight_data$Origin %in% countries_in_all_datasets,]

flight_data <- reshape2::dcast(flight_data,Origin~Destination)
flight_data <- as.matrix(flight_data[,-1])
flight_data <- flight_data - diag(diag(flight_data))
contact_data <- lapply(seq_along(contact_data), function(x) condense_age_groups(contact_data[[x]], as.numeric(demographic_data[x,3:6])))
contact_data <- lapply(contact_data, function(x) matrix(x, nrow = 1))
contact_data <- do.call(rbind, contact_data)
contact_data <- cbind(countries_in_all_datasets, as.data.frame(contact_data))

write.table(demographic_data,file = "unified/demographic_data_intersect.csv",sep = ",", row.names = FALSE)
write.table(contact_data,file = "unified/contact_data_intersect.csv",sep = ",", row.names = FALSE)
write.table(flight_data,file = "unified/flight_data_intersect.csv",sep = ",", row.names = FALSE)
write(countries_in_all_datasets, file = "unified/countries_intersect.csv", sep = ",")
