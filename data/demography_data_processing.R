# script to process demography and contact matrices csv files
# the age groups are 0-4, 5-14, 15-64, 65+

setwd("~/git_repos/vaxedemic/data")
library(xlsx)


# first pass at demographic data
demographic_data <- read.csv("demography_data_world_population_prospects.csv",
                    stringsAsFactors = FALSE)
demographic_data <- demographic_data[demographic_data$Reference.date..as.of.1.July. == 2015 &
                     demographic_data$Country.code < 900,
                   c("Region..subregion..country.or.area..","Total","X0.4","X5.14","X15.64","X65.")]
numeric_cols <- seq(2,ncol(demographic_data))
demographic_data[,numeric_cols] <- as.data.frame(lapply(demographic_data[,numeric_cols], function(x) as.numeric(gsub(" ","",x))),
                               stringsAsFactors = FALSE)
demographic_data[,numeric_cols] <- demographic_data[,numeric_cols]*1000
demographic_data[,3:6] <- lapply(demographic_data[,3:6], function(x) x/demographic_data$Total)
colnames(demographic_data) <- c("countryID", "N", "propn_age1", "propn_age2", "propn_age3", "propn_age4")
demographic_country_names <- unlist(demographic_data$countryID)
saveRDS(demographic_data,"demography_data.rds")
write.table(demographic_data,file = "demography_data_clean.csv",sep = ",", row.names = FALSE)

# get contact data names
n_files <- 2
contact_filenames <- paste0("MUestimates_all_locations_", seq_len(n_files), ".xlsx")
workbooks <- lapply(contact_filenames, loadWorkbook)
sheets <- lapply(workbooks, getSheets)
contact_country_names <- lapply(sheets, names)
first_name_in_sheet_2 <- contact_country_names[[2]][1]
contact_country_names <- unlist(contact_country_names)

# use for extrapolation of contact matrices later.  
# For now, just take intersection of country names from the demographic and
# contact data.

# locations <- read.xlsx("locations.xlsx",1,startRow = 31, header = FALSE,colIndex = c(2,6,9))
# colnames(locations) <- c("country_name", "location_code", "subregion_name")
# locations <- locations[locations$location_code == 4,c(1,3)]

countries_in_both <- intersect(demographic_country_names, contact_country_names)
demographic_data <- demographic_data[demographic_country_names %in% countries_in_both,]
countries_in_both <- sort(countries_in_both)
demographic_data <- demographic_data[order(demographic_data$countryID),]

read_format_contact_data_closure <- function(first_name_in_sheet_2, contact_filenames){
  age_group_spacing <- 5
  nrow_contact_data <- 16
  age_groups <- age_group_spacing * (seq_len(nrow_contact_data) - 1)
  age_groups_wanted <- c(0, 5, 15, 65)
  age_groups_wanted_idx <- c(which(age_groups %in% age_groups_wanted), 
                             length(age_groups) + 1)
  
  sum_over_age_groups <- function(matrix_in){
    matrix_out <- outer(age_groups_wanted, age_groups_wanted) * 0
    for(i in seq_along(age_groups_wanted)) {
      for(j in seq_along(age_groups_wanted)) {
       matrix_out[i,j] <- sum(matrix_in[seq(age_groups_wanted_idx[i],age_groups_wanted_idx[i + 1] - 1),
                                        seq(age_groups_wanted_idx[j],age_groups_wanted_idx[j + 1] - 1)]) 
      }
    }
    matrix_out
  }
  
  f <- function(country_name){
    name_in_2 <- (country_name == first_name_in_sheet_2) ||
      (!(sort(c(country_name, first_name_in_sheet_2))[1] == country_name))
    filename <- contact_filenames[as.numeric(name_in_2) + 1]
    contact_data <- read.xlsx(filename, sheetName = country_name, header = !name_in_2)
    contact_data_summed <- sum_over_age_groups(contact_data)
    contact_data_summed
  }
  f
}

read_format_contact_data <- read_format_contact_data_closure(first_name_in_sheet_2, contact_filenames)

contact_data <- lapply(countries_in_both, read_format_contact_data)
saveRDS(contact_data,"contact_data.rds")
