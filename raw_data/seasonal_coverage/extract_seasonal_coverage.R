coverage <- c("low", "high")
filename <- paste0("raw_data/seasonal_coverage/", coverage, "_coverage.xlsx")
coverage_data <- lapply(filename, function(x) read.xlsx(x, 2, startRow = 3, 
                                                        header = FALSE, 
                                                        colIndex = c(1,3)))
coverage_data <- do.call(rbind, coverage_data)
names(coverage_data) <- c("country", "dose_per_1000")
coverage_data$country <- as.character(coverage_data$country)
coverage_data <- coverage_data[order(coverage_data$country),]

data_name <- "coverage"
country_df <- data.frame(originalID = seq_len(nrow(coverage_data)), 
                         Location = coverage_data$country, 
                         Data = data_name)

country_df$Location <- vapply(country_df$Location, function(x) gsub(" ", "_", x),
                              character(1))
country_df$Location <- vapply(country_df$Location, function(x) gsub("&_", "", x),
                              character(1))

cleaning_dir <- "raw_data/JH_cleaning/"
write.table(country_df,paste0(cleaning_dir, data_name, "_countries.csv"),sep=",",row.names=FALSE)

my_df <- read.csv(paste0(cleaning_dir, data_name, "_countries.csv"),stringsAsFactors = FALSE)
correct_df <- read.csv(paste0(cleaning_dir, "correct", "_countries.csv"),stringsAsFactors = FALSE)

correct_idx <- vapply(country_df$Location, function(x) match(x, correct_df$Correct), double(1))

country_df$Correct <- NA
country_df[!is.na(correct_idx), "Correct"] <- correct_df[correct_idx[!is.na(correct_idx)], "Correct"]
# country_df[!is.na(correct_idx), "correctID"] <- correct_idx[!is.na(correct_idx)]

# replace incorrect country names by correct country names, by inspection

Location <- country_df$Location

replace_location <- function(Location, old_loc, new_loc) {
  Location[Location == old_loc] <- new_loc
  Location
}

Location <- replace_location(Location, "Brunei_Darussalam", "Brunei")
Location <- replace_location(Location, "Congo", "Republic_of_Congo")
Location <- replace_location(Location, "Côte_d’Ivoire", "Cote_dIvore")
Location <- replace_location(Location, "DPR_of_Korea", "Democratic_Peoples_Republic_of_Korea")
Location <- replace_location(Location, "DR_of_the_Congo", "Democratic_Republic_ of_the_ Congo")
Location <- replace_location(Location, "Iran_(Islamic_Republic_of)", "Iran")
Location <- replace_location(Location, "Lao_People’s_Democratic_Republic", "Laos")
Location <- replace_location(Location, "Macao_(China)", "Macao")
Location <- replace_location(Location, "Micronesia_(Federated_States_of)", "Micronesia")
Location <- replace_location(Location, "Papua_New_Guinea", "Papua New Guinea")
Location <- replace_location(Location, "Republic_of_Moldova", "Moldova")
Location <- replace_location(Location, "Russian_Federation", "Russia")
Location <- replace_location(Location, "St._Lucia", "Saint_Lucia")
Location <- replace_location(Location, "St._Vincent_The_Grenadines", "Saint_Vincent_and_the_Grenadines")
Location <- replace_location(Location, "The_former_Yugoslav_Republic_of_Macedonia", "Macedonia")
Location <- replace_location(Location, "Timor-Leste", "Timor_Leste")
Location <- replace_location(Location, "Trinidad_Tobago", "Trinidad_and_Tobago")
Location <- replace_location(Location, "United_Republic_of_Tanzania", "Tanzania")
Location <- replace_location(Location, "Venezuela_(Bolivarian_Republic_of)", "Venezuela")
Location <- replace_location(Location, "Viet_Nam", "Vietnam")
Location <- replace_location(Location, "Hong_Kong_(China)", "Hong_Kong")
Location <- replace_location(Location, "U.S.A.", "USA")

country_df$Correct <- Location

correct_idx <- as.numeric(vapply(country_df$Correct, function(x) match(x, correct_df$Correct), double(1)))
country_df$correctID <- correct_idx
# go get this from an older commit: earlier than 2d146a5 (2018-02-26)
countries <- read.table(paste0(cleaning_dir,"all_countries_FINAL.csv"), sep = ",", header = TRUE)
country_df <- country_df[,colnames(countries)]
countries <- rbind(countries, country_df)
countries$Data <- as.character(countries$Data)
countries <- countries[order(countries$correctID, countries$Data),]
write.table(countries,paste0(cleaning_dir, "all_countries_FINAL.csv"),sep=",",row.names=FALSE)

coverage_data <- data.frame(originalID = seq_len(nrow(coverage_data)), 
                         Location = coverage_data$country, 
                         dose_per_1000 = coverage_data$dose_per_1000)
write.table(coverage_data,"raw_data/seasonal_coverage/seasonal_coverage_clean.csv",sep=",",row.names=FALSE)
