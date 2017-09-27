library(xlsx)
n_files <- 2
filenames <- paste0("MUestimates_all_locations_", seq_len(n_files), ".xlsx")
workbooks <- lapply(filenames, loadWorkbook)
sheets <- lapply(workbooks, getSheets)
country_names <- unlist(lapply(sheets, names))

locations <- read.xlsx("locations.xlsx",1,startRow = 31, header = FALSE,colIndex = c(2,6,9))
colnames(locations) <- c("country_name", "location_code", "subregion_name")
locations <- locations[locations$location_code == 4,c(1,3)]