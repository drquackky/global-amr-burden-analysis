# import library
library(dplyr)

# load dataset
pop_data <- list.files(path = "dataset/population_dataset", pattern = "\\.csv$", full.names = TRUE)
pop_data <- do.call(rbind, lapply(pop_data, read.csv, stringsAsFactors = FALSE))
pop_data <- pop_data %>% 
   select(location_id, location_name, age_name, year, val) %>%
   rename(popsize = val)
