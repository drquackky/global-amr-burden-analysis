
# =======================================================================
# YOU CAN SKIP RUNNING THESE CODE BECAUSE THE FILE IS VERY LARGE
# ALL TRANSFORMED DATA ARE READY TO RUN IN 1B FILE
# =======================================================================



# =======================================================
# 1. IMPORT LIBRARY
# =======================================================
library(dplyr) # for data manipulation
library(readr) # read_csv() runs faster with larger data files, return tibble, does not convert string to column
library (ggplot2) # for data visualization
library (data.table) # rbindlist() run faster with larger data files


# =======================================================
# 2. LOAD DATASET
# =======================================================

# read GBD population estimates data (by country and age group)
# --------------------------------------------------------------------
df_popsize <- read.csv("dataset/population_estimates_gbd.csv") %>%
   select(location_id, location_name, age_name, year, val) %>%
   rename(pop_size = val, Age = age_name)



# read AMR Burden Measures data by Age Group and Country
# --------------------------------------------------------------------
# list all the csv files in the folder to an R object
# "recursive = TRUE" allows R to search through subfolders
list_amr <- list.files(path = "dataset/amr_dataset" , pattern = "\\.csv$", 
                        full.names = TRUE, recursive = TRUE)


# lapply() applies a function to each element of a list and returns the results as a list
# function (f) - tell the code to read one file during each loop
df_amr <- lapply(
   list_amr,
   function(f) read_csv(f, show_col_types = FALSE, progress = FALSE)
) %>%
   bind_rows() # row wise data merge



# ===============================================================
# 3. Merge the two datasets and transform variables for analysis
# ===============================================================

# Before merging, verify that value in categories match across datasets

# Identify mismatched age labels between AMR and population data
# -----------------------------------------------------------------------
setdiff(unique(df_amr$Age), unique(df_popsize$Age))

## Create a common age label in pop data
df_popsize <- df_popsize %>%
   mutate(
      Age = case_when(
         Age == "1-5 months" ~ "1-5 months",
         Age == "6-11 months" ~ "6-11 months",
         Age == "12-23 months" ~ "12 to 23 months",
         Age == "2-4 years" ~ "2 to 4",
         Age == "5-9 years" ~ "5 to 9",
         Age == "10-14 years" ~ "10 to 14",
         Age == "15-19 years" ~ "15 to 19",
         Age == "20-24 years" ~ "20 to 24",
         Age == "25-29 years" ~ "25 to 29",
         Age == "30-34 years" ~ "30 to 34",
         Age == "35-39 years" ~ "35 to 39",
         Age == "40-44 years" ~ "40 to 44",
         Age == "45-49 years" ~ "45 to 49",
         Age == "50-54 years" ~ "50 to 54",
         Age == "55-59 years" ~ "55 to 59",
         Age == "60-64 years" ~ "60 to 64",
         Age == "65-69 years" ~ "65 to 69",
         Age == "70-74 years" ~ "70 to 74",
         Age == "75-79 years" ~ "75 to 79",
         Age == "80-84 years" ~ "80 to 84",
         Age == "85-89 years" ~ "85 to 89",
         Age == "90-94 years" ~ "90 to 94",
         Age == "95+ years" ~ "95 plus",
         TRUE ~ NA_character_
      )
   )

# Identify mismatched Location between AMR and population data
# -------------------------------------------------------------------
# show Location in df_amr_agg that are NOT in df_pop
setdiff(unique(df_amr$Location), unique(df_popsize$location_name))

# changing a "Taiwan (Province of China)" in df_amr to "Taiwan" to match it with df_popsize
df_amr <- df_amr %>%
   mutate(Location = ifelse(Location == "Taiwan (Province of China)", "Taiwan", Location))



# Join two data on location, age-group and year
# -------------------------------------------------------------------------
df_amr_agg <- df_amr%>% left_join(
   df_popsize, by = c("Location" = "location_name", "Age", "Year" = "year")
) %>%
   # merge two column together
   mutate(bug_drug = paste(Pathogen, `Antibiotic class`, sep = "_")) %>%
   mutate(metric = paste(Counterfactual, Measure, sep = "_")) %>%
   
   # select only columns needed for analysis
   select(-c(3:8), -c(11:13)) %>%
   
   # rename value to amr_burden
   rename(amr_burden = Value) %>%
   
   # creating a priority group bug drug combination according to WHO document
   # https://www.who.int/publications/i/item/9789240093461
   mutate(priority_group = case_when(
      # Critical priority
      bug_drug %in% c(
         "Klebsiella pneumoniae_Carbapenems",
         "Escherichia coli_Third-generation cephalosporins",
         "Acinetobacter baumannii_Carbapenems",
         "Mycobacterium tuberculosis_Resistance to one or more antibiotics",
         "Mycobacterium tuberculosis_Extensive drug resistance in TB",
         "Mycobacterium tuberculosis_Multi-drug resistance excluding extensive drug resistance in TB",
         "Escherichia coli_Carbapenems",
         "Klebsiella pneumoniae_Third-generation cephalosporins",
         "Enterobacter spp._Carbapenems",
         "Citrobacter spp._Third-generation cephalosporins",
         "Proteus spp._Third-generation cephalosporins",
         "Serratia spp._Third-generation cephalosporins",
         "Morganella spp._Third-generation cephalosporins"
      ) ~ "Critical",
      
      # High priority
      bug_drug %in% c(
         "Salmonella enterica serovar Typhi_Fluoroquinolones",
         "Shigella spp._Fluoroquinolones",
         "Enterococcus faecium_Vancomycin",
         "Pseudomonas aeruginosa_Carbapenems",
         "Non-typhoidal Salmonella_Fluoroquinolones",
         "Neisseria gonorrhoeae_Fluoroquinolones",
         "Staphylococcus aureus_Methicillin",
         "Neisseria gonorrhoeae_Third-generation cephalosporins"
      ) ~ "High",
      
      # Medium priority
      bug_drug %in% c(
         "Group A Streptococcus_Macrolides",
         "Streptococcus pneumoniae_Macrolides",
         "Haemophilus influenzae_Aminopenicillin",
         "Group B Streptococcus_Penicillin"
      ) ~ "Medium",
      
      # Default if no match
      TRUE ~ "Not priority"
   ))

# write a new csv file for a transformed data
write.csv(df_amr_agg, file = "dataset/df_amrxpop.csv", row.names = FALSE)





   
   