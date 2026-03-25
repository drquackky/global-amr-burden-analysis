# =======================================================================
# skip to run these codes because the file is very large
# "[file_name]_4run.csv" is aggregated dataset ready for the analysis
# 
# each csv files contains one of the 4 burden indicator of almost 30 bug-drug combination
# in each age categories (~23) of 204 countries and territories from 1990-2021
# 4(indicators) x 30 (bug-drug) x 23 (age group) x 204 x 31 (year) = 17 mil. data rows 
# =======================================================================

# -----------------------------------------------------------------
# Load libraries
# -----------------------------------------------------------------
library(dplyr)
library(tibble)
library(readr)
library(epitools) # use for calculating direct age standardization



# -----------------------------------------------------------------
# Read datasets 
# -----------------------------------------------------------------

# read SDI (socio-demographic index) data (by country and year)
df_sdi <- read.csv("dataset/sdi_2023.csv") %>%
   select (location_id, location_name, year_id, mean_value) %>%
   rename (year = year_id, sdi_val = mean_value)


# read transformed data on global antimicrobial resistance burden
df_amr <- read_csv("dataset/df_amrxpop.csv")


# -------------------------------------------------------------------------------------------------------
# 4. Calculate age-standardized rate
# -------------------------------------------------------------------------------------------------------

# In order to have a fair comparison between countries, we are using WHO world standard population as weight
# https://seer.cancer.gov/stdpopulations/world.who.html
# create a df according to the link
who_std <- tibble(
   Age = c("0 to 4","5 to 9","10 to 14","15 to 19","20 to 24",
           "25 to 29","30 to 34","35 to 39","40 to 44","45 to 49",
           "50 to 54","55 to 59","60 to 64","65 to 69","70 to 74",
           "75 to 79","80 to 84","85 to 89","90 to 94","95 plus"),
   std_pop = c(8.860,8.690,8.600,8.470,8.220,
           7.930,7.610,7.150,6.590,6.040,
           5.370,4.550,3.720,2.960,2.210,
           1.520,0.910,0.440,0.150,0.040 + 0.005)
)


# Prepare data for calculating age standardized rate
# =======================================================================

df_amrburden <- df_amr %>%
   # Collapse infant ages (1-5 months, 6-11 months ...) into 0 to 4 (use same label as who_std)
   mutate(
      Age = case_when(
         Age %in% c("1-5 months","6-11 months","12 to 23 months","2 to 4") ~ "0 to 4",
         TRUE ~ Age)
   ) %>%
   # sum of the age specific rate and population number by each country, year, age-group, priority group and indicators
   group_by(Location, Year, metric, priority_group, Age) %>%
   summarise(
      deaths = sum(amr_burden),
      population = sum(pop_size),
      .groups = "drop"
   ) %>%
   # join data with WHO world standard population
   left_join(who_std, by = "Age")



# Merge a data with the social development index for group categorization
# =======================================================================

# check if the value of location in two datasets are different
setdiff(unique(df_amrburden$Location), unique(df_sdi$location_name)) # no difference

# categorize sdi value into group according to GBD guideline
# https://ghdx.healthdata.org/record/gbd-2023-socio-demographic-index-sdi
df_sdi <- df_sdi %>% 
   mutate(sdi_cat = case_when(
      sdi_val >= 0 & sdi_val < 0.46581580319161997 ~ "low",
      sdi_val >= 0.46581580319161997 & sdi_val < 0.6188294452454329 ~ "middle_low",
      sdi_val >= 0.6188294452454329 & sdi_val < 0.7119746219361235 ~ "middle",
      sdi_val >= 0.7119746219361235 & sdi_val < 0.8102959891918925 ~ "middle_high",
      sdi_val >= 0.8102959891918925 ~ "high"
   )) %>%
   select (-1)

# merging the two dataset: sdi and amr burden
df_amrburden_merged <- df_amrburden %>%
   left_join(df_sdi, by = c("Location" = "location_name",
                            "Year" = "year"),
             relationship = "many-to-many")



# ==================================================================================
# create an age-standardized rate dataset for mixed linear model
# to compare burden between each SDI categories using estimated marginal means
# ==================================================================================

df_asr_lmm <- df_amrburden_merged %>%
   
   # group the data by location, year, metric, priority group and sdi_cat
   # we want to calculate one ASR + CI for each of these combinations
   group_by(Location, Year, metric, priority_group, sdi_cat) %>%
   
   summarise(
      # Calculate age-standardized rate for each group
      # ageadjust.direct() returns multiple values (adj.rate, lci, uci)
      # Wrapping it in list() keeps all values together for each group
      adj = list(
         ageadjust.direct(
            count = deaths,      # vector of age-specific deaths for this group
            pop = population,    # vector of age-specific population for this group
            stdpop = std_pop     # standard population weights
         )
      ),
      .groups = "drop"        # remove grouping after summarise
   ) %>%
   
   mutate(
      # Extract each value from the list column into separate numeric columns
      # Multiply by 100,000 to get rate per 100k population
      ASR = sapply(adj, function(x) x["adj.rate"] * 100000),    # age-standardized rate
      lower_CI = sapply(adj, function(x) x["lci"] * 100000),    # lower 95% CI
      upper_CI = sapply(adj, function(x) x["uci"] * 100000)     # upper 95% CI
   ) %>%
   
   # Remove the intermediate list column; we only need numeric results
   select(-adj)

# save csv file
write.csv(df_asr_lmm, "dataset/asr_lmm.csv", row.names = FALSE)




# ==================================================================================
# create age-standardized rate dataset for joinpoint regression analysis
# countries were group into 5 SDI categories and their ASR were calculated
# ==================================================================================

df_asr_jp <- df_amrburden_merged %>%
   
   # group the data by location, year, metric, priority group, except SDI_CAT!!!
   group_by(Year, metric, priority_group, sdi_cat) %>%
   summarise(
      adj = list(
         ageadjust.direct(
            count = deaths,
            pop = population,
            stdpop = std_pop
         )
      ),
      .groups = "drop"
   ) %>%
   mutate(
      ASR = sapply(adj, function(x) x["adj.rate"] * 100000),
      lower_CI = sapply(adj, function(x) x["lci"] * 100000),
      upper_CI = sapply(adj, function(x) x["uci"] * 100000)
   ) %>%
   select(-adj)

# save csv file
write.csv(df_asr_jp, "dataset/asr_jp.csv", row.names = FALSE)

