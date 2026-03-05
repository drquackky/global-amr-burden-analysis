library(dplyr)
library(tibble)
library(readr)

# read SDI (socio-demographic index) data (by country and year)
# --------------------------------------------------------------------
df_sdi <- read.csv("dataset/sdi_2023.csv") %>%
   select (location_id, location_name, year_id, mean_value) %>%
   rename (year = year_id, sdi_val = mean_value)

# read transformed data on global antimicrobial resistance burden
# -------------------------------------------------------------------
df_amrburden <- read_csv("dataset/df_amrxpop.csv")


# =======================================================================================================
# 4. Calculate age-standardized rate for each pathogen priority group by Location and metrics indicators 
# =======================================================================================================

# In order to compare between countries, we are using WHO world standard population as weight
# https://seer.cancer.gov/stdpopulations/world.who.html
# ---------------------------------------------------------------------------------------------
# create a df according to the link
who_std <- tibble(
   Age = c("0 to 4","5 to 9","10 to 14","15 to 19","20 to 24",
           "25 to 29","30 to 34","35 to 39","40 to 44","45 to 49",
           "50 to 54","55 to 59","60 to 64","65 to 69","70 to 74",
           "75 to 79","80 to 84","85 to 89","90 to 94","95 plus"),
   pct = c(8.860,8.690,8.600,8.470,8.220,
           7.930,7.610,7.150,6.590,6.040,
           5.370,4.550,3.720,2.960,2.210,
           1.520,0.910,0.440,0.150,0.040 + 0.005)
) %>%
   # WHO age weights (percent -> proportion (sum up to 1) for standardization)
   mutate(weight = pct / sum(pct)) %>%
   select(Age, weight)


# Prepare data for calculating age standardized rate
# -----------------------------------------------------------------------

# set the rate per 100.000 population
r <- 100000

# Collapse infant ages (1-5 months, 6-11 months ...) into 0 to 4 (use same label as who_std)
df_asr <- df_amrburden %>%
   mutate(
      # case_when() - condition to set infant ages to 0 to 4, the rest "TRUE ~ Age" kept the same
      Age = case_when(
         Age %in% c("1-5 months","6-11 months","12 to 23 months","2 to 4") ~ "0 to 4",
         TRUE ~ Age)) %>%
   
   # aggregate by bug drug, metrics, location and year
   group_by(Location, Year, bug_drug, metric, priority_group, Age) %>%
   
   # collapse the infant age after grouping
   summarise(
      deaths = sum(amr_burden),
      population = sum(pop_size),
      .groups = "drop") %>%
   
   # we need to sum the number across each bug drug by their priority group
   group_by(Location, Year, metric, priority_group, Age) %>%
   summarise(
      deaths = sum(deaths),
      population = sum(population),
      .groups = "drop") %>%
   
   # calculate age-specific rate by priority group
   mutate(rate = deaths / population * r) %>%
   
   # join with the standard population table
   left_join(who_std, by = "Age") %>%
   
   # Age variable not included in group by because it is used for calculating age standardized rate
   group_by(Location, Year, metric, priority_group) %>%
   summarise(
      ASR = sum(rate * weight, na.rm = TRUE),
      .groups = "drop")


# save csv file
# -----------------------------------------------------------------
write.csv(df_asr, "dataset/asr_amrburden.csv", row.names = FALSE)


# ========================================================================
# 5. Categorize each country into an SDI group for each year.
# ========================================================================

setdiff(unique(df_asr$Location), unique(df_sdi$location_name))
df_sdi <- df_sdi %>% 
   mutate(sdi_cat = case_when(
      sdi_val >= 0 & sdi_val < 0.46581580319161997 ~ "low",
      sdi_val >= 0.46581580319161997 & sdi_val < 0.6188294452454329 ~ "middle_low",
      sdi_val >= 0.6188294452454329 & sdi_val < 0.7119746219361235 ~ "middle",
      sdi_val >= 0.7119746219361235 & sdi_val < 0.8102959891918925 ~ "middle_high",
      sdi_val >= 0.8102959891918925 ~ "high"
   )) %>%
   select (-1)


df_final <- df_asr %>%
   left_join(df_sdi, by = c("Location" = "location_name",
                            "Year" = "year"))


write.csv(df_final, "dataset/df_final.csv", row.names = FALSE)

