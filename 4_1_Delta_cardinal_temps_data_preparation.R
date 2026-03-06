library(tidyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(GGally)
library(patchwork)
library(emmeans)


# Import data ------------------------------------------------------------------

info_df <- read.csv("data/csv/serial_number_info.csv")
env = read.csv("data/csv/environment_data.csv")
bioclim = read.csv("data/csv/bioclim_data.csv")

betula = read.csv("output/csv/Thresholds_models_Betula_pubescens.csv")
alnus = read.csv("output/csv/Thresholds_models_Alnus_glutinosa.csv")
pinus = read.csv("output/csv/Thresholds_models_Pinus_sylvestris.csv")

pinus$sp = "P. sylvestris"
betula$sp = "B. pubescens"
alnus$sp = "A. glutinosa"


# Prep threshold models results df ---------------------------------------------

df = rbind(alnus, betula, pinus)
rm(alnus, betula, pinus)

df  = df %>%
  mutate(Estimate_thermal.width = Estimate_Tc - Estimate_Tb)

df <- df %>%
  rename_with(~ gsub("Std//.//.?Error", "Std.Error", .x)) %>%
  rename_with(~ gsub("t//.value", "t.value", .x)) %>%
  rename_with(~ gsub("p//.value", "p.value", .x))

# Import environmental data ----------------------------------------------------

info_df <- read.csv("data/csv/serial_number_info.csv")
env = read.csv("data/csv/environment_data.csv")
bioclim = read.csv("data/csv/bioclim_data.csv")

# Prep environmental data ------------------------------------------------------

bioclim = bioclim %>%
  rename(serial_number = Serial.No)

env = inner_join(env, bioclim)

env = env %>%
  mutate(date = as.Date(Date.Collected, "%d-%b-%y"),
         month = format(date, "%m"))

env <- info_df %>%
  dplyr::select(serial_number, Altitude) %>%
  mutate(serial_number = as.integer(serial_number)) %>%
  right_join(., env)

# Calculations environmental data ----------------------------------------------

# Get collection months
species_month_lookup_df <- env %>%
  group_by(Species) %>%
  arrange(month) %>%
  summarize(collection_months = list(unique(month)), .groups = "drop")

# For each row, calculate the mean of the climate variables for the relevant months.
env2 <- env %>% 
  left_join(species_month_lookup_df, by = "Species")

# Create columns to store mean values
env2$tas_mean    <- NA_real_
env2$tasmin_mean <- NA_real_
env2$tasmax_mean <- NA_real_

# Loop to calculate the averages for each row using only the columns corresponding to the collection months
for (i in seq_len(nrow(env2))) {
  
  months_i <- unlist(env2$collection_months[i])
  
  cols_tas    <- paste0("CHELSA_tas_",    months_i, "_1981.2010")
  cols_tasmin <- paste0("CHELSA_tasmin_", months_i, "_1981.2010")
  cols_tasmax <- paste0("CHELSA_tasmax_", months_i, "_1981.2010")
  
  env2$tas_mean[i]    <- mean(as.numeric(env2[i, cols_tas]),    na.rm = TRUE)
  env2$tasmin_mean[i] <- mean(as.numeric(env2[i, cols_tasmin]), na.rm = TRUE)
  env2$tasmax_mean[i] <- mean(as.numeric(env2[i, cols_tasmax]), na.rm = TRUE)
}

df$sp
df2 = df %>%
  rename(sp_abv = sp) %>%
  left_join(., env2, by = "serial_number")




# Calculations Deltas ----------------------------------------------------------

df_deltas <- df2 %>%
  group_by(serial_number, tas_mean, tasmin_mean, tasmax_mean, Altitude, sp, Genus, Species, across(starts_with("bio")), across(starts_with("hurs_"))) %>%
  summarise(
    delta_Estimate_Tb = Estimate_Tb[treat == "COLD"] - Estimate_Tb[treat == "NS"],
    delta_Estimate_To = Estimate_To[treat == "COLD"] - Estimate_To[treat == "NS"],
    delta_Estimate_Tc = Estimate_Tc[treat == "COLD"] - Estimate_Tc[treat == "NS"],
    delta_Estimate_ThetaT = Estimate_ThetaT[treat == "COLD"] - Estimate_ThetaT[treat == "NS"],
    delta_Estimate_thermal_width = Estimate_thermal.width[treat == "COLD"] - Estimate_thermal.width[treat == "NS"]
  )

saveRDS(df_deltas, "output/rds/df_deltas.rds")

