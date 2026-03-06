
library(dplyr)
library(tidyr)
library(lubridate)

## Function 

# Define non-linear rate function based on GRT.Ex() parameters
germ_rate <- function(temperature, Tb, Tc, k) {
  ifelse(temperature < Tb,
         0,
         (pmax(temperature, Tb) - Tb) *
           (1 - exp(k * (pmin(temperature, Tc) - Tc))) /
           (1 - exp(k * (Tb - Tc)))
  )
}


# Import data ------------------------------------------------------------------

env = read.csv("data/csv/environment_data.csv")

betula = read.csv("output/csv/Thresholds_models_Betula_pubescens.csv")
betula$Genus = "Betula"

alnus = read.csv("output/csv/Thresholds_models_Alnus_glutinosa.csv")
alnus$Genus = "Alnus"

pinus = read.csv("output/csv/Thresholds_models_Pinus_sylvestris.csv")
pinus$Genus = "Pinus"

df = rbind(alnus, betula, pinus)
rm(alnus, betula, pinus)

df = df %>% 
  # dplyr::filter(g == "50%" & mod == "modGR_Ex") %>%
  rename(Tb = Estimate_Tb,
        Tc = Estimate_Tc,
        To = Estimate_To,
        ThetaT = Estimate_ThetaT,
        k = Estimate_k)


## Averaging Tb for collections where treatments led to difference in Tb < 2C
df$Tb_calculated <- df$Tb

df = df %>% 
  group_by(Genus, serial_number, mod, g) %>%
  summarise(Tb_average = mean(Tb)) %>%
  right_join(., df)

diff <- df %>%
  group_by(Genus, serial_number, mod, g) %>%
  summarise(
    Tb_diff = abs(diff(Tb))
  ) %>%
  filter(!is.na(Tb_diff), Tb_diff < 2)

df$Tb = ifelse(df$serial_number %in% diff$serial_number, df$Tb_average, df$Tb)

df = df %>% ungroup() %>%
  dplyr::select(-Tb_average, -Genus)


output = merge(df, env, by = "serial_number")

output$Date.Collected = as.Date(output$Date.Collected, format = "%d-%b-%y")

names(output)


# ------------------------------------------------------------------------------
# CHELSA dataset 1981-2010
## calculate heat sum from air temperature by hour adjusted by max - START 01 Jan

# Create a sequence of dates from January 1st to December 31st for a non-leap year
start_date <- as.Date("2023-01-01") # any non leap year
end_date <- as.Date("2023-12-31")

datseq <- seq.Date(from = start_date, to = end_date, by = "day")
month_year <- format(datseq, "%Y-%m")
datseq <- format(datseq, format = "%Y-%m-%d")

# Initialize an empty dataframe to store the results
heatsum <- data.frame(Genus = as.character(), Species = as.character(), serial_number = as.character(), treat = as.character(), mod = as.character(), date = as.Date(character()), hour = integer(), temperature = numeric(), hourly_heat = numeric(), hourly_heat_sum = numeric(), cumulative_heat_sum =  numeric(), q50 =  numeric(), Tb = numeric())

output$date_reach_TT_hrs <- NA
output$date_reach_TT_hrs = as.Date(output$date_reach_TT_hrs)

for (i in 1:nrow(output)) {
  data <- output[i, ]
  
  hourly_temperatures <- data.frame(date = as.Date(character()), hour = integer(), temperature = numeric())
  
  fluct <- data.frame(date = datseq, month_year = month_year)
  
  # Calculate the number of days in each month
  days_in_months <- fluct %>%
    group_by(month_year) %>%
    summarise(days = n()) %>%
    ungroup()
  
  # Create a vector of environmental temperatures for each day
  fluct$env_tm <- c(
    rep(data$CHELSA_tas_01_1981.2010, each = days_in_months$days[1]),
    rep(data$CHELSA_tas_02_1981.2010, each = days_in_months$days[2]),
    rep(data$CHELSA_tas_03_1981.2010, each = days_in_months$days[3]),
    rep(data$CHELSA_tas_04_1981.2010, each = days_in_months$days[4]),
    rep(data$CHELSA_tas_05_1981.2010, each = days_in_months$days[5]),
    rep(data$CHELSA_tas_06_1981.2010, each = days_in_months$days[6]),
    rep(data$CHELSA_tas_07_1981.2010, each = days_in_months$days[7]),
    rep(data$CHELSA_tas_08_1981.2010, each = days_in_months$days[8]),
    rep(data$CHELSA_tas_09_1981.2010, each = days_in_months$days[9]),
    rep(data$CHELSA_tas_10_1981.2010, each = days_in_months$days[10]),
    rep(data$CHELSA_tas_11_1981.2010, each = days_in_months$days[11]),
    rep(data$CHELSA_tas_12_1981.2010, each = days_in_months$days[12])
  )
  
  fluct$env_tmin <- c(
    rep(data$CHELSA_tasmin_01_1981.2010, each = days_in_months$days[1]),
    rep(data$CHELSA_tasmin_02_1981.2010, each = days_in_months$days[2]),
    rep(data$CHELSA_tasmin_03_1981.2010, each = days_in_months$days[3]),
    rep(data$CHELSA_tasmin_04_1981.2010, each = days_in_months$days[4]),
    rep(data$CHELSA_tasmin_05_1981.2010, each = days_in_months$days[5]),
    rep(data$CHELSA_tasmin_06_1981.2010, each = days_in_months$days[6]),
    rep(data$CHELSA_tasmin_07_1981.2010, each = days_in_months$days[7]),
    rep(data$CHELSA_tasmin_08_1981.2010, each = days_in_months$days[8]),
    rep(data$CHELSA_tasmin_09_1981.2010, each = days_in_months$days[9]),
    rep(data$CHELSA_tasmin_10_1981.2010, each = days_in_months$days[10]),
    rep(data$CHELSA_tasmin_11_1981.2010, each = days_in_months$days[11]),
    rep(data$CHELSA_tasmin_12_1981.2010, each = days_in_months$days[12])
  )
  
  fluct$env_tmax <- c(
    rep(data$CHELSA_tasmax_01_1981.2010, each = days_in_months$days[1]),
    rep(data$CHELSA_tasmax_02_1981.2010, each = days_in_months$days[2]),
    rep(data$CHELSA_tasmax_03_1981.2010, each = days_in_months$days[3]),
    rep(data$CHELSA_tasmax_04_1981.2010, each = days_in_months$days[4]),
    rep(data$CHELSA_tasmax_05_1981.2010, each = days_in_months$days[5]),
    rep(data$CHELSA_tasmax_06_1981.2010, each = days_in_months$days[6]),
    rep(data$CHELSA_tasmax_07_1981.2010, each = days_in_months$days[7]),
    rep(data$CHELSA_tasmax_08_1981.2010, each = days_in_months$days[8]),
    rep(data$CHELSA_tasmax_09_1981.2010, each = days_in_months$days[9]),
    rep(data$CHELSA_tasmax_10_1981.2010, each = days_in_months$days[10]),
    rep(data$CHELSA_tasmax_11_1981.2010, each = days_in_months$days[11]),
    rep(data$CHELSA_tasmax_12_1981.2010, each = days_in_months$days[12])
  )
  
      # Calculate hourly temperatures for each day
      for (j in 1:nrow(fluct)) {
        date <- as.Date(fluct$date[j], format = "%Y-%m-%d")
        T_mean <- fluct$env_tm[j]
        T_min <- fluct$env_tmin[j]
        T_max <- fluct$env_tmax[j]
        t_peak <- 14  # Time of peak temperature (2 PM)
        amplitude <- (T_max - T_mean) #(T_max - T_min) / 2
        
        # Time array (24 points from 0 to 24 hours)
        time <- seq(0, 23)
        
        # Sinusoidal temperature model
        temperature <- T_mean + amplitude * cos(2 * pi * (time - t_peak) / 24)
        
        # Create a data frame for the current day's hourly temperatures
        df_hours <- data.frame(date = date, hour = time, temperature = temperature)
        
        # Append to the main dataframe
        hourly_temperatures <- rbind(hourly_temperatures, df_hours)
      }
  
  Tb = data$Tb
  TT = data$ThetaT*24 # convert to degrees per hour
  Tc = data$Tc
  k = data$k
  
  ## Function to calculate the heat sum and the date reaching Tb for each row without decaying term
  # hourly_temperatures = hourly_temperatures %>%
  #    mutate(
  #     hourly_heat = (temperature - Tb),
  #     hourly_heat_sum = ifelse(hourly_heat < 0, 0, hourly_heat),
  #     cumulative_heat_sum = cumsum(hourly_heat_sum))
  
  # Function to calculate the heat sum and the date reaching Tb for each row including decaying term
  hourly_temperatures <- hourly_temperatures %>%
    mutate(
      heat_increment  = germ_rate(temperature, Tb, Tc, k),
      cumulative_heat_sum = cumsum(heat_increment)
    )
  
  
  # Find the day when cumulative heat sum exceeds ThetaT
  date_reach_ThetaT <- hourly_temperatures[which(hourly_temperatures$cumulative_heat_sum >= TT)[1],]$date
  
  output[i,]$date_reach_TT_hrs <- date_reach_ThetaT
  
  hs = hourly_temperatures
  hs$q50 = ifelse(hourly_temperatures$cumulative_heat_sum <= TT, 0, 1)
  hs$serial_number <- data$serial_number
  hs$treat <- data$treat
  hs$mod <- data$mod
  hs$Genus <- data$Genus
  hs$Species <- data$Species   
  hs$Tb <- Tb
  
  heatsum = rbind(heatsum, hs)
  
}

heatsum_d = heatsum %>%
  group_by(serial_number, treat, mod, Genus, Species, date) %>%
  summarise(
    cumulative_heat_sum =  sum(cumulative_heat_sum/24),
    q50 = min(q50),     
    Tb = unique(Tb),
    T_min = min(temperature), 
    T_max = max(temperature), 
    T_mean = mean(temperature)
  )

write.csv(heatsum_d, "output/csv/heatsum_from_Jan01_present.csv", row.names = FALSE)


# ------------------------------------------------------------------------------
# CHELSA dataset 1981-2010
## calculate heat sum from air temperature by hour adjusted by max - START at collection date

output$date_reach_TT_hrs_datecollect <- NA
output$date_reach_TT_hrs_datecollect = as.Date(output$date_reach_TT_hrs_datecollect)

output$cumulative_heat_sum_dec31 <- NA

output$stand_start_date <- NA
output$stand_start_date <- as.Date(output$stand_start_date)

# Initialize an empty dataframe to store the results
heatsum <- data.frame(Genus = as.character(), Species = as.character(), serial_number = as.character(), treat = as.character(), mod = as.character(), date = as.Date(character()), hour = integer(), temperature = numeric(), hourly_heat = numeric(), hourly_heat_sum = numeric(), cumulative_heat_sum =  numeric(), q50 =  numeric(), Tb = numeric())


for (i in 1:nrow(output)) {
  
  data <- output[i, ]
  
  # Dates
  start_date <- as.Date(data$Date.Collected)
  
  if (month(start_date) < 5) {
    start_date <- ymd(paste0("2023-", month(start_date), "-", day(start_date)))
  } else {
    start_date <- ymd(paste0("2022-", month(start_date), "-", day(start_date)))
  }

  end_date <- as.Date("2023-12-31")
  
  # Create a sequence of dates from January 1st to December 31st for a non-leap year
  datseq <- seq.Date(from = start_date, to = end_date, by = "day")
  month_year <- format(datseq, "%Y-%m")
  datseq <- format(datseq, format = "%Y-%m-%d")
  
  
  
  hourly_temperatures <- data.frame(date = as.Date(character()), hour = integer(), temperature = numeric())
  
  fluct <- data.frame(date = datseq, month_year = month_year)
  
  # Calculate the number of days in each month
  days_in_months <- fluct %>%
    group_by(month_year) %>%
    summarise(days = n()) %>%
    ungroup()
  

  
  # Environmental temperatures for each month
  monthly_temperatures <- list(
    data$CHELSA_tas_01_1981.2010, data$CHELSA_tas_02_1981.2010, data$CHELSA_tas_03_1981.2010,
    data$CHELSA_tas_04_1981.2010, data$CHELSA_tas_05_1981.2010, data$CHELSA_tas_06_1981.2010,
    data$CHELSA_tas_07_1981.2010, data$CHELSA_tas_08_1981.2010, data$CHELSA_tas_09_1981.2010,
    data$CHELSA_tas_10_1981.2010, data$CHELSA_tas_11_1981.2010, data$CHELSA_tas_12_1981.2010
  )
  
  monthly_temperatures_min <- list(
    data$CHELSA_tasmin_01_1981.2010, data$CHELSA_tasmin_02_1981.2010, data$CHELSA_tasmin_03_1981.2010,
    data$CHELSA_tasmin_04_1981.2010, data$CHELSA_tasmin_05_1981.2010, data$CHELSA_tasmin_06_1981.2010,
    data$CHELSA_tasmin_07_1981.2010, data$CHELSA_tasmin_08_1981.2010, data$CHELSA_tasmin_09_1981.2010,
    data$CHELSA_tasmin_10_1981.2010, data$CHELSA_tasmin_11_1981.2010, data$CHELSA_tasmin_12_1981.2010
  )
  
  monthly_temperatures_max <- list(
    data$CHELSA_tasmax_01_1981.2010, data$CHELSA_tasmax_02_1981.2010, data$CHELSA_tasmax_03_1981.2010,
    data$CHELSA_tasmax_04_1981.2010, data$CHELSA_tasmax_05_1981.2010, data$CHELSA_tasmax_06_1981.2010,
    data$CHELSA_tasmax_07_1981.2010, data$CHELSA_tasmax_08_1981.2010, data$CHELSA_tasmax_09_1981.2010,
    data$CHELSA_tasmax_10_1981.2010, data$CHELSA_tasmax_11_1981.2010, data$CHELSA_tasmax_12_1981.2010
  )
  
  start_month <- month(start_date)
  
  if (start_month < 5) {
    ordered_months <- c(start_month:12)
  } else {
    ordered_months <- c(start_month:12, 1:12)
  }
  

  
  env_tm <- unlist(lapply(seq_along(ordered_months), function(x) {
    rep(monthly_temperatures[[ordered_months[x]]], each = days_in_months$days[x])
  }))
  
  env_tmin <- unlist(lapply(seq_along(ordered_months), function(x) {
    rep(monthly_temperatures_min[[ordered_months[x]]], each = days_in_months$days[x])
  }))
  
  env_tmax <- unlist(lapply(seq_along(ordered_months), function(x) {
    rep(monthly_temperatures_max[[ordered_months[x]]], each = days_in_months$days[x])
  }))
  
  fluct$env_tm <- env_tm
  fluct$env_tmin <- env_tmin
  fluct$env_tmax <- env_tmax
  

  # Calculate hourly temperatures for each day
  for (j in 1:nrow(fluct)) {
    date <- as.Date(fluct$date[j], format = "%Y-%m-%d")
    T_mean <- fluct$env_tm[j]
    T_min <- fluct$env_tmin[j]
    T_max <- fluct$env_tmax[j]
    t_peak <- 14  # Time of peak temperature (2 PM)
    amplitude <- (T_max - T_mean) # (T_max - T_min) / 2
    
    # Time array (24 points from 0 to 24 hours)
    time <- seq(0, 23)
    
    # Sinusoidal temperature model
    temperature <- T_mean + amplitude * cos(2 * pi * (time - t_peak) / 24)
    
    # Create a data frame for the current day's hourly temperatures
    df_hours <- data.frame(date = date, hour = time, temperature = temperature)
    
    # Append to the main dataframe
    hourly_temperatures <- rbind(hourly_temperatures, df_hours)
  }
  
  Tb = data$Tb
  TT = data$ThetaT*24 # convert to degrees per hour
  
  # # Function to calculate the heat sum and the date reaching ThetaT for each row (OLD)
  # hourly_temperatures = hourly_temperatures %>%
  #   mutate(
  #     hourly_heat = (temperature - Tb),
  #     hourly_heat_sum = ifelse(hourly_heat < 0, 0, hourly_heat),
  #     cumulative_heat_sum = cumsum(hourly_heat_sum))
  
  # Function to calculate the heat sum and the date reaching Tb for each row including decaying term
  hourly_temperatures <- hourly_temperatures %>%
    mutate(
      heat_increment  = germ_rate(temperature, Tb, Tc, k),
      cumulative_heat_sum = cumsum(heat_increment)
    )
  
  
  # Find the day when cumulative heat sum exceeds ThetaT
  date_reach_ThetaT <- hourly_temperatures[which(hourly_temperatures$cumulative_heat_sum >= TT)[1],]$date
  
  output[i,]$date_reach_TT_hrs_datecollect <- date_reach_ThetaT
  
  if (month(start_date) >= 4) {
    # Find the index of the rows for the date 2022-12-31
    index <- which(hourly_temperatures$date == as.Date("2022-12-31") & hourly_temperatures$hour == 23)
    
    # Extract the cumulative_heat_sum for the last hour of 2022-12-31
    cumulative_heat_sum_dec31 <- hourly_temperatures$cumulative_heat_sum[index]/24
  } else {
    cumulative_heat_sum_dec31 <- 0
  }
  
 
  output[i,]$cumulative_heat_sum_dec31 <- cumulative_heat_sum_dec31
  
  output[i,]$stand_start_date <- as.Date(start_date)
  
  hs = hourly_temperatures
  hs$q50 = ifelse(hourly_temperatures$cumulative_heat_sum <= TT, 0, 1)
  
  hs$serial_number <- data$serial_number
  hs$treat <- data$treat
  hs$mod <- data$mod
  hs$Genus <- data$Genus
  hs$Species <- data$Species   
  hs$Tb <- Tb
  
  heatsum = rbind(heatsum, hs)
  
}

heatsum_d = heatsum %>%
  group_by(serial_number, treat, mod, Genus, Species, date) %>%
  summarise(
    cumulative_heat_sum =  sum(cumulative_heat_sum/24),
    q50 = min(q50),     
    Tb = unique(Tb),
    T_min = min(temperature), 
    T_max = max(temperature), 
    T_mean = mean(temperature)
  )

write.csv(heatsum_d, "output/csv/heatsum_from_collect_present.csv", row.names = FALSE)


# ------------------------------------------------------------------------------
# Optimistic FUTURE (ssp126)
## calculate heat sum from air temperature by hour adjusted by max - START 01 Jan

# Constants
start_date <- as.Date("2023-01-01")
end_date <- as.Date("2023-12-31")

# Create a sequence of dates from January 1st to December 31st for a non-leap year
datseq <- seq.Date(from = start_date, to = end_date, by = "day")
month_year <- format(datseq, "%Y-%m")
datseq <- format(datseq, format = "%Y-%m-%d")

# Initialize an empty dataframe to store the results
heatsum <- data.frame(Genus = as.character(), Species = as.character(), serial_number = as.character(), treat = as.character(), mod = as.character(), date = as.Date(character()), hour = integer(), temperature = numeric(), hourly_heat = numeric(), hourly_heat_sum = numeric(), cumulative_heat_sum =  numeric(), q50 =  numeric(), Tb = numeric())

output$date_reach_TT_hrs_ssp126 <- NA
output$date_reach_TT_hrs_ssp126 = as.Date(output$date_reach_TT_hrs_ssp126)

for (i in 1:nrow(output)) {
  data <- output[i, ]
  
  hourly_temperatures <- data.frame(date = as.Date(character()), hour = integer(), temperature = numeric())
  
  fluct <- data.frame(date = datseq, month_year = month_year)
  
  # Calculate the number of days in each month
  days_in_months <- fluct %>%
    group_by(month_year) %>%
    summarise(days = n()) %>%
    ungroup()
  
  # Create a vector of environmental temperatures for each day
  fluct$env_tm <- c(
    rep(data$ssp126_tas_01, each = days_in_months$days[1]),
    rep(data$ssp126_tas_02, each = days_in_months$days[2]),
    rep(data$ssp126_tas_03, each = days_in_months$days[3]),
    rep(data$ssp126_tas_04, each = days_in_months$days[4]),
    rep(data$ssp126_tas_05, each = days_in_months$days[5]),
    rep(data$ssp126_tas_06, each = days_in_months$days[6]),
    rep(data$ssp126_tas_07, each = days_in_months$days[7]),
    rep(data$ssp126_tas_08, each = days_in_months$days[8]),
    rep(data$ssp126_tas_09, each = days_in_months$days[9]),
    rep(data$ssp126_tas_10, each = days_in_months$days[10]),
    rep(data$ssp126_tas_11, each = days_in_months$days[11]),
    rep(data$ssp126_tas_12, each = days_in_months$days[12])
  )
  
  fluct$env_tmin <- c(
    rep(data$ssp126_tasmin_01, each = days_in_months$days[1]),
    rep(data$ssp126_tasmin_02, each = days_in_months$days[2]),
    rep(data$ssp126_tasmin_03, each = days_in_months$days[3]),
    rep(data$ssp126_tasmin_04, each = days_in_months$days[4]),
    rep(data$ssp126_tasmin_05, each = days_in_months$days[5]),
    rep(data$ssp126_tasmin_06, each = days_in_months$days[6]),
    rep(data$ssp126_tasmin_07, each = days_in_months$days[7]),
    rep(data$ssp126_tasmin_08, each = days_in_months$days[8]),
    rep(data$ssp126_tasmin_09, each = days_in_months$days[9]),
    rep(data$ssp126_tasmin_10, each = days_in_months$days[10]),
    rep(data$ssp126_tasmin_11, each = days_in_months$days[11]),
    rep(data$ssp126_tasmin_12, each = days_in_months$days[12])
  )
  
  fluct$env_tmax <- c(
    rep(data$ssp126_tasmax_01, each = days_in_months$days[1]),
    rep(data$ssp126_tasmax_02, each = days_in_months$days[2]),
    rep(data$ssp126_tasmax_03, each = days_in_months$days[3]),
    rep(data$ssp126_tasmax_04, each = days_in_months$days[4]),
    rep(data$ssp126_tasmax_05, each = days_in_months$days[5]),
    rep(data$ssp126_tasmax_06, each = days_in_months$days[6]),
    rep(data$ssp126_tasmax_07, each = days_in_months$days[7]),
    rep(data$ssp126_tasmax_08, each = days_in_months$days[8]),
    rep(data$ssp126_tasmax_09, each = days_in_months$days[9]),
    rep(data$ssp126_tasmax_10, each = days_in_months$days[10]),
    rep(data$ssp126_tasmax_11, each = days_in_months$days[11]),
    rep(data$ssp126_tasmax_12, each = days_in_months$days[12])
  )
  
  # Calculate hourly temperatures for each day
  for (j in 1:nrow(fluct)) {
    date <- as.Date(fluct$date[j], format = "%Y-%m-%d")
    T_mean <- fluct$env_tm[j]
    T_min <- fluct$env_tmin[j]
    T_max <- fluct$env_tmax[j]
    t_peak <- 14  # Time of peak temperature (2 PM)
    amplitude <- (T_max - T_mean) #(T_max - T_min) / 2
    
    # Time array (24 points from 0 to 24 hours)
    time <- seq(0, 23)
    
    # Sinusoidal temperature model
    temperature <- T_mean + amplitude * cos(2 * pi * (time - t_peak) / 24)
    
    # Create a data frame for the current day's hourly temperatures
    df_hours <- data.frame(date = date, hour = time, temperature = temperature)
    
    # Append to the main dataframe
    hourly_temperatures <- rbind(hourly_temperatures, df_hours)
  }
  
  Tb = data$Tb
  TT = data$ThetaT*24 # convert to degrees per hour
  
  # # Function to calculate the heat sum and the date reaching ThetaT for each row (OLD)
  # hourly_temperatures = hourly_temperatures %>%
  #   mutate(
  #     hourly_heat = (temperature - Tb),
  #     hourly_heat_sum = ifelse(hourly_heat < 0, 0, hourly_heat),
  #     cumulative_heat_sum = cumsum(hourly_heat_sum))
  
  # Function to calculate the heat sum and the date reaching Tb for each row including decaying term
  hourly_temperatures <- hourly_temperatures %>%
    mutate(
      heat_increment  = germ_rate(temperature, Tb, Tc, k),
      cumulative_heat_sum = cumsum(heat_increment)
    )
  
  
  # Find the day when cumulative heat sum exceeds ThetaT
  date_reach_ThetaT <- hourly_temperatures[which(hourly_temperatures$cumulative_heat_sum >= TT)[1],]$date
  
  output[i,]$date_reach_TT_hrs_ssp126 <- date_reach_ThetaT
  
  hs = hourly_temperatures
  hs$q50 = ifelse(hourly_temperatures$cumulative_heat_sum <= TT, 0, 1)
  hs$serial_number <- data$serial_number
  hs$treat <- data$treat
  hs$mod <- data$mod
  hs$Genus <- data$Genus
  hs$Species <- data$Species   
  hs$Tb <- Tb
  
  heatsum = rbind(heatsum, hs)
  
}

heatsum_d = heatsum %>%
  group_by(serial_number, treat, mod, Genus, Species, date) %>%
  summarise(
    cumulative_heat_sum =  sum(cumulative_heat_sum/24),
    q50 = min(q50),     
    Tb = unique(Tb),
    T_min = min(temperature), 
    T_max = max(temperature), 
    T_mean = mean(temperature)
  )

write.csv(heatsum_d, "output/csv/heatsum_from_Jan01_ssp126.csv", row.names = FALSE)


# ------------------------------------------------------------------------------
# Optimistic FUTURE (ssp126)

## calculate heat sum from air temperature by hour adjusted by max - START at collection date

output$date_reach_TT_hrs_datecollect_ssp126 <- NA
output$date_reach_TT_hrs_datecollect_ssp126 = as.Date(output$date_reach_TT_hrs_datecollect_ssp126)

output$cumulative_heat_sum_dec31 <- NA

output$stand_start_date <- NA
output$stand_start_date <- as.Date(output$stand_start_date)

# Initialize an empty dataframe to store the results
heatsum <- data.frame(Genus = as.character(), Species = as.character(), serial_number = as.character(), treat = as.character(), mod = as.character(), date = as.Date(character()), hour = integer(), temperature = numeric(), hourly_heat = numeric(), hourly_heat_sum = numeric(), cumulative_heat_sum =  numeric(), q50 =  numeric(), Tb = numeric())


for (i in 1:nrow(output)) {
  
  data <- output[i, ]
  
  # Dates
  start_date <- as.Date(data$Date.Collected)
  
  if (month(start_date) < 5) {
    start_date <- ymd(paste0("2023-", month(start_date), "-", day(start_date)))
  } else {
    start_date <- ymd(paste0("2022-", month(start_date), "-", day(start_date)))
  }
  
  end_date <- as.Date("2023-12-31")
  
  # Create a sequence of dates from January 1st to December 31st for a non-leap year
  datseq <- seq.Date(from = start_date, to = end_date, by = "day")
  month_year <- format(datseq, "%Y-%m")
  datseq <- format(datseq, format = "%Y-%m-%d")
  
  
  
  hourly_temperatures <- data.frame(date = as.Date(character()), hour = integer(), temperature = numeric())
  
  fluct <- data.frame(date = datseq, month_year = month_year)
  
  # Calculate the number of days in each month
  days_in_months <- fluct %>%
    group_by(month_year) %>%
    summarise(days = n()) %>%
    ungroup()
  
  
  
  # Environmental temperatures for each month
  monthly_temperatures <- list(
    data$ssp126_tas_01, data$ssp126_tas_02, data$ssp126_tas_03,
    data$ssp126_tas_04, data$ssp126_tas_05, data$ssp126_tas_06,
    data$ssp126_tas_07, data$ssp126_tas_08, data$ssp126_tas_09,
    data$ssp126_tas_10, data$ssp126_tas_11, data$ssp126_tas_12
  )
  
  monthly_temperatures_min <- list(
    data$ssp126_tasmin_01, data$ssp126_tasmin_02, data$ssp126_tasmin_03,
    data$ssp126_tasmin_04, data$ssp126_tasmin_05, data$ssp126_tasmin_06,
    data$ssp126_tasmin_07, data$ssp126_tasmin_08, data$ssp126_tasmin_09,
    data$ssp126_tasmin_10, data$ssp126_tasmin_11, data$ssp126_tasmin_12
  )
  
  monthly_temperatures_max <- list(
    data$ssp126_tasmax_01, data$ssp126_tasmax_02, data$ssp126_tasmax_03,
    data$ssp126_tasmax_04, data$ssp126_tasmax_05, data$ssp126_tasmax_06,
    data$ssp126_tasmax_07, data$ssp126_tasmax_08, data$ssp126_tasmax_09,
    data$ssp126_tasmax_10, data$ssp126_tasmax_11, data$ssp126_tasmax_12
  )
  
  start_month <- month(start_date)
  
  if (start_month < 5) {
    ordered_months <- c(start_month:12)
  } else {
    ordered_months <- c(start_month:12, 1:12)
  }
  
  
  
  env_tm <- unlist(lapply(seq_along(ordered_months), function(x) {
    rep(monthly_temperatures[[ordered_months[x]]], each = days_in_months$days[x])
  }))
  
  env_tmin <- unlist(lapply(seq_along(ordered_months), function(x) {
    rep(monthly_temperatures_min[[ordered_months[x]]], each = days_in_months$days[x])
  }))
  
  env_tmax <- unlist(lapply(seq_along(ordered_months), function(x) {
    rep(monthly_temperatures_max[[ordered_months[x]]], each = days_in_months$days[x])
  }))
  
  fluct$env_tm <- env_tm
  fluct$env_tmin <- env_tmin
  fluct$env_tmax <- env_tmax
  
  
  # Calculate hourly temperatures for each day
  for (j in 1:nrow(fluct)) {
    date <- as.Date(fluct$date[j], format = "%Y-%m-%d")
    T_mean <- fluct$env_tm[j]
    T_min <- fluct$env_tmin[j]
    T_max <- fluct$env_tmax[j]
    t_peak <- 14  # Time of peak temperature (2 PM)
    amplitude <- (T_max - T_mean) # (T_max - T_min) / 2
    
    # Time array (24 points from 0 to 24 hours)
    time <- seq(0, 23)
    
    # Sinusoidal temperature model
    temperature <- T_mean + amplitude * cos(2 * pi * (time - t_peak) / 24)
    
    # Create a data frame for the current day's hourly temperatures
    df_hours <- data.frame(date = date, hour = time, temperature = temperature)
    
    # Append to the main dataframe
    hourly_temperatures <- rbind(hourly_temperatures, df_hours)
  }
  
  Tb = data$Tb
  TT = data$ThetaT*24 # convert to degrees per hour
  
  # Function to calculate the heat sum and the date reaching ThetaT for each row (OLD)
  # hourly_temperatures = hourly_temperatures %>%
  #   mutate(
  #     hourly_heat = (temperature - Tb),
  #     hourly_heat_sum = ifelse(hourly_heat < 0, 0, hourly_heat),
  #     cumulative_heat_sum = cumsum(hourly_heat_sum))
  
  # Function to calculate the heat sum and the date reaching Tb for each row including decaying term
  hourly_temperatures <- hourly_temperatures %>%
    mutate(
      heat_increment  = germ_rate(temperature, Tb, Tc, k),
      cumulative_heat_sum = cumsum(heat_increment)
    )
  
  
  # Find the day when cumulative heat sum exceeds ThetaT
  date_reach_ThetaT <- hourly_temperatures[which(hourly_temperatures$cumulative_heat_sum >= TT)[1],]$date
  
  output[i,]$date_reach_TT_hrs_datecollect_ssp126 <- date_reach_ThetaT
  
  if (month(start_date) >= 4) {
    # Find the index of the rows for the date 2022-12-31
    index <- which(hourly_temperatures$date == as.Date("2022-12-31") & hourly_temperatures$hour == 23)
    
    # Extract the cumulative_heat_sum for the last hour of 2022-12-31
    cumulative_heat_sum_dec31 <- hourly_temperatures$cumulative_heat_sum[index]/24
  } else {
    cumulative_heat_sum_dec31 <- 0
  }
  
  
  output[i,]$cumulative_heat_sum_dec31 <- cumulative_heat_sum_dec31
  
  output[i,]$stand_start_date <- as.Date(start_date)
  
  hs = hourly_temperatures
  hs$q50 = ifelse(hourly_temperatures$cumulative_heat_sum <= TT, 0, 1)
  
  hs$serial_number <- data$serial_number
  hs$treat <- data$treat
  hs$mod <- data$mod
  hs$Genus <- data$Genus
  hs$Species <- data$Species   
  hs$Tb <- Tb
  
  heatsum = rbind(heatsum, hs)
  
}

heatsum_d = heatsum %>%
  group_by(serial_number, treat, mod, Genus, Species, date) %>%
  summarise(
    cumulative_heat_sum =  sum(cumulative_heat_sum/24),
    q50 = min(q50),     
    Tb = unique(Tb),
    T_min = min(temperature), 
    T_max = max(temperature), 
    T_mean = mean(temperature)
  )

write.csv(heatsum_d, "output/csv/heatsum_from_collect_ssp126.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Pessimistic FUTURE (ssp585)
## calculate heat sum from air temperature by hour adjusted by max - START 01 Jan

# Constants
start_date <- as.Date("2023-01-01")
end_date <- as.Date("2023-12-31")

# Create a sequence of dates from January 1st to December 31st for a non-leap year
datseq <- seq.Date(from = start_date, to = end_date, by = "day")
month_year <- format(datseq, "%Y-%m")
datseq <- format(datseq, format = "%Y-%m-%d")

# Initialize an empty dataframe to store the results
heatsum <- data.frame(Genus = as.character(), Species = as.character(), serial_number = as.character(), treat = as.character(), mod = as.character(), date = as.Date(character()), hour = integer(), temperature = numeric(), hourly_heat = numeric(), hourly_heat_sum = numeric(), cumulative_heat_sum =  numeric(), q50 =  numeric(), Tb = numeric())

output$date_reach_TT_hrs_ssp585 <- NA
output$date_reach_TT_hrs_ssp585 = as.Date(output$date_reach_TT_hrs_ssp585)

for (i in 1:nrow(output)) {
  data <- output[i, ]
  
  hourly_temperatures <- data.frame(date = as.Date(character()), hour = integer(), temperature = numeric())
  
  fluct <- data.frame(date = datseq, month_year = month_year)
  
  # Calculate the number of days in each month
  days_in_months <- fluct %>%
    group_by(month_year) %>%
    summarise(days = n()) %>%
    ungroup()
  
  # Create a vector of environmental temperatures for each day
  fluct$env_tm <- c(
    rep(data$ssp585_tas_01, each = days_in_months$days[1]),
    rep(data$ssp585_tas_02, each = days_in_months$days[2]),
    rep(data$ssp585_tas_03, each = days_in_months$days[3]),
    rep(data$ssp585_tas_04, each = days_in_months$days[4]),
    rep(data$ssp585_tas_05, each = days_in_months$days[5]),
    rep(data$ssp585_tas_06, each = days_in_months$days[6]),
    rep(data$ssp585_tas_07, each = days_in_months$days[7]),
    rep(data$ssp585_tas_08, each = days_in_months$days[8]),
    rep(data$ssp585_tas_09, each = days_in_months$days[9]),
    rep(data$ssp585_tas_10, each = days_in_months$days[10]),
    rep(data$ssp585_tas_11, each = days_in_months$days[11]),
    rep(data$ssp585_tas_12, each = days_in_months$days[12])
  )
  
  fluct$env_tmin <- c(
    rep(data$ssp585_tasmin_01, each = days_in_months$days[1]),
    rep(data$ssp585_tasmin_02, each = days_in_months$days[2]),
    rep(data$ssp585_tasmin_03, each = days_in_months$days[3]),
    rep(data$ssp585_tasmin_04, each = days_in_months$days[4]),
    rep(data$ssp585_tasmin_05, each = days_in_months$days[5]),
    rep(data$ssp585_tasmin_06, each = days_in_months$days[6]),
    rep(data$ssp585_tasmin_07, each = days_in_months$days[7]),
    rep(data$ssp585_tasmin_08, each = days_in_months$days[8]),
    rep(data$ssp585_tasmin_09, each = days_in_months$days[9]),
    rep(data$ssp585_tasmin_10, each = days_in_months$days[10]),
    rep(data$ssp585_tasmin_11, each = days_in_months$days[11]),
    rep(data$ssp585_tasmin_12, each = days_in_months$days[12])
  )
  
  fluct$env_tmax <- c(
    rep(data$ssp585_tasmax_01, each = days_in_months$days[1]),
    rep(data$ssp585_tasmax_02, each = days_in_months$days[2]),
    rep(data$ssp585_tasmax_03, each = days_in_months$days[3]),
    rep(data$ssp585_tasmax_04, each = days_in_months$days[4]),
    rep(data$ssp585_tasmax_05, each = days_in_months$days[5]),
    rep(data$ssp585_tasmax_06, each = days_in_months$days[6]),
    rep(data$ssp585_tasmax_07, each = days_in_months$days[7]),
    rep(data$ssp585_tasmax_08, each = days_in_months$days[8]),
    rep(data$ssp585_tasmax_09, each = days_in_months$days[9]),
    rep(data$ssp585_tasmax_10, each = days_in_months$days[10]),
    rep(data$ssp585_tasmax_11, each = days_in_months$days[11]),
    rep(data$ssp585_tasmax_12, each = days_in_months$days[12])
  )
  
  # Calculate hourly temperatures for each day
  for (j in 1:nrow(fluct)) {
    date <- as.Date(fluct$date[j], format = "%Y-%m-%d")
    T_mean <- fluct$env_tm[j]
    T_min <- fluct$env_tmin[j]
    T_max <- fluct$env_tmax[j]
    t_peak <- 14  # Time of peak temperature (2 PM)
    amplitude <- (T_max - T_mean) #(T_max - T_min) / 2
    
    # Time array (24 points from 0 to 24 hours)
    time <- seq(0, 23)
    
    # Sinusoidal temperature model
    temperature <- T_mean + amplitude * cos(2 * pi * (time - t_peak) / 24)
    
    # Create a data frame for the current day's hourly temperatures
    df_hours <- data.frame(date = date, hour = time, temperature = temperature)
    
    # Append to the main dataframe
    hourly_temperatures <- rbind(hourly_temperatures, df_hours)
  }
  
  Tb = data$Tb
  TT = data$ThetaT*24 # convert to degrees per hour
  
  # # Function to calculate the heat sum and the date reaching ThetaT for each row (OLD)
  # hourly_temperatures = hourly_temperatures %>%
  #   mutate(
  #     hourly_heat = (temperature - Tb),
  #     hourly_heat_sum = ifelse(hourly_heat < 0, 0, hourly_heat),
  #     cumulative_heat_sum = cumsum(hourly_heat_sum))
  
  # Function to calculate the heat sum and the date reaching Tb for each row including decaying term
  hourly_temperatures <- hourly_temperatures %>%
    mutate(
      heat_increment  = germ_rate(temperature, Tb, Tc, k),
      cumulative_heat_sum = cumsum(heat_increment)
    )
  
  
  # Find the day when cumulative heat sum exceeds ThetaT
  date_reach_ThetaT <- hourly_temperatures[which(hourly_temperatures$cumulative_heat_sum >= TT)[1],]$date
  
  output[i,]$date_reach_TT_hrs_ssp585 <- date_reach_ThetaT
  
  hs = hourly_temperatures
  hs$q50 = ifelse(hourly_temperatures$cumulative_heat_sum <= TT, 0, 1)
  hs$serial_number <- data$serial_number
  hs$treat <- data$treat
  hs$mod <- data$mod
  hs$Genus <- data$Genus
  hs$Species <- data$Species   
  hs$Tb <- Tb
  
  heatsum = rbind(heatsum, hs)
  
}

heatsum_d = heatsum %>%
  group_by(serial_number, treat, mod, Genus, Species, date) %>%
  summarise(
    cumulative_heat_sum =  sum(cumulative_heat_sum/24),
    q50 = min(q50),     
    Tb = unique(Tb),
    T_min = min(temperature), 
    T_max = max(temperature), 
    T_mean = mean(temperature)
  )

write.csv(heatsum_d, "output/csv/heatsum_from_Jan01_ssp585.csv", row.names = FALSE)



# ------------------------------------------------------------------------------
# Pessimistic FUTURE
## calculate heat sum from air temperature by hour adjusted by max - START at collection date

output$date_reach_TT_hrs_datecollect_ssp585 <- NA
output$date_reach_TT_hrs_datecollect_ssp585 = as.Date(output$date_reach_TT_hrs_datecollect_ssp585)

output$cumulative_heat_sum_dec31 <- NA

output$stand_start_date <- NA
output$stand_start_date <- as.Date(output$stand_start_date)

# Initialize an empty dataframe to store the results
heatsum <- data.frame(Genus = as.character(), Species = as.character(), serial_number = as.character(), treat = as.character(), mod = as.character(), date = as.Date(character()), hour = integer(), temperature = numeric(), hourly_heat = numeric(), hourly_heat_sum = numeric(), cumulative_heat_sum =  numeric(), q50 =  numeric(), Tb = numeric())


for (i in 1:nrow(output)) {
  
  data <- output[i, ]
  
  # Dates
  start_date <- as.Date(data$Date.Collected)
  
  if (month(start_date) < 5) {
    start_date <- ymd(paste0("2023-", month(start_date), "-", day(start_date)))
  } else {
    start_date <- ymd(paste0("2022-", month(start_date), "-", day(start_date)))
  }
  
  end_date <- as.Date("2023-12-31")
  
  # Create a sequence of dates from January 1st to December 31st for a non-leap year
  datseq <- seq.Date(from = start_date, to = end_date, by = "day")
  month_year <- format(datseq, "%Y-%m")
  datseq <- format(datseq, format = "%Y-%m-%d")
  
  
  
  hourly_temperatures <- data.frame(date = as.Date(character()), hour = integer(), temperature = numeric())
  
  fluct <- data.frame(date = datseq, month_year = month_year)
  
  # Calculate the number of days in each month
  days_in_months <- fluct %>%
    group_by(month_year) %>%
    summarise(days = n()) %>%
    ungroup()
  
  
  
  # Environmental temperatures for each month
  monthly_temperatures <- list(
    data$ssp585_tas_01, data$ssp585_tas_02, data$ssp585_tas_03,
    data$ssp585_tas_04, data$ssp585_tas_05, data$ssp585_tas_06,
    data$ssp585_tas_07, data$ssp585_tas_08, data$ssp585_tas_09,
    data$ssp585_tas_10, data$ssp585_tas_11, data$ssp585_tas_12
  )
  
  monthly_temperatures_min <- list(
    data$ssp585_tasmin_01, data$ssp585_tasmin_02, data$ssp585_tasmin_03,
    data$ssp585_tasmin_04, data$ssp585_tasmin_05, data$ssp585_tasmin_06,
    data$ssp585_tasmin_07, data$ssp585_tasmin_08, data$ssp585_tasmin_09,
    data$ssp585_tasmin_10, data$ssp585_tasmin_11, data$ssp585_tasmin_12
  )
  
  monthly_temperatures_max <- list(
    data$ssp585_tasmax_01, data$ssp585_tasmax_02, data$ssp585_tasmax_03,
    data$ssp585_tasmax_04, data$ssp585_tasmax_05, data$ssp585_tasmax_06,
    data$ssp585_tasmax_07, data$ssp585_tasmax_08, data$ssp585_tasmax_09,
    data$ssp585_tasmax_10, data$ssp585_tasmax_11, data$ssp585_tasmax_12
  )
  
  start_month <- month(start_date)
  
  if (start_month < 5) {
    ordered_months <- c(start_month:12)
  } else {
    ordered_months <- c(start_month:12, 1:12)
  }
  
  
  
  env_tm <- unlist(lapply(seq_along(ordered_months), function(x) {
    rep(monthly_temperatures[[ordered_months[x]]], each = days_in_months$days[x])
  }))
  
  env_tmin <- unlist(lapply(seq_along(ordered_months), function(x) {
    rep(monthly_temperatures_min[[ordered_months[x]]], each = days_in_months$days[x])
  }))
  
  env_tmax <- unlist(lapply(seq_along(ordered_months), function(x) {
    rep(monthly_temperatures_max[[ordered_months[x]]], each = days_in_months$days[x])
  }))
  
  fluct$env_tm <- env_tm
  fluct$env_tmin <- env_tmin
  fluct$env_tmax <- env_tmax
  
  
  # Calculate hourly temperatures for each day
  for (j in 1:nrow(fluct)) {
    date <- as.Date(fluct$date[j], format = "%Y-%m-%d")
    T_mean <- fluct$env_tm[j]
    T_min <- fluct$env_tmin[j]
    T_max <- fluct$env_tmax[j]
    t_peak <- 14  # Time of peak temperature (2 PM)
    amplitude <- (T_max - T_mean) # (T_max - T_min) / 2
    
    # Time array (24 points from 0 to 24 hours)
    time <- seq(0, 23)
    
    # Sinusoidal temperature model
    temperature <- T_mean + amplitude * cos(2 * pi * (time - t_peak) / 24)
    
    # Create a data frame for the current day's hourly temperatures
    df_hours <- data.frame(date = date, hour = time, temperature = temperature)
    
    # Append to the main dataframe
    hourly_temperatures <- rbind(hourly_temperatures, df_hours)
  }
  
  Tb = data$Tb
  TT = data$ThetaT*24 # convert to degrees per hour
  
  # # Function to calculate the heat sum and the date reaching ThetaT for each row (OLD)
  # hourly_temperatures = hourly_temperatures %>%
  #   mutate(
  #     hourly_heat = (temperature - Tb),
  #     hourly_heat_sum = ifelse(hourly_heat < 0, 0, hourly_heat),
  #     cumulative_heat_sum = cumsum(hourly_heat_sum))
  
  # Function to calculate the heat sum and the date reaching Tb for each row including decaying term
  hourly_temperatures <- hourly_temperatures %>%
    mutate(
      heat_increment  = germ_rate(temperature, Tb, Tc, k),
      cumulative_heat_sum = cumsum(heat_increment)
    )
  
  
  # Find the day when cumulative heat sum exceeds ThetaT
  date_reach_ThetaT <- hourly_temperatures[which(hourly_temperatures$cumulative_heat_sum >= TT)[1],]$date
  
  output[i,]$date_reach_TT_hrs_datecollect_ssp585 <- date_reach_ThetaT
  
  if (month(start_date) >= 4) {
    # Find the index of the rows for the date 2022-12-31
    index <- which(hourly_temperatures$date == as.Date("2022-12-31") & hourly_temperatures$hour == 23)
    
    # Extract the cumulative_heat_sum for the last hour of 2022-12-31
    cumulative_heat_sum_dec31 <- hourly_temperatures$cumulative_heat_sum[index]/24
  } else {
    cumulative_heat_sum_dec31 <- 0
  }
  
  
  output[i,]$cumulative_heat_sum_dec31 <- cumulative_heat_sum_dec31
  
  output[i,]$stand_start_date <- as.Date(start_date)
  
  hs = hourly_temperatures
  hs$q50 = ifelse(hourly_temperatures$cumulative_heat_sum <= TT, 0, 1)
  
  hs$serial_number <- data$serial_number
  hs$treat <- data$treat
  hs$mod <- data$mod
  hs$Genus <- data$Genus
  hs$Species <- data$Species   
  hs$Tb <- Tb
  
  heatsum = rbind(heatsum, hs)
  
}


heatsum_d = heatsum %>%
  group_by(serial_number, treat, mod, Genus, Species, date) %>%
  summarise(
    cumulative_heat_sum =  sum(cumulative_heat_sum/24),
    q50 = min(q50),     
    Tb = unique(Tb),
    T_min = min(temperature), 
    T_max = max(temperature), 
    T_mean = mean(temperature)
  )

write.csv(heatsum_d, "output/csv/heatsum_from_collect_ssp585.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Save output

write.csv(output, "output/csv/date_reach_TT_present_n_future.csv", row.names = FALSE)
