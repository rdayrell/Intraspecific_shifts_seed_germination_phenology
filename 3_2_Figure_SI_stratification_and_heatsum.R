library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(cowplot)
library(patchwork)


# import environmental data ----------------------------------------------------
env1 = read.csv("data/csv/environment_data.csv")
env2 = read.csv("data/csv/environment_data_temp_n_water.csv")
env2 <- env2 %>%
  dplyr::select(-dplyr::starts_with("CHELSA_tas")) %>%
  rename(serial_number = Serial.No)  

env = env1 %>%
  left_join(., env2)

rm(env1, env2)

env <- env %>%
  rename(fcf_1981.2010_CHELSA = CHELSA_fcf_1981.2010) 


env_long <- env %>%
  pivot_longer(
    cols = starts_with("CHELSA") | starts_with("ssp126_") | starts_with("ssp585_") | starts_with("soil_temp"),
    names_to = c("variable", "month", "year"),
    names_pattern = "(CHELSA_[^_]+|ssp126_[^_]+|ssp585_[^_]+|soil_temp_0_5_cm_month)_(\\d{2})(?:_(.*))?",
    values_to = "value"
  ) %>%
  mutate(
    year = ifelse(year == "" & grepl("soil", variable), "1981.2010", 
                  ifelse(year == "", "2071.2100", year)), ## this is not the correct year range for soil
    month = as.numeric(month)
  )


env_long = env_long %>%
  dplyr::select(-year)


# Separate the environmental variable and the time component
env_wide <- env_long %>%
  pivot_wider(
    names_from = "variable",
    values_from = "value"
  )

# import heatsum results -------------------------------------------------------

output0 <- read.csv("output/csv/heatsum_from_collect_present.csv")
output0$start <- "Dispersal"
output0$type <- "1981 - 2010" # "max

output1 <- read.csv("output/csv/heatsum_from_Jan01_present.csv")
output1$start <- "Jan"
output1$type <- "1981 - 2010" #"max"

output2 <- read.csv("output/csv/heatsum_from_collect_ssp126.csv")
output2$start <- "Dispersal"
output2$type <- "2071 - 2100 (optimistic)"

output3 <- read.csv("output/csv/heatsum_from_Jan01_ssp126.csv")
output3$start <- "Jan"
output3$type <- "2071 - 2100 (optimistic)"

output4 <- read.csv("output/csv/heatsum_from_collect_ssp585.csv")
output4$start <- "Dispersal"
output4$type <- "2071 - 2100 (pessimistic)"

output5 <- read.csv("output/csv/heatsum_from_Jan01_ssp585.csv")
output5$start <- "Jan"
output5$type <- "2071 - 2100 (pessimistic)"

df = rbind(output0, output1, output2, output3, output4, output5)

df$date <- as.Date(df$date)
df$month <- month(df$date)

rm(list = paste0("output", 0:5))

dat = merge(df, env_wide, by=c("serial_number","Genus","Species", "month"))
dat$cumulative_heat_sum = ifelse(dat$q50 == 1, NA, dat$cumulative_heat_sum)

# ------------------------------------------------------------------------------
# PLOT MONTHLY TEMPERATURE AND COLD STRAT (PURPLE LINES) 
## supplementary figure

## Alnus glutinosa -------------------------------------------------------------
plot_dat = dat %>%
  filter(Genus == "Alnus" & Species == "glutinosa" & mod == "modGR_Ex" & start == "Jan") %>%
  arrange(lat_dd)

plot_dat$serial_number <- factor(plot_dat$serial_number, levels = unique(plot_dat$serial_number))

plot_dat$Geographical.Location <- sub("^UK:", "", plot_dat$Geographical.Location)
plot_dat$Geographical.Location <- sub(":", "\n", plot_dat$Geographical.Location)

plot_dat$type <- sub(" (", "\n(", plot_dat$type, fixed = T)

plot_dat <- plot_dat %>%
  mutate(tasmin = case_when(
    type == "1981 - 2010" ~ CHELSA_tasmin,
    type == "2071 - 2100\n(optimistic)" ~ ssp126_tasmin,
    type == "2071 - 2100\n(pessimistic)" ~ ssp585_tasmin
  ))

plot_dat %>%
  group_by(serial_number, month, type) %>%
  summarise(
    tasmin = unique(tasmin)  
  )

## ordering by latitude (Alnus and Betula)

loc_labels <- plot_dat %>%
  group_by(serial_number) %>%
  summarise(Geographical.Location = unique(Geographical.Location)) %>%
  distinct(serial_number, Geographical.Location) %>%
  deframe()

plot_temp_alnus = ggplot(plot_dat, aes(x = date)) +
   facet_grid(type ~ serial_number, scales = "fixed", labeller = labeller(serial_number = loc_labels)) +
   geom_line(aes(y = T_mean, color = "T_mean"), linewidth = 0.5, colour = "grey15") +
   geom_ribbon(aes(ymax = T_max, ymin = T_min), fill = "grey15", alpha = 0.3) +
   geom_hline(aes(yintercept = 5), linewidth = 1, alpha = 1, colour = "navy") + 
   geom_hline(aes(yintercept = 7), linewidth = 0.7, linetype = "dashed", alpha = 1, colour = "navy") +
   ggtitle(paste(unique(plot_dat$Genus), unique(plot_dat$Species), sep = " ")) +
   
   # Scale for the primary y-axis (cumulative_heat_sum)
   scale_y_continuous(name = "Temperature (°C)") +
   scale_x_date(name = "", date_labels = "%b", date_breaks = "1 month") +
   theme_minimal(base_size = 10) +
   theme(
     plot.title = element_text(size = 12, face = "italic"), #vjust = -8
     panel.grid = element_blank(),
     title = element_text(size = 8),
     panel.grid.minor = element_blank(),
     axis.title.y = element_text(size = 10),
     plot.background = element_rect(fill = "white", colour = NA)
   )

plot_temp_alnus = plot_temp_alnus + scale_x_date(name = "", date_breaks = "1 month", labels = function(x) substr(month.abb[as.POSIXlt(x)$mon + 1], 1, 1))

# file_name <- sprintf("output/imgs/Fig_SI_monthly_temperature_%s_%s.png", unique(plot_dat$Genus), unique(plot_dat$Species))
# ggsave(file_name, plot_temp_alnus, width = 10, height = 7, units = "in")


## Betula pubescens ------------------------------------------------------------

plot_dat = dat %>%
  filter(Genus == "Betula" & Species == "pubescens" & mod == "modGR_Ex" & start == "Jan") %>%
  arrange(lat_dd)

plot_dat$serial_number <- factor(plot_dat$serial_number, levels = unique(plot_dat$serial_number))

plot_dat$Geographical.Location <- sub("^UK:", "", plot_dat$Geographical.Location)
plot_dat$Geographical.Location <- sub(":", "\n", plot_dat$Geographical.Location)

plot_dat$type <- sub(" (", "\n(", plot_dat$type, fixed = T)

plot_dat <- plot_dat %>%
  mutate(tasmin = case_when(
    type == "1981 - 2010" ~ CHELSA_tasmin,
    type == "2071 - 2100\n(optimistic)" ~ ssp126_tasmin,
    type == "2071 - 2100\n(pessimistic)" ~ ssp585_tasmin
  ))

plot_dat %>%
  group_by(serial_number, month, type) %>%
  summarise(
    tasmin = unique(tasmin)  
  )

## ordering by latitude (Alnus and Betula)

loc_labels <- plot_dat %>%
  group_by(serial_number) %>%
  summarise(Geographical.Location = unique(Geographical.Location)) %>%
  distinct(serial_number, Geographical.Location) %>%
  deframe()

plot_temp_betula = ggplot(plot_dat, aes(x = date)) +
    facet_grid(type ~ serial_number, scales = "fixed", labeller = labeller(serial_number = loc_labels)) +
    geom_line(aes(y = T_mean, color = "T_mean"), linewidth = 0.5, colour = "grey15") +
    geom_ribbon(aes(ymax = T_max, ymin = T_min), fill = "grey15", alpha = 0.3) +
    geom_hline(aes(yintercept = 5), linewidth = 1, alpha = 1, colour = "navy") + 
    geom_hline(aes(yintercept = 7), linewidth = 0.7, linetype = "dashed", alpha = 1, colour = "navy") +
    ggtitle(paste(unique(plot_dat$Genus), unique(plot_dat$Species), sep = " ")) +
    
    # Scale for the primary y-axis (cumulative_heat_sum)
    scale_y_continuous(name = "Temperature (°C)") +
    scale_x_date(name = "", date_labels = "%b", date_breaks = "1 month") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 12, face = "italic"), #vjust = -8
      panel.grid = element_blank(),
      title = element_text(size = 8),
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(size = 10),
      plot.background = element_rect(fill = "white", colour = NA)
    )

plot_temp_betula = plot_temp_betula + scale_x_date(name = "", date_breaks = "1 month", labels = function(x) substr(month.abb[as.POSIXlt(x)$mon + 1], 1, 1))

# file_name <- sprintf("output/imgs/Fig_SI_monthly_temperature_%s_%s.png", unique(plot_dat$Genus), unique(plot_dat$Species))
# ggsave(file_name, plot_temp_betula, width = 10, height = 7, units = "in")


# Pinus sylvestris -------------------------------------------------------------

plot_dat = dat %>%
  filter(Genus == "Pinus" & Species == "sylvestris" & mod == "modGR_Ex" & start == "Jan") %>%
  arrange(Altitude)

#unique(plot_dat$Geographical.Location)
plot_dat$serial_number <- factor(plot_dat$serial_number, levels = unique(plot_dat$serial_number))

plot_dat$Geographical.Location <- sub("^UK:", "", plot_dat$Geographical.Location)
plot_dat$Geographical.Location <- sub(":", "\n", plot_dat$Geographical.Location)

plot_dat$type <- sub(" (", "\n(", plot_dat$type, fixed = T)

plot_dat <- plot_dat %>%
  mutate(tasmin = case_when(
    type == "1981 - 2010" ~ CHELSA_tasmin,
    type == "2071 - 2100\n(optimistic)" ~ ssp126_tasmin,
    type == "2071 - 2100\n(pessimistic)" ~ ssp585_tasmin
  ))

plot_dat %>%
  group_by(serial_number, month, type) %>%
  summarise(
  tasmin = unique(tasmin)  
  )

# Order by altitude

loc_labels <- plot_dat %>%
  group_by(serial_number) %>%
  summarise(Geographical.Location = paste(unique(Geographical.Location), "\n", as.character(unique(Altitude)),"m")) %>% #, unique(serial_number)
  distinct(serial_number, Geographical.Location) %>%
  deframe()

plot_temp_pinus = ggplot(plot_dat, aes(x = date)) +
    facet_grid(type ~ serial_number, scales = "fixed", labeller = labeller(serial_number = loc_labels)) +
    geom_line(aes(y = T_mean, color = "T_mean"), linewidth = 0.5, colour = "grey15") +
    geom_ribbon(aes(ymax = T_max, ymin = T_min), fill = "grey15", alpha = 0.3) +
    geom_hline(aes(yintercept = 5), linewidth = 1, alpha = 1, colour = "navy") + 
    geom_hline(aes(yintercept = 7), linewidth = 0.7, linetype = "dashed", alpha = 1, colour = "navy") +
    ggtitle(paste(unique(plot_dat$Genus), unique(plot_dat$Species), sep = " ")) +
    
    # Scale for the primary y-axis (cumulative_heat_sum)
    scale_y_continuous(name = "Temperature (°C)") +
    scale_x_date(name = "", date_labels = "%b", date_breaks = "1 month") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 12, face = "italic"), #vjust = -8
      panel.grid = element_blank(),
      title = element_text(size = 8),
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(size = 10),
      plot.background = element_rect(fill = "white", colour = NA)
    )

plot_temp_pinus = plot_temp_pinus + scale_x_date(name = "", date_breaks = "1 month", labels = function(x) substr(month.abb[as.POSIXlt(x)$mon + 1], 1, 1))

# file_name <- sprintf("output/imgs/Fig_SI_monthly_temperature_%s_%s.png", unique(plot_dat$Genus), unique(plot_dat$Species))
# ggsave(file_name, plot_temp_pinus, width = 10, height = 7, units = "in")

plot_temps = plot_temp_alnus / plot_temp_betula / plot_temp_pinus
ggsave("output/imgs/Fig_SI_monthly_temperature_stratification.png", plot_temps, width = 8, height = 12, units = "in")


# ------------------------------------------------------------------------------
# Prepare data to plot heat sum for each population at each location
## Supplementary figure

# ------------------------------------------------------------------------------
# COLD STRAT risk

names(dat)

dat$strat_temp_cutoff = 7 # 5C (+ 2C of incubator upper limit variation) 

dat1 <- dat %>%
  pivot_longer(cols = starts_with("CHELSA_tas") | starts_with("ssp126_tas") | starts_with("ssp585_tas"), 
               names_to = "scenario", 
               values_to = "tas") 

#checks whether dormancy is broken
df_filtered <- dat1 %>%
  filter(mod == "modGR_Ex"  & start == "Jan" & type == "1981 - 2010") %>%
  mutate(date = as.Date(date),
         month_year = format(date, "%Y-%m")) %>%
  group_by(serial_number, Genus, Species, month_year, lat_dd, Altitude, type, Geographical.Location, strat_temp_cutoff, scenario) %>%
  summarise(
  tas = unique(tas)) %>%
  filter(tas <= strat_temp_cutoff & (month_year == "2023-01" | month_year == "2023-02" | month_year == "2023-12") & (scenario == "ssp126_tas" | scenario == "ssp585_tas")) %>% #scenario == "CHELSA_tas" | 
  ungroup() %>%
group_by(serial_number, Genus, Species, scenario) %>%
  summarise(
    n = length(serial_number)
  )

df_filtered$serial_number <- as.factor(df_filtered$serial_number)

df_filtered$treat_time = ifelse(df_filtered$scenario == "ssp126_tas", "COLD - 2071 - 2100 (optimistic)" , NA)
df_filtered$treat_time = ifelse(df_filtered$scenario == "ssp585_tas", "COLD - 2071 - 2100 (pessimistic)" , df_filtered$treat_time)


# ------------------------------------------------------------------------------
# Prepare dataset for plot synthesis figures in script 3_3_graph_heatsum_synthesis

start_NS = "Dispersal"
start_COLD = "Jan"
mod1 = "modGR_Ex" #"modGR_BS"

plot_data <- dat %>%
  separate(Geographical.Location, into = c("Country", "Region", "Area"), sep = ":") %>%
    mutate(
      location = paste(Region, Area, sep = "\n"),
      location2 = paste(Region, Area, sep = " - "),
      serial_number = as.factor(serial_number),
      treat_time = paste(treat, type, sep = " - ")
      ) %>%
  
  # Arrange by lat_dd for Genus == "Betula" and Genus == "Alnus"
  group_by(Genus) %>%
  arrange(lat_dd) %>%
  mutate(
    location = if_else(Genus %in% c("Betula", "Alnus"), fct_reorder(location, lat_dd), location),
    serial_number = if_else(Genus %in% c("Betula", "Alnus"), fct_reorder(serial_number, lat_dd), serial_number)
  ) %>%
  ungroup() %>%
  
  # Arrange by Altitude for Genus == "Pinus"
  group_by(Genus) %>%
  arrange(Altitude) %>%
  mutate(
    location = if_else(Genus == "Pinus", fct_reorder(location, Altitude), location),
    serial_number = if_else(Genus == "Pinus", fct_reorder(serial_number, Altitude), serial_number)
  ) %>%
  ungroup() %>%
  
  # Apply the filter conditions
  filter(
    (treat == "NS" & start == start_NS & mod == mod1) |
      (treat == "COLD" & start == start_COLD & mod == mod1)
    )

plot_data$type
plot_data <- plot_data %>%
  full_join(df_filtered, by = c("serial_number", "Genus", "Species", "treat_time")) %>%
  mutate(
    show_in_plot = case_when(
      treat == "NS" ~ 1,
      type == "1981 - 2010" ~ 1,
      # this line is removing COLD in the population that NS germinates first
      # ONLY FOR GRAPH HEATSUM SYNTHESIS
      serial_number == "1010271" & treat == "COLD" & scenario == "ssp585_tas" ~ 0,
      !is.na(n) ~ 1,
      TRUE ~ 0
    )
  )

x = plot_data %>%
  group_by(serial_number, Genus, treat, type, lat_dd, Altitude, Region, Area) %>%
  summarise(
    show_in_plot = mean(show_in_plot)
  )

plot_data <- plot_data %>%
  dplyr::filter(show_in_plot == 1)


plot_data <- plot_data %>%
  mutate(serial_number = droplevels(serial_number))

total_species <- length(unique(plot_data$serial_number))
species_list <- unique(plot_data$serial_number)


y2 = ifelse(max(plot_data$T_max, na.rm = T) > max(plot_data$cumulative_heat_sum, na.rm = T)/80,
            max(plot_data$T_max, na.rm = T),
            max(plot_data$cumulative_heat_sum, na.rm = T)/80)
y2 = ifelse(max(plot_data$Tb, na.rm = T) > y2,
            max(plot_data$Tb, na.rm = T),
            y2)
x2 = min(plot_data$date)
x1 = max(plot_data$date)

unique(plot_data$treat_time)
plot_data$treat_time <- as.factor(plot_data$treat_time)
plot_data$treat_time <- factor(plot_data$treat_time, 
                               levels = c("COLD - 1981 - 2010",
                                          "COLD - 2071 - 2100 (optimistic)",
                                          "COLD - 2071 - 2100 (pessimistic)",
                                          "NS - 1981 - 2010",
                                          "NS - 2071 - 2100 (optimistic)",
                                          "NS - 2071 - 2100 (pessimistic)"))

save(plot_data,file="output/rda/plot_heatsum_data.Rda")
## Plot synthesis figures in script 3_3_graph_heatsum_synthesis

# ------------------------------------------------------------------------------
# Prepare dataset for individual figures in script 3_4_Figure_SI_heatsum

plot_data <- dat %>%
  separate(Geographical.Location, into = c("Country", "Region", "Area"), sep = ":") %>%
  mutate(
    location = paste(Region, Area, sep = "\n"),
    location2 = paste(Region, Area, sep = " - "),
    serial_number = as.factor(serial_number),
    treat_time = paste(treat, type, sep = " - ")
  ) %>%
  
  # Arrange by lat_dd for Genus == "Betula" and Genus == "Alnus"
  group_by(Genus) %>%
  arrange(lat_dd) %>%
  mutate(
    location = if_else(Genus %in% c("Betula", "Alnus"), fct_reorder(location, lat_dd), location),
    serial_number = if_else(Genus %in% c("Betula", "Alnus"), fct_reorder(serial_number, lat_dd), serial_number)
  ) %>%
  ungroup() %>%
  
  # Arrange by Altitude for Genus == "Pinus"
  group_by(Genus) %>%
  arrange(Altitude) %>%
  mutate(
    location = if_else(Genus == "Pinus", fct_reorder(location, Altitude), location),
    serial_number = if_else(Genus == "Pinus", fct_reorder(serial_number, Altitude), serial_number)
  ) %>%
  ungroup() %>%
  
  # Apply the filter conditions
  filter(
    (treat == "NS" & start == start_NS & mod == mod1) |
      (treat == "COLD" & start == start_COLD & mod == mod1)
  )

plot_data$type
plot_data <- plot_data %>%
  full_join(df_filtered, by = c("serial_number", "Genus", "Species", "treat_time")) %>%
  mutate(
    show_in_plot = case_when(
      treat == "NS" ~ 1,
      type == "1981 - 2010" ~ 1,
      # this line is removing COLD in the population that NS germinates first
      # ONLY FOR GRAPH HEATSUM SYNTHESIS
      # serial_number == "1010271" & treat == "COLD" & scenario == "ssp585_tas" ~ 0,
      !is.na(n) ~ 1,
      TRUE ~ 0
    )
  )

x = plot_data %>%
  group_by(serial_number, Genus, treat, type, lat_dd, Altitude, Region, Area) %>%
  summarise(
    show_in_plot = mean(show_in_plot)
  )

plot_data <- plot_data %>%
  dplyr::filter(show_in_plot == 1)


plot_data <- plot_data %>%
  mutate(serial_number = droplevels(serial_number))

total_species <- length(unique(plot_data$serial_number))
species_list <- unique(plot_data$serial_number)


y2 = ifelse(max(plot_data$T_max, na.rm = T) > max(plot_data$cumulative_heat_sum, na.rm = T)/80,
            max(plot_data$T_max, na.rm = T),
            max(plot_data$cumulative_heat_sum, na.rm = T)/80)
y2 = ifelse(max(plot_data$Tb, na.rm = T) > y2,
            max(plot_data$Tb, na.rm = T),
            y2)
x2 = min(plot_data$date)
x1 = max(plot_data$date)

unique(plot_data$treat_time)
plot_data$treat_time <- as.factor(plot_data$treat_time)
plot_data$treat_time <- factor(plot_data$treat_time, 
                               levels = c("COLD - 1981 - 2010",
                                          "COLD - 2071 - 2100 (optimistic)",
                                          "COLD - 2071 - 2100 (pessimistic)",
                                          "NS - 1981 - 2010",
                                          "NS - 2071 - 2100 (optimistic)",
                                          "NS - 2071 - 2100 (pessimistic)"))

save(plot_data,file="output/rda/plot_heatsum_individual_data.Rda")
## Plot individual figures in script 3_4_Figure_SI_heatsum
