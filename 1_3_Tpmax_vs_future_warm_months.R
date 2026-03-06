library(tidyr)
library(dplyr)
library(ggplot2)

# Import data ------------------------------------------------------------------
gp <- read.csv("output/csv/finalgermination_6weeks.csv")
env <- read.csv("data/csv/environment_data.csv")
bioclim <- read.csv("data/csv/bioclim_data.csv")

# Prepare datasets -------------------------------------------------------------

bioclim = bioclim %>%
  rename(serial_number = Serial.No)

env = inner_join(env, bioclim)
# glimpse(env)

env = env %>%
  mutate(date = as.Date(Date.Collected, "%d-%b-%y"),
         month = format(date, "%m"))


# Calculate predicted average environmental temperature ------------------------
# of the warmest three months in the future pessimistic scenario ---------------
# ------------------------------------------------------------------------------

# Create columns to store mean values
env$pess_tas_month    <- NA_real_
env$pess_tas_3months  <- NA_real_
env$pess_tasmax_month <- NA_real_
env$pess_tasmax_3months <- NA_real_

# Get column names for tas and tasmax
tas_cols    <- grep("^ssp585_tas_", names(env), value = TRUE)
tasmax_cols <- grep("^ssp585_tasmax_", names(env), value = TRUE)

# Loop through each row
for (i in seq_len(nrow(env))) {
  
  # Extract values for the current row
  tas_vals    <- as.numeric(env[i, tas_cols])
  tasmax_vals <- as.numeric(env[i, tasmax_cols])
  
  # Get top 3 values
  top3_tas    <- sort(tas_vals, decreasing = TRUE)[1:3]
  top3_tasmax <- sort(tasmax_vals, decreasing = TRUE)[1:3]
  
  # Assign values to new columns
  env$pess_tas_month[i]    <- max(top3_tas, na.rm = TRUE)
  env$pess_tas_3months[i]  <- mean(top3_tas, na.rm = TRUE)
  
  # env$pess_tasmax_month[i]    <- max(top3_tasmax, na.rm = TRUE)
  # env$pess_tasmax_3months[i]  <- mean(top3_tasmax, na.rm = TRUE)
}

gp = gp %>%
  dplyr::select(-lat_dd)

df = gp %>%
  left_join(., env, by = c("serial_number", "sp"))


# Plot figure ------------------------------------------------------------------

plot1 = ggplot(df) +
  
  facet_grid(. ~ group) +
  
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = pess_tas_3months, y = Topt_upper, colour = sp),
             size = 2, alpha = 0.7) +
  scale_color_viridis_d(name = "Species", option = "G", end = 0.8) +
  scale_x_continuous(name = "Predicted average environmental temperature of the\nwarmest three months in the future pessimistic scenario", limits = c(18,34)) +
  scale_y_continuous(name = expression(italic(T[pmax_upper])), limits = c(18,34)) +
  theme_bw(base_size = 11) +
  theme(legend.position   = "right",
        legend.text       = element_text(size = 9, face = "italic"),
        # strip.text.x      = element_text(size = 9, face = "italic"),
        strip.background.x= element_rect(fill = "white"),
        plot.background   = element_rect(fill = "white"))

ggsave("output/imgs/Fig_SI_pess_tas_3months_VS_Topt_upper.png", plot1, width = 5.7, height = 2.5, units = "in")

