library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(cowplot)
library(patchwork)

load("output/rda/plot_heatsum_data.Rda")

data_plot <- plot_data %>%
  dplyr::filter(cumulative_heat_sum != 0 & !is.na(cumulative_heat_sum))

d1 <- data_plot %>%
  group_by(serial_number, Genus, Species, treat, mod, type, long_dd, lat_dd, Altitude, Date.Collected, Region, Area, location, location2, treat_time) %>%
  dplyr::filter(cumulative_heat_sum == min(cumulative_heat_sum) |
                cumulative_heat_sum == max(cumulative_heat_sum))

dat_diff <- d1 %>%
  group_by(serial_number, Genus, Species, treat, mod, type, long_dd, lat_dd, Altitude, Date.Collected, Region, Area, location, location2, treat_time) %>%
  summarise(
    min_date = date[which.min(cumulative_heat_sum)],
    max_date = date[which.max(cumulative_heat_sum)],
    ndays_0_to_50 = as.numeric(max_date - min_date) # difference in days
  ) %>%
  ungroup()

dat_diff$type = as.factor(dat_diff$type)


# Separate data into different scenarios ---------------------------------------

dat_optimistic <- dat_diff %>%
  filter(type == "2071 - 2100 (optimistic)")

dat_optimistic <- dat_optimistic %>%
  group_by(serial_number) %>%
  filter(if(n() > 1) treat == "COLD" else TRUE) %>%
  ungroup() 

dat_pessimistic <- dat_diff %>%
  filter(type == "2071 - 2100 (pessimistic)")

dat_pessimistic <- dat_pessimistic %>%
  group_by(serial_number) %>%
  filter(if(n() > 1) treat == "COLD" else TRUE) %>%
  ungroup() 

dat_1981_2010 <- dat_diff %>%
  filter(type == "1981 - 2010" & treat == "COLD")



# Join optimistic with 1981 - 2010 and calculate date differences for max_date and min_date
diff_optimistic <- dat_optimistic %>%
  left_join(dat_1981_2010, by = c("serial_number", "Genus", "Species", "mod", "long_dd", "lat_dd", "Altitude", "Date.Collected", "Region", "Area", "location", "location2")) %>%
  mutate(
    max_date_diff = as.numeric(max_date.x - max_date.y),
    min_date_diff = as.numeric(min_date.x - min_date.y),
    duration = ndays_0_to_50.x - ndays_0_to_50.y
  ) 

# Join pessimistic with 1981 - 2010 and calculate date differences for max_date and min_date
diff_pessimistic <- dat_pessimistic %>%
  left_join(dat_1981_2010, by = c("serial_number", "Genus", "Species", "mod", "long_dd", "lat_dd", "Altitude", "Date.Collected", "Region", "Area", "location", "location2")) %>%
  mutate(
    max_date_diff = as.numeric(max_date.x - max_date.y),
    min_date_diff = as.numeric(min_date.x - min_date.y),
    duration = ndays_0_to_50.x - ndays_0_to_50.y
  ) %>%  
 arrange(lat_dd) %>%
mutate(location2 = fct_reorder(location2, lat_dd))

diff <- rbind(diff_optimistic, diff_pessimistic)

diff <- diff %>%
  mutate(
    location2 = if_else(Genus == "Pinus", paste0(Area, " - ", Altitude, "m"), as.character(location2)),
    sp = paste(Genus, Species)
  )

diff$location2
diff <- diff %>%  
  arrange(lat_dd) %>%
  mutate(location2 = fct_reorder(location2, lat_dd)) 


# Alnus ------------------------------------------------------------------------

diff = droplevels(diff)
(plot_diff1 <- diff %>%
    filter(Genus == "Alnus") %>%
    ggplot(aes(x = location2)) +
  facet_wrap(sp~., scales = "free") + #or "free_y"
  geom_point(aes(y=max_date_diff, colour = type.x), size = 3) +
  scale_color_manual(name = "", values = c("2071 - 2100 (optimistic)" = "#D55F00", "2071 - 2100 (pessimistic)" = "darkred")) +
  geom_hline(yintercept = 0, linetype = "dotted")+
  labs(x = "",
       y = expression(Delta~"day-of-year to reach 50% germination (days)")) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 10),
        #axis.text.x = element_text(angle = 90), # hjust = -1, 
        strip.background = element_blank(),
        strip.text = element_text(face ="italic"),
        legend.position = "none",
        plot.background = element_rect(fill = "white")) +
  coord_flip())

# file_name <- sprintf("output/plot_synthesis_diff_days_latitude.png")
# ggsave(file_name, plot_diff1, width = 8, height = 3, units = "in")


(plot_spread1 <- diff %>%
    filter(Genus == "Alnus") %>%
    ggplot(aes(x = location2)) +
    facet_wrap(sp~., scales = "free") + #or "free_y"
    geom_point(aes(y=duration, colour = type.x), size = 3) +
    scale_color_manual(name = "", values = c("2071 - 2100 (optimistic)" = "#D55F00", "2071 - 2100 (pessimistic)" = "darkred")) +
    geom_hline(yintercept = 0, linetype = "dotted")+
    labs(x = "",
         y = expression(Delta~"days from first germination to 50% germination")) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 10),
          #axis.text.x = element_text(angle = 90), # hjust = -1,
          strip.background = element_blank(),
          strip.text = element_text(face ="italic"),
          legend.position = "none",
          plot.background = element_rect(fill = "white")) +
    coord_flip())

# file_name <- sprintf("output/plot_synthesis_spread_days_latitude.png")
# ggsave(file_name, plot_spread1, width = 8, height = 3, units = "in")



# Betula -----------------------------------------------------------------------

(plot_diff2 <- diff %>%
    filter(Genus == "Betula") %>%
    ggplot(aes(x = location2)) +
    facet_wrap(sp~., scales = "free") + #or "free_y"
    geom_point(aes(y=max_date_diff, colour = type.x), size = 3) +
    scale_color_manual(name = "", values = c("2071 - 2100 (optimistic)" = "#D55F00", "2071 - 2100 (pessimistic)" = "darkred")) +
    geom_hline(yintercept = 0, linetype = "dotted")+
    labs(x = "",
         y = expression(Delta~"day-of-year to reach 50% germination (days)")) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 10),
          #axis.text.x = element_text(angle = 90), # hjust = -1, 
          strip.background = element_blank(),
          strip.text = element_text(face ="italic"),
          legend.position = "none",
          plot.background = element_rect(fill = "white")) +
    coord_flip())

# file_name <- sprintf("output/plot_synthesis_diff_days_latitude.png")
# ggsave(file_name, plot_diff1, width = 8, height = 3, units = "in")


(plot_spread2 <- diff %>%
    filter(Genus == "Betula") %>%
    ggplot(aes(x = location2)) +
    facet_wrap(sp~., scales = "free") + #or "free_y"
    geom_point(aes(y=duration, colour = type.x), size = 3) +
    scale_color_manual(name = "", values = c("2071 - 2100 (optimistic)" = "#D55F00", "2071 - 2100 (pessimistic)" = "darkred")) +
    geom_hline(yintercept = 0, linetype = "dotted")+
    labs(x = "",
         y = expression(Delta~"days from first germination to 50% germination")) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 10),
          #axis.text.x = element_text(angle = 90), # hjust = -1,
          strip.background = element_blank(),
          strip.text = element_text(face ="italic"),
          legend.position = "none",
          plot.background = element_rect(fill = "white")) +
    coord_flip())

# file_name <- sprintf("output/plot_synthesis_spread_days_latitude.png")
# ggsave(file_name, plot_spread1, width = 8, height = 3, units = "in")


#### PINUS

diff$location2
diff <- diff %>%
  arrange(Altitude) %>%
  mutate(
  location2 = if_else(Genus == "Pinus", fct_reorder(factor(location2), Altitude), fct_reorder(factor(location2), lat_dd))
) 

levels(as.factor(diff$type.x))
diff = droplevels(diff)
(plot_diff3 <- diff %>%
    filter(Genus == "Pinus") %>%
    ggplot(aes(x = location2)) +
    facet_wrap(sp~., scales = "free") + #or "free_y"
    geom_point(aes(y=max_date_diff, colour = type.x), size = 3) +
    scale_color_manual(name = "", values = c("2071 - 2100 (optimistic)" = "#D55F00", "2071 - 2100 (pessimistic)" = "darkred")) +
    geom_hline(yintercept = 0, linetype = "dotted")+
    labs(x = "",
         y = expression(Delta~"day-of-year to reach 50% germination (days)")) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 10),
          #axis.text.x = element_text(angle = 90), # hjust = -1, 
          strip.background = element_blank(),
          strip.text = element_text(face ="italic"),
          legend.position = "right",
          plot.background = element_rect(fill = "white")) +
    coord_flip())

# file_name <- sprintf("output/plot_synthesis_diff_days_altitude.png")
# ggsave(file_name, plot_diff2, width = 5.4, height = 3, units = "in")



(plot_spread3 <- diff %>%
    filter(Genus == "Pinus") %>%
    ggplot(aes(x = location2)) +
    facet_wrap(sp~., scales = "free") + #or "free_y"
    geom_point(aes(y=duration, colour = type.x), size = 3) +
    scale_color_manual(name = "", values = c("2071 - 2100 (optimistic)" = "#D55F00", "2071 - 2100 (pessimistic)" = "darkred")) +
    geom_hline(yintercept = 0, linetype = "dotted")+
    labs(x = "",
         y = expression(Delta~"days from first germination to 50% germination")) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 10),
          #axis.text.x = element_text(angle = 90), # hjust = -1, 
          strip.background = element_blank(),
          strip.text = element_text(face ="italic"),
          legend.position = "right",
          plot.background = element_rect(fill = "white")) +
    coord_flip())


# file_name <- sprintf("output/plot_synthesis_spread_days_altitude.png")
# ggsave(file_name, plot_spread2, width = 5.4, height = 3, units = "in")

# write.csv(diff, "output/days_difference_germination.csv")

# (combined_plot <- plot_diff1 + plot_diff2 +
#   plot_layout(axis_titles = "collect_x", guides = "collect", heights = 1, widths = c(3, 1.2)))
# 
# Remove guides = "collect" from the subplots
combined_diff <- (plot_diff1 + plot_diff2 + plot_diff3) +
  plot_layout(axis_titles = "collect_x") &
  theme(legend.position = "top")

combined_spread <- (plot_spread1 + plot_spread2 + plot_spread3) +
  plot_layout(axis_titles = "collect_x") &
  theme(legend.position = "top")

# Apply guides = "collect" only at the final combination level
combined_plot <- (combined_diff / combined_spread) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

file_name <- sprintf("output/imgs/Fig_4_heatsum_synthesis.png")
ggsave(file_name, combined_plot, width = 9, height = 5, units = "in")

# ggsave("output/imgs/production/Fig_4_heatsum_synthesis.pdf", combined_plot, width = 9, height = 5, units = "in")
# 
# ggsave("output/imgs/production/Fig_4_heatsum_synthesis.tif", combined_plot, width = 9, height = 5, units = "in", dpi = 600, compression = "lzw")
