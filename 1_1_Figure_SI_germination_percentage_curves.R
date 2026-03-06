library(ggplot2)
library(tidyverse)
library(abind)
library(patchwork)
library(dplyr)
library(ggh4x)

options(scipen=999)


# Import data ------------------------------------------------------------------

df1 <- read.csv("data/csv/germination_timetoevent_Alnus_glutinosa.csv")
df1$species = "Alnus glunitosa"

df2 <- read.csv("data/csv/germination_timetoevent_Betula_pubescens.csv")
df2$species = "Betula pubescens"

df3 <- read.csv("data/csv/germination_timetoevent_Pinus_sylvestris.csv")
df3$species = "Pinus sylvestris"

df = rbind(df1, df2, df3)
rm(df1, df2, df3)


df$treat <- as.factor(df$treat)
df$treat <- factor(df$treat, levels = c("NS", "COLD"))

plot = df %>%
  group_by(species, serial_number, treat, temp, timeAf, group) %>%
  summarise(propCum = mean(propCum), .groups = "drop") %>%
  ggplot() +
  facet_nested(rows = vars(species, serial_number), cols = vars(treat)) +  # <- use facet_nested
  geom_line(aes(x = timeAf / (24 * 60), y = propCum * 100, 
                color = as.factor(temp), group = group)) +
  scale_x_continuous(
    name = "Time (days)", 
    breaks = seq(0, 42, by = 7)
  ) +
  coord_cartesian(xlim = c(0, 42)) +
  scale_y_continuous(
    name = "Cumulative proportion of germinated seeds",
    breaks = seq(0, 100, by = 20)
  ) +
  scale_colour_viridis_d(
    option = "H", 
    direction = 1, 
    name = "Temperature (\u00B0C)"
  ) +
  theme_bw(base_size = 8) +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    strip.background = element_rect(fill = "white"),
    strip.placement = "outside",
    strip.text = element_text(size = 8),
    strip.text.y.right = element_text(face = "italic"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

ggsave("output/imgs/Germination_all.png", plot, width = 7, height = 10, units = "in")

