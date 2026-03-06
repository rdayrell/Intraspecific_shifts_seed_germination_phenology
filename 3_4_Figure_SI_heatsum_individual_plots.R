
library(tidyverse)
library(forcats)

load("output/rda/plot_heatsum_individual_data.Rda")


# Sequence for Alnus and Betula by lat_dd
ab_seq <- plot_data %>%
  dplyr::select(Genus, lat_dd, serial_number) %>%
  unique() %>%
  filter(Genus %in% c("Alnus", "Betula")) %>%
  arrange(Genus, lat_dd) %>%
  mutate(seq = row_number()) %>%
  dplyr::select(-lat_dd)

# Sequence for Pinus by Altitude, continuing from last seq
pinus_seq <- plot_data %>%
  dplyr::select(Genus, Altitude, serial_number) %>%
  unique() %>%
  filter(Genus == "Pinus") %>%
  arrange(Altitude) %>%
  mutate(seq = row_number() + max(ab_seq$seq, na.rm = TRUE)) %>%
  dplyr::select(-Altitude)

all_seq <- rbind(ab_seq, pinus_seq)

# Combine sequences
plot_data <- merge(all_seq, plot_data, by = c("serial_number", "Genus"))


# ------------------------------------------------------------------------------
# PRESENT

plot_data_1 <- plot_data %>%
  dplyr::filter(treat_time == "COLD - 1981 - 2010" | treat_time == "NS - 1981 - 2010")

total_species <- length(unique(plot_data$seq))
species_list <- unique(plot_data$seq)


y2 = ifelse(max(plot_data$T_max, na.rm = T) > max(plot_data$cumulative_heat_sum, na.rm = T)/80,
            max(plot_data$T_max, na.rm = T),
            max(plot_data$cumulative_heat_sum, na.rm = T)/80)
y2 = ifelse(max(plot_data$Tb, na.rm = T) > y2,
            max(plot_data$Tb, na.rm = T),
            y2)
x2 = min(plot_data$date)
x1 = max(plot_data$date) - 20

plots <- list()
plots2 <- list()
plots3 <- list()

for (i in seq_along(species_list)) {
  
  s <- species_list[i] 
  
  dp <- filter(plot_data_1, seq == s) 
  dp$sp <- paste(dp$Genus, dp$Species, sep = " ")
  
  lab = unique(dp$location2)
  lab2 = paste(unique(dp$sp), round(unique(dp$lat_dd),2), unique(dp$Altitude))
  Tb = unique(dp$Tb)
  (st = paste(unique(dp$start), collapse = " - "))
  
  plot = ggplot(dp, aes(x = date)) +
      
      geom_area(aes(y = cumulative_heat_sum/80, fill = treat_time, colour = treat_time), alpha = 0.7, position = "identity")+
      scale_fill_manual(name = "",
                        values = c("NS - 1981 - 2010" = "yellowgreen",
                                   "COLD - 1981 - 2010" = "darkgreen"
                        ))+
      scale_colour_manual(name = "",
                          values = c("NS - 1981 - 2010" = "yellowgreen",
                                     "COLD - 1981 - 2010" = "darkgreen"
                          ))+
      
      geom_line(aes(y = CHELSA_pr / 30), colour = "deepskyblue3", size = 0.6, alpha = 0.7) +
      
      geom_line(aes(y = CHELSA_tas), size = 0.3, color = "grey30") +
      geom_ribbon(aes(ymin = CHELSA_tasmin, ymax = CHELSA_tasmax), fill = "grey50", alpha = 0.1, position = "identity") +
      # geom_line(aes(y = CHELSA_tasmin), size = 0.1, linetype = "dashed", color = "grey50") +
      # geom_line(aes(y = CHELSA_tasmax), size = 0.1, linetype = "dashed", color = "grey50") +
      
      geom_hline(aes(yintercept = Tb, colour = treat_time), size = 0.5, linetype = "dotted") +
      annotate(geom = "text", y = Tb + 1.2, x = x1, label = "Tb", size = 2, fontface = "italic", hjust = 0) + 
      annotate(geom = "text", x = x2, y = y2-2, label = lab, size = 2.2, hjust = 0) + 
      # ggtitle(if (i == 1) bquote(italic(.(lab2)) ~ .(unique(dp$mod)) ~ .(unique(dp$type)) ~ .(st)) else "") +
      # ggtitle(bquote(italic(.(lab2)))) +
      
      scale_y_continuous(
        name = "Temperature (°C)",
        sec.axis = sec_axis(~ .  , name = expression("Precipitation (kg"~ m^{-2} ~ d^{-1} * 10^{-2}~")")),
        limits = c(-1, y2)) +
      scale_x_date(name = "", date_labels = "%b", date_breaks = "1 month", limits = c(as.Date(x2), as.Date(x1))) + 
      theme_bw(base_size =9) +
      theme(
        strip.background.y = element_rect(fill = "white"),
        legend.position = "none", #if (i == 1) "bottom" else "none",
        title = element_text(size = 6),
        strip.text.y = if (i == total_species) element_text() else element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        axis.title.y.left = element_text(size = 6),
        axis.title.y.right = element_blank() 
      )
  plots[[s]] <- plot
  
  
  p2 = ggplot(dp, aes(x = date)) +
      geom_area(aes(y = cumulative_heat_sum/80, fill = treat_time, colour = treat_time), alpha = 0.7, position = "identity")+
      scale_fill_manual(name = "",
                        values = c("NS - 1981 - 2010" = "yellowgreen",
                                   "COLD - 1981 - 2010" = "darkgreen"
                        ))+
      scale_colour_manual(name = "",
                          values = c("NS - 1981 - 2010" = "yellowgreen",
                                     "COLD - 1981 - 2010" = "darkgreen"
                          ))+
      scale_y_continuous(
        labels = function(x) x * 80,
        name = "Cumulative Heat Sum (°Cd)",
        limits = c(0, y2)) +
      scale_x_date(date_labels = "%b", date_breaks = "1 month", limits = c(as.Date(x2), as.Date(x1))) +
      theme_bw(base_size =6) +
      theme(
        legend.position = "none"
      )
  
  plots2[[s]] <- p2
  
  
  p3 = wrap_elements(get_plot_component(p2, "ylab-l")) + 
      wrap_elements(get_y_axis(p2)) + 
      plot +
      plot_layout(widths = c(2, 2, 40))
  
  plots3[[s]] <- p3
  
}

    
plot3a = plots3[[1]] / plots3[[2]] / plots3[[3]] / plots3[[4]] / plots3[[5]]/ plots3[[6]] +
          plot_layout(axis_titles = "collect", axes = "collect", guides = "collect") & 
          scale_x_date(name = "", date_breaks = "1 month", labels = function(x) substr(month.abb[as.POSIXlt(x)$mon + 1], 1, 1)) &
          theme(plot.margin = margin(t=-1, r=1, b=-1, l=4),
                legend.position = "none")

# save_name = paste0("output/heatsum_Alnus.png")
# ggsave(save_name, plot3a, width = 4, height = 9, units = "in")


plot3b = plots3[[7]] / plots3[[8]] / plots3[[9]] / plots3[[10]] / plots3[[11]]/ plots3[[12]] +
          plot_layout(axis_titles = "collect", axes = "collect", guides = "collect") & 
          scale_x_date(name = "", date_breaks = "1 month", labels = function(x) substr(month.abb[as.POSIXlt(x)$mon + 1], 1, 1)) &
          theme(plot.margin = margin(t=-1, r=1, b=-1, l=4),
                legend.position = "none")

# save_name = paste0("output/heatsum_Betula.png")
# ggsave(save_name, plot3b, width = 4, height = 9, units = "in")


plot3p = plots3[[13]] / plots3[[14]] / plots3[[15]] / plots3[[16]] / plots3[[17]]/ plots3[[18]] +    
          plot_layout(axis_titles = "collect", axes = "collect", guides = "collect") &
          scale_x_date(name = "", date_breaks = "1 month", labels = function(x) substr(month.abb[as.POSIXlt(x)$mon + 1], 1, 1)) &
          theme(plot.margin = margin(t=-1, r=1, b=-1, l=4),
                legend.position = "none")

# save_name = paste0("output/heatsum_Pinus.png")
# ggsave(save_name, plot3p, width = 4, height = 9, units = "in")


# ------------------------------------------------------------------------------
# OPTIMISTIC FUTURE

plots <- list()
plots2 <- list()
plots3 <- list()

plot_data_2 <- plot_data %>%
  dplyr::filter(treat_time == "COLD - 2071 - 2100 (optimistic)" | treat_time == "NS - 2071 - 2100 (optimistic)") %>%
  mutate(treat_time = gsub(" \\(", "\n(", treat_time))



for (i in seq_along(species_list)) {
  
  s <- species_list[i] 
  
  dp <- filter(plot_data_2, seq == s) 
  dp$sp <- paste(dp$Genus, dp$Species, sep = " ")
  
  lab = unique(dp$location2)
  lab2 = paste(unique(dp$sp), round(unique(dp$lat_dd),2), unique(dp$Altitude))
  Tb = unique(dp$Tb)
  (st = paste(unique(dp$start), collapse = " - "))
  
  plot = ggplot(dp, aes(x = date)) +
      
      #facet_grid(location ~ ., scales = "fixed") +
      
      geom_area(aes(y = cumulative_heat_sum/80, fill = treat_time, colour = treat_time), alpha = 0.7, position = "identity")+
      scale_fill_manual(name = "",
                        values = c("NS - 2071 - 2100\n(optimistic)" = "goldenrod1",
                                   "COLD - 2071 - 2100\n(optimistic)" = "goldenrod4"
                        )) +
      scale_colour_manual(name = "",
                          values = c("NS - 2071 - 2100\n(optimistic)" = "goldenrod1",
                                     "COLD - 2071 - 2100\n(optimistic)" = "goldenrod4"
                                     )) +

      geom_line(aes(y = ssp126_pr / 30), colour = "deepskyblue3", size = 0.6, alpha = 0.7) +
      
      geom_line(aes(y = ssp126_tas), size = 0.3, color = "grey30") +
      # geom_line(aes(y = ssp126_tasmin), size = 0.2, linetype = "dashed", color = "grey50") +
      # geom_line(aes(y = ssp126_tasmax), size = 0.2, linetype = "dashed", color = "grey50") +
      geom_ribbon(aes(ymin = ssp126_tasmin, ymax = ssp126_tasmax), fill = "grey50", alpha = 0.1, position = "identity") +

      geom_hline(aes(yintercept = Tb, colour = treat_time), size = 0.5, linetype = "dotted") +
      annotate(geom = "text", y = Tb + 1.2, x = x1, label = "Tb", size = 2, fontface = "italic", hjust = 0) + 
      annotate(geom = "text", x = x2, y = y2-2, label = lab, size = 2.2, hjust = 0) + 
 
      # ggtitle(bquote(italic(.(lab2)))) +
      
      scale_y_continuous(
        name = "Temperature (°C)",
        sec.axis = sec_axis(~ .  , name = expression("Precipitation (kg"~ m^{-2} ~ d^{-1} * 10^{-2}~")")),
        limits = c(-1, y2)) +
      scale_x_date(name = "", date_labels = "%b", date_breaks = "1 month", limits = c(as.Date(x2), as.Date(x1))) + 
      theme_bw(base_size =9) +
      theme(
        strip.background.y = element_rect(fill = "white"),
        legend.position = "none",
        title = element_text(size = 6),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        strip.text.y = if (i == total_species) element_text() else element_blank(), 
        #axis.text.x = if (i == total_species) element_text() else element_blank(),
        axis.title.y.left = element_blank(), # element_text(size = 6),
        axis.title.y.right = element_blank() # element_text(size = 6)
      )
  plots[[s]] <- plot
  
  
  p2 = ggplot(dp, aes(x = date)) +
      #facet_grid(location ~ ., scales = "fixed") +
      
      #annotate("segment", x=x2,xend=x2, y=0,yend = y2) +
      
      # geom_area(aes(y = cumulative_heat_sum/80, fill = treat, colour = treat), alpha = 0.6)+
      # scale_fill_manual(name = "Stratification", values = c("NS" = "seagreen4", "COLD" = "limegreen")) +
      # scale_colour_manual(name = "Stratification", values = c("NS" = "seagreen4", "COLD" = "limegreen")) +
      
      geom_area(aes(y = cumulative_heat_sum/80, fill = treat_time, colour = treat_time), alpha = 0.7, position = "identity")+
      scale_fill_manual(name = "",
                        values = c("NS - 2071 - 2100\n(optimistic)" = "goldenrod1",
                                   "COLD - 2071 - 2100\n(optimistic)" = "goldenrod4"
                        )) +
      scale_colour_manual(name = "",
                          values = c("NS - 2071 - 2100\n(optimistic)" = "goldenrod1",
                                     "COLD - 2071 - 2100\n(optimistic)" = "goldenrod4"
                                     )) +
      scale_y_continuous(
        labels = function(x) x * 80,
        name = "Cumulative Heat Sum (°Cd)",
        limits = c(0, y2)) +
      scale_x_date(date_labels = "%b", date_breaks = "1 month", limits = c(as.Date(x2), as.Date(x1))) +
      guides(fill = guide_legend(order =1),
             colour = guide_legend(order =1)) +
      theme_bw(base_size =6) +
      theme(
        axis.title.y.left = element_blank(),
        legend.position = "none"
      )
  
  plots2[[s]] <- p2
  
  
  p3 = wrap_elements(get_plot_component(p2, "ylab-l")) + 
      wrap_elements(get_y_axis(p2)) + 
      plot +
      plot_layout(widths = c(1, 2, 40))
  
  plots3[[s]] <- p3
  
}

plot4a = plots3[[1]] / plots3[[2]] / plots3[[3]] / plots3[[4]] / plots3[[5]] / plots3[[6]] +
          plot_layout(axis_titles = "collect", axes = "collect", guides = "collect") &
          scale_x_date(name = "", date_breaks = "1 month", labels = function(x) substr(month.abb[as.POSIXlt(x)$mon + 1], 1, 1)) &
          theme(plot.margin = margin(t=-1, r=1, b=-1, l=0),
                legend.position = "none")

plot4b = plots3[[7]] / plots3[[8]] / plots3[[9]] / plots3[[10]] / plots3[[11]] / plots3[[12]] +
          plot_layout(axis_titles = "collect", axes = "collect", guides = "auto") &
          scale_x_date(name = "", date_breaks = "1 month", labels = function(x) substr(month.abb[as.POSIXlt(x)$mon + 1], 1, 1)) &
          theme(plot.margin = margin(t=-1, r=1, b=-1, l=0))
  
plot4p = plots3[[13]] / plots3[[14]] / plots3[[15]] / plots3[[16]] / plots3[[17]] / plots3[[18]] +
    plot_layout(axis_titles = "collect", axes = "collect", guides = "collect") &
    scale_x_date(name = "", date_breaks = "1 month", labels = function(x) substr(month.abb[as.POSIXlt(x)$mon + 1], 1, 1)) &
    theme(plot.margin = margin(t=-1, r=1, b=-1, l=0),
          legend.position = "none")



# ------------------------------------------------------------------------------
# PESSIMISTIC FUTURE

plots <- list()
plots2 <- list()
plots3 <- list()

plot_data_3 <- plot_data %>%
  dplyr::filter(treat_time == "COLD - 2071 - 2100 (pessimistic)" | treat_time == "NS - 2071 - 2100 (pessimistic)") %>%
  mutate(treat_time = gsub(" \\(", "\n(", treat_time))

unique(plot_data_3$treat_time)

for (i in seq_along(species_list)) {
  
  s <- species_list[i] 
  
  dp <- filter(plot_data_3, seq == s) 
  dp$sp <- paste(dp$Genus, dp$Species, sep = " ")
  
  lab = unique(dp$location2)
  lab2 = paste(unique(dp$sp), round(unique(dp$lat_dd),2), unique(dp$Altitude))
  Tb = unique(dp$Tb)
  (st = paste(unique(dp$start), collapse = " - "))
  
  plot = ggplot(dp, aes(x = date)) +
      
      #facet_grid(location ~ ., scales = "fixed") +
      
      geom_area(aes(y = cumulative_heat_sum/80, fill = treat_time, colour = treat_time), alpha = 0.7, position = "identity")+
      scale_fill_manual(name = "",
                        values = c("NS - 2071 - 2100\n(pessimistic)" = "coral1",
                                   "COLD - 2071 - 2100\n(pessimistic)" = "red4")) +
      scale_colour_manual(name = "",
                          values = c("NS - 2071 - 2100\n(pessimistic)" = "coral1",
                                     "COLD - 2071 - 2100\n(pessimistic)" = "red4")) +
      # geom_area(aes(y = cumulative_heat_sum/80, fill = treat, colour = treat), alpha = 0.7)+
      #   scale_fill_manual(name = "Stratification", values = c("NS" = "seagreen4", "COLD" = "limegreen")) +
      #   scale_colour_manual(name = "Stratification", values = c("NS" = "seagreen4", "COLD" = "limegreen")) +
      
      geom_line(aes(y = ssp585_pr / 30), colour = "deepskyblue3", size = 0.6, alpha = 0.7) +
      
      geom_line(aes(y = ssp585_tas), size = 0.3, color = "grey30") +
      # geom_line(aes(y = ssp585_tasmin), size = 0.2, linetype = "dashed", color = "grey50") +
      # geom_line(aes(y = ssp585_tasmax), size = 0.2, linetype = "dashed", color = "grey50") +
      
      geom_ribbon(aes(ymin = ssp585_tasmin, ymax = ssp585_tasmax), fill = "grey50", alpha = 0.1, position = "identity") +
      
      geom_hline(aes(yintercept = Tb, colour = treat_time), size = 0.5, linetype = "dotted") +
      annotate(geom = "text", y = Tb + 1.2, x = x1, label = "Tb", size = 2, fontface = "italic", hjust = 0) + 
      annotate(geom = "text", x = x2, y = y2-2, label = lab, size = 2.2, hjust = 0) + 
      # ggtitle(if (i == 1) bquote(italic(.(lab2)) ~ .(unique(dp$mod)) ~ .(unique(dp$type)) ~ .(st)) else "") +
      
      scale_y_continuous(
        name = "Temperature (°C)",
        sec.axis = sec_axis(~ .  , name = expression("Precipitation (kg"~ m^{-2} ~ d^{-1} * 10^{-2}~")")),
        limits = c(-1, y2)) +
      scale_x_date(name = "", date_labels = "%b", date_breaks = "1 month", limits = c(as.Date(x2), as.Date(x1))) + 
      theme_bw(base_size =9) +
      theme(
        strip.background.y = element_rect(fill = "white"),
        legend.position = "none", # if (i == total_species) "bottom" else "none",
        title = element_text(size = 6),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        strip.text.y = if (i == total_species) element_text() else element_blank(), 
        #axis.text.x = if (i == total_species) element_text() else element_blank(),
        axis.title.y.left = element_blank(), #element_text(size = 6),
        axis.title.y.right = element_text(size = 6)
      )
  plots[[s]] <- plot
  
  
  p2 = ggplot(dp, aes(x = date)) +
      #facet_grid(location ~ ., scales = "fixed") +
      
      #annotate("segment", x=x2,xend=x2, y=0,yend = y2) +
      
      # geom_area(aes(y = cumulative_heat_sum/80, fill = treat, colour = treat), alpha = 0.6)+
      # scale_fill_manual(name = "Stratification", values = c("NS" = "seagreen4", "COLD" = "limegreen")) +
      # scale_colour_manual(name = "Stratification", values = c("NS" = "seagreen4", "COLD" = "limegreen")) +
      
      geom_area(aes(y = cumulative_heat_sum/80, fill = treat_time, colour = treat_time), alpha = 0.7, position = "identity")+
      scale_fill_manual(name = "",
                        values = c("NS - 2071 - 2100\n(pessimistic)" = "coral1",
                                   "COLD - 2071 - 2100\n(pessimistic)" = "red4")) +
      scale_colour_manual(name = "",
                          values = c("NS - 2071 - 2100\n(pessimistic)" = "coral1",
                                     "COLD - 2071 - 2100\n(pessimistic)" = "red4")) +
      scale_y_continuous(
        labels = function(x) x * 80,
        name = "Cumulative Heat Sum (°Cd)",
        limits = c(0, y2)) +
      scale_x_date(date_labels = "%b", date_breaks = "1 month", limits = c(as.Date(x2), as.Date(x1))) +
      theme_bw(base_size =6) +
      theme(
        axis.title.y.left = element_blank(),
        legend.position = "none"
      )
  
  plots2[[s]] <- p2
  
  
  p3 = wrap_elements(get_plot_component(p2, "ylab-l")) + 
      wrap_elements(get_y_axis(p2)) + 
      plot +
      plot_layout(widths = c(1, 2, 40))
  
  plots3[[s]] <- p3
  
}


plot5p = plots3[[13]] / plots3[[14]] / plots3[[15]] / plots3[[16]] / plots3[[17]] / plots3[[18]] +
    plot_layout(axis_titles = "collect", axes = "collect", guides = "auto") &
    scale_x_date(name = "", date_breaks = "1 month", labels = function(x) substr(month.abb[as.POSIXlt(x)$mon + 1], 1, 1)) &
    theme(plot.margin = margin(t=-1, r=2, b=-1, l=0),
          legend.position = "none")

plot5a = plots3[[1]] / plots3[[2]] / plots3[[3]] / plots3[[4]] / plots3[[5]] / plots3[[6]] + plot_layout(axis_titles = "collect", axes = "collect", guides = "auto") & scale_x_date(name = "", date_breaks = "1 month", labels = function(x) substr(month.abb[as.POSIXlt(x)$mon + 1], 1, 1)) &
    theme(plot.margin = margin(t=-1, r=2, b=-1, l=0),
          legend.position = "none")

plot5b = plots3[[7]] / plots3[[8]] / plots3[[9]] / plots3[[10]] / plots3[[11]] / plots3[[12]] +
    plot_layout(axis_titles = "collect", axes = "collect", guides = "auto") &
    scale_x_date(name = "", date_breaks = "1 month", labels = function(x) substr(month.abb[as.POSIXlt(x)$mon + 1], 1, 1)) &
    theme(plot.margin = margin(t=-1, r=2, b=-1, l=0),
          legend.position = "none")


# ------------------------------------------------------------------------------

plot_legend <- plot_data %>%
  filter(is.finite(cumulative_heat_sum), !is.na(date)) %>%
  mutate(treat_time = gsub(" \\(", "\n(", treat_time))

plot_legend$treat_time <- factor(plot_legend$treat_time, levels = c(
  "NS - 1981 - 2010",
  "COLD - 1981 - 2010",
  "NS - 2071 - 2100\n(optimistic)",
  "COLD - 2071 - 2100\n(optimistic)",
  "NS - 2071 - 2100\n(pessimistic)",
  "COLD - 2071 - 2100\n(pessimistic)"
))

p <- ggplot(plot_legend, aes(x = date)) +
  geom_area(aes(y = cumulative_heat_sum, fill = treat_time), alpha = 0.7) +
  scale_fill_manual(name = "", values = c(
    "NS - 1981 - 2010" = "yellowgreen",
    "COLD - 1981 - 2010" = "darkgreen",
    "NS - 2071 - 2100\n(optimistic)" = "goldenrod1",
    "COLD - 2071 - 2100\n(optimistic)" = "goldenrod4",
    "NS - 2071 - 2100\n(pessimistic)" = "coral1",
    "COLD - 2071 - 2100\n(pessimistic)" = "red4"
  )) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        guide_legend(nrow =1)
        )

legend = cowplot::get_plot_component(p, 'guide-box-bottom', return_all = TRUE)
cowplot::ggdraw(legend)

## account for different margins
total_margin <- c(5, 1, 2)
inverse_margin <- 1 / total_margin
normalized_widths <- inverse_margin / sum(inverse_margin)


plot_vert_a = plot3a | plot4a | plot5a + plot_layout(axis_titles = "collect", axes = "collect", guides = "auto", widths = normalized_widths)
plot_vert_b = plot3b | plot4b | plot5b + plot_layout(axis_titles = "collect", axes = "collect", guides = "auto", widths = normalized_widths)
plot_vert_p = plot3p | plot4p | plot5p + plot_layout(axis_titles = "collect", axes = "collect", guides = "auto", widths = normalized_widths)

plot_vert_al <- (plot_vert_a / legend) +  plot_layout(heights = c(16, 1))
plot_vert_bl <- (plot_vert_b / legend) +  plot_layout(heights = c(16, 1))
plot_vert_pl <- (plot_vert_p / legend) +  plot_layout(heights = c(16, 1))

save_name = paste0("output/imgs/Fig_SI_heatsum_Alnus.png")
ggsave(save_name, plot_vert_al, width = 10, height = 9, units = "in")

save_name = paste0("output/imgs/Fig_SI_heatsum_Betula.png")
ggsave(save_name, plot_vert_bl, width = 10, height = 9, units = "in")

save_name = paste0("output/imgs/Fig_SI_heatsum_Pinus.png")
ggsave(save_name, plot_vert_pl, width = 10, height = 9, units = "in")

