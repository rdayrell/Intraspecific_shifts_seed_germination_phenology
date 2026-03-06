library(car)
library(ggplot2)
library(tidyverse)
library(abind)
library(patchwork)
library(dplyr)
library(ggpp)


# Import data ------------------------------------------------------------------

info_df <- read.csv("data/csv/serial_number_info.csv")

df1 <- read.csv("data/csv/germination_timetoevent_Alnus_glutinosa.csv")
df1$species = "Alnus glutinosa"

df2 <- read.csv("data/csv/germination_timetoevent_Betula_pubescens.csv")
df2$species = "Betula pubescens"

df3 <- read.csv("data/csv/germination_timetoevent_Pinus_sylvestris.csv")
df3$species = "Pinus sylvestris"

df = rbind(df1, df2, df3)
rm(df1, df2, df3)


# Prepare datasets -------------------------------------------------------------

df <- df %>%
  mutate(
    ID    = paste(serial_number, treat, temp, Units, sep = "_"),
    group = paste(serial_number, treat, temp, sep = "_"),
    treat = factor(treat, levels = c("NS", "COLD")),
    sp = species
  )

df <- df %>%
  mutate(across(c(treat, serial_number, ID, group, Units), as.factor),
         temp   = as.numeric(temp),
         count  = as.integer(count),
         nCum   = as.integer(nCum),
         timeBef = timeBef / (24 * 60),
         timeAf  = timeAf  / (24 * 60)) %>%
  filter(! ID %in% '945534_NS_15_4',
         !(serial_number == "912563" & temp == 15))



## Get ungerminated seeds
df_ungerminated <- df %>%
  dplyr::filter(is.infinite(timeAf)) %>%
  dplyr::select(ID, serial_number, treat, temp, Units, sp, group, total_ungerminated = count)

## Retain one representative row per ID for metadata
df_germinated <- df %>%
  dplyr::filter(is.finite(timeAf) & timeAf < 43) %>%
  group_by(ID) %>%
  arrange(timeBef) %>%
  slice_tail() %>%
  ungroup() %>%
  dplyr::select(ID, serial_number, treat, temp, Units, sp, group, total_germinated = nCum, propCum)

df_final <- df_ungerminated %>%
  left_join(df_germinated, by = c("ID", "serial_number", "treat", "temp", "Units", "sp", "group"))

info_df$serial_number <- as.factor(info_df$serial_number) 
df_group <- df_final %>%
  group_by(temp, treat, serial_number, Units, sp) %>%
  left_join(info_df, by = c("serial_number", "sp")) %>%
  arrange(lat_dd) %>%
  mutate(serial_number = factor(serial_number, levels = unique(serial_number))
         )

df_sum <- df_group %>%
  group_by(sp, temp, treat, serial_number, lat_dd, Geographical.Location, Altitude) %>%
  summarise(
    med       = median(propCum, na.rm = TRUE) * 100,
    error_top = quantile(propCum, 0.75, na.rm = TRUE) * 100,
    error_bot = quantile(propCum, 0.25, na.rm = TRUE) * 100,
    .groups   = "drop"
  )  %>% 
  arrange(lat_dd)

# Define inverse logit
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Logistic thresholds
logit_5 <- qlogis(0.05)  # logit(0.05)

# Fit quadratic model
smooth_data_all <- do.call(rbind, lapply(unique(df_group$sp), function(species) {
  df_sp <- df_group[df_group$sp == species, ]
  
  do.call(rbind, lapply(unique(df_sp$serial_number), function(serial) {
    df_serial <- df_sp[df_sp$serial_number == serial, ]
    
    do.call(rbind, lapply(unique(df_serial$treat), function(tr) {
      df_tmp <- df_serial[df_serial$treat == tr, ]
      
      fit <- glm(cbind(total_germinated, total_ungerminated) ~ temp + I(temp^2), data = df_tmp, family = binomial(link = "logit"))
      
      new_temp <- seq(min(df_tmp$temp), max(df_tmp$temp), length.out = 200)
      y_pred <- predict(fit, newdata = data.frame(temp = new_temp), type = "response")
      x_at_max <- new_temp[which.max(y_pred)]
      
      # Assume coefficients from the glm object (GLM version, but can be extended to GLMM fixed effects)
      b0 <- coef(fit)[1]  # intercept
      b1 <- coef(fit)[2]  # temp
      b2 <- coef(fit)[3]  # temp^2
      
      # Calculate Tmin and Tmax using the formula from the quadratic curve
      sqrt_term_5 <- sqrt(b1^2 - 4 * b2 * (b0 - logit_5))
      Tmin <- (-b1 - sqrt_term_5) / (2 * b2)
      Tmax <- (-b1 + sqrt_term_5) / (2 * b2)
      
      # Optimal temperature (vertex of the parabola)
      Topt <- -b1 / (2 * b2)
      
      # Calculate maximum germination (at Topt)
      max_germ <- inv_logit(b0 + b1 * Topt + b2 * Topt^2)
      logit_95 <- qlogis(0.95 * max_germ)
      
      # Calculate Topt lower and upper (95% of maximum)
      sqrt_term_95 <- sqrt(b1^2 - 4 * b2 * (b0 - logit_95))
      Topt_lower <- (-b1 - sqrt_term_95) / (2 * b2)
      Topt_upper <- (-b1 + sqrt_term_95) / (2 * b2)
      
      # Ensure Tmin < Tmax
      if (Tmin > Tmax) {
        tmp <- Tmin
        Tmin <- Tmax
        Tmax <- tmp
      }
      
      # Ensure Topt_lower < Topt_upper
      if (Topt_lower > Topt_upper) {
        tmp <- Topt_lower
        Topt_lower <- Topt_upper
        Topt_upper <- tmp
      }
      
      breadth <- Tmax - Tmin
      opt_breadth <- Topt_upper - Topt_lower
      
      # Return results
      data.frame(
        temp = new_temp,
        y = y_pred * 100, 
        Tmin = Tmin,
        Tmax = Tmax,
        x_at_max = x_at_max,
        Topt = Topt,
        Topt_lower = Topt_lower,
        Topt_upper = Topt_upper,
        Breadth = breadth,
        Optimum_Breadth = opt_breadth,
        group = tr,
        serial_number = serial,
        sp = species,
        lat_dd = unique(df_tmp$lat_dd),
        stringsAsFactors = FALSE
      )
    }))
  }))
}))

names(smooth_data_all)
results_GLM = smooth_data_all %>% 
  dplyr::select(-temp, -y) %>%
  unique() 

write.csv(results_GLM, "output/csv/finalgermination_6weeks.csv", row.names = F)

loc_labels <- df_group %>%
  group_by(serial_number) %>%
  summarise(Geographical.Location = unique(Geographical.Location),
            sp = unique(sp),
            Altitude = unique(Altitude)) %>%
  mutate(Geographical.Location = sub(" - ", "\n", Geographical.Location),
         Geographical.Location = ifelse(sp == "Pinus sylvestris", 
                                        paste0(Geographical.Location, "\n", Altitude, " m"), 
                                        Geographical.Location)) %>%
  dplyr::select(-Altitude) %>%
  distinct(serial_number, Geographical.Location) %>%
  deframe()


# function to make one species‐level facet plot
make_sp_plot <- function(species_name, show_legend = FALSE) {
  
  summary_opt = smooth_data_all %>% dplyr::filter(sp == species_name) %>%
    dplyr::select(-temp, -y) %>%
    unique() %>%
    mutate(
      yplot = ifelse(group == "NS", 109, 115)
    )
  # write.csv(summary_opt, paste0("output/finalgermination_6months_", species_name,".csv"), row.names = FALSE)
  
  sd_sp  <- smooth_data_all %>% dplyr::filter(sp == species_name)
  sum_sp <- df_sum %>% dplyr::filter(sp == species_name)
  
  if (species_name == "Pinus sylvestris") {
    sum_sp <- sum_sp %>% arrange(Altitude)
    sum_sp$serial_number <- factor(sum_sp$serial_number,
                                   levels = unique(sum_sp$serial_number))
    
    # reorder smooth data to match
    sd_sp$serial_number <- factor(sd_sp$serial_number,
                                  levels = levels(sum_sp$serial_number))
  }
  
  else {
    sum_sp <- sum_sp %>% arrange(lat_dd)
    sum_sp$serial_number <- factor(sum_sp$serial_number,
                                   levels = unique(sum_sp$serial_number))
    
    # reorder smooth data to match
    sd_sp$serial_number <- factor(sd_sp$serial_number,
                                  levels = levels(sum_sp$serial_number))
  }
  
  p <- ggplot() +
    

    geom_line(data = sd_sp,
              aes(x = temp, y = pmax(y, 0), colour = group, group = group),
              linewidth = 0.6) +
    geom_rect(data = summary_opt,
              aes(xmin = Topt_lower, xmax = Topt_upper,
                  ymin = yplot, ymax = yplot + 2.5,
                  fill = group)) +
    geom_rect(data = summary_opt,
             aes(xmin = Topt_lower, xmax = Topt_upper,
                 ymin = yplot, ymax = yplot + 2.5,
                 fill = group)) +
    geom_vline(data = summary_opt,
               aes(xintercept = Topt, colour = group),
               linewidth = 0.4, linetype = "dashed") +
    geom_errorbar(data = sum_sp,
                  aes(x = temp, ymin = error_bot, ymax = error_top, colour = treat),
                  width = 1.5, alpha = 0.7) +
    geom_point(data = sum_sp,
               aes(x = temp, y = med, colour = treat),
               size = 2.5, alpha = 0.7) +
    facet_wrap(. ~ serial_number, ncol = 6, 
               labeller = labeller(serial_number = loc_labels)) +
    ggtitle(species_name) +
    scale_fill_manual(name   = "Treatment", guide = "none",
                        values = c("NS" = "deepskyblue3", "COLD" = "navy")) +
    scale_colour_manual(name   = "Treatment",
                        values = c("NS" = "deepskyblue3", "COLD" = "navy")) +
    scale_x_continuous(name = "Temperature (\u00B0 C)", limits = c(0,45),
                       breaks = seq(0, 45, by = 10)) +
    scale_y_continuous("Germination % (6 weeks)", limits = c(0, 119),
                       breaks = seq(0, 100, by = 25)) +
    theme_bw() +
    theme(
      legend.position  = if (show_legend) "top" else "none",
      plot.title = element_text(face = "italic", size = 12),
      strip.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank()
    )
  
  p
}

plot_Alnus  <- make_sp_plot("Alnus glutinosa",   show_legend = TRUE)
plot_Betula <- make_sp_plot("Betula pubescens",  show_legend = FALSE)
plot_Pinus  <- make_sp_plot("Pinus sylvestris",  show_legend = FALSE)

final_plot <- (plot_Alnus / plot_Betula / plot_Pinus) +
  plot_layout(guides = "collect", axis_titles = "collect_x") &
  theme(legend.position = "bottom")


final_plot

ggsave("output/imgs/Fig_2_Germination_percentage_at_6_months.png", final_plot, width = 9.5, height = 7.5, units = "in")

# ggsave("output/imgs/production/Fig_2_Germination_percentage_at_6_months.pdf", final_plot, width = 9.5, height = 7.5, units = "in")
# 
# ggsave("output/imgs/production/Fig_2_Germination_percentage_at_6_months.tif", final_plot, width = 9.5, height = 7.5, units = "in", dpi = 600, compression = "lzw")
