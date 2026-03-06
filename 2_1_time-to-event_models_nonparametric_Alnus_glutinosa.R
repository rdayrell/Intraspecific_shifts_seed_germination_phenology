library(tidyverse)

library(drcSeedGerm)
library(drcte)

library(lme4)
library(lmtest)
library(sandwich)

library(ggplot2)
library(abind)
library(patchwork)

options(scipen=999)

# Import data ------------------------------------------------------------------

df <- read.csv("data/csv/germination_timetoevent_Alnus_glutinosa.csv")

## Check dataset

ggplot(df, aes(timeAf/(24*60), propCum)) +
  facet_grid(serial_number ~ treat) +
  geom_line(aes(group=ID, colour=as.factor(temp))) +
  scale_x_continuous(name = "Time (days)") +
  scale_y_continuous(name = "Cumulative proportion of germinated seeds") +
  theme_bw()


df %>%
  group_by(serial_number,treat,temp,timeAf,group) %>%
  summarise(propCum=mean(propCum)) %>%
  ggplot() +
  facet_grid(serial_number ~ treat) +
  geom_line(aes(x= timeAf/(24*60), y=propCum*100, color=as.factor(temp), group=group)) +
  scale_x_continuous(name = "Time (days)") +
  scale_y_continuous(name = "Cumulative proportion of germinated seeds") +
  theme_bw()


## Data prep

df <- df %>% 
  mutate(across(c(treat, serial_number, ID, group, Units), .fns = factor),
         temp = as.numeric(temp),
         count = as.integer(count), 
         nCum = as.integer(nCum),
         timeBef = timeBef/(24*60), # converting from minutes to days
         timeAf = timeAf/(24*60)) # converting from minutes to days


# Fit model --------------------------------------------------------------------
## Germination progress in each dish (Units) 

mod1 <- drmte(count ~ timeBef + timeAf,
              curveid = ID,
              fct = NPMLE(), 
              data = df,
              separate = T)
 

## check output
tab <- plotData(mod1)

tab2 <- tab$plotPoints %>% separate(col = "ID", into = c("serial_number","treat", "temp","dish"),
                                    sep = "_", remove = FALSE) #%>% 

ggplot() +
  geom_line(data = tab2, mapping = aes(x = timeAf, y = CDF, group = ID, colour = as.factor(temp))) +
  facet_wrap(serial_number ~ treat) +
  scale_x_continuous(name = "Time (d)") +
  scale_y_continuous(name = "Cumulative proportion of germinated seeds") +
  theme_bw()


# Calculate GR50 ---------------------------------------------------------------

g <- 0.5

set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")

## This takes a very long time to run!!! File can be directly loaded below under fit threshold model
GR <- quantile(mod1, probs = g, robust = T, interval = T, units = ID, rate = T,  restricted = F,  display = F)

GRlist <- tibble(serial_number = row.names(GR), GR, row.names = NULL) %>% 
  separate(col = "serial_number", into = c("serial_number","treat", "temp","dish"),
           sep = "_", remove = FALSE) %>% 
  mutate(temp = as.numeric(temp)) 

GRlist$g <- "50%"

saveRDS(GRlist,file="output/rds/GR_mod1_Alnus_glutinosa.rds")

# Plot the GR values against temperature:
ggplot(GRlist,aes(x= temp, y=Mean, color=g)) +
  facet_grid(serial_number ~ treat) +
  geom_point()+
  labs(x= "Temperature (\u00B0 C)" , y=expression(paste("Germination rate (", d^-1, ")", sep = "")), color="Percentile" ) +
  theme_bw()


# Fit threshold model ----------------------------------------------------------

# GRlist <- readRDS("output/rds/GR_mod1_Alnus_glutinosa.rds")

GRlist$ID <- paste(GRlist$serial_number, GRlist$treat, GRlist$temp, GRlist$dish, sep = "_")
GRlist_group = droplevels(GRlist)

estfun.drc <- drc::estfun.drc # fixing bug sandwich

Out<-NULL
Predictions<-NULL
mod_names <- c("modGR_Ex") # "modGR_BS"

error_list <- list()
warning_list <- list()
error_index <- 1
warning_index <- 1


for (s in unique(GRlist_group$serial_number)) {
  for (t in unique(GRlist_group$treat)) {
    for (i in unique(GRlist_group$g)) {
      
        dat = subset(GRlist_group, serial_number == s & treat == t & g == i)
        dat = dat %>%
          group_by_at(vars(temp, g, serial_number, treat, dish)) %>%
          summarise(
            Median = median(Median))
        
        #https://www.statforbiology.com/2023/stat_drcte_12-htt2step/#grt.bs
        models <- list(
          # modGR_BS = drm(Median ~ temp, data = dat, fct = GRT.BS()),
          modGR_Ex = drm(Median ~ temp, data = dat, fct = GRT.Ex())
          #modGR_RF = drm(Median ~ temp, data = dat, fct = GRT.RF())
          #modGR_M = drm(Median ~ temp, data = dat, fct = GRT.M()),
          #modGR_YL = drm(Median ~ temp, data = dat, fct = GRT.YL())
        )
        
        for (m in mod_names) {
          
          tryCatch({
            withCallingHandlers({
              
          
          model <- models[[m]]
          sm2 <- summary(model)
          sm2_df <- as.data.frame(sm2$coefficients)

          # I should do this according to Onofri, but it is giving me unrealistically small std errors          
          # coef_mod <- coeftest(model, vcov. = sandwich)
          # sm2_df <- coef_mod[complete.cases(coef_mod), ]
          # sm2_df <- as.data.frame(sm2_df)
          
        }, warning = function(w) {
          warning_list[[warning_index]] <<- paste("serial_number =", s, "treat =", t, "g =", i, "m =", m, warning = conditionMessage(w))
          warning_index <<- warning_index + 1
          invokeRestart("muffleWarning")
        })
        }, error = function(e) {
          error_list[[error_index]] <<- paste("serial_number =", s, "treat =", t, "g =", i, "m =", m, error = conditionMessage(e))
          error_index <<- error_index + 1
        })
      
    
          # out_temp = tibble(coefficients_mod2 = row.names(sm2_df), sm2_df, row.names = NULL) %>% 
          #  separate_wider_delim("coefficients_mod2", delim = ":", names = c("coefficients_mod2", "factor")) %>%
          #   remove_rownames()
          out_temp <- sm2_df %>%
            rownames_to_column(var = "coefficients_mod2") %>%
            separate_wider_delim(coefficients_mod2, delim = ":", 
                                 names = c("coefficients_mod2", "factor"))
          
          out_temp$serial_number <- paste(s)
          out_temp$treat <- paste(t)
          out_temp$g <- paste(i)
          out_temp$mod <- paste(m)
      
          #t_grid = expand.grid(temp=seq(0,40, by=0.1),GR = g*100)
          t_grid = data.frame(temp=seq(0,43, by=0.1))
          t_pred = as.data.frame(cbind(t_grid, predict(models[[m]], t_grid)))
          colnames(t_pred)[2] <-"pred"
          t_pred$serial_number <- paste(s)
          t_pred$treat <- paste(t)
          t_pred$g <- paste(i)
          t_pred$mod <- paste(m)
          
          out_temp2 = out_temp[1,c(1:length(colnames(out_temp)))]
          out_temp2[1,c(1:6)] = NA
          out_temp2$coefficients_mod2 <- "To"
          out_temp2$Estimate = t_pred$temp[c(which.max(t_pred$pred))]
          
          # Find the first positive number
          first_positive_index <- which(t_pred$pred > 0)[1]
          
          # Find the first zero after the positive numbers start
          first_zero_after_positive_index <- which(t_pred$pred == 0 & seq_along(t_pred$pred) > first_positive_index)[1]
          
          out_temp3 = out_temp[1,c(1:length(colnames(out_temp)))]
          out_temp3[1,c(1:6)] = NA
          out_temp3$coefficients_mod2 <- "Tc"
          out_temp3$Estimate = t_pred$temp[first_zero_after_positive_index]
          
          out_temp = rbind(out_temp, out_temp2, out_temp3)
          out_temp$max_GR = t_pred$pred[which.max(t_pred$pred)]
          
          Out<-rbind(Out, out_temp)
          Predictions<-rbind(Predictions, t_pred)
        }
        
        # ## save individual plots    
        # i_name = gsub("%", "", i)
        #         
        # pred_plot = subset(Predictions, serial_number == s & treat == t & g == i)
        # 
        # (plot <- ggplot() +
        #   geom_line(data = pred_plot,aes(x = temp, y=pred, colour = mod, linetype = mod)) +
        #   geom_point(data = dat, aes(x = temp, y = Median), alpha = 0.6, size = 2)+
        #   scale_x_continuous(name = "Temperature (\u00B0 C)") +
        #   scale_y_continuous(name=expression(paste("Germination rate (", d^-1, ")", sep = "")),limits=c(0,1))+
        #   theme_bw() +
        #   ggtitle(paste(s, t, i, sep = " ")))
        # 
        # # Save the plot
        # ggsave(paste("Thresh_mod_", s, t, i_name, ".png", sep = "_"), plot, width = 2000, height = 1200, units = "px")
    
        rm(dat, models)
    }
  }
}

error_list
warning_list

Out |>
  dplyr::summarise(n = dplyr::n(), .by = c(serial_number, treat, g, mod, max_GR, coefficients_mod2)) |>
  dplyr::filter(n > 1L)

## Keeping just one Tc value (loop calculates is twice for exponential model)
Out  = Out %>%
  dplyr::filter(
    # mod == "modGR_BS" & coefficients_mod2 != "To" |
    # mod == "modGR_BS" & coefficients_mod2 == "To" & !is.na(factor) |
    mod == "modGR_Ex" & coefficients_mod2 != "Tc"| 
    mod == "modGR_Ex" & coefficients_mod2 == "Tc" & !is.na(factor)
    ) 

Out_w  = Out %>%
  dplyr::select(-factor) %>%
  group_by(serial_number, treat,	g, mod) %>%
  pivot_wider(names_from = coefficients_mod2,
              values_from = c(Estimate, `Std. Error`, `t-value`, `p-value`))

write.csv(Out_w, "output/csv/Thresholds_models_Alnus_glutinosa.csv", row.names = FALSE)
write.csv(Out, "output/csv/results_models_Alnus_glutinosa.csv", row.names = FALSE)
write.csv(Predictions, "output/csv/predictions_Alnus_glutinosa.csv", row.names = FALSE)



# plot -------------------------------------------------------------------------

df_To = subset(Out, coefficients_mod2 == "To")

GRlist_group_med = GRlist_group %>%
  group_by_at(vars(temp, g, treat, serial_number)) %>%
  summarise(
    sd = sd(Median),
    error_top = quantile(Median, probs = 0.75, na.rm = T),
    error_bot = quantile(Median, probs = 0.25, na.rm = T),
    Median = median(Median)
  )

## load serial_number information
info_df <- read.csv("data/csv/serial_number_info.csv")

GRlist_group_med = merge(GRlist_group_med, info_df, all.x = T, all.y =F)

# Arranging the dataframe by lat_dd in ascending order
GRlist_group_med <- GRlist_group_med %>% arrange(lat_dd)
GRlist_group_med$serial_number = as.factor(GRlist_group_med$serial_number)
levels(GRlist_group_med$serial_number)
GRlist_group_med$serial_number <- factor(GRlist_group_med$serial_number, levels = unique(GRlist_group_med$serial_number))

df_To$treat <- factor(df_To$treat, levels = c("NS", "COLD"))
Predictions$treat <- factor(Predictions$treat, levels = c("NS", "COLD"))
GRlist_group_med$treat <- factor(GRlist_group_med$treat, levels = c("NS", "COLD"))

Predictions = merge(Predictions, info_df, all.x = T, all.y =F)

# Arranging the dataframe by lat_dd in ascending order
Predictions <- Predictions %>% arrange(lat_dd)
Predictions$serial_number = as.factor(Predictions$serial_number)
levels(Predictions$serial_number)
Predictions$serial_number <- factor(Predictions$serial_number, levels = unique(Predictions$serial_number))

loc_labels <- Predictions %>%
  group_by(serial_number) %>%
  summarise(Geographical.Location = unique(Geographical.Location)) %>%
  mutate(Geographical.Location = sub(" - ", "\n", Geographical.Location)) %>%
  distinct(serial_number, Geographical.Location) %>%
  deframe()

#select_mod = "modGR_BS"
select_mod = "modGR_Ex"
Predictions1 = subset(Predictions, mod == select_mod)

df_To1 = subset(df_To, mod == select_mod & g == "50%")
df_To1 = merge(df_To1, info_df, all.x = T, all.y =F)

# Arranging the dataframe by lat_dd in ascending order
df_To1 <- df_To1 %>% arrange(lat_dd)
df_To1$serial_number = as.factor(df_To1$serial_number)
levels(df_To1$serial_number)
df_To1$serial_number <- factor(df_To1$serial_number, levels = unique(df_To1$serial_number))

GRlist_group_med1 = subset(GRlist_group_med, g == "50%")

(plot <- ggplot() +
    facet_grid(. ~ serial_number, labeller = labeller(serial_number = loc_labels)) +
    geom_vline(data = df_To1, aes(xintercept=Estimate, colour=treat), linetype="dashed", linewidth=0.3)+
    geom_line(data = subset(Predictions1, Predictions1$g == "50%"), aes(x = temp, y=pred, colour=treat, group = treat), linewidth = 0.5) + #linetype = treat
    #geom_line(data = Predictions, aes(x = temp, y=pred, colour=mod, group = mod, linetype = mod), linewidth = 0.9) +
    geom_errorbar(data = GRlist_group_med1, aes(x = temp, ymin = Median-sd, ymax = Median+sd, colour = treat), width = 1.7, linewidth = 0.3, alpha = 0.6) +
    geom_point(data = GRlist_group_med1, aes(x = temp, y = Median, colour = treat, fill = treat), size = 2.3, alpha = 0.5, shape = 21)+ 
    scale_color_manual(name = "Treatment", values = c("NS" = "deepskyblue3", "COLD" = "navy")) +
    scale_fill_manual(name = "Treatment", values = c("NS" = "deepskyblue3", "COLD" = "navy")) +
    # scale_linetype_discrete(name = "Treatment") +
    scale_shape_discrete(name = "Treatment") +
    coord_cartesian(ylim = c(-0.03, 0.75)) +
    scale_x_continuous(name = "Temperature (\u00B0 C)") +
    scale_y_continuous(name=expression(paste("Germination rate (", d^-1, ")", sep = ""))) + #,limits=c(0,max_pred))+ 
    theme_bw() +
    theme(legend.position = "top",
          strip.background.x = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white")))


save(plot, file = "output/rda/plot_alnus.rda")
