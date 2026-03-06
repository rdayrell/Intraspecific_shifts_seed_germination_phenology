library(dplyr)
library(tidyr)
library(lubridate)

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


df = df %>% dplyr::select(-Genus)

env = read.csv("data/csv/environment_data.csv")
output = merge(df, env, by = "serial_number")

output$Date.Collected = as.Date(output$Date.Collected, format = "%d-%b-%y")

names(output)

table_thresh <- output %>%
  mutate(
    Species_Name = paste(Genus, Species, sep = " "),
    Tb_SE = paste0(round(Tb, 2), " (", signif(Std..Error_Tb, 2), ")"),
    Tc_SE = paste0(round(Tc, 2), " (", signif(Std..Error_Tc, 2), ")"),
    ThetaT_SE = paste0(round(ThetaT, 2), " (", signif(Std..Error_ThetaT, 2), ")"),
    k_SE = paste0(round(k, 3), " (", round(Std..Error_k, 3), ")"),
    max_GR = round(max_GR, 2)
  ) %>%
  dplyr::select(Species_Name, Geographical.Location, serial_number, treat, Tb_SE, Tc_SE, To, ThetaT_SE, k_SE, max_GR)

table_thresh_sorted <- table_thresh %>%
  arrange(desc(treat)) %>%
  arrange(Species_Name, serial_number)

print(table_thresh_sorted)

write.csv(table_thresh_sorted, "output/tables/Table_threshold_models_all_collections.csv", row.names = FALSE)
