# Load libraries
library(patchwork)
library(dplyr)
library(emmeans)
library(tidyr)
library(purrr)
library(caret)
library(factoextra)
library(ggrepel)
library(tidyverse)
library(AICcmodavg)
library(MuMIn)


options(scipen = 999)

df <- readRDS("output/rds/df_deltas.rds")


# CLIMATE NICHE ----------------------------------------------------------------

# # PCA ------------------------------------------------------------------------

# Select bioclim + temperature at time of dispersal
env_data <- df %>% 
  ungroup() %>%
  dplyr::select(starts_with("bio"), tas_mean) %>%
  rename(
    T_disp = tas_mean
  )

# Run PCA
pca_res <- prcomp(env_data, scale. = TRUE)

summary(pca_res) # Check how many PCs explain >80% of variance

df_sp <- cbind(df, pca_res$x[, 1:2]) # Keep only first 2 PCs


pal = viridis::viridis(3, end = 0.7, option = "D")

plot_pca1 <- fviz_pca_biplot(pca_res,
                # Individuals (Points)
                geom.ind = "point",
                col.ind = df_sp$sp,      # Color points by Species
                palette = pal,          
                addEllipses = TRUE,      # Add concentration ellipses
                legend.title = "Species",
                
                # Variables (Arrows/Loadings)
                col.var = "grey40",      # Color arrows black
                repel = TRUE,            # Avoid text overlapping
                alpha.var = 0.5,         
                
                # Labels
                title = ""
)


var_cor <- as.data.frame(get_pca_var(pca_res)$cor)[,c(1:2)]
var_cor$term <- row.names(var_cor)

bioclim <- read.csv("data/csv/bioclim_names.csv")
bioclim <- bioclim %>%
  mutate(term = str_remove(term, "CHELSA_")) %>%
  full_join(., var_cor, by = "term")
bioclim$term_name <- ifelse(bioclim$term == "T_disp", "mean temperature during the seed dispersal months (°C)", bioclim$term_name)
bioclim$term_name2 <- paste(bioclim$term_name, bioclim$term, sep = " - ") 


# calculate the Euclidean distance between variables based on their Dim.1 and Dim.2 scores
# To reorder variables
loadings_matrix <- bioclim %>% dplyr::select(Dim.1, Dim.2)
row.names(loadings_matrix) <- bioclim$term_name2

cluster_order <- hclust(dist(loadings_matrix))$order

# apply cluster order
bioclim$term_name2 <- factor(bioclim$term_name2, 
                            levels = bioclim$term_name2[cluster_order])

bioclim_long <- bioclim %>%
  # 1. Rename columns first (New = Old)
  rename(
    PC1 = Dim.1,
    PC2 = Dim.2
  ) %>%
  
  # 2. Pivot (names are already correct)
  pivot_longer(
    cols = c("PC1", "PC2"),
    names_to = "PC",
    values_to = "Loading"
  )

plot_pca2 <- 
    ggplot(bioclim_long, aes(x = PC, y = term_name2, fill = Loading)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      "Correlations",
      low = "#00AFBB", 
      mid = "white", 
      high = "#FC4E07", 
      midpoint = 0, 
      limit = c(-1, 1)
    ) +
    geom_text(aes(label = round(Loading, 2)), 
            size = 3) +
    scale_x_discrete(position = "top") + 
      labs(
      title = "",
      x = NULL, 
      y = NULL
    ) +
    theme_minimal() +
    theme(
      # panel.border = element_blank(),
      legend.position = "right",
      # legend.margin = margin(0, 0, 0, 0),
      # legend.justification = c(5, 0),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    ) #+
    #coord_fixed(ratio = 0.5)


plot_pca1 <- plot_pca1 + labs(tag = "a") +  theme(plot.tag.position = c(0.02, 0.94))
plot_pca2 <- plot_pca2 + labs(tag = "b") +  theme(plot.tag.position = c(0.02, 0.94))

layout <- c(
  patchwork::area(t = 0, l = 0, b = 5, r = 12),
  patchwork::area(t = 6, l = 0, b = 11, r = 12)
)

(plot_pca = wrap_plots(free(plot_pca1), plot_pca2) + 
     plot_layout(design = layout))

# # Figure for Supplementary Material
# ggsave("output/imgs/Fig_SI_PCA.png", plot_pca, width = 8, height = 10, units = "in")


# ------------------------------------------------------------------------------
# CLIMATIC MODEL CONFIGURATION

pred <- paste0("PC", 1:2)
predictors_sp <- paste0("PC", 1:2, collapse = "+")

# Define response variables
responses  <- c(
  "delta_Estimate_Tb",
  "delta_Estimate_To",
  "delta_Estimate_Tc",
  "delta_Estimate_ThetaT",
  "delta_Estimate_thermal_width"
)


# test below for reference only. Code continues below.

# # # -----------------------------------------------------------------------------
# # # TEST HOW THE MODEL PERFORMS WITH INTERACTIONS
# # # -----------------------------------------------------------------------------
# 
# # Manually fit the suspicious model (most clear relationship from previous tests) with interactions
# suspect_mod <- lm(delta_Estimate_To ~ PC2 * sp, data = df_sp)
# 
# summary(suspect_mod)
# 
# # -------------------------------------------------------
# # Visual Inspection
# # Does the interaction look true, or does a single point pull a line?
# # -------------------------------------------------------
# p1 <- ggplot(df_sp, aes(x = PC2, y = delta_Estimate_To, color = sp, fill = sp)) +
#   geom_point(size = 3, alpha = 0.8) +
#   geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
#   theme_minimal()
# 
# print(p1)
# # P. sylvestris clearly has one outlier pulling the line
# 
# # -------------------------------------------------------
# # The "Influence Check" (Cook's Distance)
# # Rule of Thumb: D > 4/N (0.22) is suspicious. D > 1 is FATAL.
# # -------------------------------------------------------
# # Calculate Cook's Distance
# cooksD <- cooks.distance(suspect_mod)
# threshold <- 4 / nrow(df_sp)
# 
# # Plot it
# plot(cooksD, pch = 19, main = "Cook's Distance per Observation", ylab = "Cook's D")
# abline(h = threshold, col = "red", lty = 2)
# text(x = 1:length(cooksD), y = cooksD, labels = ifelse(cooksD > threshold, names(cooksD), ""), pos = 3)
# 
# print(paste("Threshold for suspicion (4/N):", round(threshold, 3)))
# print(paste("Max Cook's Distance found:", round(max(cooksD), 3)))
# 
# # -------------------------------------------------------
# # The "Stability Test" (Jackknife)
# # We re-run the model 18 times, dropping 1 point each time.
# # If the Interaction P-value jumps > 0.05, the result is unstable.
# # -------------------------------------------------------
# p_values_interaction <- numeric(nrow(df_sp))
# 
# for(i in 1:nrow(df_sp)) {
# 
#   # Drop observation 'i'
#   sub_data <- df_sp[-i, ]
# 
#   # Refit model
#   mod_iter <- lm(delta_Estimate_To ~ PC2 * sp, data = sub_data)
# 
#   # Get P-value of the Interaction Term specifically
#   # We use drop1() to test the significance of the interaction term as a whole
#   # (This is safer than checking single coefficients for 3-level factors)
#   test_res <- drop1(mod_iter, test = "F")
# 
#   # Extract P-value for the interaction row (usually the 2nd row of results)
#   p_val <- test_res$`Pr(>F)`[grep(":", rownames(test_res))]
# 
#   p_values_interaction[i] <- p_val
# }
# 
# results_stability <- data.frame(
#   Dropped_Row = 1:nrow(df_sp),
#   Interaction_P_Value = round(p_values_interaction, 4),
#   Status = ifelse(p_values_interaction < 0.05, "Significant", "NOT SIGNIFICANT (Fragile)")
# )
# 
# print(results_stability)
# 
# if(any(p_values_interaction > 0.05)) {
#   print("WARNING: The interaction effect disappears when removing certain points. RESULT IS UNSTABLE.")
# } else {
#   print("PASS: Interaction remains significant regardless of which point is dropped.")
# }
# 
# # Identify which drops are non-significant (p >= 0.05)
# is_outlier <- p_values_interaction >= 0.05
# 
# # X positions are the dropped row indices
# x <- seq_along(p_values_interaction)
# y <- p_values_interaction
# 
# # Plot all points in grey
# plot(x, y, pch = 19, col = "black", cex = 1.1,
#      xlab = "Dropped row", ylab = "Interaction p-value",
#      main = "Leave-one-out sensitivity: interaction term")
# abline(h = 0.05, col = "red", lty = 2)
# 
# # Overplot outliers (non-significant) in red and larger
# points(x[is_outlier], y[is_outlier], pch = 19, col = "red", cex = 1.1)
# 
# # Label outliers with their row numbers (e.g., 11)
# text(x[is_outlier], y[is_outlier], labels = x[is_outlier],
#      pos = 3, cex = 0.9, col = "red")
# 
# # Optional: connect points with a line
# lines(x, y, col = "grey70")
# 
# 
# # # -----------------------------------------------------------------------------
# # # END OF TEST
# # # -----------------------------------------------------------------------------


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Plot Model Diagnostics
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

make_resid_plot <- function(mod) {
  df <- data.frame(
    fitted = fitted(mod),
    resid  = resid(mod)
  )
  ggplot(df, aes(fitted, resid)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    ggtitle(label = resp) +
    labs(x = "Fitted values", y = "Residuals") +
    theme_bw()
}


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Analysis Loop
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# comparing model without interactions as N=18 is too small for interactions (* sp)
# Results of interactions are unstable and driven by outliers - see tests above

# Storage Lists
mod_diag <- list()
plots_top_mod <- list()
table_mod_selection <- list()
table_res <- list()

for (resp in responses) {
  
  #---------------------------------------------------
  # Fit additive model without interactions
  # --------------------------------------------------
  form_full <- as.formula(paste(resp, "~ ",  predictors_sp, " + sp"))
  mod_full <- lm(form_full, data = df_sp, na.action = "na.fail")
  
  pdiag <- make_resid_plot(mod_full)
  
  # Run dredge for all subsets of the model
  model_set <- dredge(mod_full, fixed = "sp", evaluate = TRUE)
  
  # --------------------------------------
  # Get the models that meet the criteria
  # --------------------------------------
  
  # # Keep models with cumulative sum of AIC weights is 0.95
  cum_w <- cumsum(model_set$weight)
  idx <- cum_w <= 0.95
  idx[which.max(cum_w >= 0.95)] <- TRUE  # include the first model that pushes over
  top_models <- get.models(model_set, subset = idx)
  
  
  if (length(top_models) > 1) {
    # CASE A: Multiple models - Perform Model Averaging
    selected_mod <- model.avg(model_set, subset = idx, fit = TRUE)
    
    conf_int  <- confint(selected_mod, full = TRUE)
    sum_mod  <- summary(selected_mod)$coefmat.full
    
  } else {
    # CASE B: Only one model - Extract stats from that single model
    selected_mod  <- top_models[[1]]
    
    conf_int  <- confint(selected_mod)
    sum_mod <- summary(selected_mod)$coefficients
  }
  
  # # Save results
  
  sw(selected_mod)
  
  res_table <- as.data.frame(cbind(sum_mod, conf_int))
  res_table[] <- lapply(res_table, function(x) if(is.numeric(x)) round(x, 4) else x)
  colnames(res_table)[colnames(res_table) %in% c("z value", "t value")] <- "z or t value"
  colnames(res_table)[colnames(res_table) %in% c("Pr(>|z|)", "Pr(>|t|)")] <- "p value"
  
  model_df <- as.data.frame(model_set)
  model_df$Status <- ifelse(idx, "Top Model", "Rejected")
  
  # Add metadata
  res_table$Response <- resp
  res_table$Full_Model <- paste(form_full)[3]
  
  model_df$Response <- resp
  model_df$Full_Model <- paste(form_full)[3]
  
  # # Capture the current model_set and labels
  # plot_fn <- (function(ms, title_txt = NULL) {
  #   force(ms)
  #   force(title_txt)
  #   function() {
  #     plot(ms, labAsExpr = TRUE, main = NULL)
  #     if (!is.null(title_txt)) {
  #       mtext(
  #         title_txt,
  #         side = 3,
  #         line = 0.5,
  #         cex = 0.6,
  #         font = 2
  #       )
  #     }
  #   }
  # })(model_set, title_txt = paste(resp))
  
  # Store
  key <- paste(resp)
  table_mod_selection[[key]] <- model_df
  table_res[[key]] <- res_table
  mod_diag[[key]] <- pdiag
  # plots_top_mod[[key]] <- plot_fn
}

# # End of Loop ----------------------------------------------------------------

# # Reporting results

# # COMBINE DIAGNOSTIC PLOTS FOR VISUAL CHECK
p_row1 <- wrap_plots(mod_diag[1:3], nrow = 1, ncol = 3)    
p_row2 <- wrap_plots(mod_diag[4:5], nrow = 1, ncol = 2) 
p_row1 / p_row2


# Summary table with estimates
final_table_res <- bind_rows(table_res)
final_table_res$pred_label <- gsub("\\.\\.\\..*$", "", row.names(final_table_res))
row.names(final_table_res) <- NULL
final_table_res$pred_label <- gsub("^sp", "", final_table_res$pred_label)
final_table_res

# Summary table model selection
final_table_mod_selection <- bind_rows(table_mod_selection)
names(final_table_mod_selection)

# Desired column order
new_order <- c("Response", "(Intercept)", "sp", "PC1", "PC2", "df", 
               "logLik", "AICc", "delta", "weight", "Status", "Full_Model")

# Reorder columns
final_table_mod_selection <- final_table_mod_selection[, new_order[new_order %in% names(final_table_mod_selection)]]
final_table_mod_selection$mod_number <- gsub("\\.\\.\\..*$", "", row.names(final_table_mod_selection))

final_table_mod_selection


# # COMBINE ALL PLOTS:
# plot_list_on_page <- function(plot_fn_list, nrow = 3, ncol = 5,
#                               mar = c(3, 3, 2, 1)) {
#   op <- par(mfrow = c(nrow, ncol), mar = mar)
#   on.exit(par(op), add = TRUE)
#   for (fn in plot_fn_list) fn()
# }
# plot_list_on_page(plots_top_mod, nrow = 2, ncol = 3)
  

# # SAVE RESULTS 

write.csv(final_table_mod_selection, "output/tables/Table_Model_Selection_Climate.csv", row.names = FALSE)
write.csv(final_table_res, "output/tables/Table_Estimates_Climate.csv", row.names = FALSE)

# -------------------------------------------------------
# FIGURE main text - Effect size
# -------------------------------------------------------

# Filter out Intercepts to focus on effects
plot_data <- final_table_res %>%
  filter(pred_label %in% c("PC1", "PC2")) %>%
  mutate(Response = dplyr::recode(Response,
                           "delta_Estimate_Tb" = "Delta~T[b]",
                           "delta_Estimate_To" = "Delta~T[o]",
                           "delta_Estimate_Tc" = "Delta~T[c]",
                           "delta_Estimate_ThetaT" = "Delta~Theta[T]",
                           "delta_Estimate_thermal_width" = "atop(Delta~Thermal, 'Width')"
                            ),
         pred_label = case_when(
           pred_label == "P. sylvestris" ~ "italic('P. sylvestris')",
           pred_label == "B. pubescens"  ~ "italic('B. pubescens')",
           pred_label == "PC1" ~ "(a)~PC1",
           pred_label == "PC2" ~ "(b)~PC2",
           TRUE ~ pred_label
         )
  )

plot1_main <- ggplot(plot_data, aes(x = Estimate, y = pred_label)) +
  # Vertical line at 0 (No Effect)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray20") +
  # Error Bars (95% CI)
  geom_errorbar(aes(xmin = `2.5 %`, xmax = `97.5 %`), width = 0.5) +
  # Point Estimates
  geom_point(size = 2.5, color = "black") +
  # Facet by Response Variable to show all models at once
  facet_grid(Response ~ pred_label, scales = "free_x", labeller = label_parsed, switch   = "y") +
  scale_y_discrete(labels = scales::label_parse()) +
  labs(x = "Estimate (Effect Size)",
       y = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_blank(),
        strip.text.x = element_text(size = 12),
        # strip.text.x = element_text(size = 11, face = "italic"),
        strip.text.y.left = element_text(size = 11, angle = 0, lineheight = 0.55),
        strip.background = element_blank(),
        
        plot.title = element_text(
          hjust = 0.5,   # 0 = left, 0.5 = center, 1 = right
          size = 12,     # optional: adjust font size
          face = "bold"  # optional: make it bold
        ),
        plot.background = element_rect(fill = "white"),
        plot.margin = margin(t = 0, l = 0, b = 0, r = 10)
  ) 


plot1_main


# # Prepare data for plotting 


df_long <- df_sp %>%
  pivot_longer(
    cols = starts_with("delta_"),
    names_to = "delta_label",
    values_to = "resp_value"
  )


df_longer <- df_long %>%
  pivot_longer(
    cols = c(starts_with("PC"), starts_with("bio"), starts_with("tas"), starts_with("hurs")),
    names_to = "pred_label",
    values_to = "pred_value"
  )


df_longer <- df_longer %>%
  filter(str_starts(pred_label, "PC")) %>%
  rename(Response = delta_label,
         Predictor = pred_label)   


df_longer <- df_longer %>%
  mutate(Response = dplyr::recode(Response,
                              "delta_Estimate_Tb" = "Delta~T[b]",
                              "delta_Estimate_To" = "Delta~T[o]",
                              "delta_Estimate_Tc" = "Delta~T[c]",
                              "delta_Estimate_ThetaT" = "Delta~Theta[T]",
                              "delta_Estimate_thermal_width" = "Delta~'Thermal Width'"
                           )
  )

# Fit the additive model per facet (Response x Predictor), no interaction term
#    so slopes are shared across species; species only shifts intercepts.
models_by_facet <- df_longer %>%
  group_by(Response, Predictor) %>%
  tidyr::nest() %>%
  mutate(
    fit = map(data, ~ lm(resp_value ~ pred_value + sp, data = .x))
  )

# Build a prediction grid (sequence of x for each facet x species)
pred_grid <- models_by_facet %>%
  mutate(
    newdata = map2(
      data, fit,
      ~ {
        sp_levels <- unique(.x$sp)
        x_seq <- seq(min(.x$pred_value, na.rm = TRUE),
                     max(.x$pred_value, na.rm = TRUE),
                     length.out = 200)
        expand.grid(pred_value = x_seq, sp = sp_levels) %>%
          as_tibble()
      }
    ),
    preds = map2(fit, newdata, ~ mutate(.y, .fitted = predict(.x, newdata = .y)))
  ) %>%
  dplyr::select(Response, Predictor, preds) %>%
  tidyr::unnest(preds)


# # Plot

(plot1 <- 
    ggplot(df_longer, aes(x = pred_value, y = resp_value)) +
    facet_grid(Response ~ Predictor,
               scales = "free",
               labeller = label_parsed) +
    
    # Additive model lines: shared slope, species-specific intercept
    geom_line(
      data = pred_grid,
      aes(x = pred_value, y = .fitted, colour = sp),
      linewidth = 0.4
    ) +
    
    # Raw per-species lm fits: these show different slopes (interaction).
    geom_smooth(aes(colour = sp),
                method = "lm", se = FALSE, alpha = 0.15,
                show.legend = FALSE, linetype = "dashed", linewidth = 0.6) +
    
    geom_point(aes(colour = sp), size = 1.7, alpha = 0.7) +
    
    scale_alpha_identity(name = "", guide = "none") +
    scale_linetype_identity(name = "", guide = "none") +
    scale_color_viridis_d(name = "", option = "D", end = 0.7) +
    scale_x_continuous(name = "score") +
    scale_y_continuous(name = "value (°C)") +
    theme_bw(base_size = 11) +
    theme(
      legend.position    = "top",
      legend.text        = element_text(size = 9, face = "italic"),
      strip.text.y       = element_text(size = 9, face = "italic"),
      strip.background   = element_blank(),
      plot.background    = element_rect(fill = "white"),
      plot.margin        = margin(t = 0, l = 0, b = 0, r = 0)
    )
)

# ggsave("output/deltas_pca_per_sp.png", plot1, width = 6.5, height = 3, units = "in")


# ------------------------------------------------------------------------------
# SEED MASS --------------------------------------------------------------------
# ------------------------------------------------------------------------------

wgt = read.csv("data/csv/Seed_mass_ALL.csv")
names(wgt)

wgt_sum = wgt %>%
  dplyr::select(serial_number, mass_single_seed_mg) %>%
  group_by(serial_number) %>%
  summarise(seed_mass = median(mass_single_seed_mg))

df_sp_mass = df_sp %>%
  inner_join(., wgt_sum, by = "serial_number") %>%
  mutate(
    log_seed_mass = log10(seed_mass)
  )


# # # -----------------------------------------------------------------------------
# # # TEST HOW THE MODEL PERFORMS WITH INTERACTIONS
# # # -----------------------------------------------------------------------------
# 
# # Manually fit the suspicious model (most clear relationship from previous tests) with interactions
# suspect_mod <- lm(delta_Estimate_To ~ log_seed_mass * sp, data = df_sp_mass)
# 
# summary(suspect_mod)
# 
# # -------------------------------------------------------
# # Visual Inspection
# # Does the interaction look true, or does a single point pull a line?
# # -------------------------------------------------------
# p1 <- ggplot(df_sp_mass, aes(x = log_seed_mass, y = delta_Estimate_To, color = sp, fill = sp)) +
#   geom_point(size = 3, alpha = 0.8) +
#   geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
#   theme_minimal()
# 
# print(p1)
# 
# # -------------------------------------------------------
# # The "Influence Check" (Cook's Distance)
# # Rule of Thumb: D > 4/N (0.22) is suspicious. D > 1 is FATAL.
# # -------------------------------------------------------
# # Calculate Cook's Distance
# cooksD <- cooks.distance(suspect_mod)
# threshold <- 4 / nrow(df_sp_mass)
# 
# # Plot it
# plot(cooksD, pch = 19, main = "Cook's Distance per Observation", ylab = "Cook's D")
# abline(h = threshold, col = "red", lty = 2)
# text(x = 1:length(cooksD), y = cooksD, labels = ifelse(cooksD > threshold, names(cooksD), ""), pos = 3)
# 
# print(paste("Threshold for suspicion (4/N):", round(threshold, 3)))
# print(paste("Max Cook's Distance found:", round(max(cooksD), 3)))
# 
# # -------------------------------------------------------
# # The "Stability Test" (Jackknife)
# # We re-run the model 18 times, dropping 1 point each time.
# # If the Interaction P-value jumps > 0.05, the result is unstable.
# # -------------------------------------------------------
# p_values_interaction <- numeric(nrow(df_sp_mass))
# 
# for(i in 1:nrow(df_sp_mass)) {
#   
#   # Drop observation 'i'
#   sub_data <- df_sp_mass[-i, ]
#   
#   # Refit model
#   mod_iter <- lm(delta_Estimate_To ~ log_seed_mass * sp, data = sub_data)
#   
#   # Get P-value of the Interaction Term specifically
#   # We use drop1() to test the significance of the interaction term as a whole
#   # (This is safer than checking single coefficients for 3-level factors)
#   test_res <- drop1(mod_iter, test = "F")
#   
#   # Extract P-value for the interaction row (usually the 2nd row of results)
#   p_val <- test_res$`Pr(>F)`[grep(":", rownames(test_res))]
#   
#   p_values_interaction[i] <- p_val
# }
# 
# results_stability <- data.frame(
#   Dropped_Row = 1:nrow(df_sp_mass),
#   Interaction_P_Value = round(p_values_interaction, 4),
#   Status = ifelse(p_values_interaction < 0.05, "Significant", "NOT SIGNIFICANT (Fragile)")
# )
# 
# print(results_stability)
# 
# if(any(p_values_interaction > 0.05)) {
#   print("WARNING: The interaction effect disappears when removing certain points. RESULT IS UNSTABLE.")
# } else {
#   print("PASS: Interaction remains significant regardless of which point is dropped.")
# }
# 
# # Identify which drops are non-significant (p >= 0.05)
# is_outlier <- p_values_interaction >= 0.05
# 
# # X positions are the dropped row indices
# x <- seq_along(p_values_interaction)
# y <- p_values_interaction
# 
# # Plot all points in grey
# plot(x, y, pch = 19, col = "black", cex = 1.1,
#      xlab = "Dropped row", ylab = "Interaction p-value",
#      main = "Leave-one-out sensitivity: interaction term")
# abline(h = 0.05, col = "red", lty = 2)
# 
# # Overplot outliers (non-significant) in red and larger
# points(x[is_outlier], y[is_outlier], pch = 19, col = "red", cex = 1.1)
# 
# # Label outliers with their row numbers (e.g., 11)
# text(x[is_outlier], y[is_outlier], labels = x[is_outlier],
#      pos = 3, cex = 0.9, col = "red")
# 
# # Optional: connect points with a line
# lines(x, y, col = "grey70")
# 
# 
# # -----------------------------------------------------------------------------
# # END OF TEST
# # -----------------------------------------------------------------------------
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# CONFIGURATION
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
pred <- "log_seed_mass"
predictors_sp <- "log_seed_mass"

# Define response variables
responses  <- c(
  "delta_Estimate_Tb",
  "delta_Estimate_To",
  "delta_Estimate_Tc",
  "delta_Estimate_ThetaT",
  "delta_Estimate_thermal_width"
)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Analysis Loop
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# comparing model without interactions as N=15 is too small for interactions (* sp)
# Results of interactions are unstable and driven by outliers as for PC1 and PC2

# Storage Lists
mod_diag_mass <- list()
plots_top_mod_mass <- list()
table_mod_selection_mass <- list()
table_res_mass <- list()

for (resp in responses) {
  
  #---------------------------------------------------
  # Fit additive model without interactions
  # --------------------------------------------------
  form_full <- as.formula(paste(resp, "~ ",  predictors_sp, " + sp"))
  mod_full <- lm(form_full, data = df_sp_mass, na.action = "na.fail")
  
  pdiag <- make_resid_plot(mod_full)
  
  # Run dredge for all subsets of the model
  model_set <- dredge(mod_full, fixed = "sp", evaluate = TRUE)
  
  # --------------------------------------
  # Get the models that meet the criteria
  # --------------------------------------
  
  # # Keep models with cumulative sum of AIC weights is 0.95
  cum_w <- cumsum(model_set$weight)
  idx <- cum_w <= 0.95
  idx[which.max(cum_w >= 0.95)] <- TRUE  # include the first model that pushes over
  top_models <- get.models(model_set, subset = idx)

  
  if (length(top_models) > 1) {
    # CASE A: Multiple models - Perform Model Averaging
    selected_mod <- model.avg(model_set, subset = idx, fit = TRUE)
    
    conf_int  <- confint(selected_mod, full = TRUE)
    sum_mod  <- summary(selected_mod)$coefmat.full
  
    
  } else {
    # CASE B: Only one model - Extract stats from that single model
    selected_mod  <- top_models[[1]]
    
    conf_int  <- confint(selected_mod)
    sum_mod <- summary(selected_mod)$coefficients
  }
  
  # # Save results ------------------------------------------------------------
  
  res_table <- as.data.frame(cbind(sum_mod, conf_int))
  res_table[] <- lapply(res_table, function(x) if(is.numeric(x)) round(x, 4) else x)
  colnames(res_table)[colnames(res_table) %in% c("z value", "t value")] <- "z or t value"
  colnames(res_table)[colnames(res_table) %in% c("Pr(>|z|)", "Pr(>|t|)")] <- "p value"
  
  model_df <- as.data.frame(model_set)
  model_df$Status <- ifelse(idx, "Top Model", "Rejected")
  
  # Add metadata
  res_table$Response <- resp
  res_table$Full_Model <- paste(form_full)[3]
  
  model_df$Response <- resp
  model_df$Full_Model <- paste(form_full)[3]
  
  # # Capture the current model_set and labels
  # plot_fn <- (function(ms, title_txt = NULL) {
  #   force(ms)
  #   force(title_txt)
  #   function() {
  #     plot(ms, labAsExpr = TRUE, main = NULL)
  #     if (!is.null(title_txt)) {
  #       mtext(
  #         title_txt,
  #         side = 3,
  #         line = 0.5,
  #         cex = 0.6,
  #         font = 2
  #       )
  #     }
  #   }
  # })(model_set, title_txt = paste(resp))
  
  # Store
  key <- paste(resp)
  mod_diag_mass[[key]] <- pdiag
  table_mod_selection_mass[[key]] <- model_df
  table_res_mass[[key]] <- res_table
  # plots_top_mod_mass[[key]] <- plot_fn
}

# # End of Loop -----------------------------------------------------------------

# #–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# # Reporting results
# #–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# # COMBINE DIAGNOSTIC PLOTS FOR VISUAL CHECK:
p_row1 <- wrap_plots(mod_diag_mass[1:3], nrow = 1, ncol = 3)    
p_row2 <- wrap_plots(mod_diag_mass[4:5], nrow = 1, ncol = 2) 

p_row1 / p_row2


# Summary table with estimates
final_table_res_mass <- bind_rows(table_res_mass)
final_table_res_mass$pred_label <- gsub("\\.\\.\\..*$", "", row.names(final_table_res_mass))
row.names(final_table_res_mass) <- NULL
final_table_res_mass$pred_label <- gsub("^sp", "", final_table_res_mass$pred_label)
subset(final_table_res_mass, pred_label == "log_seed_mass")

# Summary table model selection
final_table_mod_selection_mass <- bind_rows(table_mod_selection_mass)
names(final_table_mod_selection_mass)

# Desired column order
new_order <- c("Response", "(Intercept)", "sp", "log_seed_mass", "df", 
               "logLik", "AICc", "delta", "weight", "Status", "Full_Model")

# Reorder columns
final_table_mod_selection_mass <- final_table_mod_selection_mass[, new_order[new_order %in% names(final_table_mod_selection_mass)]]
final_table_mod_selection_mass$mod_number <- gsub("\\.\\.\\..*$", "", row.names(final_table_mod_selection_mass))

final_table_mod_selection_mass


# # SAVE RESULTS ----------------------------------------------------------

write.csv(final_table_mod_selection_mass, "output/tables/Table_Model_Selection_Seed_Mass.csv", row.names = FALSE)
write.csv(final_table_res_mass, "output/tables/Table_Estimates_Seed_Mass.csv", row.names = FALSE)

# -------------------------------------------------------
# FIGURE main text - Effect size
# -------------------------------------------------------

# Filter out Intercepts to focus on effects
plot_data_mass <- final_table_res_mass %>%
  filter(pred_label == "log_seed_mass") %>%
  mutate(Response = dplyr::recode(Response,
                           "delta_Estimate_Tb" = "Delta~T[b]",
                           "delta_Estimate_To" = "Delta~T[o]",
                           "delta_Estimate_Tc" = "Delta~T[c]",
                           "delta_Estimate_ThetaT" = "Delta~Theta[T]",
                           "delta_Estimate_thermal_width" = "atop(Delta~Thermal, 'Width')"
  ),
  pred_label = case_when(
    pred_label == "P. sylvestris" ~ "italic('P. sylvestris')",
    pred_label == "B. pubescens"  ~ "italic('B. pubescens')",
    pred_label == "log_seed_mass" ~ paste0("(c)~log","[10]"," (Seed~Mass)"),
    TRUE ~ pred_label
  )
  )

plot2_main <- ggplot(plot_data_mass, aes(x = Estimate, y = pred_label)) +
  # Vertical line at 0 (No Effect)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray20") +
  # Error Bars (95% CI)
  geom_errorbar(aes(xmin = `2.5 %`, xmax = `97.5 %`), width = 0.5) +
  # Point Estimates
  geom_point(size = 2.5, color = "black") +
  # Facet by Response Variable to show all models at once
  facet_grid(Response ~ pred_label, scales = "free_x", labeller = label_parsed, switch   = "y") +
  # ggtitle(expression(log[10]~"(Seed"~Mass*")")) +
  labs(x = "Estimate (Effect Size)",
       y = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_blank(),
        strip.text.x = element_text(size = 12),
        strip.text.y.left = element_text(size = 11, angle = 0, lineheight = 0.55),
        strip.background = element_blank(),
        
        plot.title = element_text(
          hjust = 0.5,   # 0 = left, 0.5 = center, 1 = right
          size = 12,     # optional: adjust font size
          face = "bold"  # optional: make it bold
        ),
        plot.background = element_rect(fill = "white"),
        plot.margin = margin(t = 0, l = 0, b = 0, r = 0)
  ) 


plot2_main


# #–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# # Prepare data for plotting 
# #–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# head(final_table_res_mass)
# head(final_table_mod_selection_mass)

df_long_mass <- df_sp_mass %>%
  pivot_longer(
    cols = starts_with("delta_"),
    names_to = "delta_label",
    values_to = "resp_value"
  )


df_longer_mass <- df_long_mass %>%
  pivot_longer(
    cols = c("log_seed_mass", starts_with("PC"), starts_with("bio"), starts_with("tas"), starts_with("hurs")),
    names_to = "pred_label",
    values_to = "pred_value"
  )


df_longer_mass <- df_longer_mass %>%
  filter(str_starts(pred_label, "log_seed_mass")) %>%
  rename(Response = delta_label,
         Predictor = pred_label)   


df_longer_mass <- df_longer_mass %>%
  mutate(Response = dplyr::recode(Response,
                                  "delta_Estimate_Tb" = "Delta~T[b]",
                                  "delta_Estimate_To" = "Delta~T[o]",
                                  "delta_Estimate_Tc" = "Delta~T[c]",
                                  "delta_Estimate_ThetaT" = "Delta~Theta[T]",
                                  "delta_Estimate_thermal_width" = "Delta~'Thermal Width'"
         ),
         Predictor = dplyr::recode(Predictor,
                                   "log_seed_mass" = "'Seed mass'" 
         )
  )

# Fit the additive model per facet (Response x Predictor), no interaction term
#    so slopes are shared across species; species only shifts intercepts.
models_by_facet_mass <- df_longer_mass %>%
  group_by(Response, Predictor) %>%
  tidyr::nest() %>%
  mutate(
    fit = map(data, ~ lm(resp_value ~ pred_value + sp, data = .x))
  )

# Build a prediction grid (sequence of x for each facet x species)
pred_grid_mass <- models_by_facet_mass %>%
  mutate(
    newdata = map2(
      data, fit,
      ~ {
        sp_levels <- unique(.x$sp)
        x_seq <- seq(min(.x$pred_value, na.rm = TRUE),
                     max(.x$pred_value, na.rm = TRUE),
                     length.out = 200)
        expand.grid(pred_value = x_seq, sp = sp_levels) %>%
          as_tibble()
      }
    ),
    preds = map2(fit, newdata, ~ mutate(.y, .fitted = predict(.x, newdata = .y)))
  ) %>%
  dplyr::select(Response, Predictor, preds) %>%
  tidyr::unnest(preds)

# #–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# # Plot
# #–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

(plot2 <- 
    ggplot(df_longer_mass, aes(x = pred_value, y = resp_value)) +
    facet_grid(Response ~ Predictor,
               scales = "free",
               labeller = label_parsed) +
    
    # Additive model lines: shared slope, species-specific intercept
    geom_line(
      data = pred_grid_mass,
      aes(x = pred_value, y = .fitted, colour = sp),
      linewidth = 0.4
    ) +
    
    # Raw per-species lm fits: these show different slopes (interaction).
    geom_smooth(aes(colour = sp),
                method = "lm", se = FALSE, alpha = 0.15,
                show.legend = FALSE, linetype = "dashed", linewidth = 0.6) +
    
    geom_point(aes(colour = sp), size = 1.7, alpha = 0.7) +
    
    scale_alpha_identity(name = "", guide = "none") +
    scale_linetype_identity(name = "", guide = "none") +
    scale_color_viridis_d(name = "", option = "D", end = 0.7) +
    scale_x_continuous(name = "score") +
    scale_y_continuous(name = "value (°C)") +
    theme_bw(base_size = 11) +
    theme(
      legend.position    = "top",
      legend.text        = element_text(size = 9, face = "italic"),
      strip.text.y       = element_text(size = 9, face = "italic"),
      strip.background   = element_blank(),
      plot.background    = element_rect(fill = "white"),
      plot.margin        = margin(t = 0, l = 0, b = 0, r = 0)
    )
)

# ggsave("output/deltas_pca_per_sp.png", plot1, width = 6.5, height = 3, units = "in")




# #–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# # Combine Plots
# #–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

values_y_scale <- df_longer %>%
  ungroup() %>%
  group_by(Response) %>%
  summarise(
    min = floor(min(resp_value)*1.15),
    max = ceiling(max(resp_value)*1.15)
    )


# Create a dummy dataframe with the min/max limits for each Response
dummy_limits <- values_y_scale %>%
  pivot_longer(cols = c(min, max), values_to = "resp_value") %>%
  mutate(pred_value = 0) 

# Add dummy data to BOTH plots to standardise scale
(plot1 <- plot1 + geom_blank(data = dummy_limits) +
    theme(
      strip.text.y = element_blank()  # Remove the facet labels on the right side of plot2 (e.g., 'Thermal Width')
    )
)
(plot2 <- plot2 + geom_blank(data = dummy_limits) +
  theme(
    axis.title.y = element_blank(), # Remove the main 'value (°C)' title
    axis.text.y  = element_blank(), # Remove the numerical tick labels (0, 4, 8, etc.)
    axis.ticks.y = element_blank(), # Remove the tick marks
  )
)

(combined = plot1 + plot2 + 
    plot_layout(guides = "collect", axis_titles = "collect_y", widths = c(2, 1)) 
  & theme(legend.position = 'top')
)

ggsave("output/imgs/Fig_SI_delta_pca_and_mass_supplementary.png", combined, width = 5.8, height = 7, units = "in")

# PLOT MAIN TEXT

# Filter out Intercepts to focus on effects
combined_plot_data <- rbind(plot_data, plot_data_mass)
combined_plot_data$y_dummy <- "Dummy"


(plot_main <- ggplot(combined_plot_data, aes(x = Estimate, y = y_dummy)) +
  # Vertical line at 0 (No Effect)
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  # Error Bars (95% CI)
  geom_errorbar(aes(xmin = `2.5 %`, xmax = `97.5 %`), width = 0.4, colour = "grey20") +
  # Point Estimates
  geom_point(size = 3, colour = "grey20") +
  # Facet by Response Variable to show all models at once
  facet_grid(Response ~ pred_label, scales = "free", labeller = label_parsed, switch = "y", drop = TRUE) +
  labs(x = "Estimate (Effect Size)",
       y = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.text.y.left = element_text(size = 11, angle = 0, lineheight = 0.55),
        strip.background = element_blank(),
        plot.background = element_rect(fill = "white"),
        plot.margin = margin(t = 0, l = 0, b = 0, r = 0)
  ) 
) 


ggsave("output/imgs/Fig_5_delta_pca_and_mass.png", plot_main, width = 5.5, height = 3, units = "in")

# ggsave("output/imgs/production/Fig_5_delta_pca_and_mass.pdf", plot_main, width = 5.6, height = 3, units = "in")
# 
# ggsave("output/imgs/production/Fig_5_delta_pca_and_mass.tif", plot_main, width = 5.6, height = 3, units = "in", dpi = 600, compression = "lzw")
