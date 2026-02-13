# ==============================================================================
# PROJECT:   Global Amenable Tuberculosis Mortality & Economic Burden (GBD 2021)
# SCRIPT:    visualization.R
# AUTHOR:    Yueyang Liu
# LICENSE:   MIT License
#
# DESCRIPTION:
#   This script generates the visualization figures for the manuscript.
#   
#   FIGURE LIST:
#   1. Mortality Maps : Global distribution of:
#      - Total TB Mortality Rate
#      - Amenable TB Mortality Rate
#      - Amenable Mortality Percentage (Amenable/Total)
#      (Stratified by Total TB, DS-TB, MDR-TB, XDR-TB)
#
#   2. Economic Maps : Global distribution of:
#      - Absolute Value of Lost Welfare (VLW)
#      - VLW as a percentage of GDP (VLW/GDP ratio)
#
#   3. Regional Analysis : Stacked bar charts showing VLW/GDP ratios 
#      stratified by age groups across 21 GBD regions.
#
#   4. Correlation Analysis: Scatter plots with segmented regression 
#      analyzing the relationship between:
#      - the percentage of TB mortality is amenable vs. SDI or Health Expenditure
#      - VLW/GDP vs. SDI or Health Expenditure
#
# REQUIREMENTS:
#   - Shapefile: 'data/shapefiles/word-map-with-correct-china-serbia.shp'
#   - Results data: 
#     - 'results/result_am_all_name1.xlsx' (Mortality data)
#     - 'results/vlw_country.xlsx' (Economic data)
# ==============================================================================

# ==============================================================================
# 0. SETUP & LIBRARIES
# ==============================================================================

rm(list = ls())

# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, sf, openxlsx, ggplot2, ggsci, patchwork, ggpubr, segmented, ggtext)

# Set Global Theme
theme_set(theme_classic() + 
            theme(
              plot.title = element_text(face = "bold", size = 12),
              plot.subtitle = element_text(hjust = 0.5, size = 11),
              legend.position = c(0.1, 0.2),
              legend.background = element_rect(fill = NA, color = NA),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank()
            ))

# Output directory
if(!dir.exists("figures")) dir.create("figures")

# ==============================================================================
# SECTION 1: HELPER FUNCTIONS
# ==============================================================================

# --- Function A: Data Imputation for China Territories ---
# Applies Mainland China estimates to Taiwan, Hong Kong, and Macao
impute_china_regions <- function(df, target_cols) {
  china_row <- df %>% filter(worldname == "China") %>% slice(1)
  
  if(nrow(china_row) > 0) {
    regions <- c("Taiwan", "Hong Kong", "Macao")
    for (region in regions) {
      if (region %in% df$worldname) {
        for (col in target_cols) {
          # Only fill if NA or explicitly required
          df[df$worldname == region, col] <- china_row[[col]]
        }
      }
    }
  }
  return(df)
}

# --- Function B: Categorizer ---
make_categories <- function(data_vec, breaks) {
  cut(data_vec, breaks = breaks, include.lowest = TRUE, right = FALSE,
      labels = paste0("[", round(head(breaks, -1), 3), ", ", round(tail(breaks, -1), 3), ")"))
}

# --- Function C: Map Plotter Engine ---
plot_map_layer <- function(map_data, fill_col, title, legend_title, palette) {
  ggplot(map_data) +
    geom_sf(aes(fill = !!sym(fill_col)), color = "#65605E", size = 0.1) +
    scale_fill_manual(name = legend_title, values = palette, na.translate = TRUE, na.value = "grey90") +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    labs(title = title) +
    theme(legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7))
}

# Define Color Palettes (Based on original code)
pal_blue  <- c("#B7E4EA", "#acf1e4", "#bafbcd", "#e2ffae", "#f0dc8b", "#f7b87a", "#f3957b")
pal_purp  <- c("#f1eef6", "#d4b9da", "#c994c7", "#df65b0", "#e7298a", "#ce1256", "#91003f")

# ==============================================================================
# SECTION 2: MORTALITY MAPS (Rate & Percentage)
# ==============================================================================
cat("Generating Mortality Maps...\n")

# Load Shapefile
map_base <- st_read("data/shapefiles/word-map-with-correct-china-serbia.shp", quiet = TRUE) %>%
  filter(fid != 342, !is.na(NAME), NAME != "Antarctica") %>%
  rename(worldname = NAME)

# Wrapper for Mortality Plots
plot_mortality_maps <- function(sheet_name, breaks_rate, breaks_pct, title_prefix) {
  
  # 1. Prepare Data
  df <- read.xlsx("data/result_am_all_name1.xlsx", sheet = sheet_name)
  df_map <- map_base %>% left_join(df, by = "worldname")
  
  # 2. Impute & Clean
  df_map <- impute_china_regions(df_map, c("weighted_mortality_rate", "tol_rate", "p_n"))
  
  # 3. Categorize
  df_map$cat_tol <- make_categories(df_map$tol_rate, breaks_rate)
  df_map$cat_am  <- make_categories(df_map$weighted_mortality_rate, breaks_rate)
  df_map$cat_pct <- make_categories(df_map$p_n, breaks_pct)
  
  # 4. Plot
  p1 <- plot_map_layer(df_map, "cat_tol", paste("A. Total", title_prefix, "Mortality"), "Rate (per 10^5)", pal_blue)
  p2 <- plot_map_layer(df_map, "cat_am",  paste("B. Amenable", title_prefix, "Mortality"), "Rate (per 10^5)", pal_blue)
  p3 <- plot_map_layer(df_map, "cat_pct", paste("C. Amenable Percentage"), "Percentage (%)", pal_purp)
  
  return(p1 / p2 / p3)
}

# --- Example Usage (Uncomment to run) ---
# breaks_tb_rate <- c(0, 0.27, 0.71, 1.64, 3.63, 11.2, 24, 146)
# breaks_tb_pct  <- quantile(c(0, 81.23), probs = 0:7/7) 
# p_mortality <- plot_mortality_maps("297", breaks_tb_rate, breaks_tb_pct, "TB")
# ggsave("figures/Fig1_Mortality_Maps.png", p_mortality, width = 8, height = 12)


# ==============================================================================
# SECTION 3: ECONOMIC MAPS (VLW & VLW/GDP)
# ==============================================================================
cat("Generating Economic Maps (VLW)...\n")

# Wrapper for Economic Plots
plot_economic_maps <- function(sheet_name, breaks_vlw, breaks_ratio, title_prefix) {
  
  # 1. Prepare Data
  df <- read.xlsx("data/vlw_country.xlsx", sheet = sheet_name)
  df_map <- map_base %>% left_join(df, by = "worldname")
  
  # 2. Impute (Specific to Economic Data)
  df_map <- impute_china_regions(df_map, c("total_vlw", "GDP", "P")) # P is VLW/GDP ratio
  
  # 3. Categorize
  df_map$cat_vlw   <- make_categories(df_map$total_vlw, breaks_vlw)
  df_map$cat_ratio <- make_categories(df_map$P, breaks_ratio)
  
  # 4. Plot
  # Map A: Absolute VLW
  p1 <- plot_map_layer(df_map, "cat_vlw", 
                       paste("A.", title_prefix, "Value of Lost Welfare (VLW)"), 
                       "VLW (Millions)", pal_blue)
  
  # Map B: VLW / GDP Ratio
  p2 <- plot_map_layer(df_map, "cat_ratio", 
                       paste("B.", title_prefix, "VLW / GDP Ratio"), 
                       "Percentage (%)", pal_purp)
  
  return(p1 / p2)
}

# --- Example Usage (Uncomment to run) ---
# breaks_vlw_total <- c(0, 7.8, 36.3, 87.5, 316, 705, 2050, 199000)
# breaks_ratio_total <- quantile(c(0, 10.25), probs = 0:7/7)
# p_economic <- plot_economic_maps("297", breaks_vlw_total, breaks_ratio_total, "Total TB")
# ggsave("figures/Fig2_Economic_Maps.png", p_economic, width = 10, height = 10)


# ==============================================================================
# SECTION 4: REGIONAL BAR CHARTS (Age-Stratified VLW)
# ==============================================================================
cat("Generating Regional Bar Charts...\n")

plot_regional_bars <- function(data_path, sheet, title_text, y_max, bold_label = NULL) {
  df <- read.xlsx(data_path, sheet = sheet)
  df$age_id_a <- factor(df$age_id_a)
  
  # Sorting logic
  df$worldname <- factor(df$worldname, levels = unique(df$worldname[order(df$P, decreasing = FALSE)]))
  
  p <- ggplot(df, aes(x = worldname, y = P)) +
    geom_col(aes(fill = age_id_a), width = 0.7) +
    scale_fill_npg(name = "Age Group", labels = c("0-14", "15-49", "50-64", "65-74")) +
    coord_flip() +
    scale_y_continuous(expand = c(0,0), limits = c(0, y_max)) +
    labs(title = title_text, y = "VLW/GDP (%)", x = "") +
    theme(axis.line = element_line(color = "black"), 
          axis.text = element_text(color = "black", size = 8))
  
  if(!is.null(bold_label)) {
    p <- p + scale_x_discrete(labels = function(x) ifelse(x == bold_label, paste0("**", x, "**"), x)) +
      theme(axis.text.y = element_markdown())
  }
  return(p)
}

# ==============================================================================
# SECTION 5: SCATTER PLOTS (Segmented Regression)
# ==============================================================================
cat("Generating Scatter Plots...\n")

plot_segmented_regression <- function(sheet_name, x_col, y_col, psi, x_lab, y_lab) {
  df <- read.xlsx("data/result_am_all_name1.xlsx", sheet = sheet_name) %>%
    mutate(across(c(all_of(x_col), all_of(y_col)), as.numeric)) %>%
    filter(!is.na(!!sym(x_col)))
  
  # Fit Model
  fit_lm <- lm(as.formula(paste(y_col, "~", x_col)), data = df)
  fit_seg <- segmented(fit_lm, seg.Z = as.formula(paste("~", x_col)), psi = psi)
  bkpt <- fit_seg$psi[2]
  
  # Predict
  pred_df <- data.frame(x = df[[x_col]]); names(pred_df) <- x_col
  preds <- predict(fit_seg, newdata = pred_df, interval = "confidence")
  df_plot <- cbind(df, preds)
  
  ggplot(df_plot, aes_string(x = x_col, y = y_col)) +
    geom_point(color = "orange3", alpha = 0.5) +
    geom_line(aes(y = fit), color = "steelblue", size = 1) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "steelblue", alpha = 0.2) +
    geom_vline(xintercept = bkpt, linetype = "dashed") +
    annotate("text", x = bkpt, y = max(df[[y_col]]), label = paste("Breakpoint:", round(bkpt)), hjust = -0.2) +
    labs(x = x_lab, y = y_lab)
}

cat("Script loaded. Use plot functions to generate specific figures.\n")