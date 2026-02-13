# ==============================================================================
# PROJECT:   Global Amenable Tuberculosis Mortality & Economic Burden (GBD 2021)
# SCRIPT:    main_analysis.R
# AUTHOR:    Yueyang Liu
# LICENSE:   MIT License
#
# DESCRIPTION:
#   This script performs the core statistical analysis for estimating amenable 
#   TB mortality and the associated Value of Lost Welfare (VLW).
#
# METHODOLOGY SUMMARY:
#   1. Uncertainty Propagation: Two-stage hybrid Monte Carlo simulation.
#      - Stage 1: Gamma distribution for raw GBD mortality rates (Right-skewed).
#      - Stage 2: Triangular distribution for derived parameters (CFR, UM, AM) 
#        to strictly enforce interval bounds propagated from Stage 1.
#   2. Benchmark Identification: Regional minimum Case Fatality Rate (CFR).
#   3. Economic Valuation: 
#      - Country-specific VSL Peak adjusted by GDP per capita (Eq. 7).
#      - VSLY calculated using polynomial age-weighting and discounting.
#   4. Aggregation: Delta Method for variance synthesis at global/regional levels.
#
# DATA REQUIREMENTS:
#   Files must be in 'data/' folder:
#   - death1.csv, risk_deathall.csv, incidence1.csv (GBD Estimates)
#   - population1.csv, le.xlsx (Demographics)
#   - gdp_data.csv (contain 'location_id' and 'gdp_per_capita')
# ==============================================================================

# ==============================================================================
# 0. ENVIRONMENT SETUP
# ==============================================================================

rm(list = ls())

# Auto-install and load required packages
required_packages <- c("tidyverse", "readxl", "writexl", "triangle", "yll")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

cat("Step 0: Environment setup complete.\n")

# ==============================================================================
# SECTION 1: HELPER FUNCTIONS
# ==============================================================================

# --- Func 1: Fit Gamma Parameters (Method of Moments) ---
# Used for: Simulating raw mortality rates from GBD data
get_gamma_params <- function(mean, lower, upper){
  # Approximate SD from 95% Uncertainty Interval
  sd <- (upper - lower) / 3.92
  sd <- max(sd, 1e-12) # Prevent zero division
  mean <- max(mean, 1e-12)
  
  var <- sd^2
  scale <- var / mean
  shape <- mean^2 / var
  return(list(shape = shape, scale = scale))
}

# --- Func 2: Calculate VSLY (Value of Statistical Life Year) ---
# Logic: Polynomial inverted U-shape curve for age-specific valuation
calculate_VSLY <- function(VSL_peak, a, LE_i, LE_remaining_a, r = 0.03) {
  # Polynomial function f(a) for age weighting
  f_a <- function(a, LE_i) {
    ratio <- a / LE_i
    return(-0.36456 * ratio^4 - 5.65959 * ratio^3 + 4.3871 * ratio^2 + 1.1505 * ratio + 0.0008447)
  }
  
  f_val <- f_a(a, LE_i)
  exp_part <- exp(-r * LE_remaining_a)
  denom <- exp_part - 1
  if(abs(denom) < 1e-9) denom <- -1e-9 
  
  # Final calculation
  VSLY <- VSL_peak * f_val * (-r / denom)
  return(max(0, VSLY)) 
}

# --- Func 3: Discounted YLL Integral ---
calculate_YLL_integral <- function(a, LE, r = 0.03) {
  integrand <- function(x) exp(-r * (x - a))
  tryCatch({
    res <- integrate(integrand, lower = a, upper = LE)
    return(res$value)
  }, error = function(e) return(NA))
}

# ==============================================================================
# SECTION 2: DATA LOADING & PRE-PROCESSING
# ==============================================================================
cat("Step 1: Loading and processing raw data...\n")

# Load Datasets (Relative Paths)
# Note: Ensure these files are present in the 'data/' folder
if(!dir.exists("data")) stop("Error: 'data' folder not found. Please create it and add data files.")

df_death <- read.csv("data/death1.csv")
df_risk  <- read.csv("data/risk_deathall.csv")
df_inc   <- read.csv("data/incidence1.csv")
df_pop   <- read.csv("data/population1.csv")
df_le    <- read_excel("data/le.xlsx") 
# IMPORTANT: You need a file with GDP per capita for each location
# df_gdp <- read.csv("data/gdp_data.csv") 

# Helper to reshape GBD data (Long -> Wide)
pivot_gbd <- function(df, metric_prefix) {
  df %>%
    pivot_wider(
      names_from = metric,
      values_from = c(contains("val"), contains("upper"), contains("lower")),
      names_glue = "{metric}_{.value}"
    )
}

# Reshape and Merge
df_death_w <- pivot_gbd(df_death, "Rate")
df_risk_w  <- pivot_gbd(df_risk, "Rate")
df_inc_w   <- pivot_gbd(df_inc, "in")

data_merged <- df_death_w %>%
  left_join(df_risk_w, by = c("location_id", "sex_id", "age_id", "cause_id")) %>%
  left_join(df_inc_w, by = c("location_id", "sex_id", "age_id", "cause_id")) %>%
  select(-matches("name")) 

# Clean Data
data_merged <- data_merged %>%
  mutate(across(starts_with("Rate_r"), ~replace_na(., 0))) %>%
  mutate(across(starts_with("Rate_r"), ~pmax(0, .))) %>%
  filter(!(age_id %in% c(20, 21, 22, 27))) # Exclude aggregates/elderly

# ==============================================================================
# SECTION 3: HYBRID MONTE CARLO SIMULATION
# ==============================================================================
cat("Step 2: Running Monte Carlo Simulation (Gamma + Triangle)...\n")

set.seed(123)
n_sim <- 10000 

# --- Phase 1: Gamma Simulation for Remaining Mortality (RM) ---
# Logic: RM = Total Mortality (Gamma) - Preventable Mortality (Gamma)

vec_rm_mean <- numeric(nrow(data_merged))
vec_rm_low  <- numeric(nrow(data_merged))
vec_rm_high <- numeric(nrow(data_merged))
data_merged$rm_direct <- data_merged$Rate_tol_val - data_merged$Rate_r_val

for(i in 1:nrow(data_merged)){
  # Extract Params
  A_est <- data_merged$Rate_tol_val[i]; A_lo <- data_merged$Rate_tol_lower[i]; A_hi <- data_merged$Rate_tol_upper[i]
  B_est <- data_merged$Rate_r_val[i];   B_lo <- data_merged$Rate_r_lower[i];   B_hi <- data_merged$Rate_r_upper[i]
  
  if (any(is.na(c(A_est, B_est)))) next
  
  # Gamma Sampling
  if (A_hi <= A_lo) A_sim <- rep(A_est, n_sim) else {
    pA <- get_gamma_params(A_est, A_lo, A_hi)
    A_sim <- rgamma(n_sim, shape = pA$shape, scale = pA$scale)
  }
  if (B_hi <= B_lo) B_sim <- rep(B_est, n_sim) else {
    pB <- get_gamma_params(B_est, B_lo, B_hi)
    B_sim <- rgamma(n_sim, shape = pB$shape, scale = pB$scale)
  }
  
  # Constraint & Calculation
  B_sim <- pmin(B_sim, A_sim)
  C_sim <- A_sim - B_sim 
  
  vec_rm_mean[i] <- mean(C_sim)
  qs <- quantile(C_sim, c(0.025, 0.975), na.rm=TRUE)
  vec_rm_low[i]  <- qs[1]; vec_rm_high[i] <- qs[2]
}

data_merged$RM_mean  <- vec_rm_mean
data_merged$RM_lower <- vec_rm_low
data_merged$RM_upper <- vec_rm_high

# --- Phase 2: Identify Benchmark CFR (with Uncertainty) ---
# Logic: Find region-specific minimum CFR and preserve its uncertainty bounds

data_merged$CFR_est <- data_merged$rm_direct / data_merged$Rate_in_val
data_merged$CFR_est <- pmin(pmax(data_merged$CFR_est, 0), 1)

# Propagate bounds for CFR (Approximate for benchmark selection context)
data_merged$CFR_lower <- data_merged$RM_lower / data_merged$Rate_in_val
data_merged$CFR_upper <- data_merged$RM_upper / data_merged$Rate_in_val

data_merged <- data_merged %>%
  mutate(level_m = case_when(
    location_id == 1 ~ 0, # Global
    # NOTE: Insert your full location mapping logic here
    TRUE ~ 3
  ))

# Identify Benchmark (Min CFR per region)
benchmark_cfr <- data_merged %>%
  filter(level_m != 0) %>%
  group_by(age_id, sex_id, cause_id) %>%
  slice(which.min(CFR_est)) %>%
  summarise(
    min_CFR_val = CFR_est,
    min_CFR_low = CFR_lower,
    min_CFR_high = CFR_upper,
    .groups = "drop"
  )

data_final <- data_merged %>%
  left_join(benchmark_cfr, by = c("age_id", "sex_id", "cause_id")) %>%
  mutate(UM_direct = Rate_in_val * min_CFR_val)

# --- Phase 3: Triangular Simulation for Amenable Mortality (AM) ---
# Logic: Unavoidable (UM) = min_CFR (Triangle) * Incidence (Triangle)
#        Amenable (AM)    = Remaining (Triangle) - UM

vec_am_mean <- numeric(nrow(data_final))
vec_am_low  <- numeric(nrow(data_final))
vec_am_high <- numeric(nrow(data_final))

for(i in 1:nrow(data_final)){
  # 1. Prepare RM (Remaining)
  rm_mode <- data_final$rm_direct[i]; rm_min <- data_final$RM_lower[i]; rm_max <- data_final$RM_upper[i]
  
  # 2. Prepare Benchmark CFR
  cfr_mode <- data_final$min_CFR_val[i]; cfr_min <- data_final$min_CFR_low[i]; cfr_max <- data_final$min_CFR_high[i]
  
  # 3. Prepare Incidence
  inc_mode <- data_final$Rate_in_val[i]; inc_min <- data_final$Rate_in_lower[i]; inc_max <- data_final$Rate_in_upper[i]
  
  if (any(is.na(c(rm_mode, cfr_mode, inc_mode)))) next
  
  # Simulate Components (Triangle)
  rm_sim  <- rtriangle(n_sim, a=max(0, rm_min), b=max(rm_min, rm_max), c=rm_mode)
  cfr_sim <- rtriangle(n_sim, a=max(0, cfr_min), b=max(cfr_min, cfr_max), c=cfr_mode)
  inc_sim <- rtriangle(n_sim, a=max(0, inc_min), b=max(inc_min, inc_max), c=inc_mode)
  
  # Calculate Unavoidable (UM)
  um_sim <- cfr_sim * inc_sim
  um_sim <- pmin(um_sim, rm_sim) # Constraint
  
  # Calculate Amenable (AM)
  am_sim <- rm_sim - um_sim
  
  vec_am_mean[i] <- mean(am_sim)
  qs <- quantile(am_sim, c(0.025, 0.975), na.rm=TRUE)
  vec_am_low[i]  <- qs[1]; vec_am_high[i] <- qs[2]
}

data_final$AM_val <- vec_am_mean; data_final$AM_lower <- vec_am_low; data_final$AM_upper <- vec_am_high

# ==============================================================================
# SECTION 4: ECONOMIC EVALUATION (VLW)
# ==============================================================================
cat("Step 3: Calculating VLW (Adjusted by GDP per capita)...\n")

data_econ <- data_final %>%
  left_join(df_pop %>% select(location_id, age_id, sex_id, pop_val, pop_lower, pop_upper), 
            by = c("location_id", "sex_id", "age_id")) %>%
  left_join(df_le, by = c("location_id", "sex_id", "age_id"))

# Constants
r_rate <- 0.03
VSL_US_2021 <- 11.8 * 1e6
GDP_US_CAPITA_2021 <- 70248 # USD (World Bank)

# --- VSL Peak Calculation (Eq 7) ---
# Note: You should merge your actual GDP per capita data here.
# data_econ <- data_econ %>% left_join(df_gdp_per_capita, by="location_id")

# Placeholder: If gdp_per_capita is missing, warn user
if(!"gdp_per_capita" %in% names(data_econ)) {
  cat("WARNING: 'gdp_per_capita' column not found. Using placeholder.\n")
  data_econ$gdp_per_capita <- 10000 
}

data_econ <- data_econ %>%
  mutate(
    # Income Elasticity (IE): Baseline = 1.0 (As per methods)
    # Sensitivity analysis can use: ifelse(gdp_per_capita >= GDP_US_CAPITA_2021, 0.8, 1.5)
    ie = 1.0, 
    
    # Calculate Country-Specific Peak VSL using GDP per capita
    vsl_peak_country = VSL_US_2021 * (gdp_per_capita / GDP_US_CAPITA_2021)^ie
  )

# Calculate VSLY & Final VLW
data_econ <- data_econ %>%
  mutate(
    age_start = case_when(age_id == 6 ~ 5, age_id == 7 ~ 10, TRUE ~ 0), 
    YLL_per_death = mapply(calculate_YLL_integral, a=age_start, LE=LE, MoreArgs=list(r=r_rate)),
    
    # Use calculated country peak VSL here
    VSLY = mapply(calculate_VSLY, VSL_peak=vsl_peak_country, a=age_start, LE_i=LE, LE_remaining_a=LE, MoreArgs=list(r=r_rate)),
    
    AM_Deaths_Num = (AM_val * pop_val) / 100000,
    Total_VLW = AM_Deaths_Num * YLL_per_death * VSLY
  )

# ==============================================================================
# SECTION 5: AGGREGATION (DELTA METHOD) & EXPORT
# ==============================================================================
cat("Step 4: Aggregating variances and exporting results...\n")

# Variance calculation for Delta Method
data_econ <- data_econ %>%
  mutate(
    am_var  = ((AM_upper - AM_lower) / 3.92)^2,
    pop_var = ((pop_upper - pop_lower) / 3.92)^2,
    
    # Delta Method: Var(Rate * Pop) ≈ Pop^2*Var(Rate) + Rate^2*Var(Pop)
    deaths_count_var = (pop_val/100000)^2 * am_var + (AM_val/100000)^2 * pop_var
  )

# Global Aggregation
agg_global <- data_econ %>%
  group_by(sex_id, cause_id) %>%
  summarise(
    Sum_Deaths = sum(AM_Deaths_Num, na.rm = TRUE),
    Sum_VLW    = sum(Total_VLW, na.rm = TRUE),
    Sum_Pop    = sum(pop_val, na.rm = TRUE),
    Sum_Deaths_Var = sum(deaths_count_var, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Deaths_SE    = sqrt(Sum_Deaths_Var),
    Deaths_Lower = pmax(0, Sum_Deaths - 1.96 * Deaths_SE),
    Deaths_Upper = Sum_Deaths + 1.96 * Deaths_SE
  )

# Export
if(!dir.exists("results")) dir.create("results")
write_xlsx(agg_global, "results/Global_Summary_2021.xlsx")

cat("Analysis Complete. Results saved in 'results/'.\n")