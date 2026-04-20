################################################################################
#       HETEROGENEOUS TRIGGERING BY FATALITY STATUS
#
#       Tests whether protests with fatalities have higher triggering rates
#       than non-fatal protests using decomposed exposure terms
#
#       H0: α_fatal = α_nonfatal (equal triggering effects)
#       H1: α_fatal ≠ α_nonfatal (differential triggering effects)
#
#       Uses optimal parameters from main analysis: θ = 200 km, κ = 0.85
################################################################################

library(tidyverse)
library(readxl)
library(lubridate)
library(sf)
library(geosphere)
library(MASS)
library(mgcv)

select <- dplyr::select

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   HETEROGENEOUS TRIGGERING BY FATALITY STATUS                           ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# OPTIMAL PARAMETERS FROM MAIN ANALYSIS
# =============================================================================

THETA_OPTIMAL <- 200    # Spatial decay distance (km)
KAPPA_OPTIMAL <- 0.85   # Temporal decay parameter
MAX_LAG <- 30           # Temporal window (days)

cat(sprintf("Using optimal parameters from main analysis:\n"))
cat(sprintf("  θ (spatial decay) = %d km\n", THETA_OPTIMAL))
cat(sprintf("  κ (temporal decay) = %.2f\n", KAPPA_OPTIMAL))
cat(sprintf("  Temporal window = %d days\n\n", MAX_LAG))

# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================

cat("=== LOADING DATA ===\n\n")

# Load ACLED data
acled_raw <- read_excel("ACLED Data_2025-10-27_Indonesia20142024_OFFICIAL.xlsx")

col_names <- as.character(acled_raw[1, ])
acled <- acled_raw[-1, ]
names(acled) <- col_names

acled <- acled %>%
  mutate(
    event_date = if (is.numeric(event_date)) {
      as.Date(event_date, origin = "1899-12-30")
    } else {
      lubridate::ymd(event_date)
    },
    year = as.integer(year),
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude),
    fatalities = as.integer(fatalities)
  )

# Filter to protest events
protests <- acled %>%
  filter(event_type %in% c("Protests", "Riots")) %>%
  mutate(
    is_violent = sub_event_type %in% c("Violent demonstration",
                                       "Protest with intervention",
                                       "Mob violence",
                                       "Looting/property destruction",
                                       "Excessive force against protesters"),
    is_peaceful = sub_event_type == "Peaceful protest",
    has_fatalities = fatalities > 0,
    is_fatal = fatalities > 0  # Binary indicator for stratification
  ) %>%
  arrange(event_date)

cat(sprintf("Total protest events: %d\n", nrow(protests)))
cat(sprintf("  Fatal protests (≥1 death): %d (%.1f%%)\n",
            sum(protests$is_fatal), 100 * mean(protests$is_fatal)))
cat(sprintf("  Non-fatal protests: %d (%.1f%%)\n",
            sum(!protests$is_fatal), 100 * mean(!protests$is_fatal)))
cat(sprintf("  Total fatalities: %d\n\n", sum(protests$fatalities)))

# =============================================================================
# LOAD GADM SHAPEFILE AND CREATE DISTRICT COORDINATES
# =============================================================================

cat("=== LOADING SHAPEFILE ===\n\n")

gadm <- st_read("Indonesian district-level shapefile copy 2/gadm41_IDN_2.shp",
                quiet = TRUE)

# Create district covariates from GADM - use admin1+admin2 for unique identification
# but keep admin2 as the panel key (matches original analysis)
gadm_to_pop_name <- function(name_2, engtype_2) {
  if (is.na(engtype_2) || is.na(name_2)) return(NA_character_)
  if (engtype_2 == "City") {
    base_name <- sub("^Kota\\s+", "", name_2)
    return(paste0(base_name, ", Kota"))
  } else {
    return(paste0(name_2, ", Kab."))
  }
}

# Create unique district identifier (admin1_admin2) to handle duplicate names
# Then use this as the panel key instead of just admin2
district_covariates <- gadm %>%
  mutate(
    centroid = st_point_on_surface(geometry),
    lon = st_coordinates(centroid)[, 1],
    lat = st_coordinates(centroid)[, 2]
  ) %>%
  st_drop_geometry() %>%
  transmute(
    admin2 = NAME_2,
    admin1 = NAME_1,
    # Create unique ID combining province and district
    admin2_unique = paste(NAME_1, NAME_2, sep = "_"),
    district = mapply(gadm_to_pop_name, NAME_2, ENGTYPE_2),
    lon = lon,
    lat = lat
  )

cat(sprintf("Districts from GADM: %d\n", nrow(district_covariates)))
cat(sprintf("Unique admin2 names: %d\n", n_distinct(district_covariates$admin2)))
cat(sprintf("Unique admin2_unique: %d\n", n_distinct(district_covariates$admin2_unique)))

# =============================================================================
# LOAD POPULATION DATA
# =============================================================================

pop_raw <- read.csv("district_level_population_2014_2020.csv",
                    stringsAsFactors = FALSE)

pop_data <- pop_raw %>%
  filter(Series.Name == "Total Population (in number of people)") %>%
  select(district = Provinces.Name,
         `2014` = X2014..YR2014.,
         `2015` = X2015..YR2015.,
         `2016` = X2016..YR2016.,
         `2017` = X2017..YR2017.,
         `2018` = X2018..YR2018.,
         `2019` = X2019..YR2019.,
         `2020` = X2020..YR2020.)

pop_data$district <- gsub('^"', '', pop_data$district)

# Reshape to long format
pop_long <- pop_data %>%
  pivot_longer(cols = starts_with("20"),
               names_to = "year",
               values_to = "population") %>%
  mutate(year = as.integer(year),
         population = as.numeric(population)) %>%
  filter(!is.na(population))

# Extrapolate population to 2021-2024
districts_to_extrapolate <- unique(pop_long$district)

extrapolated_pop <- lapply(districts_to_extrapolate, function(dist) {
  hist_data <- pop_long %>% filter(district == dist)

  if (nrow(hist_data) >= 3) {
    model <- lm(population ~ year, data = hist_data)
    future_years <- data.frame(year = 2021:2024)
    predictions <- predict(model, newdata = future_years)

    future_data <- data.frame(
      district = dist,
      year = 2021:2024,
      population = pmax(predictions, min(hist_data$population)),
      extrapolated = TRUE
    )

    hist_data$extrapolated <- FALSE
    return(rbind(hist_data, future_data))
  } else {
    hist_data$extrapolated <- FALSE
    return(hist_data)
  }
})

pop_complete <- bind_rows(extrapolated_pop) %>%
  mutate(log_pop = log(population))

# =============================================================================
# LOAD POVERTY DATA
# =============================================================================

poverty_raw <- read.csv("poverty_rate_2015_2024.csv", stringsAsFactors = FALSE)

poverty_clean <- poverty_raw %>%
  filter(!is.na(poverty_rate)) %>%
  mutate(year = as.integer(year))

# Extrapolate poverty to cover all years
poverty_trends <- poverty_clean %>%
  filter(year >= 2018, year <= 2020) %>%
  group_by(district) %>%
  summarise(trend = if (n() >= 2) coef(lm(poverty_rate ~ year))[2] else 0,
            base_2020 = poverty_rate[year == max(year)],
            .groups = "drop")

poverty_extended <- lapply(2021:2024, function(yr) {
  poverty_trends %>%
    mutate(
      year = yr,
      poverty_rate = pmax(base_2020 + trend * (yr - 2020), 0.1),
      extrapolated = TRUE
    ) %>%
    select(district, year, poverty_rate, extrapolated)
})

poverty_full <- bind_rows(
  poverty_clean %>% mutate(extrapolated = FALSE),
  bind_rows(poverty_extended)
)

# =============================================================================
# LOAD CPI DATA
# =============================================================================

cpi_data <- read.csv("indonesia_cpi_processed.csv", stringsAsFactors = FALSE)

cpi_monthly <- cpi_data %>%
  select(year_month, cpi) %>%
  mutate(year_month = as.character(year_month))

# =============================================================================
# JOIN COVARIATES TO DISTRICT DATA
# =============================================================================

# Join population (use 2019 as reference)
district_covariates <- district_covariates %>%
  left_join(
    pop_complete %>% filter(year == 2019) %>% select(district, log_pop),
    by = "district"
  ) %>%
  left_join(
    poverty_full %>% filter(year == 2019) %>% select(district, poverty_rate),
    by = "district"
  ) %>%
  filter(!is.na(lon), !is.na(lat), !is.na(log_pop), !is.na(poverty_rate))

cat(sprintf("Districts with complete covariates: %d\n\n", nrow(district_covariates)))

# =============================================================================
# SPATIAL MATCHING OF PROTESTS
# =============================================================================

cat("=== SPATIAL MATCHING ===\n\n")

# Spatial join protests to GADM
protests_with_coords <- protests %>%
  filter(!is.na(latitude) & !is.na(longitude))

protests_sf <- st_as_sf(protests_with_coords,
                        coords = c("longitude", "latitude"),
                        crs = 4326)

protests_matched <- st_join(protests_sf, gadm, join = st_intersects)

# Use unique district ID for matching
protests_final <- protests_matched %>%
  st_drop_geometry() %>%
  mutate(
    admin2 = NAME_2,
    admin1 = NAME_1,
    admin2_unique = paste(NAME_1, NAME_2, sep = "_")
  ) %>%
  filter(admin2_unique %in% district_covariates$admin2_unique) %>%
  select(event_date, admin2, admin1, admin2_unique, is_fatal, fatalities, is_violent, is_peaceful)

cat(sprintf("Protests matched to panel districts: %d\n", nrow(protests_final)))
cat(sprintf("  Fatal: %d\n", sum(protests_final$is_fatal)))
cat(sprintf("  Non-fatal: %d\n\n", sum(!protests_final$is_fatal)))

# =============================================================================
# CREATE DISTRICT-DAY PANEL WITH STRATIFIED COUNTS
# =============================================================================

cat("=== CREATING STRATIFIED PANEL ===\n\n")

# Use admin2_unique for panel construction to avoid duplicate name issues
all_districts <- unique(district_covariates$admin2_unique)
all_dates <- seq(min(protests_final$event_date),
                 max(protests_final$event_date), by = "day")

# Aggregate counts by district-day, stratified by fatality status
daily_counts_stratified <- protests_final %>%
  group_by(admin2_unique, event_date) %>%
  summarise(
    count = n(),
    count_fatal = sum(is_fatal, na.rm = TRUE),
    count_nonfatal = sum(!is_fatal, na.rm = TRUE),
    .groups = "drop"
  )

# Create balanced panel
panel <- expand.grid(
  admin2_unique = all_districts,
  date = all_dates,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  left_join(daily_counts_stratified %>%
              rename(date = event_date),
            by = c("admin2_unique", "date")) %>%
  mutate(
    count = replace_na(count, 0),
    count_fatal = replace_na(count_fatal, 0),
    count_nonfatal = replace_na(count_nonfatal, 0),
    year_month = format(date, "%Y-%m"),
    year = year(date)
  ) %>%
  left_join(district_covariates %>% select(admin2_unique, admin2, lon, lat, log_pop, poverty_rate),
            by = "admin2_unique") %>%
  left_join(cpi_monthly, by = "year_month") %>%
  arrange(admin2_unique, date)

cat(sprintf("Panel dimensions: %d district-days\n", nrow(panel)))
cat(sprintf("Districts: %d\n", length(unique(panel$admin2_unique))))
cat(sprintf("Date range: %s to %s\n\n", min(panel$date), max(panel$date)))

# Verification: stratified counts sum to total
total_count <- sum(panel$count)
total_fatal <- sum(panel$count_fatal)
total_nonfatal <- sum(panel$count_nonfatal)

cat("Verification of stratified counts:\n")
cat(sprintf("  Total count: %d\n", total_count))
cat(sprintf("  Fatal count: %d\n", total_fatal))
cat(sprintf("  Non-fatal count: %d\n", total_nonfatal))
cat(sprintf("  Sum of stratified: %d\n", total_fatal + total_nonfatal))
cat(sprintf("  Match: %s\n\n", ifelse(total_count == total_fatal + total_nonfatal,
                                       "YES ✓", "NO ✗")))

# =============================================================================
# COMPUTE DISTANCE MATRIX
# =============================================================================

cat("=== COMPUTING DISTANCE MATRIX ===\n\n")

coords <- as.matrix(district_covariates[, c("lon", "lat")])
dist_matrix <- distm(coords, coords, fun = distHaversine) / 1000
rownames(dist_matrix) <- district_covariates$admin2_unique
colnames(dist_matrix) <- district_covariates$admin2_unique
diag(dist_matrix) <- 0

n_districts <- nrow(district_covariates)

cat(sprintf("Distance matrix: %d x %d\n", n_districts, n_districts))
cat(sprintf("Mean distance: %.0f km\n", mean(dist_matrix[dist_matrix > 0])))
cat(sprintf("Median distance: %.0f km\n\n", median(dist_matrix[dist_matrix > 0])))

# =============================================================================
# CREATE COUNT MATRICES
# =============================================================================

# Total counts - use admin2_unique for unique columns
panel_wide_total <- panel %>%
  select(admin2_unique, date, count) %>%
  pivot_wider(names_from = admin2_unique, values_from = count, values_fill = 0)

# Fatal counts
panel_wide_fatal <- panel %>%
  select(admin2_unique, date, count_fatal) %>%
  pivot_wider(names_from = admin2_unique, values_from = count_fatal, values_fill = 0)

# Non-fatal counts
panel_wide_nonfatal <- panel %>%
  select(admin2_unique, date, count_nonfatal) %>%
  pivot_wider(names_from = admin2_unique, values_from = count_nonfatal, values_fill = 0)

dates_vec <- panel_wide_total$date

count_matrix_total <- as.matrix(panel_wide_total[, -1])
count_matrix_total <- count_matrix_total[, rownames(dist_matrix)]

count_matrix_fatal <- as.matrix(panel_wide_fatal[, -1])
count_matrix_fatal <- count_matrix_fatal[, rownames(dist_matrix)]

count_matrix_nonfatal <- as.matrix(panel_wide_nonfatal[, -1])
count_matrix_nonfatal <- count_matrix_nonfatal[, rownames(dist_matrix)]

n_days <- nrow(count_matrix_total)
n_dist <- ncol(count_matrix_total)

cat(sprintf("Count matrices: %d days x %d districts\n\n", n_days, n_dist))

# =============================================================================
# COMPUTE SPATIAL WEIGHTS
# =============================================================================

compute_spatial_weights <- function(theta, dist_matrix) {
  if (is.infinite(theta)) {
    W <- (dist_matrix > 0) * 1
  } else {
    W <- exp(-dist_matrix / theta)
  }
  diag(W) <- 0
  sweep(W, 1, pmax(rowSums(W), 1e-10), FUN = "/")
}

# =============================================================================
# COMPUTE CROSS-DISTRICT EXPOSURE TERMS
# =============================================================================

compute_cross_exposure <- function(W_norm, kappa, count_matrix, n_days, n_dist, max_lag = 30) {
  temporal_weights <- kappa^(1:max_lag)
  temporal_weights <- temporal_weights / sum(temporal_weights)

  cross_exposure <- matrix(0, n_days, n_dist)
  for (lag in 1:max_lag) {
    lagged <- rbind(matrix(0, lag, n_dist), count_matrix[1:(n_days - lag), ])
    cross_exposure <- cross_exposure + temporal_weights[lag] * t(W_norm %*% t(lagged))
  }
  cross_exposure
}

cat("=== COMPUTING EXPOSURE TERMS ===\n\n")

# Spatial weights matrix
W_norm <- compute_spatial_weights(THETA_OPTIMAL, dist_matrix)

# Compute separate exposure terms
cross_exp_total <- compute_cross_exposure(W_norm, KAPPA_OPTIMAL, count_matrix_total,
                                          n_days, n_dist, MAX_LAG)
cross_exp_fatal <- compute_cross_exposure(W_norm, KAPPA_OPTIMAL, count_matrix_fatal,
                                          n_days, n_dist, MAX_LAG)
cross_exp_nonfatal <- compute_cross_exposure(W_norm, KAPPA_OPTIMAL, count_matrix_nonfatal,
                                             n_days, n_dist, MAX_LAG)

# Compute national counts (for fixed effects)
national_count <- rowSums(count_matrix_total)

cat("Exposure terms computed:\n")
cat(sprintf("  Cross_total: mean = %.4f, sd = %.4f\n",
            mean(cross_exp_total), sd(cross_exp_total)))
cat(sprintf("  Cross_fatal: mean = %.4f, sd = %.4f\n",
            mean(cross_exp_fatal), sd(cross_exp_fatal)))
cat(sprintf("  Cross_nonfatal: mean = %.4f, sd = %.4f\n\n",
            mean(cross_exp_nonfatal), sd(cross_exp_nonfatal)))

# Check for collinearity
cor_fatal_nonfatal <- cor(as.vector(cross_exp_fatal), as.vector(cross_exp_nonfatal))
cat(sprintf("Correlation between fatal and non-fatal exposure: %.3f\n", cor_fatal_nonfatal))
cat(sprintf("Collinearity check: %s\n\n",
            ifelse(abs(cor_fatal_nonfatal) < 0.95, "PASSED ✓", "WARNING: High collinearity!")))

# =============================================================================
# MERGE EXPOSURE TERMS INTO PANEL
# =============================================================================

exp_long <- data.frame(
  date = rep(dates_vec, n_dist),
  admin2_unique = rep(colnames(count_matrix_total), each = n_days),
  cross_total = as.vector(cross_exp_total),
  cross_fatal = as.vector(cross_exp_fatal),
  cross_nonfatal = as.vector(cross_exp_nonfatal),
  national_count = rep(national_count, n_dist)
)

# Merge with panel
panel_fit <- panel %>%
  left_join(exp_long, by = c("admin2_unique", "date")) %>%
  mutate(admin2_factor = factor(admin2_unique)) %>%
  filter(!is.na(log_pop), !is.na(cpi))

cat(sprintf("Panel for fitting: %d observations\n\n", nrow(panel_fit)))

# =============================================================================
# MODEL 1: BASELINE (COMBINED EXPOSURE)
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL 1: BASELINE (COMBINED CROSS-DISTRICT EXPOSURE)       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

m_baseline <- bam(count ~ cpi + national_count + cross_total +
                    s(admin2_factor, bs = "re") +
                    offset(log_pop),
                  data = panel_fit,
                  family = quasipoisson(link = "log"),
                  method = "fREML",
                  discrete = TRUE)

coef_baseline <- coef(m_baseline)
se_baseline <- sqrt(diag(vcov(m_baseline)))

cat("Baseline model (combined exposure):\n")
cat(sprintf("  α_cross (combined): %.4f (SE: %.4f)\n",
            coef_baseline["cross_total"], se_baseline["cross_total"]))
cat(sprintf("  Deviance: %.2f\n", deviance(m_baseline)))
cat(sprintf("  Dispersion: %.3f\n\n", summary(m_baseline)$dispersion))

# =============================================================================
# MODEL 2: DECOMPOSED EXPOSURE (FATAL vs NON-FATAL)
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL 2: DECOMPOSED EXPOSURE (FATAL vs NON-FATAL)          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

m_decomposed <- bam(count ~ cpi + national_count + cross_fatal + cross_nonfatal +
                      s(admin2_factor, bs = "re") +
                      offset(log_pop),
                    data = panel_fit,
                    family = quasipoisson(link = "log"),
                    method = "fREML",
                    discrete = TRUE)

coef_decomposed <- coef(m_decomposed)
se_decomposed <- sqrt(diag(vcov(m_decomposed)))

cat("Decomposed model results:\n")
cat(sprintf("  α_fatal: %.4f (SE: %.4f)\n",
            coef_decomposed["cross_fatal"], se_decomposed["cross_fatal"]))
cat(sprintf("  α_nonfatal: %.4f (SE: %.4f)\n",
            coef_decomposed["cross_nonfatal"], se_decomposed["cross_nonfatal"]))
cat(sprintf("  Deviance: %.2f\n", deviance(m_decomposed)))
cat(sprintf("  Dispersion: %.3f\n\n", summary(m_decomposed)$dispersion))

# Multiplicative interpretation
cat("Multiplicative interpretation:\n")
cat(sprintf("  exp(α_fatal) = %.3f (%.1f%% increase per unit fatal exposure)\n",
            exp(coef_decomposed["cross_fatal"]),
            100 * (exp(coef_decomposed["cross_fatal"]) - 1)))
cat(sprintf("  exp(α_nonfatal) = %.3f (%.1f%% increase per unit non-fatal exposure)\n\n",
            exp(coef_decomposed["cross_nonfatal"]),
            100 * (exp(coef_decomposed["cross_nonfatal"]) - 1)))

# =============================================================================
# WALD TEST FOR COEFFICIENT EQUALITY
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  WALD TEST: α_fatal = α_nonfatal                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Extract variance-covariance matrix
vcov_decomposed <- vcov(m_decomposed)

alpha_fatal <- coef_decomposed["cross_fatal"]
alpha_nonfatal <- coef_decomposed["cross_nonfatal"]

var_fatal <- vcov_decomposed["cross_fatal", "cross_fatal"]
var_nonfatal <- vcov_decomposed["cross_nonfatal", "cross_nonfatal"]
cov_fatal_nonfatal <- vcov_decomposed["cross_fatal", "cross_nonfatal"]

# Wald statistic: z = (α_fatal - α_nonfatal) / SE(α_fatal - α_nonfatal)
diff <- alpha_fatal - alpha_nonfatal
se_diff <- sqrt(var_fatal + var_nonfatal - 2 * cov_fatal_nonfatal)
z_stat <- diff / se_diff

# Two-sided p-value
p_value_twosided <- 2 * pnorm(-abs(z_stat))

# One-sided p-value (H1: α_fatal > α_nonfatal)
p_value_onesided <- pnorm(-z_stat)

cat("H0: α_fatal = α_nonfatal (equal triggering effects)\n")
cat("H1: α_fatal ≠ α_nonfatal (differential triggering effects)\n\n")

cat("Wald test results:\n")
cat(sprintf("  α_fatal - α_nonfatal = %.4f\n", diff))
cat(sprintf("  SE(difference) = %.4f\n", se_diff))
cat(sprintf("  z-statistic = %.3f\n", z_stat))
cat(sprintf("  p-value (two-sided) = %.4f\n", p_value_twosided))
cat(sprintf("  p-value (one-sided, H1: fatal > nonfatal) = %.4f\n\n", p_value_onesided))

if (p_value_twosided < 0.05) {
  cat("*** RESULT: Reject H0 at α = 0.05 ***\n")
  cat("    Fatal and non-fatal protests have DIFFERENT triggering effects.\n")
  if (diff > 0) {
    cat("    Fatal protests have HIGHER triggering effect than non-fatal.\n\n")
  } else {
    cat("    Non-fatal protests have HIGHER triggering effect than fatal.\n\n")
  }
} else {
  cat("RESULT: Fail to reject H0 at α = 0.05\n")
  cat("  No significant difference in triggering effects.\n\n")
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SUMMARY TABLE                                               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Create summary table
summary_table <- data.frame(
  Parameter = c("α_fatal", "α_nonfatal", "Difference (α_fatal - α_nonfatal)"),
  Estimate = c(alpha_fatal, alpha_nonfatal, diff),
  SE = c(sqrt(var_fatal), sqrt(var_nonfatal), se_diff),
  z_value = c(alpha_fatal / sqrt(var_fatal),
              alpha_nonfatal / sqrt(var_nonfatal),
              z_stat),
  p_value = c(2 * pnorm(-abs(alpha_fatal / sqrt(var_fatal))),
              2 * pnorm(-abs(alpha_nonfatal / sqrt(var_nonfatal))),
              p_value_twosided)
)

print(summary_table, digits = 4, row.names = FALSE)

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("Notes:\n")
cat("  - Two-sided p-values reported\n")
cat("  - Coefficients are on log scale (Poisson link)\n")
cat(sprintf("  - Spatial decay: θ = %d km\n", THETA_OPTIMAL))
cat(sprintf("  - Temporal decay: κ = %.2f (half-life = %.1f days)\n",
            KAPPA_OPTIMAL, log(0.5) / log(KAPPA_OPTIMAL)))
cat("═══════════════════════════════════════════════════════════════\n")

# =============================================================================
# MODEL COMPARISON: BASELINE vs DECOMPOSED
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL COMPARISON: BASELINE vs DECOMPOSED                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Deviance comparison (quasi-LRT)
dev_baseline <- deviance(m_baseline)
dev_decomposed <- deviance(m_decomposed)
dispersion <- summary(m_decomposed)$dispersion

# Additional parameter in decomposed model
df_diff <- 1  # split one coefficient into two

F_stat <- ((dev_baseline - dev_decomposed) / df_diff) / dispersion
p_value_F <- pf(F_stat, df_diff, df.residual(m_decomposed), lower.tail = FALSE)

cat("Quasi-likelihood ratio test:\n")
cat(sprintf("  Deviance (baseline): %.2f\n", dev_baseline))
cat(sprintf("  Deviance (decomposed): %.2f\n", dev_decomposed))
cat(sprintf("  Deviance reduction: %.2f\n", dev_baseline - dev_decomposed))
cat(sprintf("  Dispersion: %.3f\n", dispersion))
cat(sprintf("  F-statistic: %.3f\n", F_stat))
cat(sprintf("  p-value: %.4f\n\n", p_value_F))

if (p_value_F < 0.05) {
  cat("*** Decomposed model fits significantly better ***\n")
} else {
  cat("No significant improvement from decomposition.\n")
}

# =============================================================================
# EFFECT SIZE COMPARISON
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  EFFECT SIZE COMPARISON                                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Weighted average should approximate baseline coefficient
mean_fatal <- mean(panel_fit$cross_fatal)
mean_nonfatal <- mean(panel_fit$cross_nonfatal)
total_mean <- mean_fatal + mean_nonfatal

if (total_mean > 0) {
  weight_fatal <- mean_fatal / total_mean
  weight_nonfatal <- mean_nonfatal / total_mean

  weighted_avg <- weight_fatal * alpha_fatal + weight_nonfatal * alpha_nonfatal

  cat("Consistency check with baseline model:\n")
  cat(sprintf("  Baseline α_cross: %.4f\n", coef_baseline["cross_total"]))
  cat(sprintf("  Weighted average (fatal: %.1f%%, nonfatal: %.1f%%): %.4f\n",
              100 * weight_fatal, 100 * weight_nonfatal, weighted_avg))
  cat(sprintf("  Difference: %.4f\n\n", coef_baseline["cross_total"] - weighted_avg))
}

# Ratio of effects
if (alpha_nonfatal != 0) {
  ratio <- alpha_fatal / alpha_nonfatal
  cat(sprintf("Ratio of effects: α_fatal / α_nonfatal = %.2f\n", ratio))
  if (ratio > 1) {
    cat(sprintf("  → Fatal protests are %.0f%% more effective at triggering\n",
                100 * (ratio - 1)))
  } else if (ratio < 1) {
    cat(sprintf("  → Non-fatal protests are %.0f%% more effective at triggering\n",
                100 * (1/ratio - 1)))
  }
}

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("\n")
cat("=== SAVING RESULTS ===\n\n")

results <- list(
  optimal_params = list(
    theta = THETA_OPTIMAL,
    kappa = KAPPA_OPTIMAL,
    max_lag = MAX_LAG
  ),
  baseline_model = list(
    coefficients = coef_baseline,
    se = se_baseline,
    deviance = dev_baseline
  ),
  decomposed_model = list(
    alpha_fatal = alpha_fatal,
    alpha_nonfatal = alpha_nonfatal,
    se_fatal = sqrt(var_fatal),
    se_nonfatal = sqrt(var_nonfatal),
    covariance = cov_fatal_nonfatal,
    deviance = dev_decomposed,
    dispersion = dispersion
  ),
  wald_test = list(
    difference = diff,
    se_difference = se_diff,
    z_statistic = z_stat,
    p_value_twosided = p_value_twosided,
    p_value_onesided = p_value_onesided
  ),
  model_comparison = list(
    F_statistic = F_stat,
    p_value = p_value_F,
    deviance_reduction = dev_baseline - dev_decomposed
  ),
  summary_table = summary_table
)

saveRDS(results, "fatality_heterogeneity_results.rds")
cat("Results saved to: fatality_heterogeneity_results.rds\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FINAL SUMMARY                                               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("RESEARCH QUESTION:\n")
cat("  Do protests with fatalities trigger more subsequent protests\n")
cat("  than non-fatal protests?\n\n")

cat("METHOD:\n")
cat("  Decomposed cross-district exposure into fatal and non-fatal\n")
cat("  components, using optimal spatial-temporal parameters from\n")
cat("  the main analysis.\n\n")

cat("KEY FINDINGS:\n")
cat(sprintf("  1. Fatal protest triggering effect (α_fatal): %.4f\n", alpha_fatal))
cat(sprintf("  2. Non-fatal protest triggering effect (α_nonfatal): %.4f\n", alpha_nonfatal))
cat(sprintf("  3. Difference: %.4f (SE: %.4f)\n", diff, se_diff))
cat(sprintf("  4. z-statistic: %.3f\n", z_stat))
cat(sprintf("  5. p-value (two-sided): %.4f\n\n", p_value_twosided))

if (p_value_twosided < 0.05) {
  cat("CONCLUSION: Significant difference in triggering effects\n")
  if (diff > 0) {
    cat("  → Fatal protests trigger MORE subsequent protests than non-fatal.\n")
  } else {
    cat("  → Non-fatal protests trigger MORE subsequent protests than fatal.\n")
  }
} else {
  cat("CONCLUSION: No significant difference in triggering effects\n")
  cat("  → Fatal and non-fatal protests have similar contagion effects.\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
