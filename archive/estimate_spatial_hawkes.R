################################################################################
#       DISCRETE SPATIAL HAWKES MODEL: CROSS-DISTRICT CONTAGION
#
#       Tests whether protests in other districts trigger local protests,
#       and whether distance matters for this cross-district contagion.
#
#       Hypothesis: Cross-district contagion exists but distance doesn't matter
#                   (i.e., national-level mechanisms rather than local spillovers)
################################################################################

library(tidyverse)
library(geosphere)
library(MASS)

select <- dplyr::select  # Avoid MASS conflict

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   DISCRETE SPATIAL HAWKES MODEL: CROSS-DISTRICT CONTAGION               ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

protests <- readRDS("protests_daily.rds")
cpi_data <- read_csv("indonesia_cpi_processed.csv", show_col_types = FALSE)

cat(sprintf("Total events: %d\n", nrow(protests)))
cat(sprintf("Date range: %s to %s\n", min(protests$date), max(protests$date)))
cat(sprintf("Districts: %d\n", n_distinct(protests$admin2)))

# =============================================================================
# PREPARE DISTRICT-LEVEL COVARIATES
# =============================================================================

# Get district-level covariates (time-invariant: use mean across observations)
district_covariates <- protests %>%
  group_by(admin2) %>%
  summarise(
    log_pop = mean(log_pop, na.rm = TRUE),
    poverty_rate = mean(poverty_rate, na.rm = TRUE),
    lon = mean(longitude, na.rm = TRUE),
    lat = mean(latitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(admin2), !is.na(lon), !is.na(lat))

cat(sprintf("\nDistricts with complete covariates: %d\n", nrow(district_covariates)))

# Prepare CPI data (monthly, merge by year-month)
cpi_monthly <- cpi_data %>%
  select(year_month, cpi, inflation_yoy) %>%
  mutate(year_month = as.character(year_month))

# =============================================================================
# CREATE DISTRICT-DAY PANEL
# =============================================================================

cat("\nCreating district-day panel...\n")

all_districts <- district_covariates$admin2
all_dates <- seq(min(protests$date), max(protests$date), by = "day")

# Count protests per district-day
daily_counts <- protests %>%
  filter(admin2 %in% all_districts) %>%
  group_by(admin2, date) %>%
  summarise(count = n(), .groups = "drop")

# Create full panel
panel <- expand.grid(
  admin2 = all_districts,
  date = all_dates,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  left_join(daily_counts, by = c("admin2", "date")) %>%
  mutate(
    count = replace_na(count, 0),
    year_month = format(date, "%Y-%m"),
    year = year(date)
  ) %>%
  left_join(district_covariates, by = "admin2") %>%
  left_join(cpi_monthly, by = "year_month") %>%
  arrange(admin2, date)

cat(sprintf("Panel size: %d district-days\n", nrow(panel)))
cat(sprintf("Mean daily count per district: %.4f\n", mean(panel$count)))
cat(sprintf("Districts with protests: %.1f%%\n", 100 * mean(panel$count > 0)))

# =============================================================================
# COMPUTE DISTANCE MATRIX
# =============================================================================

cat("\nComputing inter-district distances...\n")

n_districts <- nrow(district_covariates)
dist_matrix <- matrix(0, n_districts, n_districts)
rownames(dist_matrix) <- district_covariates$admin2
colnames(dist_matrix) <- district_covariates$admin2

for (i in 1:n_districts) {
  for (j in 1:n_districts) {
    if (i != j) {
      dist_matrix[i, j] <- distHaversine(
        c(district_covariates$lon[i], district_covariates$lat[i]),
        c(district_covariates$lon[j], district_covariates$lat[j])
      ) / 1000  # km
    }
  }
}

cat(sprintf("Mean inter-district distance: %.0f km\n", mean(dist_matrix[dist_matrix > 0])))
cat(sprintf("Median inter-district distance: %.0f km\n", median(dist_matrix[dist_matrix > 0])))

# Distance band thresholds
NEAR_THRESHOLD <- 100    # 0-100 km
MED_THRESHOLD <- 300     # 100-300 km
# Far = 300+ km

# Count neighbors in each band
n_near <- rowSums(dist_matrix > 0 & dist_matrix <= NEAR_THRESHOLD)
n_med <- rowSums(dist_matrix > NEAR_THRESHOLD & dist_matrix <= MED_THRESHOLD)
n_far <- rowSums(dist_matrix > MED_THRESHOLD)

cat(sprintf("\nNeighbors per district:\n"))
cat(sprintf("  Near (0-%d km):      mean=%.1f, median=%.0f\n", NEAR_THRESHOLD, mean(n_near), median(n_near)))
cat(sprintf("  Medium (%d-%d km):  mean=%.1f, median=%.0f\n", NEAR_THRESHOLD, MED_THRESHOLD, mean(n_med), median(n_med)))
cat(sprintf("  Far (%d+ km):       mean=%.1f, median=%.0f\n", MED_THRESHOLD, mean(n_far), median(n_far)))

# =============================================================================
# COMPUTE SPATIAL EXPOSURE VARIABLES
# =============================================================================

cat("\nComputing spatial exposure variables (this may take a few minutes)...\n")

# Create binary weight matrices for each distance band
W_near <- (dist_matrix > 0 & dist_matrix <= NEAR_THRESHOLD) * 1
W_med <- (dist_matrix > NEAR_THRESHOLD & dist_matrix <= MED_THRESHOLD) * 1
W_far <- (dist_matrix > MED_THRESHOLD) * 1

# Row-normalize (average of neighbors, not sum)
W_near_norm <- W_near / pmax(rowSums(W_near), 1)
W_med_norm <- W_med / pmax(rowSums(W_med), 1)
W_far_norm <- W_far / pmax(rowSums(W_far), 1)

# Create wide format: rows = dates, cols = districts
panel_wide <- panel %>%
  select(admin2, date, count) %>%
  pivot_wider(names_from = admin2, values_from = count, values_fill = 0)

dates_vec <- panel_wide$date
count_matrix <- as.matrix(panel_wide[, -1])  # Remove date column
colnames(count_matrix) <- names(panel_wide)[-1]

# Reorder columns to match distance matrix
count_matrix <- count_matrix[, rownames(dist_matrix)]

n_days <- nrow(count_matrix)
n_dist <- ncol(count_matrix)

# Initialize exposure matrices
own_exposure <- matrix(0, n_days, n_dist)
near_exposure <- matrix(0, n_days, n_dist)
med_exposure <- matrix(0, n_days, n_dist)
far_exposure <- matrix(0, n_days, n_dist)

# Compute lagged sums (lags 1-7)
cat("Computing 7-day lagged exposures...\n")

for (lag in 1:7) {
  if (lag <= n_days) {
    # Lagged count matrix
    lagged_counts <- rbind(
      matrix(0, lag, n_dist),
      count_matrix[1:(n_days - lag), ]
    )

    # Own-district exposure (just the lagged count itself)
    own_exposure <- own_exposure + lagged_counts

    # Cross-district exposure by distance band
    # Each row of lagged_counts gets multiplied by the weight matrix
    near_exposure <- near_exposure + t(W_near_norm %*% t(lagged_counts))
    med_exposure <- med_exposure + t(W_med_norm %*% t(lagged_counts))
    far_exposure <- far_exposure + t(W_far_norm %*% t(lagged_counts))
  }
}

cat("Done.\n")

# Convert back to long format and merge with panel
exposure_long <- data.frame(
  date = rep(dates_vec, n_dist),
  admin2 = rep(colnames(count_matrix), each = n_days),
  own_lag = as.vector(own_exposure),
  near_lag = as.vector(near_exposure),
  med_lag = as.vector(med_exposure),
  far_lag = as.vector(far_exposure)
)

panel <- panel %>%
  left_join(exposure_long, by = c("admin2", "date"))

# Remove rows with missing covariates
panel_clean <- panel %>%
  filter(!is.na(log_pop), !is.na(poverty_rate), !is.na(cpi), !is.na(own_lag))

cat(sprintf("\nClean panel size: %d observations\n", nrow(panel_clean)))

# =============================================================================
# MODEL ESTIMATION
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL ESTIMATION                                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# -----------------------------------------------------------------------------
# M0: Baseline with covariates only (no contagion)
# -----------------------------------------------------------------------------

cat("--- M0: Baseline (Covariates Only, No Contagion) ---\n\n")

# Use log link for numerical stability with sparse counts
m0 <- glm(count ~ log_pop + poverty_rate + cpi + factor(year),
          data = panel_clean, family = quasipoisson(link = "log"))

cat("Covariate effects:\n")
print(round(summary(m0)$coefficients[1:4, ], 5))
cat(sprintf("\nDispersion: %.3f\n", summary(m0)$dispersion))

# -----------------------------------------------------------------------------
# M1: Add own-district lagged effect (within-district contagion)
# -----------------------------------------------------------------------------

cat("\n--- M1: Within-District Contagion Only ---\n\n")

m1 <- glm(count ~ log_pop + poverty_rate + cpi + factor(year) + own_lag,
          data = panel_clean, family = quasipoisson(link = "log"))

own_coef <- coef(m1)["own_lag"]
own_se <- summary(m1)$coefficients["own_lag", "Std. Error"]

cat(sprintf("Own-district effect (α_own): %.5f (SE: %.5f)\n", own_coef, own_se))
cat(sprintf("  t-value: %.2f, p < 0.001\n", own_coef / own_se))

# F-test: M1 vs M0
f_test_m1 <- anova(m0, m1, test = "F")
cat(sprintf("\nF-test (M1 vs M0): F = %.2f, p = %.2e\n",
            f_test_m1$F[2], f_test_m1$`Pr(>F)`[2]))

# -----------------------------------------------------------------------------
# M2: Add cross-district effects (single coefficient for all other districts)
# -----------------------------------------------------------------------------

cat("\n--- M2: Cross-District Contagion (Pooled) ---\n\n")

# Create pooled cross-district exposure (average of near, med, far)
panel_clean <- panel_clean %>%
  mutate(other_lag = (near_lag + med_lag + far_lag) / 3)

m2 <- glm(count ~ log_pop + poverty_rate + cpi + factor(year) + own_lag + other_lag,
          data = panel_clean, family = quasipoisson(link = "log"))

other_coef <- coef(m2)["other_lag"]
other_se <- summary(m2)$coefficients["other_lag", "Std. Error"]

cat(sprintf("Own-district effect:    %.5f (SE: %.5f)\n", coef(m2)["own_lag"],
            summary(m2)$coefficients["own_lag", "Std. Error"]))
cat(sprintf("Other-district effect:  %.5f (SE: %.5f)\n", other_coef, other_se))

# F-test: M2 vs M1
f_test_m2 <- anova(m1, m2, test = "F")
cat(sprintf("\nF-test (M2 vs M1): F = %.2f, p = %.2e\n",
            f_test_m2$F[2], f_test_m2$`Pr(>F)`[2]))

if (f_test_m2$`Pr(>F)`[2] < 0.05) {
  cat("  *** Cross-district contagion is significant ***\n")
} else {
  cat("  Cross-district contagion is NOT significant\n")
}

# -----------------------------------------------------------------------------
# M3: Decompose cross-district by distance bands
# -----------------------------------------------------------------------------

cat("\n--- M3: Cross-District by Distance Bands ---\n\n")

m3 <- glm(count ~ log_pop + poverty_rate + cpi + factor(year) +
            own_lag + near_lag + med_lag + far_lag,
          data = panel_clean, family = quasipoisson(link = "log"))

spatial_coefs <- summary(m3)$coefficients[c("own_lag", "near_lag", "med_lag", "far_lag"), ]

cat("Spatial coefficients:\n")
print(round(spatial_coefs, 6))

cat("\n─────────────────────────────────────────────────────────────────\n")
cat("CROSS-DISTRICT EFFECTS BY DISTANCE:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Own district:         %.5f (SE: %.5f)\n", spatial_coefs["own_lag", 1], spatial_coefs["own_lag", 2]))
cat(sprintf("  Near (0-%d km):       %.5f (SE: %.5f)\n", NEAR_THRESHOLD, spatial_coefs["near_lag", 1], spatial_coefs["near_lag", 2]))
cat(sprintf("  Medium (%d-%d km):   %.5f (SE: %.5f)\n", NEAR_THRESHOLD, MED_THRESHOLD, spatial_coefs["med_lag", 1], spatial_coefs["med_lag", 2]))
cat(sprintf("  Far (%d+ km):        %.5f (SE: %.5f)\n", MED_THRESHOLD, spatial_coefs["far_lag", 1], spatial_coefs["far_lag", 2]))

# -----------------------------------------------------------------------------
# TEST: Does distance matter? (H0: α_near = α_med = α_far)
# -----------------------------------------------------------------------------

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  KEY TEST: DOES DISTANCE MATTER FOR CROSS-DISTRICT CONTAGION?║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("H0: α_near = α_medium = α_far (no spatial decay)\n")
cat("H1: Distance affects triggering (spatial decay exists)\n\n")

# F-test: M3 vs M2 (do the separate distance bands improve fit over pooled?)
f_test_distance <- anova(m2, m3, test = "F")
f_stat <- f_test_distance$F[2]
p_val <- f_test_distance$`Pr(>F)`[2]

cat(sprintf("F-test (separate bands vs pooled): F = %.3f, p = %.4f\n\n", f_stat, p_val))

if (p_val < 0.05) {
  cat("RESULT: Distance DOES matter (reject H0)\n")
  cat("  → Spatial decay exists in cross-district contagion\n")

  # Check direction
  if (spatial_coefs["near_lag", 1] > spatial_coefs["far_lag", 1]) {
    cat("  → Near protests trigger MORE than far protests (local diffusion)\n")
  } else {
    cat("  → Far protests trigger MORE than near protests (unexpected pattern)\n")
  }
} else {
  cat("RESULT: Distance does NOT matter (fail to reject H0)\n")
  cat("  → No spatial decay: near and far protests trigger equally\n")
  cat("  → Supports national-level mechanism (media, coordination)\n")
}

# Pairwise comparisons
cat("\n--- Pairwise Coefficient Comparisons ---\n")

# Near vs Far
diff_near_far <- spatial_coefs["near_lag", 1] - spatial_coefs["far_lag", 1]
se_diff <- sqrt(spatial_coefs["near_lag", 2]^2 + spatial_coefs["far_lag", 2]^2)
z_near_far <- diff_near_far / se_diff
p_near_far <- 2 * pnorm(-abs(z_near_far))

cat(sprintf("Near vs Far: diff = %.5f, z = %.2f, p = %.4f\n", diff_near_far, z_near_far, p_near_far))

# Near vs Medium
diff_near_med <- spatial_coefs["near_lag", 1] - spatial_coefs["med_lag", 1]
se_diff_nm <- sqrt(spatial_coefs["near_lag", 2]^2 + spatial_coefs["med_lag", 2]^2)
z_near_med <- diff_near_med / se_diff_nm
p_near_med <- 2 * pnorm(-abs(z_near_med))

cat(sprintf("Near vs Medium: diff = %.5f, z = %.2f, p = %.4f\n", diff_near_med, z_near_med, p_near_med))

# =============================================================================
# PERMUTATION TEST FOR CROSS-DISTRICT EFFECTS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PERMUTATION TEST: IS CROSS-DISTRICT CONTAGION REAL?        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

N_PERM <- 500
cat(sprintf("Running %d permutations...\n", N_PERM))

observed_other <- other_coef

null_other <- numeric(N_PERM)

set.seed(42)
pb <- txtProgressBar(min = 0, max = N_PERM, style = 3)

for (i in 1:N_PERM) {
  # Shuffle dates within each district (breaks temporal structure)
  panel_perm <- panel_clean %>%
    group_by(admin2) %>%
    mutate(count = sample(count)) %>%
    ungroup()

  # Re-fit model with pooled other-district effect
  m_perm <- glm(count ~ log_pop + poverty_rate + cpi + factor(year) + own_lag + other_lag,
                data = panel_perm, family = quasipoisson(link = "log"))

  null_other[i] <- coef(m_perm)["other_lag"]

  setTxtProgressBar(pb, i)
}
close(pb)

# Calculate p-value
p_perm <- mean(null_other >= observed_other)

cat("\n\nPermutation Test Results:\n")
cat(sprintf("  Observed cross-district effect: %.5f\n", observed_other))
cat(sprintf("  Null mean: %.5f, sd: %.5f\n", mean(null_other), sd(null_other)))
cat(sprintf("  Null 99th percentile: %.5f\n", quantile(null_other, 0.99)))
cat(sprintf("  p-value: %.4f\n", p_perm))

if (p_perm < 0.05) {
  cat("\n  *** Cross-district contagion is REAL (not an artifact) ***\n")
} else {
  cat("\n  Cross-district effect may be spurious (artifact of common shocks)\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SUMMARY                                                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("BACKGROUND RATE CONTROLS:\n")
cat(sprintf("  Log population:  %.5f (SE: %.5f)\n", coef(m3)["log_pop"], summary(m3)$coefficients["log_pop", 2]))
cat(sprintf("  Poverty rate:    %.5f (SE: %.5f)\n", coef(m3)["poverty_rate"], summary(m3)$coefficients["poverty_rate", 2]))
cat(sprintf("  CPI:             %.5f (SE: %.5f)\n", coef(m3)["cpi"], summary(m3)$coefficients["cpi", 2]))

cat("\nSPATIAL CONTAGION EFFECTS:\n")
cat(sprintf("  Within-district (own):  %.5f ***\n", spatial_coefs["own_lag", 1]))
cat(sprintf("  Near (0-%d km):         %.5f\n", NEAR_THRESHOLD, spatial_coefs["near_lag", 1]))
cat(sprintf("  Medium (%d-%d km):     %.5f\n", NEAR_THRESHOLD, MED_THRESHOLD, spatial_coefs["med_lag", 1]))
cat(sprintf("  Far (%d+ km):          %.5f\n", MED_THRESHOLD, spatial_coefs["far_lag", 1]))

cat("\nKEY FINDINGS:\n")
cat(sprintf("  1. Cross-district contagion: %s (p = %.2e)\n",
            ifelse(f_test_m2$`Pr(>F)`[2] < 0.05, "EXISTS", "NOT SIGNIFICANT"),
            f_test_m2$`Pr(>F)`[2]))
cat(sprintf("  2. Spatial decay (distance matters): %s (p = %.4f)\n",
            ifelse(p_val < 0.05, "YES", "NO"),
            p_val))
cat(sprintf("  3. Permutation validation: %s (p = %.4f)\n",
            ifelse(p_perm < 0.05, "CONFIRMED", "INCONCLUSIVE"),
            p_perm))

# =============================================================================
# SAVE RESULTS
# =============================================================================

results <- list(
  m0_baseline = list(
    dispersion = summary(m0)$dispersion
  ),
  m1_within = list(
    own_effect = own_coef,
    own_se = own_se,
    f_test = f_test_m1
  ),
  m2_pooled_cross = list(
    own_effect = coef(m2)["own_lag"],
    other_effect = other_coef,
    f_test = f_test_m2
  ),
  m3_distance_bands = list(
    coefficients = spatial_coefs,
    own_effect = spatial_coefs["own_lag", 1],
    near_effect = spatial_coefs["near_lag", 1],
    med_effect = spatial_coefs["med_lag", 1],
    far_effect = spatial_coefs["far_lag", 1],
    f_test_distance = f_test_distance
  ),
  distance_test = list(
    f_stat = f_stat,
    p_value = p_val,
    distance_matters = p_val < 0.05
  ),
  permutation = list(
    observed = observed_other,
    null_distribution = null_other,
    null_mean = mean(null_other),
    null_sd = sd(null_other),
    p_value = p_perm
  ),
  covariate_effects = list(
    log_pop = coef(m3)["log_pop"],
    poverty = coef(m3)["poverty_rate"],
    cpi = coef(m3)["cpi"]
  )
)

saveRDS(results, "spatial_hawkes_results.rds")
cat("\nResults saved to: spatial_hawkes_results.rds\n")
