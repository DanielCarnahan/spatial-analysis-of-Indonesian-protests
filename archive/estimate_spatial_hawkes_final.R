################################################################################
#       SPATIAL HAWKES MODEL - CROSS-DISTRICT CONTAGION ONLY
#
#       Model: E[Y_tr] = exp(X'β + α × Cross_tr(θ))
#
#       Tests:
#       1. Does cross-district contagion exist? (M0 vs M1)
#       2. Does distance matter? (θ = ∞ vs θ < ∞)
################################################################################

library(tidyverse)
library(geosphere)
library(MASS)

select <- dplyr::select

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   SPATIAL HAWKES MODEL - CROSS-DISTRICT CONTAGION                       ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

protests <- readRDS("protests_daily.rds")
cpi_data <- read_csv("indonesia_cpi_processed.csv", show_col_types = FALSE)

cat(sprintf("Total events: %d\n", nrow(protests)))
cat(sprintf("Date range: %s to %s\n", min(protests$date), max(protests$date)))

# =============================================================================
# PREPARE DATA
# =============================================================================

district_covariates <- protests %>%
  group_by(admin2) %>%
  summarise(
    log_pop = mean(log_pop, na.rm = TRUE),
    poverty_rate = mean(poverty_rate, na.rm = TRUE),
    lon = mean(longitude, na.rm = TRUE),
    lat = mean(latitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(admin2), !is.na(lon), !is.na(lat), !is.na(log_pop))

cat(sprintf("Districts: %d\n", nrow(district_covariates)))

cpi_monthly <- cpi_data %>%
  select(year_month, cpi) %>%
  mutate(year_month = as.character(year_month))

# =============================================================================
# CREATE PANEL
# =============================================================================

cat("\nCreating district-day panel...\n")

all_districts <- district_covariates$admin2
all_dates <- seq(min(protests$date), max(protests$date), by = "day")

daily_counts <- protests %>%
  filter(admin2 %in% all_districts) %>%
  group_by(admin2, date) %>%
  summarise(count = n(), .groups = "drop")

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

cat(sprintf("Panel: %d district-days\n", nrow(panel)))
cat(sprintf("Mean count: %.4f\n", mean(panel$count)))

# =============================================================================
# COMPUTE DISTANCE MATRIX
# =============================================================================

cat("\nComputing distances...\n")

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
      ) / 1000
    }
  }
}

cat(sprintf("Mean distance: %.0f km\n", mean(dist_matrix[dist_matrix > 0])))

# =============================================================================
# CREATE COUNT MATRIX
# =============================================================================

panel_wide <- panel %>%
  select(admin2, date, count) %>%
  pivot_wider(names_from = admin2, values_from = count, values_fill = 0)

dates_vec <- panel_wide$date
count_matrix <- as.matrix(panel_wide[, -1])
count_matrix <- count_matrix[, rownames(dist_matrix)]

n_days <- nrow(count_matrix)
n_dist <- ncol(count_matrix)

# =============================================================================
# FUNCTION: Compute cross-district exposure
# =============================================================================

compute_cross_exposure <- function(theta, dist_matrix, count_matrix, n_days, n_dist) {
  if (is.infinite(theta)) {
    W <- (dist_matrix > 0) * 1
  } else {
    W <- exp(-dist_matrix / theta)
    diag(W) <- 0
  }

  # Row-normalize
  W_norm <- W / pmax(rowSums(W), 1e-10)

  # Sum over lags 1-30
  cross_exposure <- matrix(0, n_days, n_dist)
  for (lag in 1:30) {
    lagged <- rbind(matrix(0, lag, n_dist), count_matrix[1:(n_days - lag), ])
    cross_exposure <- cross_exposure + t(W_norm %*% t(lagged))
  }

  return(cross_exposure)
}

# =============================================================================
# M0: BASELINE MODEL (NO CONTAGION)
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  M0: BASELINE MODEL (Background Rate Only)                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

panel_clean <- panel %>%
  filter(!is.na(log_pop), !is.na(poverty_rate), !is.na(cpi))

m0 <- glm(count ~ log_pop + poverty_rate + cpi + factor(year),
          data = panel_clean, family = quasipoisson(link = "log"))

cat("Coefficients:\n")
print(round(summary(m0)$coefficients[1:4, ], 5))
cat(sprintf("\nDeviance: %.2f\n", deviance(m0)))
cat(sprintf("Dispersion: %.3f\n", summary(m0)$dispersion))

dev_m0 <- deviance(m0)

# =============================================================================
# M1: CROSS-DISTRICT CONTAGION (with profile likelihood over θ)
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  M1: CROSS-DISTRICT CONTAGION MODEL                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

theta_grid <- c(50, 100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000, Inf)

results_grid <- data.frame(
  theta = theta_grid,
  deviance = NA,
  alpha = NA,
  se = NA,
  t_value = NA
)

cat("Profile likelihood over θ:\n\n")

for (i in seq_along(theta_grid)) {
  theta <- theta_grid[i]
  theta_label <- if (is.infinite(theta)) "∞" else as.character(theta)

  cross_exp <- compute_cross_exposure(theta, dist_matrix, count_matrix, n_days, n_dist)

  exp_long <- data.frame(
    date = rep(dates_vec, n_dist),
    admin2 = rep(colnames(count_matrix), each = n_days),
    cross_lag = as.vector(cross_exp)
  )

  panel_fit <- panel_clean %>%
    left_join(exp_long, by = c("admin2", "date"))

  m <- glm(count ~ log_pop + poverty_rate + cpi + factor(year) + cross_lag,
           data = panel_fit, family = quasipoisson(link = "log"))

  results_grid$deviance[i] <- deviance(m)
  results_grid$alpha[i] <- coef(m)["cross_lag"]
  results_grid$se[i] <- summary(m)$coefficients["cross_lag", "Std. Error"]
  results_grid$t_value[i] <- coef(m)["cross_lag"] / summary(m)$coefficients["cross_lag", "Std. Error"]

  cat(sprintf("  θ = %5s km: deviance = %.1f, α = %.4f (t = %.1f)\n",
              theta_label, deviance(m), coef(m)["cross_lag"], results_grid$t_value[i]))
}

# Find optimal
best_idx <- which.min(results_grid$deviance)
theta_hat <- results_grid$theta[best_idx]
inf_idx <- which(is.infinite(results_grid$theta))

cat(sprintf("\nOptimal θ: %s\n", ifelse(is.infinite(theta_hat), "∞ (no spatial decay)", theta_hat)))

# =============================================================================
# MODEL COMPARISON: M0 vs M1
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  TEST 1: DOES CROSS-DISTRICT CONTAGION EXIST?               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

dev_m1 <- results_grid$deviance[best_idx]

# Fit best model to get dispersion
cross_exp_best <- compute_cross_exposure(theta_hat, dist_matrix, count_matrix, n_days, n_dist)
exp_long_best <- data.frame(
  date = rep(dates_vec, n_dist),
  admin2 = rep(colnames(count_matrix), each = n_days),
  cross_lag = as.vector(cross_exp_best)
)
panel_best <- panel_clean %>%
  left_join(exp_long_best, by = c("admin2", "date"))

m1_best <- glm(count ~ log_pop + poverty_rate + cpi + factor(year) + cross_lag,
               data = panel_best, family = quasipoisson(link = "log"))

dispersion <- summary(m1_best)$dispersion

# F-test: M0 vs M1
F_stat_contagion <- ((dev_m0 - dev_m1) / 1) / dispersion
p_contagion <- pf(F_stat_contagion, 1, nrow(panel_clean) - length(coef(m1_best)), lower.tail = FALSE)

cat("Comparing M0 (no contagion) vs M1 (cross-district contagion):\n\n")
cat(sprintf("  M0 deviance: %.2f\n", dev_m0))
cat(sprintf("  M1 deviance: %.2f\n", dev_m1))
cat(sprintf("  Deviance reduction: %.2f (%.2f%%)\n", dev_m0 - dev_m1, 100 * (dev_m0 - dev_m1) / dev_m0))
cat(sprintf("  Dispersion: %.3f\n", dispersion))
cat(sprintf("  F-statistic: %.2f\n", F_stat_contagion))
cat(sprintf("  p-value: %.2e\n", p_contagion))

if (p_contagion < 0.001) {
  cat("\n  *** RESULT: Cross-district contagion EXISTS (p < 0.001) ***\n")
} else if (p_contagion < 0.05) {
  cat("\n  * RESULT: Cross-district contagion exists (p < 0.05) *\n")
} else {
  cat("\n  RESULT: No evidence of cross-district contagion\n")
}

# =============================================================================
# TEST 2: DOES DISTANCE MATTER?
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  TEST 2: DOES DISTANCE MATTER FOR CROSS-DISTRICT CONTAGION? ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("H0: θ = ∞ (no spatial decay)\n")
cat("H1: θ < ∞ (spatial decay exists)\n\n")

dev_inf <- results_grid$deviance[inf_idx]
dev_best <- results_grid$deviance[best_idx]

if (is.infinite(theta_hat)) {
  cat("θ = ∞ has lowest deviance → Distance does NOT matter\n")
  cat(sprintf("\nDeviance at θ = ∞: %.2f (optimal)\n", dev_inf))
  distance_matters <- FALSE
  p_distance <- 1.0
} else {
  LR_stat <- (dev_inf - dev_best) / dispersion
  p_distance <- pchisq(LR_stat, df = 1, lower.tail = FALSE)

  cat(sprintf("Deviance at θ = ∞:   %.2f\n", dev_inf))
  cat(sprintf("Deviance at θ = %d: %.2f\n", theta_hat, dev_best))
  cat(sprintf("LR statistic: %.2f\n", LR_stat))
  cat(sprintf("p-value: %.4f\n", p_distance))

  distance_matters <- p_distance < 0.05
}

if (distance_matters) {
  cat(sprintf("\n  RESULT: Distance MATTERS (θ = %d km)\n", theta_hat))
} else {
  cat("\n  *** RESULT: Distance does NOT matter (θ = ∞) ***\n")
  cat("  → Cross-district contagion operates uniformly across space\n")
}

# =============================================================================
# FINAL MODEL SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FINAL MODEL SUMMARY                                        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Model: E[Y_tr] = exp(X'β + α × Cross_tr)\n\n")

cat("Background Rate Coefficients:\n")
print(round(summary(m1_best)$coefficients[c("(Intercept)", "log_pop", "poverty_rate", "cpi"), ], 5))

cat("\nCross-District Contagion:\n")
cat(sprintf("  α = %.4f (SE: %.4f, t = %.2f)\n",
            coef(m1_best)["cross_lag"],
            summary(m1_best)$coefficients["cross_lag", "Std. Error"],
            coef(m1_best)["cross_lag"] / summary(m1_best)$coefficients["cross_lag", "Std. Error"]))
cat(sprintf("  exp(α) = %.2f\n", exp(coef(m1_best)["cross_lag"])))
cat(sprintf("  Interpretation: A unit increase in cross-district exposure\n"))
cat(sprintf("                  multiplies expected protest count by %.1f\n", exp(coef(m1_best)["cross_lag"])))

cat(sprintf("\nOptimal spatial decay: θ = %s\n", ifelse(is.infinite(theta_hat), "∞ (uniform)", theta_hat)))

# =============================================================================
# LAG STRUCTURE ANALYSIS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  LAG STRUCTURE: WHICH DAYS MATTER?                          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Compute weekly lag exposures (4 weeks)
W_uniform <- (dist_matrix > 0) * 1
W_uniform_norm <- W_uniform / pmax(rowSums(W_uniform), 1e-10)

# Week 1: lags 1-7, Week 2: lags 8-14, Week 3: lags 15-21, Week 4: lags 22-30
week_exposures <- list(
  week1 = matrix(0, n_days, n_dist),
  week2 = matrix(0, n_days, n_dist),
  week3 = matrix(0, n_days, n_dist),
  week4 = matrix(0, n_days, n_dist)
)

for (lag in 1:30) {
  lagged <- rbind(matrix(0, lag, n_dist), count_matrix[1:(n_days - lag), ])
  lag_exp <- t(W_uniform_norm %*% t(lagged))
  if (lag <= 7) {
    week_exposures$week1 <- week_exposures$week1 + lag_exp
  } else if (lag <= 14) {
    week_exposures$week2 <- week_exposures$week2 + lag_exp
  } else if (lag <= 21) {
    week_exposures$week3 <- week_exposures$week3 + lag_exp
  } else {
    week_exposures$week4 <- week_exposures$week4 + lag_exp
  }
}

lag_df <- data.frame(
  date = rep(dates_vec, n_dist),
  admin2 = rep(colnames(count_matrix), each = n_days),
  week1 = as.vector(week_exposures$week1),
  week2 = as.vector(week_exposures$week2),
  week3 = as.vector(week_exposures$week3),
  week4 = as.vector(week_exposures$week4)
)

panel_lags <- panel_clean %>%
  left_join(lag_df, by = c("admin2", "date"))

m_lags <- glm(count ~ log_pop + poverty_rate + cpi + factor(year) +
                week1 + week2 + week3 + week4,
              data = panel_lags, family = quasipoisson(link = "log"))

cat("Cross-district effect by week:\n\n")
lag_coefs <- summary(m_lags)$coefficients[paste0("week", 1:4), ]
print(round(lag_coefs, 4))

cat(sprintf("\nSum of weekly coefficients (total effect): %.4f\n", sum(lag_coefs[, 1])))

# =============================================================================
# SAVE RESULTS
# =============================================================================

results <- list(
  m0 = list(
    deviance = dev_m0,
    coefficients = coef(m0)
  ),
  profile_likelihood = results_grid,
  optimal_theta = theta_hat,
  test_contagion = list(
    dev_m0 = dev_m0,
    dev_m1 = dev_m1,
    F_stat = F_stat_contagion,
    p_value = p_contagion,
    contagion_exists = p_contagion < 0.05
  ),
  test_distance = list(
    dev_inf = dev_inf,
    dev_best = dev_best,
    p_value = p_distance,
    distance_matters = distance_matters
  ),
  final_model = list(
    coefficients = coef(m1_best),
    se = summary(m1_best)$coefficients[, "Std. Error"],
    dispersion = dispersion
  ),
  lag_structure = list(
    coefficients = lag_coefs[, 1],
    se = lag_coefs[, 2],
    n_lags = 30
  )
)

saveRDS(results, "spatial_hawkes_final_results.rds")
cat("\nResults saved to: spatial_hawkes_final_results.rds\n")

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SUMMARY                                                    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("┌─────────────────────────────────────────────────────────────┐\n")
cat("│ Question                          │ Result                 │\n")
cat("├─────────────────────────────────────────────────────────────┤\n")
cat(sprintf("│ Cross-district contagion exists?  │ %-22s │\n",
            ifelse(p_contagion < 0.05, "YES (p < 0.001)", "NO")))
cat(sprintf("│ Distance matters?                 │ %-22s │\n",
            ifelse(distance_matters, paste0("YES (θ = ", theta_hat, " km)"), "NO (θ = ∞)")))
cat(sprintf("│ Cross-district effect (α)         │ %-22s │\n",
            sprintf("%.3f (SE: %.3f)", coef(m1_best)["cross_lag"],
                    summary(m1_best)$coefficients["cross_lag", "Std. Error"])))
cat(sprintf("│ Multiplicative effect exp(α)      │ %-22s │\n",
            sprintf("%.1fx", exp(coef(m1_best)["cross_lag"]))))
cat("└─────────────────────────────────────────────────────────────┘\n")
