################################################################################
#       SPATIAL HAWKES MODEL WITH CONTINUOUS DISTANCE KERNEL
#
#       Uses exponential decay kernel: f(d) = exp(-d/θ)
#       Estimates optimal θ via profile likelihood
#
#       Test: H0: θ = ∞ (no spatial decay) vs H1: θ < ∞ (decay exists)
################################################################################

library(tidyverse)
library(geosphere)
library(MASS)

select <- dplyr::select

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   SPATIAL HAWKES MODEL WITH EXPONENTIAL KERNEL                          ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

protests <- readRDS("protests_daily.rds")
cpi_data <- read_csv("indonesia_cpi_processed.csv", show_col_types = FALSE)

cat(sprintf("Total events: %d\n", nrow(protests)))
cat(sprintf("Date range: %s to %s\n", min(protests$date), max(protests$date)))

# =============================================================================
# PREPARE DISTRICT-LEVEL DATA
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

cat(sprintf("Districts with complete data: %d\n", nrow(district_covariates)))

cpi_monthly <- cpi_data %>%
  select(year_month, cpi) %>%
  mutate(year_month = as.character(year_month))

# =============================================================================
# CREATE DISTRICT-DAY PANEL
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

# =============================================================================
# COMPUTE DISTANCE MATRIX
# =============================================================================

cat("\nComputing distance matrix...\n")

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

cat(sprintf("Mean distance: %.0f km, Median: %.0f km\n",
            mean(dist_matrix[dist_matrix > 0]),
            median(dist_matrix[dist_matrix > 0])))

# =============================================================================
# CREATE COUNT MATRIX (dates x districts)
# =============================================================================

panel_wide <- panel %>%
  select(admin2, date, count) %>%
  pivot_wider(names_from = admin2, values_from = count, values_fill = 0)

dates_vec <- panel_wide$date
count_matrix <- as.matrix(panel_wide[, -1])
count_matrix <- count_matrix[, rownames(dist_matrix)]  # Align with distance matrix

n_days <- nrow(count_matrix)
n_dist <- ncol(count_matrix)

# =============================================================================
# COMPUTE OWN-DISTRICT LAGGED EXPOSURE (always the same)
# =============================================================================

cat("Computing own-district lagged exposure...\n")

own_exposure <- matrix(0, n_days, n_dist)
for (lag in 1:7) {
  lagged <- rbind(matrix(0, lag, n_dist), count_matrix[1:(n_days - lag), ])
  own_exposure <- own_exposure + lagged
}

# =============================================================================
# FUNCTION: Compute cross-district exposure for given theta
# =============================================================================

compute_cross_exposure <- function(theta, dist_matrix, count_matrix, n_days, n_dist) {
  # Compute kernel weights: W[i,j] = exp(-d[i,j] / theta)
  # Row-normalize so weights sum to 1 for each district

  if (is.infinite(theta)) {
    # theta = Inf means uniform weights (no decay)
    W <- (dist_matrix > 0) * 1  # Binary: 1 if different district, 0 otherwise
  } else {
    W <- exp(-dist_matrix / theta)
    diag(W) <- 0  # No self-influence in cross-district
  }

  # Row-normalize
  row_sums <- rowSums(W)
  W_norm <- W / pmax(row_sums, 1e-10)

  # Compute cross-district exposure: sum over lags 1-7 of weighted neighbor counts
  cross_exposure <- matrix(0, n_days, n_dist)

  for (lag in 1:7) {
    lagged <- rbind(matrix(0, lag, n_dist), count_matrix[1:(n_days - lag), ])
    # For each day, multiply lagged counts by weight matrix
    cross_exposure <- cross_exposure + t(W_norm %*% t(lagged))
  }

  return(cross_exposure)
}

# =============================================================================
# PROFILE LIKELIHOOD OVER THETA
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PROFILE LIKELIHOOD: ESTIMATING OPTIMAL THETA               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Grid of theta values (km)
theta_grid <- c(50, 100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000, Inf)

results_grid <- data.frame(
  theta = theta_grid,
  loglik = NA,
  deviance = NA,
  aic = NA,
  alpha_own = NA,
  alpha_cross = NA,
  se_cross = NA
)

cat("Fitting models for each theta...\n\n")

for (i in seq_along(theta_grid)) {
  theta <- theta_grid[i]
  theta_label <- if (is.infinite(theta)) "Inf" else as.character(theta)

  cat(sprintf("  theta = %s km: ", theta_label))

  # Compute cross-district exposure
  cross_exposure <- compute_cross_exposure(theta, dist_matrix, count_matrix, n_days, n_dist)

  # Convert to long format
  exposure_long <- data.frame(
    date = rep(dates_vec, n_dist),
    admin2 = rep(colnames(count_matrix), each = n_days),
    own_lag = as.vector(own_exposure),
    cross_lag = as.vector(cross_exposure)
  )

  # Merge with panel
  panel_fit <- panel %>%
    left_join(exposure_long, by = c("admin2", "date")) %>%
    filter(!is.na(log_pop), !is.na(poverty_rate), !is.na(cpi))

  # Fit model
  m <- glm(count ~ log_pop + poverty_rate + cpi + factor(year) + own_lag + cross_lag,
           data = panel_fit, family = quasipoisson(link = "log"))

  # Store results
  results_grid$loglik[i] <- logLik(m)
  results_grid$deviance[i] <- deviance(m)
  results_grid$aic[i] <- AIC(m)
  results_grid$alpha_own[i] <- coef(m)["own_lag"]
  results_grid$alpha_cross[i] <- coef(m)["cross_lag"]
  results_grid$se_cross[i] <- summary(m)$coefficients["cross_lag", "Std. Error"]

  cat(sprintf("deviance = %.1f, α_cross = %.4f\n", deviance(m), coef(m)["cross_lag"]))
}

# Find optimal theta (minimum deviance)
best_idx <- which.min(results_grid$deviance)
theta_hat <- results_grid$theta[best_idx]
inf_idx <- which(is.infinite(results_grid$theta))

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("PROFILE LIKELIHOOD RESULTS:\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

print(results_grid %>% mutate(theta = ifelse(is.infinite(theta), "Inf", theta)))

cat(sprintf("\nOptimal theta: %s km\n",
            ifelse(is.infinite(theta_hat), "Inf (no decay)", as.character(theta_hat))))

# =============================================================================
# LIKELIHOOD RATIO TEST: theta = Inf vs theta = theta_hat
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  TEST: DOES DISTANCE MATTER?                                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("H0: θ = ∞ (no spatial decay, distance doesn't matter)\n")
cat("H1: θ < ∞ (spatial decay exists)\n\n")

dev_null <- results_grid$deviance[inf_idx]  # theta = Inf
dev_alt <- results_grid$deviance[best_idx]  # theta = optimal

# Quasi-likelihood ratio test
# Under H0, 2*(logL_alt - logL_null) / dispersion ~ chi-squared(1)
# For deviance: (dev_null - dev_alt) / dispersion ~ chi-squared(1)

# Get dispersion from best model
cross_exp_best <- compute_cross_exposure(theta_hat, dist_matrix, count_matrix, n_days, n_dist)
exp_long_best <- data.frame(
  date = rep(dates_vec, n_dist),
  admin2 = rep(colnames(count_matrix), each = n_days),
  own_lag = as.vector(own_exposure),
  cross_lag = as.vector(cross_exp_best)
)
panel_best <- panel %>%
  left_join(exp_long_best, by = c("admin2", "date")) %>%
  filter(!is.na(log_pop), !is.na(poverty_rate), !is.na(cpi))
m_best <- glm(count ~ log_pop + poverty_rate + cpi + factor(year) + own_lag + cross_lag,
              data = panel_best, family = quasipoisson(link = "log"))
dispersion <- summary(m_best)$dispersion

# LR test statistic
LR_stat <- (dev_null - dev_alt) / dispersion
p_value <- pchisq(LR_stat, df = 1, lower.tail = FALSE)

cat(sprintf("Deviance (θ = ∞):        %.2f\n", dev_null))
cat(sprintf("Deviance (θ = %s):    %.2f\n",
            ifelse(is.infinite(theta_hat), "Inf", theta_hat), dev_alt))
cat(sprintf("Dispersion:              %.3f\n", dispersion))
cat(sprintf("LR statistic:            %.2f\n", LR_stat))
cat(sprintf("p-value:                 %.4f\n", p_value))

if (p_value < 0.05) {
  if (!is.infinite(theta_hat)) {
    cat(sprintf("\nRESULT: Distance MATTERS (reject H0)\n"))
    cat(sprintf("  → Optimal decay distance: %d km\n", theta_hat))
    cat(sprintf("  → Spatial contagion decays with distance\n"))
  }
} else {
  cat(sprintf("\nRESULT: Distance does NOT matter (fail to reject H0)\n"))
  cat(sprintf("  → No evidence of spatial decay\n"))
  cat(sprintf("  → Cross-district contagion operates uniformly\n"))
}

# =============================================================================
# INTERPRETATION OF OPTIMAL MODEL
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  OPTIMAL MODEL SUMMARY                                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Coefficients:\n")
print(round(summary(m_best)$coefficients[c("(Intercept)", "log_pop", "poverty_rate",
                                            "cpi", "own_lag", "cross_lag"), ], 5))

cat(sprintf("\nInterpretation (log link, so coefficients are multiplicative):\n"))
cat(sprintf("  Own-district lag:   exp(%.4f) = %.3f (%.1f%% increase per lagged protest)\n",
            coef(m_best)["own_lag"], exp(coef(m_best)["own_lag"]),
            100 * (exp(coef(m_best)["own_lag"]) - 1)))
cat(sprintf("  Cross-district lag: exp(%.4f) = %.3f (%.1f%% increase per unit exposure)\n",
            coef(m_best)["cross_lag"], exp(coef(m_best)["cross_lag"]),
            100 * (exp(coef(m_best)["cross_lag"]) - 1)))

# =============================================================================
# PLOT: DEVIANCE VS THETA
# =============================================================================

cat("\n")
cat("Generating deviance profile plot...\n")

# For plotting, replace Inf with a large number
plot_data <- results_grid %>%
  mutate(theta_plot = ifelse(is.infinite(theta), 10000, theta),
         theta_label = ifelse(is.infinite(theta), "∞", as.character(theta)))

p <- ggplot(plot_data, aes(x = theta_plot, y = deviance)) +
  geom_line(color = "#3498db", linewidth = 1) +
  geom_point(color = "#3498db", size = 3) +
  geom_point(data = plot_data[best_idx, ], color = "#e74c3c", size = 5) +
  geom_vline(xintercept = plot_data$theta_plot[best_idx], color = "#e74c3c",
             linetype = "dashed", linewidth = 0.8) +
  scale_x_log10(breaks = c(50, 100, 200, 500, 1000, 2000, 5000, 10000),
                labels = c("50", "100", "200", "500", "1000", "2000", "5000", "∞")) +
  labs(
    x = "Decay Distance θ (km)",
    y = "Deviance",
    title = "Profile Likelihood for Spatial Decay Parameter",
    subtitle = sprintf("Optimal θ = %s km (red point)",
                       ifelse(is.infinite(theta_hat), "∞", theta_hat))
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank())

ggsave("figures/fig_profile_likelihood.pdf", p, width = 8, height = 5)
cat("Saved: figures/fig_profile_likelihood.pdf\n")

# =============================================================================
# SAVE RESULTS
# =============================================================================

results <- list(
  profile_likelihood = results_grid,
  optimal_theta = theta_hat,
  lr_test = list(
    dev_null = dev_null,
    dev_alt = dev_alt,
    dispersion = dispersion,
    lr_stat = LR_stat,
    p_value = p_value,
    distance_matters = p_value < 0.05
  ),
  optimal_model = list(
    coefficients = coef(m_best),
    se = summary(m_best)$coefficients[, "Std. Error"],
    dispersion = dispersion
  )
)

saveRDS(results, "spatial_kernel_results.rds")
cat("\nResults saved to: spatial_kernel_results.rds\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FINAL SUMMARY                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("MODEL: E[Y_tr] = exp(Xβ + α_own·Own_tr + α_cross·Cross_tr(θ))\n\n")

cat("Where Cross_tr(θ) = Σ_s [exp(-d_rs/θ) / Σ_j exp(-d_rj/θ)] × (7-day lag of Y_s)\n\n")

cat("FINDINGS:\n")
cat(sprintf("  1. Optimal decay distance: θ = %s km\n",
            ifelse(is.infinite(theta_hat), "∞ (no decay)", theta_hat)))
cat(sprintf("  2. Distance matters? %s (p = %.4f)\n",
            ifelse(p_value < 0.05, "YES", "NO"), p_value))
cat(sprintf("  3. Own-district effect: α = %.4f (SE: %.4f)\n",
            coef(m_best)["own_lag"], summary(m_best)$coefficients["own_lag", 2]))
cat(sprintf("  4. Cross-district effect: α = %.4f (SE: %.4f)\n",
            coef(m_best)["cross_lag"], summary(m_best)$coefficients["cross_lag", 2]))
