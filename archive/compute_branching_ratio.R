################################################################################
#       COMPUTE BRANCHING RATIO AND ENDOGENOUS FRACTION
#
#       For log-link Hawkes model: E[Y] = exp(X'β + α × Cross)
################################################################################

library(tidyverse)
library(geosphere)

select <- dplyr::select

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   BRANCHING RATIO & ENDOGENOUS FRACTION                                 ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

protests <- readRDS("protests_daily.rds")
cpi_data <- read_csv("indonesia_cpi_processed.csv", show_col_types = FALSE)
results <- readRDS("spatial_hawkes_final_results.rds")

theta_hat <- results$optimal_theta
alpha <- results$final_model$coefficients["cross_lag"]

cat(sprintf("Optimal θ: %s km\n", ifelse(is.infinite(theta_hat), "∞", theta_hat)))
cat(sprintf("Cross-district effect α: %.4f\n", alpha))
cat(sprintf("exp(α) = %.3f\n\n", exp(alpha)))

# =============================================================================
# REBUILD PANEL AND EXPOSURE
# =============================================================================

cat("Building panel data...\n")

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

cpi_monthly <- cpi_data %>%
  select(year_month, cpi) %>%
  mutate(year_month = as.character(year_month))

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

panel_clean <- panel %>%
  filter(!is.na(log_pop), !is.na(poverty_rate), !is.na(cpi))

cat("Computing distance matrix...\n")

# Distance matrix
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

# Count matrix
panel_wide <- panel %>%
  select(admin2, date, count) %>%
  pivot_wider(names_from = admin2, values_from = count, values_fill = 0)

dates_vec <- panel_wide$date
count_matrix <- as.matrix(panel_wide[, -1])
count_matrix <- count_matrix[, rownames(dist_matrix)]

n_days <- nrow(count_matrix)
n_dist <- ncol(count_matrix)

# Spatial weight matrix
if (is.infinite(theta_hat)) {
  W <- (dist_matrix > 0) * 1
} else {
  W <- exp(-dist_matrix / theta_hat)
  diag(W) <- 0
}
W_norm <- W / pmax(rowSums(W), 1e-10)

# Cross-district exposure (sum over lags 1-30)
cat("Computing cross-district exposure...\n")

cross_exposure <- matrix(0, n_days, n_dist)
for (lag in 1:30) {
  lagged <- rbind(matrix(0, lag, n_dist), count_matrix[1:(n_days - lag), ])
  cross_exposure <- cross_exposure + t(W_norm %*% t(lagged))
}

exp_long <- data.frame(
  date = rep(dates_vec, n_dist),
  admin2 = rep(colnames(count_matrix), each = n_days),
  cross_lag = as.vector(cross_exposure)
)

panel_fit <- panel_clean %>%
  left_join(exp_long, by = c("admin2", "date"))

# =============================================================================
# FIT MODELS
# =============================================================================

cat("Fitting models...\n\n")

# Full model with contagion
m1 <- glm(count ~ log_pop + poverty_rate + cpi + factor(year) + cross_lag,
          data = panel_fit, family = quasipoisson(link = "log"))

# =============================================================================
# ENDOGENOUS FRACTION VIA COUNTERFACTUAL
# =============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("  ENDOGENOUS FRACTION: Counterfactual Decomposition\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# For log-link: λ = exp(X'β + α×Cross) = exp(X'β) × exp(α×Cross)
# Background intensity: λ_bg = exp(X'β)
# Triggered multiplier: exp(α×Cross)
# Triggered intensity: λ × (1 - exp(-α×Cross)) ≈ λ - λ_bg for small α×Cross

# Get linear predictor WITHOUT cross_lag
X_bg <- model.matrix(~ log_pop + poverty_rate + cpi + factor(year), data = panel_fit)
coefs_bg <- coef(m1)[colnames(X_bg)]
eta_bg <- X_bg %*% coefs_bg
lambda_bg <- exp(eta_bg)

# Full fitted values
lambda_full <- predict(m1, type = "response")

# Triggered intensity at each district-day
# λ_triggered = λ_full - λ_bg = exp(X'β) × (exp(α×Cross) - 1)
lambda_triggered <- lambda_full - lambda_bg

# Sum to get totals
total_observed <- sum(panel_fit$count)
total_expected <- sum(lambda_full)
total_background <- sum(lambda_bg)
total_triggered <- sum(lambda_triggered)

endogenous_fraction <- total_triggered / total_expected

cat(sprintf("Total observed events: %d\n", total_observed))
cat(sprintf("Total expected (full model): %.0f\n", total_expected))
cat(sprintf("Background component: %.0f\n", total_background))
cat(sprintf("Triggered component: %.0f\n", total_triggered))
cat(sprintf("\nENDOGENOUS FRACTION: %.1f%%\n", 100 * endogenous_fraction))

# =============================================================================
# BRANCHING RATIO
# =============================================================================

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  BRANCHING RATIO: Expected Offspring per Event\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# For one event in district r at time t:
# - Contributes to Cross_{t+ℓ,s} for all s≠r and ℓ ∈ {1,...,7}
# - Contribution amount: W_norm[s,r] (weight from r to s)
# - This multiplies expected count at (t+ℓ,s) by exp(α × W_norm[s,r])
# - Extra events at (t+ℓ,s): λ_{t+ℓ,s} × (exp(α × W_norm[s,r]) - 1)

# To compute BR, we average over all events:
# BR = mean over events of [Σ_s Σ_ℓ λ_s × (exp(α × W_norm[s,r]) - 1)]

# For efficiency, use the average λ per district
avg_lambda_by_district <- panel_fit %>%
  group_by(admin2) %>%
  summarise(avg_lambda = mean(predict(m1, newdata = cur_data(), type = "response")),
            .groups = "drop")

# Match to column order
avg_lambda_vec <- avg_lambda_by_district$avg_lambda[match(colnames(count_matrix),
                                                          avg_lambda_by_district$admin2)]
names(avg_lambda_vec) <- colnames(count_matrix)

# For each origin district r, compute total offspring
offspring_by_origin <- numeric(n_dist)

for (r in 1:n_dist) {
  total <- 0
  for (s in 1:n_dist) {
    if (s == r) next
    w <- W_norm[s, r]
    # Offspring at s from one event at r (summed over 30 lags)
    offspring_s <- 30 * avg_lambda_vec[s] * (exp(alpha * w) - 1)
    total <- total + offspring_s
  }
  offspring_by_origin[r] <- total
}

# Weight by how many events come from each district
events_by_district <- panel_fit %>%
  group_by(admin2) %>%
  summarise(n_events = sum(count), .groups = "drop")

event_weights <- events_by_district$n_events[match(colnames(count_matrix),
                                                    events_by_district$admin2)]
event_weights <- event_weights / sum(event_weights)

# Weighted average branching ratio
branching_ratio <- sum(offspring_by_origin * event_weights)

cat(sprintf("Mean offspring per event: %.4f\n", branching_ratio))
cat(sprintf("Min offspring (by district): %.4f\n", min(offspring_by_origin)))
cat(sprintf("Max offspring (by district): %.4f\n", max(offspring_by_origin)))
cat(sprintf("\nBRANCHING RATIO: %.4f\n", branching_ratio))

# Alternative: from endogenous fraction
# In a branching process: endogenous_frac = BR / (1 + BR)
# So: BR = endogenous_frac / (1 - endogenous_frac)
br_from_endog <- endogenous_fraction / (1 - endogenous_fraction)
cat(sprintf("\nBR implied by endogenous fraction: %.4f\n", br_from_endog))

# =============================================================================
# INTERPRETATION
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  INTERPRETATION\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("The cross-district contagion effect:\n")
cat(sprintf("  α = %.3f means a unit increase in spatial exposure\n", alpha))
cat(sprintf("  multiplies expected protests by %.1f×\n\n", exp(alpha)))

cat("However, the actual exposure values are small:\n")
cat(sprintf("  Mean cross-district exposure: %.4f\n", mean(panel_fit$cross_lag)))
cat(sprintf("  Max cross-district exposure: %.4f\n", max(panel_fit$cross_lag)))
cat(sprintf("  Typical multiplicative effect: exp(%.3f × %.3f) = %.3f×\n",
            alpha, mean(panel_fit$cross_lag[panel_fit$cross_lag > 0]),
            exp(alpha * mean(panel_fit$cross_lag[panel_fit$cross_lag > 0]))))

cat("\nThis explains why:\n")
cat(sprintf("  - α is large (%.2f) but BR is small (%.4f)\n", alpha, branching_ratio))
cat("  - The spatial weights are row-normalized, so one event contributes\n")
cat(sprintf("    only ~1/%d = %.4f to each neighbor's exposure\n", n_dist - 1, 1/(n_dist-1)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SUMMARY                                                    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("┌─────────────────────────────────────────────────────────────┐\n")
cat("│ Quantity                          │ Value                  │\n")
cat("├─────────────────────────────────────────────────────────────┤\n")
cat(sprintf("│ Cross-district effect (α)         │ %.4f                 │\n", alpha))
cat(sprintf("│ Multiplicative effect exp(α)      │ %.2fx                  │\n", exp(alpha)))
cat(sprintf("│ Optimal decay distance (θ)        │ %s km                │\n",
            ifelse(is.infinite(theta_hat), "∞", theta_hat)))
cat(sprintf("│ Endogenous fraction               │ %.1f%%                  │\n", 100 * endogenous_fraction))
cat(sprintf("│ Branching ratio                   │ %.4f                 │\n", branching_ratio))
cat(sprintf("│ Process type                      │ %-22s │\n",
            ifelse(branching_ratio < 1, "SUBCRITICAL", "CRITICAL")))
cat("└─────────────────────────────────────────────────────────────┘\n")

cat("\n")
cat("KEY FINDINGS:\n")
cat(sprintf("  • %.1f%% of expected protests are triggered by cross-district contagion\n",
            100 * endogenous_fraction))
cat(sprintf("  • Each protest triggers an average of %.3f additional protests elsewhere\n",
            branching_ratio))
cat(sprintf("  • The process is SUBCRITICAL (BR = %.3f < 1): protests spread but don't explode\n",
            branching_ratio))

# Save results
br_results <- list(
  alpha = alpha,
  exp_alpha = exp(alpha),
  theta = theta_hat,
  endogenous_fraction = endogenous_fraction,
  branching_ratio = branching_ratio,
  total_observed = total_observed,
  total_triggered = total_triggered
)

saveRDS(br_results, "branching_ratio_results.rds")
cat("\nResults saved to: branching_ratio_results.rds\n")
