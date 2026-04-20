################################################################################
#                    CROSS-DAY HAWKES MODEL ESTIMATION
#
#   Alternative to jittered univariate Hawkes that eliminates same-day artifacts
#
#   Key differences from estimate_all_models.R:
#     1. Uses discrete daily time (no jitter)
#     2. Events can only trigger events on FUTURE days (not same day)
#     3. Eliminates the jitter artifact that caused branching ratio ~ 1.0
#
#   Models:
#     Model 0: Poisson (Background Only)
#     Model 1: Cross-Day Hawkes (previous days trigger, same day excluded)
#
################################################################################

library(tidyverse)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║     CROSS-DAY HAWKES: DISCRETE DAILY TRIGGERING (NO JITTER)             ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

cat("KEY CHANGE: Events on day d can only be triggered by events on days < d\n")
cat("           Same-day events CANNOT trigger each other\n")
cat("           This eliminates the jitter artifact\n\n")

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

N_STARTS <- 10       # Number of starting points per model
MAX_ITER <- 2000     # Max iterations for optimization

# =============================================================================
# 2. LOAD AND PREPARE DATA (DISCRETE DAILY)
# =============================================================================

cat("=== LOADING DATA ===\n\n")

protests <- readRDS("protests_daily.rds") %>%
  filter(!is.na(poverty_decimal) & !is.na(log_pop) & !is.na(log_cpi)) %>%
  mutate(
    type_numeric = ifelse(is_severe, 1L, 2L)  # 1=Severe, 2=Peaceful
  ) %>%
  arrange(date)

# Compute day index (integer days from start) - NO JITTER
start_date <- min(protests$date)
protests$day <- as.numeric(protests$date - start_date)

# Sort by day (within-day order doesn't matter since same-day can't trigger)
protests <- protests %>% arrange(day)

D_max <- max(protests$day) + 1  # Total days

n_total <- nrow(protests)
n_severe <- sum(protests$is_severe)
n_peaceful <- sum(!protests$is_severe)
n_days_with_events <- length(unique(protests$day))

cat(sprintf("Events: %d total (%d severe, %d peaceful)\n", n_total, n_severe, n_peaceful))
cat(sprintf("Days: %d total, %d with events (%.1f%%)\n",
            D_max, n_days_with_events, 100 * n_days_with_events / D_max))
cat(sprintf("Date range: %s to %s\n\n", min(protests$date), max(protests$date)))

# Prepare data vectors
days <- protests$day
log_pop <- protests$log_pop
poverty <- protests$poverty_decimal
log_cpi <- protests$log_cpi
years <- protests$year

# Compute day differences (for decay computation)
day_diffs <- c(0, diff(days))

# =============================================================================
# 3. BUILD FULL DISTRICT-MONTH GRID FOR INTEGRAL
# =============================================================================

cat("=== BUILDING FULL INTEGRATION GRID ===\n\n")

# Load population data (all 514 districts)
pop_raw <- read.csv("district_level_population_2014_2020.csv", stringsAsFactors = FALSE)
pop_data <- pop_raw %>%
  filter(Series.Name == "Total Population (in number of people)") %>%
  select(district = Provinces.Name,
         `2014` = X2014..YR2014., `2015` = X2015..YR2015., `2016` = X2016..YR2016.,
         `2017` = X2017..YR2017., `2018` = X2018..YR2018., `2019` = X2019..YR2019.,
         `2020` = X2020..YR2020.) %>%
  mutate(district = gsub('^"', '', district)) %>%
  pivot_longer(-district, names_to = "year", values_to = "population") %>%
  mutate(year = as.integer(year), population = as.numeric(population)) %>%
  filter(!is.na(population))

# Extrapolate to 2021-2024
pop_extrapolated <- pop_data %>%
  group_by(district) %>%
  do({
    df <- .
    if (nrow(df) >= 3) {
      model <- lm(population ~ year, data = df)
      future <- data.frame(year = 2021:2024)
      future$population <- pmax(predict(model, future), min(df$population))
      future$district <- df$district[1]
      rbind(df, future[, c("district", "year", "population")])
    } else {
      df
    }
  }) %>%
  ungroup() %>%
  mutate(log_pop = log(population))

# Load poverty data
poverty_data <- readRDS("poverty_rate_2015_2024.rds") %>%
  mutate(poverty_decimal = poverty_rate / 100)

# Load CPI data
cpi_data <- readRDS("indonesia_cpi.rds") %>%
  mutate(log_cpi = log(cpi)) %>%
  select(year_month, log_cpi)

# Create full grid
all_districts <- unique(pop_extrapolated$district)
all_months <- seq.Date(min(protests$date), max(protests$date), by = "month")
all_year_months <- format(all_months, "%Y-%m")

full_grid <- expand.grid(
  district = all_districts,
  year_month = all_year_months,
  stringsAsFactors = FALSE
) %>%
  mutate(
    year = as.integer(substr(year_month, 1, 4)),
    month = as.integer(substr(year_month, 6, 7)),
    days_in_month = case_when(
      month %in% c(1,3,5,7,8,10,12) ~ 31,
      month %in% c(4,6,9,11) ~ 30,
      month == 2 & year %% 4 == 0 ~ 29,
      TRUE ~ 28
    )
  )

# Merge covariates
district_year_months <- full_grid %>%
  left_join(pop_extrapolated %>% select(district, year, log_pop),
            by = c("district", "year")) %>%
  left_join(poverty_data %>% select(district, year, poverty_decimal),
            by = c("district", "year")) %>%
  left_join(cpi_data, by = "year_month") %>%
  filter(!is.na(log_pop) & !is.na(poverty_decimal) & !is.na(log_cpi))

cat(sprintf("Full integration grid: %d district-months\n", nrow(district_year_months)))
cat(sprintf("  Districts: %d\n", length(unique(district_year_months$district))))
cat(sprintf("  Months: %d\n\n", length(unique(district_year_months$year_month))))

# =============================================================================
# 4. PRECOMPUTE VECTORIZED DATA FOR FAST LIKELIHOOD
# =============================================================================

# Precompute year effects lookup for events
year_effect_idx_events <- pmax(0, years - 2015)

# Precompute year effects lookup for grid
year_effect_idx_grid <- pmax(0, district_year_months$year - 2015)

# Extract grid vectors
grid_log_pop <- district_year_months$log_pop
grid_poverty <- district_year_months$poverty_decimal
grid_log_cpi <- district_year_months$log_cpi
grid_days <- district_year_months$days_in_month

cat("Precomputed vectorized data for fast likelihood evaluation.\n\n")

# =============================================================================
# 5. LIKELIHOOD FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# Model 0: Poisson (Background Only) - Same as before
# -----------------------------------------------------------------------------

poisson_loglik <- function(params) {
  beta_0 <- params[1]
  gamma <- params[2]
  delta <- params[3]
  zeta <- params[4]
  beta_years <- params[5:13]

  # Year effects (vectorized)
  year_effects <- c(0, beta_years)[year_effect_idx_events + 1]

  # Background rate at each event
  mu <- exp(beta_0 + gamma * log_pop + delta * poverty + zeta * log_cpi + year_effects)

  # Log-likelihood: sum of log(mu) at event times
  loglik <- sum(log(pmax(mu, 1e-10)))

  # Background integral (vectorized)
  year_effects_grid <- c(0, beta_years)[year_effect_idx_grid + 1]
  mu_grid <- exp(beta_0 + gamma * grid_log_pop + delta * grid_poverty + zeta * grid_log_cpi + year_effects_grid)
  integral <- sum(mu_grid * grid_days)

  loglik <- loglik - integral
  if (!is.finite(loglik)) return(-1e10)
  return(loglik)
}

# -----------------------------------------------------------------------------
# Model 1: Cross-Day Hawkes (Discrete Daily - No Same-Day Triggering)
# Parameters: β₀, γ, δ, ζ, β_2016...β_2024, log_α, log_ω (15 total)
# where ω = exp(-β) is the daily decay factor
# -----------------------------------------------------------------------------

crossday_hawkes_loglik <- function(params) {
  beta_0 <- params[1]
  gamma <- params[2]
  delta <- params[3]
  zeta <- params[4]
  beta_years <- params[5:13]
  alpha <- exp(params[14])       # Triggering intensity
  omega <- exp(params[15])       # Daily decay factor ω = exp(-β_decay)

  # Constraint: ω must be in (0, 1) for proper decay
  if (omega <= 0 || omega >= 1) return(-1e10)

  # Convert to decay rate for interpretation
  beta_decay <- -log(omega)  # β = -log(ω)

  n <- length(days)

  # Year effects (vectorized)
  year_effects <- c(0, beta_years)[year_effect_idx_events + 1]

  # Background rate at each event
  mu <- exp(beta_0 + gamma * log_pop + delta * poverty + zeta * log_cpi + year_effects)

  # Cross-day recursive triggering sum
  # R[i] = triggering pressure from events on PREVIOUS days only
  R <- numeric(n)
  for (i in 2:n) {
    if (day_diffs[i] > 0) {
      # New day: previous event(s) can now contribute
      # Apply decay for the day gap
      decay <- omega^day_diffs[i]
      R[i] <- decay * (1 + R[i-1])
    } else {
      # Same day: just carry forward (same-day events don't trigger each other)
      R[i] <- R[i-1]
    }
  }

  # Log-likelihood at event times (vectorized)
  lambda <- mu + alpha * R
  if (any(lambda <= 1e-10)) return(-1e10)
  loglik <- sum(log(lambda))

  # Background integral (same as Poisson)
  year_effects_grid <- c(0, beta_years)[year_effect_idx_grid + 1]
  mu_grid <- exp(beta_0 + gamma * grid_log_pop + delta * grid_poverty + zeta * grid_log_cpi + year_effects_grid)
  integral_bg <- sum(mu_grid * grid_days)

  # Triggering integral (discrete daily)
  # Each event on day d contributes to days d+1, d+2, ..., D_max
  # Integral = Σ_j α × Σ_{s=1}^{D_max - d_j} ω^s
  #          = Σ_j α × ω × (1 - ω^{D_max - d_j}) / (1 - ω)
  days_remaining <- D_max - days
  integral_trig <- sum(alpha * omega * (1 - omega^days_remaining) / (1 - omega))

  loglik <- loglik - integral_bg - integral_trig
  if (!is.finite(loglik)) return(-1e10)
  return(loglik)
}

# =============================================================================
# 6. FIT MODEL 0: POISSON
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL 0: POISSON (Background Only)                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Empirical rate
total_district_days <- sum(district_year_months$days_in_month)
emp_rate_per_district_day <- n_total / total_district_days
cat(sprintf("Empirical rate: %.6f per district-day\n\n", emp_rate_per_district_day))

# Mean covariates from grid
mean_log_pop_grid <- mean(district_year_months$log_pop)
mean_poverty_grid <- mean(district_year_months$poverty_decimal)
mean_log_cpi_grid <- mean(district_year_months$log_cpi)

# Starting values
beta0_start <- log(emp_rate_per_district_day)

start_m0 <- list(
  c(beta0_start - 0.1 * mean_log_pop_grid, 0.1, 0.0, 0.0, rep(0, 9)),
  c(beta0_start - 0.2 * mean_log_pop_grid, 0.2, -1.0, 0.0, rep(0, 9)),
  c(beta0_start - 0.05 * mean_log_pop_grid, 0.05, 1.0, -0.5, rep(0, 9))
)

# Bounds
lower_m0 <- c(-20, -3, -10, -10, rep(-5, 9))
upper_m0 <- c(10, 3, 10, 10, rep(5, 9))

best_ll_m0 <- -Inf
best_result_m0 <- NULL

cat("Running multi-start optimization...\n")
for (s in seq_along(start_m0)) {
  cat(sprintf("  Start %d/%d... ", s, length(start_m0)))

  result <- tryCatch({
    optim(
      par = start_m0[[s]],
      fn = function(p) -poisson_loglik(p),
      method = "L-BFGS-B",
      lower = lower_m0,
      upper = upper_m0,
      control = list(maxit = MAX_ITER)
    )
  }, error = function(e) NULL)

  if (!is.null(result)) {
    ll <- -result$value
    cat(sprintf("LL=%.1f, conv=%d\n", ll, result$convergence))
    if (ll > best_ll_m0) {
      best_ll_m0 <- ll
      best_result_m0 <- result
    }
  } else {
    cat("FAILED\n")
  }
}

params_m0 <- best_result_m0$par
ll_m0 <- -best_result_m0$value

cat(sprintf("\nModel 0 Results:\n"))
cat(sprintf("  Log-likelihood: %.2f\n", ll_m0))
cat(sprintf("  β₀ = %.4f\n", params_m0[1]))
cat(sprintf("  γ (population) = %.4f\n", params_m0[2]))
cat(sprintf("  δ (poverty) = %.4f\n", params_m0[3]))
cat(sprintf("  ζ (log_cpi) = %.4f\n\n", params_m0[4]))

# =============================================================================
# 7. FIT MODEL 1: CROSS-DAY HAWKES
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL 1: CROSS-DAY HAWKES (Discrete Daily, No Same-Day)    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Starting points (chain from M0)
# params[14] = log(α), params[15] = log(ω) where ω ∈ (0,1)
# For half-life of ~3 days: ω = exp(-log(2)/3) ≈ 0.79
# For half-life of ~7 days: ω = exp(-log(2)/7) ≈ 0.91

start_m1 <- list(
  c(params_m0[1:13], log(0.1), log(0.8)),    # half-life ~3 days
  c(params_m0[1:13], log(0.05), log(0.7)),   # half-life ~2 days
  c(params_m0[1:13], log(0.2), log(0.9)),    # half-life ~7 days
  c(params_m0[1:13], log(0.02), log(0.85)),  # weak, ~5 days
  c(params_m0[1:13], log(0.15), log(0.75)),  # stronger, ~2.5 days
  c(params_m0[1:13], log(0.3), log(0.6)),    # strong, ~1.5 days
  c(params_m0[1:13], log(0.08), log(0.88)),  # variation
  c(params_m0[1:13], log(0.12), log(0.82)),  # variation
  c(params_m0[1:13], log(0.25), log(0.95)),  # slow decay
  c(params_m0[1:13], log(0.04), log(0.5))    # fast decay
)

# Bounds:
# log(α) ∈ [-10, 3] → α ∈ [0.00005, 20]
# log(ω) ∈ [-3, -0.001] → ω ∈ [0.05, 0.999] → half-life ∈ [0.23, 693] days
lower_m1 <- c(lower_m0, -10, -3)
upper_m1 <- c(upper_m0, 3, -0.001)

best_ll_m1 <- -Inf
best_result_m1 <- NULL

cat("Running multi-start optimization...\n")
for (s in seq_along(start_m1)) {
  cat(sprintf("  Start %d/%d... ", s, length(start_m1)))

  result <- tryCatch({
    optim(
      par = start_m1[[s]],
      fn = function(p) -crossday_hawkes_loglik(p),
      method = "L-BFGS-B",
      lower = lower_m1,
      upper = upper_m1,
      control = list(maxit = MAX_ITER)
    )
  }, error = function(e) NULL)

  if (!is.null(result)) {
    ll <- -result$value
    cat(sprintf("LL=%.1f, conv=%d\n", ll, result$convergence))
    if (ll > best_ll_m1) {
      best_ll_m1 <- ll
      best_result_m1 <- result
    }
  } else {
    cat("FAILED\n")
  }
}

params_m1 <- best_result_m1$par
ll_m1 <- -best_result_m1$value
alpha_m1 <- exp(params_m1[14])
omega_m1 <- exp(params_m1[15])
beta_decay_m1 <- -log(omega_m1)
halflife_m1 <- log(2) / beta_decay_m1
branching_m1 <- alpha_m1 * omega_m1 / (1 - omega_m1)  # E[offspring] = α × ω / (1-ω)

cat(sprintf("\nModel 1 Results (Cross-Day Hawkes):\n"))
cat(sprintf("  Log-likelihood: %.2f\n", ll_m1))
cat(sprintf("  Convergence code: %d\n", best_result_m1$convergence))
cat(sprintf("  γ (population) = %.4f\n", params_m1[2]))
cat(sprintf("  δ (poverty) = %.4f\n", params_m1[3]))
cat(sprintf("  ζ (log_cpi) = %.4f\n", params_m1[4]))
cat(sprintf("  α (triggering intensity) = %.4f\n", alpha_m1))
cat(sprintf("  ω (daily decay factor) = %.4f\n", omega_m1))
cat(sprintf("  β (decay rate) = %.4f per day\n", beta_decay_m1))
cat(sprintf("  Half-life = %.2f days\n", halflife_m1))
cat(sprintf("  Branching ratio = %.4f\n", branching_m1))

# Diagnostics
if (branching_m1 > 0.9) {
  cat("  ⚠️  WARNING: Branching ratio > 0.9 (near-critical process)\n")
}
if (sign(params_m1[3]) != sign(params_m0[3])) {
  cat("  ⚠️  WARNING: Poverty coefficient sign differs from Poisson model!\n")
}
cat("\n")

# =============================================================================
# 8. MODEL COMPARISON
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL COMPARISON                                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

n_params_m0 <- 13
n_params_m1 <- 15

# AIC and BIC
aic_m0 <- -2 * ll_m0 + 2 * n_params_m0
aic_m1 <- -2 * ll_m1 + 2 * n_params_m1

bic_m0 <- -2 * ll_m0 + n_params_m0 * log(n_total)
bic_m1 <- -2 * ll_m1 + n_params_m1 * log(n_total)

# LR test
lr_stat <- 2 * (ll_m1 - ll_m0)
df <- n_params_m1 - n_params_m0
p_value <- pchisq(lr_stat, df = df, lower.tail = FALSE)

cat("Model fit summary:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-30s %10s %8s %12s %12s\n", "Model", "Log-Lik", "Params", "AIC", "BIC"))
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-30s %10.1f %8d %12.1f %12.1f\n", "M0: Poisson", ll_m0, n_params_m0, aic_m0, bic_m0))
cat(sprintf("%-30s %10.1f %8d %12.1f %12.1f\n", "M1: Cross-Day Hawkes", ll_m1, n_params_m1, aic_m1, bic_m1))
cat("─────────────────────────────────────────────────────────────────\n\n")

cat("Likelihood Ratio Test:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("M1 vs M0 (Does cross-day triggering exist?)\n"))
cat(sprintf("  LR statistic: %.2f, df: %d, p-value: %.2e\n", lr_stat, df, p_value))
if (p_value < 0.001) {
  cat("  *** SIGNIFICANT: Cross-day triggering improves fit ***\n\n")
} else {
  cat("  Not significant at α = 0.001\n\n")
}

# =============================================================================
# 9. SAVE RESULTS
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SAVING RESULTS                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

results <- list(
  model_type = "Cross-Day Hawkes (Discrete Daily)",
  n_events = n_total,
  n_days = D_max,

  model0 = list(
    name = "Poisson (Background Only)",
    params = params_m0,
    loglik = ll_m0,
    aic = aic_m0,
    bic = bic_m0,
    gamma = params_m0[2],
    delta = params_m0[3],
    zeta = params_m0[4]
  ),

  model1 = list(
    name = "Cross-Day Hawkes",
    params = params_m1,
    loglik = ll_m1,
    aic = aic_m1,
    bic = bic_m1,
    gamma = params_m1[2],
    delta = params_m1[3],
    zeta = params_m1[4],
    alpha = alpha_m1,
    omega = omega_m1,
    beta_decay = beta_decay_m1,
    halflife = halflife_m1,
    branching_ratio = branching_m1
  ),

  comparison = list(
    lr_statistic = lr_stat,
    df = df,
    p_value = p_value
  )
)

saveRDS(results, "model_results_crossday.rds")
cat("Saved: model_results_crossday.rds\n\n")

# =============================================================================
# 10. SUMMARY
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SUMMARY: CROSS-DAY VS JITTERED HAWKES                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("CROSS-DAY HAWKES RESULTS:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Half-life: %.2f days\n", halflife_m1))
cat(sprintf("  Branching ratio: %.4f\n", branching_m1))
cat(sprintf("  γ (population): %.4f\n", params_m1[2]))
cat(sprintf("  δ (poverty): %.4f\n", params_m1[3]))
cat(sprintf("  ζ (log_cpi): %.4f\n\n", params_m1[4]))

cat("COMPARISON TO JITTERED RESULTS (from sensitivity analysis):\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("  Jittered (1.0): Half-life=1.14 days, BR=1.0002\n")
cat(sprintf("  Cross-Day:      Half-life=%.2f days, BR=%.4f\n\n", halflife_m1, branching_m1))

if (abs(branching_m1 - 1.0) > 0.1) {
  cat("*** FINDING: Branching ratio differs from jittered model! ***\n")
  cat("    This suggests the jitter artifact was indeed biasing results.\n\n")
} else {
  cat("*** FINDING: Branching ratio still near 1.0 ***\n")
  cat("    This suggests triggering dynamics are genuinely near-critical.\n\n")
}

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  ESTIMATION COMPLETE                                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
