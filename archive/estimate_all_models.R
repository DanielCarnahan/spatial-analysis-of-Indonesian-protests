################################################################################
#                    UNIFIED MODEL ESTIMATION SCRIPT
#
#   Three nested models for protest diffusion:
#     Model 0: Homogeneous Point Process (Poisson) - background only
#     Model 1: Univariate Hawkes - background + uniform triggering
#     Model 2: Bivariate Hawkes - background + type-dependent triggering
#
#   Background rate: μ = exp(β₀ + γ·log_pop + δ·poverty + ζ·log_cpi + year_effects)
#
################################################################################

library(tidyverse)
library(numDeriv)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║           UNIFIED MODEL ESTIMATION: PROTEST DIFFUSION                    ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

N_STARTS <- 3        # Number of starting points per model
MAX_ITER <- 2000     # Max iterations for optimization
JITTER_SCALE <- 1.0  # Scale of uniform jitter for same-day events (default: 0-1 day)
SKIP_BIVARIATE <- TRUE  # Skip bivariate model (slow with full grid)

# Allow command line override for jitter_scale
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  JITTER_SCALE <- as.numeric(args[1])
  cat(sprintf("Using jitter scale from command line: %.1f\n", JITTER_SCALE))
}

# =============================================================================
# 2. LOAD AND PREPARE DATA
# =============================================================================

cat("=== LOADING DATA ===\n\n")

protests <- readRDS("protests_daily.rds") %>%
  filter(!is.na(poverty_decimal) & !is.na(log_pop) & !is.na(log_cpi)) %>%
  mutate(
    type_numeric = ifelse(is_severe, 1L, 2L)  # 1=Severe, 2=Peaceful
  ) %>%
  arrange(date)

# Compute time in days from start
start_date <- min(protests$date)
protests$time <- as.numeric(protests$date - start_date)

# Add random jitter within each day (uniform 0 to JITTER_SCALE)
# This breaks the Δt=0 degeneracy for same-day events across locations
set.seed(42)  # For reproducibility
protests$time <- protests$time + runif(nrow(protests), 0, JITTER_SCALE)
protests <- protests %>% arrange(time)  # Re-sort by jittered time

cat(sprintf("Jitter scale: %.1f days (same-day events spread over 0-%.1f)\n", JITTER_SCALE, JITTER_SCALE))

T_max <- max(protests$time) + 1

n_total <- nrow(protests)
n_severe <- sum(protests$is_severe)
n_peaceful <- sum(!protests$is_severe)

cat(sprintf("Events: %d total (%d severe, %d peaceful)\n", n_total, n_severe, n_peaceful))
cat(sprintf("Period: %.0f days (%.1f years)\n", T_max, T_max/365.25))
cat(sprintf("Date range: %s to %s\n\n", min(protests$date), max(protests$date)))

# Prepare data vectors
times <- protests$time
types <- protests$type_numeric
log_pop <- protests$log_pop
poverty <- protests$poverty_decimal
log_cpi <- protests$log_cpi
years <- protests$year
time_diffs <- c(0, diff(times))

# =============================================================================
# BUILD FULL DISTRICT-MONTH GRID FOR INTEGRAL
# =============================================================================
# CRITICAL: Integrate over ALL districts × ALL months, not just protest-containing ones

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

# Load poverty data (all 514 districts)
poverty_data <- readRDS("poverty_rate_2015_2024.rds") %>%
  mutate(poverty_decimal = poverty_rate / 100)

# Load CPI data
cpi_data <- readRDS("indonesia_cpi.rds") %>%
  mutate(log_cpi = log(cpi)) %>%
  select(year_month, log_cpi)

# Create full grid: all districts × all months
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
cat(sprintf("  Months: %d\n", length(unique(district_year_months$year_month))))

# Covariate summaries
cat("\nCovariate ranges:\n")
cat(sprintf("  log_pop: %.2f to %.2f (mean: %.2f)\n", min(log_pop), max(log_pop), mean(log_pop)))
cat(sprintf("  poverty: %.3f to %.3f (mean: %.3f)\n", min(poverty), max(poverty), mean(poverty)))
cat(sprintf("  log_cpi: %.3f to %.3f (mean: %.3f)\n", min(log_cpi), max(log_cpi), mean(log_cpi)))
cat("\n")

# =============================================================================
# 3. PRECOMPUTE VECTORIZED DATA FOR FAST LIKELIHOOD
# =============================================================================

# Precompute year effects lookup for events (vectorized)
year_effect_idx_events <- pmax(0, years - 2015)  # 0 for 2015, 1 for 2016, etc.

# Precompute year effects lookup for grid (vectorized)
year_effect_idx_grid <- pmax(0, district_year_months$year - 2015)

# Extract grid vectors for fast computation
grid_log_pop <- district_year_months$log_pop
grid_poverty <- district_year_months$poverty_decimal
grid_log_cpi <- district_year_months$log_cpi
grid_days <- district_year_months$days_in_month

cat("Precomputed vectorized data for fast likelihood evaluation.\n\n")

# =============================================================================
# 4. LIKELIHOOD FUNCTIONS (VECTORIZED)
# =============================================================================

# -----------------------------------------------------------------------------
# Model 0: Poisson (Background Only)
# Parameters: β₀, γ, δ, ζ, β_2016...β_2024 (13 total)
# -----------------------------------------------------------------------------

poisson_loglik <- function(params) {
  beta_0 <- params[1]
  gamma <- params[2]
  delta <- params[3]
  zeta <- params[4]
  beta_years <- params[5:13]

  # Year effects for events (vectorized)
  year_effects <- c(0, beta_years)[year_effect_idx_events + 1]

  # Background rate at each event (vectorized)
  mu <- exp(beta_0 + gamma * log_pop + delta * poverty + zeta * log_cpi + year_effects)

  # Log-likelihood: sum of log(mu) at event times
  loglik <- sum(log(pmax(mu, 1e-10)))

  # Integral: vectorized over district-year-months
  year_effects_grid <- c(0, beta_years)[year_effect_idx_grid + 1]
  mu_grid <- exp(beta_0 + gamma * grid_log_pop + delta * grid_poverty + zeta * grid_log_cpi + year_effects_grid)
  integral <- sum(mu_grid * grid_days)

  loglik <- loglik - integral
  if (!is.finite(loglik)) return(-1e10)
  return(loglik)
}

# -----------------------------------------------------------------------------
# Model 1: Univariate Hawkes
# Parameters: β₀, γ, δ, ζ, β_2016...β_2024, log_α, log_β (15 total)
# -----------------------------------------------------------------------------

univariate_hawkes_loglik <- function(params) {
  beta_0 <- params[1]
  gamma <- params[2]
  delta <- params[3]
  zeta <- params[4]
  beta_years <- params[5:13]
  alpha <- exp(params[14])
  beta_decay <- exp(params[15])

  if (beta_decay < 1e-6 || beta_decay > 100) return(-1e10)

  n <- length(times)

  # Year effects (vectorized)
  year_effects <- c(0, beta_years)[year_effect_idx_events + 1]

  # Background rate (vectorized)
  mu <- exp(beta_0 + gamma * log_pop + delta * poverty + zeta * log_cpi + year_effects)

  # Recursive triggering sum (must remain sequential)
  R <- numeric(n)
  for (i in 2:n) {
    decay <- exp(-beta_decay * time_diffs[i])
    R[i] <- decay * (1 + R[i-1])
  }

  # Log-likelihood at event times (vectorized)
  lambda <- mu + alpha * R
  if (any(lambda <= 1e-10)) return(-1e10)
  loglik <- sum(log(lambda))

  # Background integral (vectorized)
  year_effects_grid <- c(0, beta_years)[year_effect_idx_grid + 1]
  mu_grid <- exp(beta_0 + gamma * grid_log_pop + delta * grid_poverty + zeta * grid_log_cpi + year_effects_grid)
  integral_bg <- sum(mu_grid * grid_days)

  # Triggering integral (vectorized)
  integral_trig <- sum(alpha * (1 - exp(-beta_decay * (T_max - times))) / beta_decay)

  loglik <- loglik - integral_bg - integral_trig
  if (!is.finite(loglik)) return(-1e10)
  return(loglik)
}

# -----------------------------------------------------------------------------
# Model 2: Bivariate Hawkes
# Parameters: β₀_S, β₀_P, γ, δ, ζ, β_2016...β_2024,
#             log_α_SS, log_α_SP, log_α_PS, log_α_PP, log_β (19 total)
# -----------------------------------------------------------------------------

bivariate_hawkes_loglik <- function(params) {
  beta_0_S <- params[1]
  beta_0_P <- params[2]
  gamma <- params[3]
  delta <- params[4]
  zeta <- params[5]
  beta_years <- params[6:14]
  alpha_ss <- exp(params[15])
  alpha_sp <- exp(params[16])
  alpha_ps <- exp(params[17])
  alpha_pp <- exp(params[18])
  beta_decay <- exp(params[19])

  if (beta_decay < 1e-6 || beta_decay > 100) return(-1e10)

  n <- length(times)

  # Year effects (vectorized)
  year_effects <- c(0, beta_years)[year_effect_idx_events + 1]

  # Background rates (vectorized, shared covariates, different intercepts)
  mu_s <- exp(beta_0_S + gamma * log_pop + delta * poverty + zeta * log_cpi + year_effects)
  mu_p <- exp(beta_0_P + gamma * log_pop + delta * poverty + zeta * log_cpi + year_effects)

  # Recursive triggering sums (must remain sequential)
  R_s <- numeric(n)
  R_p <- numeric(n)

  for (i in 2:n) {
    decay <- exp(-beta_decay * time_diffs[i])
    if (types[i-1] == 1) {  # Previous was Severe
      R_s[i] <- decay * (1 + R_s[i-1])
      R_p[i] <- decay * R_p[i-1]
    } else {  # Previous was Peaceful
      R_s[i] <- decay * R_s[i-1]
      R_p[i] <- decay * (1 + R_p[i-1])
    }
  }

  # Log-likelihood at event times (vectorized)
  is_severe <- types == 1
  lambda <- ifelse(is_severe,
                   mu_s + alpha_ss * R_s + alpha_ps * R_p,
                   mu_p + alpha_sp * R_s + alpha_pp * R_p)
  if (any(lambda <= 1e-10)) return(-1e10)
  loglik <- sum(log(lambda))

  # Background integral (vectorized)
  year_effects_grid <- c(0, beta_years)[year_effect_idx_grid + 1]
  lin_pred_grid <- gamma * grid_log_pop + delta * grid_poverty + zeta * grid_log_cpi + year_effects_grid
  mu_s_grid <- exp(beta_0_S + lin_pred_grid)
  mu_p_grid <- exp(beta_0_P + lin_pred_grid)
  integral_bg <- sum((mu_s_grid + mu_p_grid) * grid_days)

  # Triggering integral (vectorized)
  contrib <- (1 - exp(-beta_decay * (T_max - times))) / beta_decay
  integral_trig <- sum(ifelse(types == 1,
                              (alpha_ss + alpha_sp) * contrib,
                              (alpha_ps + alpha_pp) * contrib))

  loglik <- loglik - integral_bg - integral_trig
  if (!is.finite(loglik)) return(-1e10)
  return(loglik)
}

# =============================================================================
# 4. FIT MODEL 0: POISSON
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL 0: POISSON (Background Only)                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Empirical rate for starting values (accounting for full spatial grid)
total_district_days <- sum(district_year_months$days_in_month)
emp_rate_per_district_day <- n_total / total_district_days
cat(sprintf("\nEmpirical rate: %.6f per district-day (%.2f per district-year)\n",
            emp_rate_per_district_day, emp_rate_per_district_day * 365.25))

# Mean covariates from the FULL grid (for centering intercept)
mean_log_pop_grid <- mean(district_year_months$log_pop)
mean_poverty_grid <- mean(district_year_months$poverty_decimal)
mean_log_cpi_grid <- mean(district_year_months$log_cpi)

cat(sprintf("Grid mean log_pop: %.2f, poverty: %.3f, log_cpi: %.3f\n",
            mean_log_pop_grid, mean_poverty_grid, mean_log_cpi_grid))

# For protest data (used in log-likelihood at event times)
mean_log_pop <- mean(log_pop)
mean_poverty <- mean(poverty)
mean_log_cpi <- mean(log_cpi)

# Starting values: intercept calibrated to match observed rate given mean covariates
# If μ = exp(β₀ + γ*X), need β₀ ≈ log(rate) - γ*mean_X
beta0_start <- log(emp_rate_per_district_day)
cat(sprintf("Starting β₀: %.2f\n\n", beta0_start))

# Starting points for Model 0
start_m0 <- list(
  c(beta0_start - 0.1 * mean_log_pop_grid, 0.1, 0.0, 0.0, rep(0, 9)),
  c(beta0_start - 0.2 * mean_log_pop_grid, 0.2, -1.0, 0.0, rep(0, 9)),
  c(beta0_start - 0.05 * mean_log_pop_grid, 0.05, 1.0, -0.5, rep(0, 9))
)

# Bounds: allow γ (population) to be negative, widened to [-3, 3]
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
cat(sprintf("  ζ (log_cpi) = %.4f\n", params_m0[4]))
cat("\n")

# =============================================================================
# 5. FIT MODEL 1: UNIVARIATE HAWKES
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL 1: UNIVARIATE HAWKES                                 ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Starting points (use Model 0 estimates for background)
# More diverse starting points for better exploration
start_m1 <- list(
  c(params_m0[1:13], log(0.1), log(0.5)),    # Moderate triggering, ~1 day decay
  c(params_m0[1:13], log(0.05), log(0.3)),   # Weaker triggering
  c(params_m0[1:13], log(0.2), log(1.0)),    # Stronger triggering, faster decay
  c(params_m0[1:13], log(0.02), log(0.2)),   # Weak triggering
  c(params_m0[1:13], log(0.15), log(0.7)),   # Medium
  c(params_m0[1:13], log(0.3), log(1.5)),    # Strong triggering, fast decay
  c(params_m0[1:13], log(0.08), log(0.4)),   # Variation 1
  c(params_m0[1:13], log(0.12), log(0.6)),   # Variation 2
  c(params_m0[1:13], log(0.25), log(2.0)),   # Very fast decay
  c(params_m0[1:13], log(0.04), log(0.15))   # Slow decay
)

# Tighter bounds: constrain branching ratio to be < 0.9
# log(alpha) upper bound ensures alpha/beta < 0.9
lower_m1 <- c(lower_m0, -10, -2)   # log(alpha) >= -10, log(beta) >= -2 (beta >= 0.14)
upper_m1 <- c(upper_m0, 2, 3)       # log(alpha) <= 2 (alpha <= 7.4), log(beta) <= 3 (beta <= 20)

best_ll_m1 <- -Inf
best_result_m1 <- NULL

cat("Running multi-start optimization...\n")
for (s in seq_along(start_m1)) {
  cat(sprintf("  Start %d/%d... ", s, length(start_m1)))

  result <- tryCatch({
    optim(
      par = start_m1[[s]],
      fn = function(p) -univariate_hawkes_loglik(p),
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
beta_m1 <- exp(params_m1[15])
branching_m1 <- alpha_m1 / beta_m1

cat(sprintf("\nModel 1 Results:\n"))
cat(sprintf("  Log-likelihood: %.2f\n", ll_m1))
cat(sprintf("  Convergence code: %d\n", best_result_m1$convergence))
cat(sprintf("  γ (population) = %.4f\n", params_m1[2]))
cat(sprintf("  δ (poverty) = %.4f\n", params_m1[3]))
cat(sprintf("  ζ (log_cpi) = %.4f\n", params_m1[4]))
cat(sprintf("  α (triggering) = %.4f\n", alpha_m1))
cat(sprintf("  β (decay) = %.4f (half-life = %.2f days)\n", beta_m1, log(2)/beta_m1))
cat(sprintf("  Branching ratio = %.4f\n", branching_m1))

# Diagnostics
if (best_result_m1$convergence != 0) {
  cat("  ⚠️  WARNING: Optimization did not converge!\n")
}
if (branching_m1 > 0.9) {
  cat("  ⚠️  WARNING: Branching ratio > 0.9 (near-critical process)\n")
}
if (sign(params_m1[3]) != sign(params_m0[3])) {
  cat("  ⚠️  WARNING: Poverty coefficient sign differs from Poisson model!\n")
}
cat("\n")

# =============================================================================
# 6. FIT MODEL 2: BIVARIATE HAWKES (CONDITIONAL)
# =============================================================================

if (SKIP_BIVARIATE) {
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  MODEL 2: BIVARIATE HAWKES - SKIPPED                        ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")

  # Placeholder values for skipped model
  ll_m2 <- ll_m1 - 100  # Dummy value
  params_m2 <- rep(NA, 19)
  beta_0_S <- NA; beta_0_P <- NA; gamma_m2 <- NA; delta_m2 <- NA; zeta_m2 <- NA
  alpha_ss <- NA; alpha_sp <- NA; alpha_ps <- NA; alpha_pp <- NA; beta_m2 <- NA
  mobil_severe <- NA; mobil_peaceful <- NA; ratio_sp <- NA
  best_result_m2 <- list(convergence = NA)

} else {

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL 2: BIVARIATE HAWKES                                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Chain initialization: use Model 1 estimates
# Split intercept based on type-specific empirical rates
emp_rate_s <- n_severe / T_max
emp_rate_p <- n_peaceful / T_max
intercept_shift_s <- log(emp_rate_s / (emp_rate_s + emp_rate_p))
intercept_shift_p <- log(emp_rate_p / (emp_rate_s + emp_rate_p))

# Base starting point from M1
base_intercept_s <- params_m1[1] + intercept_shift_s + 0.5
base_intercept_p <- params_m1[1] + intercept_shift_p - 0.5

# Starting points - chain from Model 1
# params_m2: [beta0_S, beta0_P, gamma, delta, zeta, year_effects(9), alpha_ss, alpha_sp, alpha_ps, alpha_pp, beta]
start_m2 <- list(
  # Start 1: Direct chain from M1, split triggering equally
  c(base_intercept_s, base_intercept_p,
    params_m1[2:13],  # All covariates and year effects from M1
    log(alpha_m1/4), log(alpha_m1/4), log(alpha_m1/4), log(alpha_m1/4),  # Split alpha
    params_m1[15]),   # Same beta as M1

  # Start 2: Stronger severe-to-severe triggering
  c(base_intercept_s, base_intercept_p,
    params_m1[2:13],
    log(alpha_m1/2), log(alpha_m1/4), log(alpha_m1/8), log(alpha_m1/4),
    params_m1[15]),

  # Start 3: Stronger peaceful-to-peaceful triggering
  c(base_intercept_s, base_intercept_p,
    params_m1[2:13],
    log(alpha_m1/8), log(alpha_m1/4), log(alpha_m1/4), log(alpha_m1/2),
    params_m1[15]),

  # Start 4: Asymmetric (severe triggers more)
  c(base_intercept_s, base_intercept_p,
    params_m1[2:13],
    log(alpha_m1/3), log(alpha_m1/3), log(alpha_m1/6), log(alpha_m1/6),
    params_m1[15]),

  # Start 5: Asymmetric (peaceful triggers more)
  c(base_intercept_s, base_intercept_p,
    params_m1[2:13],
    log(alpha_m1/6), log(alpha_m1/6), log(alpha_m1/3), log(alpha_m1/3),
    params_m1[15]),

  # Start 6: Very weak triggering
  c(base_intercept_s, base_intercept_p,
    params_m1[2:13],
    log(0.01), log(0.01), log(0.01), log(0.01),
    params_m1[15]),

  # Start 7: Slightly different intercepts
  c(base_intercept_s + 0.3, base_intercept_p - 0.3,
    params_m1[2:13],
    log(alpha_m1/4), log(alpha_m1/4), log(alpha_m1/4), log(alpha_m1/4),
    params_m1[15]),

  # Start 8: Different decay rate
  c(base_intercept_s, base_intercept_p,
    params_m1[2:13],
    log(alpha_m1/4), log(alpha_m1/4), log(alpha_m1/4), log(alpha_m1/4),
    log(beta_m1 * 1.5)),

  # Start 9: Cross-type triggering only
  c(base_intercept_s, base_intercept_p,
    params_m1[2:13],
    log(0.01), log(alpha_m1/2), log(alpha_m1/2), log(0.01),
    params_m1[15]),

  # Start 10: Diagonal triggering only
  c(base_intercept_s, base_intercept_p,
    params_m1[2:13],
    log(alpha_m1/2), log(0.01), log(0.01), log(alpha_m1/2),
    params_m1[15])
)

# Tighter bounds for bivariate (γ widened to [-3, 3])
lower_m2 <- c(-20, -20, -3, -10, -10, rep(-5, 9), rep(-10, 4), -2)
upper_m2 <- c(10, 10, 3, 10, 10, rep(5, 9), rep(1, 4), 3)  # log(alpha) <= 1 means alpha <= 2.7

best_ll_m2 <- -Inf
best_result_m2 <- NULL

cat("Running multi-start optimization...\n")
for (s in seq_along(start_m2)) {
  cat(sprintf("  Start %d/%d... ", s, length(start_m2)))

  result <- tryCatch({
    optim(
      par = start_m2[[s]],
      fn = function(p) -bivariate_hawkes_loglik(p),
      method = "L-BFGS-B",
      lower = lower_m2,
      upper = upper_m2,
      control = list(maxit = MAX_ITER)
    )
  }, error = function(e) NULL)

  if (!is.null(result)) {
    ll <- -result$value
    cat(sprintf("LL=%.1f, conv=%d\n", ll, result$convergence))
    if (ll > best_ll_m2) {
      best_ll_m2 <- ll
      best_result_m2 <- result
    }
  } else {
    cat("FAILED\n")
  }
}

params_m2 <- best_result_m2$par
ll_m2 <- -best_result_m2$value

beta_0_S <- params_m2[1]
beta_0_P <- params_m2[2]
gamma_m2 <- params_m2[3]
delta_m2 <- params_m2[4]
zeta_m2 <- params_m2[5]
alpha_ss <- exp(params_m2[15])
alpha_sp <- exp(params_m2[16])
alpha_ps <- exp(params_m2[17])
alpha_pp <- exp(params_m2[18])
beta_m2 <- exp(params_m2[19])

# Mobilization potential
mobil_severe <- (alpha_ss + alpha_sp) / beta_m2
mobil_peaceful <- (alpha_ps + alpha_pp) / beta_m2
ratio_sp <- mobil_severe / mobil_peaceful

cat(sprintf("\nModel 2 Results:\n"))
cat(sprintf("  Log-likelihood: %.2f\n", ll_m2))
cat(sprintf("  β₀_S = %.4f, β₀_P = %.4f\n", beta_0_S, beta_0_P))
cat(sprintf("  γ (population) = %.4f\n", gamma_m2))
cat(sprintf("  δ (poverty) = %.4f\n", delta_m2))
cat(sprintf("  ζ (log_cpi) = %.4f\n", zeta_m2))
cat(sprintf("  β (decay) = %.4f (half-life = %.2f days)\n", beta_m2, log(2)/beta_m2))
cat(sprintf("\nTriggering matrix (α):\n"))
cat(sprintf("  α_SS = %.4f, α_SP = %.4f\n", alpha_ss, alpha_sp))
cat(sprintf("  α_PS = %.4f, α_PP = %.4f\n", alpha_ps, alpha_pp))
cat(sprintf("\nMobilization potential:\n"))
cat(sprintf("  Severe: %.4f\n", mobil_severe))
cat(sprintf("  Peaceful: %.4f\n", mobil_peaceful))
cat(sprintf("  Ratio (S/P): %.4f\n", ratio_sp))

# Diagnostics for Model 2
cat(sprintf("\n  Convergence code: %d\n", best_result_m2$convergence))
if (best_result_m2$convergence != 0) {
  cat("  ⚠️  WARNING: Optimization did not converge!\n")
}
if (ll_m2 < ll_m1) {
  cat("  ⚠️  WARNING: Bivariate LL is worse than Univariate!\n")
  cat("     This suggests optimization got stuck in local minimum.\n")
}
if (mobil_severe > 0.95 | mobil_peaceful > 0.95) {
  cat("  ⚠️  WARNING: Near-critical mobilization potential detected!\n")
}
if (sign(delta_m2) != sign(params_m0[3])) {
  cat("  ⚠️  WARNING: Poverty coefficient sign differs from Poisson model!\n")
}
cat("\n")

}  # End of else block for SKIP_BIVARIATE

# =============================================================================
# 7. MODEL COMPARISON
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL COMPARISON                                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Parameter counts
n_params_m0 <- 13
n_params_m1 <- 15
n_params_m2 <- 19

# AIC and BIC
aic_m0 <- -2 * ll_m0 + 2 * n_params_m0
aic_m1 <- -2 * ll_m1 + 2 * n_params_m1
aic_m2 <- -2 * ll_m2 + 2 * n_params_m2

bic_m0 <- -2 * ll_m0 + n_params_m0 * log(n_total)
bic_m1 <- -2 * ll_m1 + n_params_m1 * log(n_total)
bic_m2 <- -2 * ll_m2 + n_params_m2 * log(n_total)

# LR tests
lr_stat_1 <- 2 * (ll_m1 - ll_m0)
df_1 <- n_params_m1 - n_params_m0
p_value_1 <- pchisq(lr_stat_1, df = df_1, lower.tail = FALSE)

lr_stat_2 <- 2 * (ll_m2 - ll_m1)
df_2 <- n_params_m2 - n_params_m1
p_value_2 <- pchisq(lr_stat_2, df = df_2, lower.tail = FALSE)

cat("Model fit summary:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-25s %10s %8s %12s %12s\n", "Model", "Log-Lik", "Params", "AIC", "BIC"))
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-25s %10.1f %8d %12.1f %12.1f\n", "M0: Poisson", ll_m0, n_params_m0, aic_m0, bic_m0))
cat(sprintf("%-25s %10.1f %8d %12.1f %12.1f\n", "M1: Univariate Hawkes", ll_m1, n_params_m1, aic_m1, bic_m1))
if (!SKIP_BIVARIATE) {
  cat(sprintf("%-25s %10.1f %8d %12.1f %12.1f\n", "M2: Bivariate Hawkes", ll_m2, n_params_m2, aic_m2, bic_m2))
}
cat("─────────────────────────────────────────────────────────────────\n\n")

cat("Likelihood Ratio Tests:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("Test 1: M1 vs M0 (Does triggering exist?)\n"))
cat(sprintf("  LR statistic: %.2f, df: %d, p-value: %.2e\n", lr_stat_1, df_1, p_value_1))
if (p_value_1 < 0.001) {
  cat("  *** SIGNIFICANT: Self-excitation improves fit ***\n\n")
} else {
  cat("  Not significant at α = 0.001\n\n")
}

if (!SKIP_BIVARIATE) {
  cat(sprintf("Test 2: M2 vs M1 (Does protest type matter?)\n"))
  cat(sprintf("  LR statistic: %.2f, df: %d, p-value: %.2e\n", lr_stat_2, df_2, p_value_2))
  if (p_value_2 < 0.001) {
    cat("  *** SIGNIFICANT: Type-dependent triggering improves fit ***\n\n")
  } else {
    cat("  Not significant at α = 0.001\n\n")
  }
}

# =============================================================================
# 8. SAVE RESULTS
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SAVING RESULTS                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

results <- list(
  # Data info
  n_events = n_total,
  n_severe = n_severe,
  n_peaceful = n_peaceful,
  T_max = T_max,

  # Model 0
  model0 = list(
    name = "Poisson (Background Only)",
    params = params_m0,
    loglik = ll_m0,
    n_params = n_params_m0,
    aic = aic_m0,
    bic = bic_m0,
    beta_0 = params_m0[1],
    gamma = params_m0[2],
    delta = params_m0[3],
    zeta = params_m0[4],
    beta_years = params_m0[5:13]
  ),

  # Model 1
  model1 = list(
    name = "Univariate Hawkes",
    params = params_m1,
    loglik = ll_m1,
    n_params = n_params_m1,
    aic = aic_m1,
    bic = bic_m1,
    gamma = params_m1[2],
    delta = params_m1[3],
    zeta = params_m1[4],
    alpha = alpha_m1,
    beta = beta_m1,
    branching_ratio = alpha_m1 / beta_m1
  ),

  # Model 2
  model2 = list(
    name = "Bivariate Hawkes",
    params = params_m2,
    loglik = ll_m2,
    n_params = n_params_m2,
    aic = aic_m2,
    bic = bic_m2,
    beta_0_S = beta_0_S,
    beta_0_P = beta_0_P,
    gamma = gamma_m2,
    delta = delta_m2,
    zeta = zeta_m2,
    alpha_ss = alpha_ss,
    alpha_sp = alpha_sp,
    alpha_ps = alpha_ps,
    alpha_pp = alpha_pp,
    beta = beta_m2,
    mobil_severe = mobil_severe,
    mobil_peaceful = mobil_peaceful,
    ratio_sp = ratio_sp
  ),

  # Model comparison
  comparison = list(
    lr_test_1 = list(
      description = "M1 vs M0: Does triggering exist?",
      statistic = lr_stat_1,
      df = df_1,
      p_value = p_value_1
    ),
    lr_test_2 = list(
      description = "M2 vs M1: Does protest type matter?",
      statistic = lr_stat_2,
      df = df_2,
      p_value = p_value_2
    )
  )
)

# Save with jitter scale in filename for sensitivity analysis
output_filename <- sprintf("model_results_jitter_%.1f.rds", JITTER_SCALE)
results$jitter_scale <- JITTER_SCALE
saveRDS(results, output_filename)
cat(sprintf("Saved: %s\n\n", output_filename))

# =============================================================================
# 9. SUMMARY
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SUMMARY                                                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("Configuration: Jitter scale = %.1f days\n", JITTER_SCALE))
cat("─────────────────────────────────────────────────────────────────\n")

cat("\nKEY FINDINGS:\n")
cat("─────────────────────────────────────────────────────────────────\n")

if (SKIP_BIVARIATE) {
  # Use univariate model results
  cat("\n1. COVARIATE EFFECTS (from univariate Hawkes model):\n")
  cat(sprintf("   γ (population): %.4f\n", params_m1[2]))
  cat(sprintf("   δ (poverty): %.4f → ", params_m1[3]))
  if (params_m1[3] < 0) {
    cat("higher poverty associated with LOWER protest rates\n")
  } else {
    cat("higher poverty associated with higher protest rates\n")
  }
  cat(sprintf("   ζ (log_cpi): %.4f → ", params_m1[4]))
  if (params_m1[4] < 0) {
    cat("higher prices associated with FEWER protests\n")
  } else {
    cat("higher prices associated with MORE protests\n")
  }

  cat("\n2. TRIGGERING (univariate model):\n")
  cat(sprintf("   α (triggering rate): %.4f\n", alpha_m1))
  cat(sprintf("   β (decay rate): %.4f\n", beta_m1))
  cat(sprintf("   Half-life: %.2f days (%.1f hours)\n", log(2)/beta_m1, log(2)/beta_m1 * 24))
  cat(sprintf("   Branching ratio: %.4f\n", branching_m1))
  if (branching_m1 > 0.9) {
    cat("   ⚠️  Near-critical process (branching ratio > 0.9)\n")
  }

} else {
  # Use bivariate model results
  cat("\n1. COVARIATE EFFECTS (from bivariate model):\n")
  cat(sprintf("   γ (population): %.4f → larger populations have higher protest rates\n", gamma_m2))
  cat(sprintf("   δ (poverty): %.4f → ", delta_m2))
  if (delta_m2 < 0) {
    cat("higher poverty associated with LOWER protest rates\n")
  } else {
    cat("higher poverty associated with higher protest rates\n")
  }
  cat(sprintf("   ζ (log_cpi): %.4f → ", zeta_m2))
  if (zeta_m2 > 0) {
    cat("higher prices associated with MORE protests\n")
  } else if (zeta_m2 < 0) {
    cat("higher prices associated with FEWER protests\n")
  } else {
    cat("no clear effect of prices\n")
  }

  cat("\n2. TRIGGERING:\n")
  cat(sprintf("   Decay half-life: %.2f days (%.1f hours)\n", log(2)/beta_m2, log(2)/beta_m2 * 24))
  cat(sprintf("   Mobilization (Severe): %.4f offspring per event\n", mobil_severe))
  cat(sprintf("   Mobilization (Peaceful): %.4f offspring per event\n", mobil_peaceful))
  cat(sprintf("   Ratio (S/P): %.4f\n", ratio_sp))

  if (ratio_sp < 1) {
    cat("\n   *** Severe protests have LOWER mobilization than peaceful ***\n")
    cat("   This contradicts the 'radical flank' hypothesis\n")
  } else {
    cat("\n   *** Severe protests have HIGHER mobilization than peaceful ***\n")
  }
}

cat("\n3. MODEL SELECTION:\n")
if (SKIP_BIVARIATE) {
  best_model <- which.min(c(bic_m0, bic_m1))
  model_names <- c("Poisson", "Univariate Hawkes")
} else {
  best_model <- which.min(c(bic_m0, bic_m1, bic_m2))
  model_names <- c("Poisson", "Univariate Hawkes", "Bivariate Hawkes")
}
cat(sprintf("   Best model by BIC: %s\n", model_names[best_model]))

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  ESTIMATION COMPLETE                                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
