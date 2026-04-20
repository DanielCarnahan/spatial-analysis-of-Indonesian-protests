################################################################################
#                 SPATIO-TEMPORAL HAWKES MODEL ESTIMATION
#
#   Combines cross-day triggering with spatial distance decay
#
#   Key features:
#     1. Cross-day triggering only (no same-day)
#     2. Gaussian spatial kernel: exp(-d²/2σ²)
#     3. Only nearby events (< cutoff) can trigger
#     4. Background covariates: population, poverty, CPI, year effects
#
#   This addresses the identification problem in the univariate cross-day model
#   by requiring triggering events to be spatially proximate.
#
################################################################################

library(tidyverse)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   SPATIO-TEMPORAL HAWKES: CROSS-DAY + DISTANCE DECAY                    ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

cat("MODEL: λ(s,t) = μ(s) + α × Σ_{j: t_j < t} exp(-β(t-t_j)) × exp(-d(s,s_j)²/2σ²)\n")
cat("  - Cross-day: events on day d triggered only by days < d\n")
cat("  - Spatial: Gaussian kernel with bandwidth σ\n")
cat("  - Only events within spatial cutoff contribute\n\n")

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

SPATIAL_CUTOFF_KM <- 500   # Only events within 500km can trigger
TEMPORAL_CUTOFF_DAYS <- 60 # Only events within 60 days can trigger
N_STARTS <- 5              # Multi-start optimization
MAX_ITER <- 1000           # Max iterations

cat(sprintf("Configuration:\n"))
cat(sprintf("  Spatial cutoff: %d km\n", SPATIAL_CUTOFF_KM))
cat(sprintf("  Temporal cutoff: %d days\n", TEMPORAL_CUTOFF_DAYS))
cat(sprintf("  Optimization starts: %d\n\n", N_STARTS))

# =============================================================================
# 2. HAVERSINE DISTANCE FUNCTION
# =============================================================================

haversine_km <- function(lon1, lat1, lon2, lat2) {
  # Convert to radians
  lon1_rad <- lon1 * pi / 180
  lat1_rad <- lat1 * pi / 180
  lon2_rad <- lon2 * pi / 180
  lat2_rad <- lat2 * pi / 180

  # Haversine formula
  dlon <- lon2_rad - lon1_rad
  dlat <- lat2_rad - lat1_rad

  a <- sin(dlat/2)^2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))

  # Earth radius in km
  R <- 6371
  return(R * c)
}

# =============================================================================
# 3. LOAD AND PREPARE DATA
# =============================================================================

cat("=== LOADING DATA ===\n\n")

protests <- readRDS("protests_daily.rds") %>%
  filter(!is.na(poverty_decimal) & !is.na(log_pop) & !is.na(log_cpi) &
         !is.na(longitude) & !is.na(latitude)) %>%
  arrange(date)

# Compute day index (integer days from start) - NO JITTER
start_date <- min(protests$date)
protests$day <- as.numeric(protests$date - start_date)

n_total <- nrow(protests)
D_max <- max(protests$day) + 1

cat(sprintf("Events with coordinates: %d\n", n_total))
cat(sprintf("Days: %d total\n", D_max))
cat(sprintf("Date range: %s to %s\n\n", min(protests$date), max(protests$date)))

# Extract vectors
days <- protests$day
lons <- protests$longitude
lats <- protests$latitude
log_pop <- protests$log_pop
poverty <- protests$poverty_decimal
log_cpi <- protests$log_cpi
years <- protests$year

# =============================================================================
# 4. BUILD FULL INTEGRATION GRID (same as cross-day model)
# =============================================================================

cat("=== BUILDING FULL INTEGRATION GRID ===\n\n")

# Load population data
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

# Load poverty and CPI
poverty_data <- readRDS("poverty_rate_2015_2024.rds") %>%
  mutate(poverty_decimal = poverty_rate / 100)

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

district_year_months <- full_grid %>%
  left_join(pop_extrapolated %>% select(district, year, log_pop), by = c("district", "year")) %>%
  left_join(poverty_data %>% select(district, year, poverty_decimal), by = c("district", "year")) %>%
  left_join(cpi_data, by = "year_month") %>%
  filter(!is.na(log_pop) & !is.na(poverty_decimal) & !is.na(log_cpi))

cat(sprintf("Full integration grid: %d district-months\n\n", nrow(district_year_months)))

# Precompute year effect indices
year_effect_idx_events <- pmax(0, years - 2015)
year_effect_idx_grid <- pmax(0, district_year_months$year - 2015)
grid_log_pop <- district_year_months$log_pop
grid_poverty <- district_year_months$poverty_decimal
grid_log_cpi <- district_year_months$log_cpi
grid_days <- district_year_months$days_in_month

# =============================================================================
# 5. PRECOMPUTE SPATIOTEMPORAL NEIGHBORS
# =============================================================================

cat("=== PRECOMPUTING SPATIOTEMPORAL NEIGHBORS ===\n\n")
cat("This may take a few minutes...\n")

n <- length(days)

# For each event, find events that:
#   1. Occurred on a PREVIOUS day (cross-day only)
#   2. Are within SPATIAL_CUTOFF_KM
#   3. Are within TEMPORAL_CUTOFF_DAYS

neighbors <- vector("list", n)
neighbor_distances <- vector("list", n)
neighbor_day_diffs <- vector("list", n)

start_time <- Sys.time()

for (i in 2:n) {
  # Candidate events: all events before event i
  candidates <- 1:(i-1)

  # Filter: different day (cross-day only)
  day_diff <- days[i] - days[candidates]
  valid_temporal <- day_diff > 0 & day_diff <= TEMPORAL_CUTOFF_DAYS
  candidates <- candidates[valid_temporal]
  day_diff <- day_diff[valid_temporal]

  if (length(candidates) > 0) {
    # Compute distances
    dists <- haversine_km(lons[candidates], lats[candidates], lons[i], lats[i])

    # Filter: within spatial cutoff
    valid_spatial <- dists <= SPATIAL_CUTOFF_KM

    neighbors[[i]] <- candidates[valid_spatial]
    neighbor_distances[[i]] <- dists[valid_spatial]
    neighbor_day_diffs[[i]] <- day_diff[valid_spatial]
  } else {
    neighbors[[i]] <- integer(0)
    neighbor_distances[[i]] <- numeric(0)
    neighbor_day_diffs[[i]] <- numeric(0)
  }

  if (i %% 1000 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat(sprintf("  Progress: %d/%d (%.1f%%) - %.1f min elapsed\r",
                i, n, 100*i/n, elapsed))
  }
}

cat("\n\n")

# Summary statistics
n_with_neighbors <- sum(sapply(neighbors, length) > 0)
avg_neighbors <- mean(sapply(neighbors, length))
max_neighbors <- max(sapply(neighbors, length))

cat(sprintf("Neighbor statistics:\n"))
cat(sprintf("  Events with ≥1 neighbor: %d (%.1f%%)\n",
            n_with_neighbors, 100 * n_with_neighbors / n))
cat(sprintf("  Average neighbors per event: %.1f\n", avg_neighbors))
cat(sprintf("  Maximum neighbors: %d\n\n", max_neighbors))

# =============================================================================
# 6. LIKELIHOOD FUNCTIONS
# =============================================================================

# Poisson baseline (same as before)
poisson_loglik <- function(params) {
  beta_0 <- params[1]
  gamma <- params[2]
  delta <- params[3]
  zeta <- params[4]
  beta_years <- params[5:13]

  year_effects <- c(0, beta_years)[year_effect_idx_events + 1]
  mu <- exp(beta_0 + gamma * log_pop + delta * poverty + zeta * log_cpi + year_effects)
  loglik <- sum(log(pmax(mu, 1e-10)))

  year_effects_grid <- c(0, beta_years)[year_effect_idx_grid + 1]
  mu_grid <- exp(beta_0 + gamma * grid_log_pop + delta * grid_poverty + zeta * grid_log_cpi + year_effects_grid)
  integral <- sum(mu_grid * grid_days)

  loglik <- loglik - integral
  if (!is.finite(loglik)) return(-1e10)
  return(loglik)
}

# Spatio-temporal Hawkes with cross-day triggering
# Parameters: β₀, γ, δ, ζ, β_2016...β_2024, log_α, log_β, log_σ (16 total)
spatiotemporal_hawkes_loglik <- function(params) {
  beta_0 <- params[1]
  gamma <- params[2]
  delta <- params[3]
  zeta <- params[4]
  beta_years <- params[5:13]
  alpha <- exp(params[14])      # Triggering intensity
  beta_decay <- exp(params[15]) # Temporal decay rate
  sigma <- exp(params[16])      # Spatial bandwidth (km)

  # Background rate at each event
  year_effects <- c(0, beta_years)[year_effect_idx_events + 1]
  mu <- exp(beta_0 + gamma * log_pop + delta * poverty + zeta * log_cpi + year_effects)

  # Compute triggered intensity at each event
  triggered <- numeric(n)

  for (i in 2:n) {
    if (length(neighbors[[i]]) > 0) {
      # Temporal kernel: exp(-β × day_diff)
      h_temporal <- exp(-beta_decay * neighbor_day_diffs[[i]])

      # Spatial kernel: exp(-d²/(2σ²))
      g_spatial <- exp(-neighbor_distances[[i]]^2 / (2 * sigma^2))

      # Combined triggering
      triggered[i] <- sum(h_temporal * g_spatial)
    }
  }

  # Log-likelihood at event times
  lambda <- mu + alpha * triggered
  if (any(lambda <= 1e-10)) return(-1e10)
  loglik <- sum(log(lambda))

  # Background integral (same as Poisson)
  year_effects_grid <- c(0, beta_years)[year_effect_idx_grid + 1]
  mu_grid <- exp(beta_0 + gamma * grid_log_pop + delta * grid_poverty + zeta * grid_log_cpi + year_effects_grid)
  integral_bg <- sum(mu_grid * grid_days)

  # Triggering integral
  # Each event j contributes: α × ∫∫ exp(-β(t-t_j)) × exp(-d²/(2σ²)) dt ds
  # Temporal: ∫_{t_j+1}^{T_max} exp(-β(t-t_j)) dt ≈ (1/β) × (exp(-β) - exp(-β×(T_max-t_j)))
  # Spatial: ∫ exp(-d²/(2σ²)) ds ≈ 2πσ² (Gaussian integral over plane)
  # But we have a cutoff, so: 2πσ² × (1 - exp(-R²/(2σ²))) where R = SPATIAL_CUTOFF_KM

  spatial_integral <- 2 * pi * sigma^2 * (1 - exp(-SPATIAL_CUTOFF_KM^2 / (2 * sigma^2)))

  # Temporal integral for each event (cross-day: starts from day+1)
  days_remaining <- pmax(D_max - days - 1, 0)  # Days after the event's day
  temporal_integral <- (1 / beta_decay) * (exp(-beta_decay) - exp(-beta_decay * (days_remaining + 1)))
  temporal_integral <- pmax(temporal_integral, 0)

  # But we also have temporal cutoff
  effective_days <- pmin(days_remaining, TEMPORAL_CUTOFF_DAYS)
  temporal_integral <- (1 / beta_decay) * (exp(-beta_decay) - exp(-beta_decay * (effective_days + 1)))
  temporal_integral <- pmax(temporal_integral, 0)

  integral_trig <- sum(alpha * spatial_integral * temporal_integral)

  loglik <- loglik - integral_bg - integral_trig
  if (!is.finite(loglik)) return(-1e10)
  return(loglik)
}

# =============================================================================
# 7. FIT MODEL 0: POISSON
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL 0: POISSON (Background Only)                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Empirical rate
total_district_days <- sum(district_year_months$days_in_month)
emp_rate <- n_total / total_district_days
cat(sprintf("Empirical rate: %.6f per district-day\n\n", emp_rate))

# Starting values and bounds
mean_log_pop_grid <- mean(district_year_months$log_pop)
beta0_start <- log(emp_rate)

start_m0 <- list(
  c(beta0_start - 0.1 * mean_log_pop_grid, 0.1, 0.0, 0.0, rep(0, 9)),
  c(beta0_start - 0.2 * mean_log_pop_grid, 0.2, -1.0, 0.0, rep(0, 9)),
  c(beta0_start - 0.05 * mean_log_pop_grid, 0.05, 1.0, -0.5, rep(0, 9))
)

lower_m0 <- c(-20, -3, -10, -10, rep(-5, 9))
upper_m0 <- c(10, 3, 10, 10, rep(5, 9))

best_ll_m0 <- -Inf
best_result_m0 <- NULL

cat("Running optimization...\n")
for (s in seq_along(start_m0)) {
  cat(sprintf("  Start %d/%d... ", s, length(start_m0)))

  result <- tryCatch({
    optim(par = start_m0[[s]], fn = function(p) -poisson_loglik(p),
          method = "L-BFGS-B", lower = lower_m0, upper = upper_m0,
          control = list(maxit = MAX_ITER))
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
cat(sprintf("  γ (population) = %.4f\n", params_m0[2]))
cat(sprintf("  δ (poverty) = %.4f\n", params_m0[3]))
cat(sprintf("  ζ (log_cpi) = %.4f\n\n", params_m0[4]))

# =============================================================================
# 8. FIT MODEL 1: SPATIO-TEMPORAL HAWKES
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL 1: SPATIO-TEMPORAL HAWKES (Cross-Day + Distance)     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Starting points
# params[14] = log(α), params[15] = log(β), params[16] = log(σ)
# Reasonable starting guesses:
#   α ~ 0.1-0.5 (triggering intensity)
#   β ~ 0.1-0.5 (half-life 1-7 days)
#   σ ~ 50-200 km (spatial bandwidth)

start_m1 <- list(
  c(params_m0[1:13], log(0.1), log(0.2), log(100)),   # half-life ~3.5 days, 100km
  c(params_m0[1:13], log(0.2), log(0.1), log(50)),    # half-life ~7 days, 50km
  c(params_m0[1:13], log(0.05), log(0.3), log(150)),  # half-life ~2.3 days, 150km
  c(params_m0[1:13], log(0.15), log(0.15), log(75)),  # half-life ~4.6 days, 75km
  c(params_m0[1:13], log(0.3), log(0.5), log(100))    # half-life ~1.4 days, 100km
)

# Bounds:
# log(α) ∈ [-5, 3] → α ∈ [0.007, 20]
# log(β) ∈ [-3, 2] → β ∈ [0.05, 7.4] → half-life ∈ [0.09, 14] days
# log(σ) ∈ [2, 6] → σ ∈ [7, 403] km
lower_m1 <- c(lower_m0, -5, -3, 2)
upper_m1 <- c(upper_m0, 3, 2, 6)

best_ll_m1 <- -Inf
best_result_m1 <- NULL

cat("Running multi-start optimization...\n")
cat("(This may take several minutes per start)\n\n")

for (s in seq_along(start_m1)) {
  cat(sprintf("  Start %d/%d... ", s, length(start_m1)))
  start_time <- Sys.time()

  result <- tryCatch({
    optim(par = start_m1[[s]], fn = function(p) -spatiotemporal_hawkes_loglik(p),
          method = "L-BFGS-B", lower = lower_m1, upper = upper_m1,
          control = list(maxit = MAX_ITER))
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    NULL
  })

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

  if (!is.null(result)) {
    ll <- -result$value
    cat(sprintf("LL=%.1f, conv=%d (%.1f min)\n", ll, result$convergence, elapsed))
    if (ll > best_ll_m1) {
      best_ll_m1 <- ll
      best_result_m1 <- result
    }
  } else {
    cat(sprintf("FAILED (%.1f min)\n", elapsed))
  }
}

params_m1 <- best_result_m1$par
ll_m1 <- -best_result_m1$value

alpha_m1 <- exp(params_m1[14])
beta_decay_m1 <- exp(params_m1[15])
sigma_m1 <- exp(params_m1[16])
halflife_m1 <- log(2) / beta_decay_m1

cat(sprintf("\nModel 1 Results (Spatio-Temporal Hawkes):\n"))
cat(sprintf("  Log-likelihood: %.2f\n", ll_m1))
cat(sprintf("  Convergence code: %d\n", best_result_m1$convergence))
cat(sprintf("  γ (population) = %.4f\n", params_m1[2]))
cat(sprintf("  δ (poverty) = %.4f\n", params_m1[3]))
cat(sprintf("  ζ (log_cpi) = %.4f\n", params_m1[4]))
cat(sprintf("  α (triggering intensity) = %.4f\n", alpha_m1))
cat(sprintf("  β (temporal decay) = %.4f per day\n", beta_decay_m1))
cat(sprintf("  Half-life = %.2f days\n", halflife_m1))
cat(sprintf("  σ (spatial bandwidth) = %.1f km\n", sigma_m1))

# Compute approximate branching ratio
# E[offspring] ≈ α × (spatial integral) × (temporal integral per event)
# This is approximate because it depends on the spatial distribution of events
spatial_integral <- 2 * pi * sigma_m1^2 * (1 - exp(-SPATIAL_CUTOFF_KM^2 / (2 * sigma_m1^2)))
temporal_integral <- 1 / beta_decay_m1  # For infinite horizon
branching_approx <- alpha_m1 * spatial_integral * temporal_integral / (pi * SPATIAL_CUTOFF_KM^2)
# Normalize by area to get per-unit-area branching
cat(sprintf("  Approximate branching (normalized): %.4f\n", branching_approx))

if (halflife_m1 > 50) {
  cat("  ⚠️  WARNING: Half-life > 50 days (may indicate model issues)\n")
}
cat("\n")

# =============================================================================
# 9. MODEL COMPARISON
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL COMPARISON                                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

n_params_m0 <- 13
n_params_m1 <- 16

aic_m0 <- -2 * ll_m0 + 2 * n_params_m0
aic_m1 <- -2 * ll_m1 + 2 * n_params_m1
bic_m0 <- -2 * ll_m0 + n_params_m0 * log(n_total)
bic_m1 <- -2 * ll_m1 + n_params_m1 * log(n_total)

lr_stat <- 2 * (ll_m1 - ll_m0)
df <- n_params_m1 - n_params_m0
p_value <- pchisq(lr_stat, df = df, lower.tail = FALSE)

cat("Model fit summary:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-35s %10s %8s %12s %12s\n", "Model", "Log-Lik", "Params", "AIC", "BIC"))
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-35s %10.1f %8d %12.1f %12.1f\n", "M0: Poisson", ll_m0, n_params_m0, aic_m0, bic_m0))
cat(sprintf("%-35s %10.1f %8d %12.1f %12.1f\n", "M1: Spatio-Temporal Hawkes", ll_m1, n_params_m1, aic_m1, bic_m1))
cat("─────────────────────────────────────────────────────────────────\n\n")

cat("Likelihood Ratio Test:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("M1 vs M0 (Does spatio-temporal triggering exist?)\n"))
cat(sprintf("  LR statistic: %.2f, df: %d, p-value: %.2e\n", lr_stat, df, p_value))
if (p_value < 0.001) {
  cat("  *** SIGNIFICANT: Spatio-temporal triggering improves fit ***\n\n")
} else {
  cat("  Not significant at α = 0.001\n\n")
}

# =============================================================================
# 10. INTERPRETATION
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  INTERPRETATION                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("SPATIAL KERNEL INTERPRETATION:\n")
cat(sprintf("  σ = %.1f km (spatial bandwidth)\n", sigma_m1))
cat(sprintf("  At distance d = σ: kernel = exp(-0.5) ≈ 0.61 (61%% of max)\n"))
cat(sprintf("  At distance d = 2σ (%.0f km): kernel ≈ 0.14 (14%% of max)\n", 2*sigma_m1))
cat(sprintf("  At distance d = 3σ (%.0f km): kernel ≈ 0.01 (1%% of max)\n\n", 3*sigma_m1))

cat("TEMPORAL KERNEL INTERPRETATION:\n")
cat(sprintf("  Half-life = %.2f days\n", halflife_m1))
cat(sprintf("  After 1 day: %.1f%% of triggering potential remains\n", 100 * exp(-beta_decay_m1)))
cat(sprintf("  After 7 days: %.1f%% of triggering potential remains\n", 100 * exp(-beta_decay_m1 * 7)))
cat(sprintf("  After 30 days: %.1f%% of triggering potential remains\n\n", 100 * exp(-beta_decay_m1 * 30)))

# =============================================================================
# 11. SAVE RESULTS
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SAVING RESULTS                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

results <- list(
  model_type = "Spatio-Temporal Hawkes (Cross-Day + Distance Decay)",
  config = list(
    spatial_cutoff_km = SPATIAL_CUTOFF_KM,
    temporal_cutoff_days = TEMPORAL_CUTOFF_DAYS
  ),
  n_events = n_total,
  n_days = D_max,
  neighbor_stats = list(
    n_with_neighbors = n_with_neighbors,
    avg_neighbors = avg_neighbors,
    max_neighbors = max_neighbors
  ),

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
    name = "Spatio-Temporal Hawkes",
    params = params_m1,
    loglik = ll_m1,
    aic = aic_m1,
    bic = bic_m1,
    gamma = params_m1[2],
    delta = params_m1[3],
    zeta = params_m1[4],
    alpha = alpha_m1,
    beta_decay = beta_decay_m1,
    halflife = halflife_m1,
    sigma = sigma_m1
  ),

  comparison = list(
    lr_statistic = lr_stat,
    df = df,
    p_value = p_value
  )
)

saveRDS(results, "model_results_spatiotemporal.rds")
cat("Saved: model_results_spatiotemporal.rds\n\n")

# =============================================================================
# 12. COMPARISON TO OTHER MODELS
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  COMPARISON ACROSS ALL MODEL SPECIFICATIONS                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Model               | Half-life | σ (km) | γ (pop) | δ (poverty)\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("Jittered (1.0 day)  | %7.2f   |   —    | %7.2f | %7.2f\n", 1.14, 2.56, -9.36))
cat(sprintf("Cross-Day Only      | %7.2f   |   —    | %7.2f | %7.2f\n", 693.15, -0.74, -1.21))
cat(sprintf("Spatio-Temporal     | %7.2f   | %5.1f  | %7.2f | %7.2f\n",
            halflife_m1, sigma_m1, params_m1[2], params_m1[3]))
cat("─────────────────────────────────────────────────────────────────\n\n")

cat("KEY FINDINGS:\n")
if (halflife_m1 < 30 && halflife_m1 > 0.5) {
  cat(sprintf("  ✓ Half-life (%.1f days) is now in a reasonable range\n", halflife_m1))
  cat("    This suggests spatial structure was key to identifying temporal dynamics\n\n")
} else if (halflife_m1 >= 30) {
  cat(sprintf("  ⚠️  Half-life (%.1f days) is still quite long\n", halflife_m1))
  cat("    May indicate persistent regional effects or remaining model issues\n\n")
}

if (sigma_m1 < 200) {
  cat(sprintf("  ✓ Spatial bandwidth (%.0f km) suggests LOCAL diffusion\n", sigma_m1))
  cat("    Protests trigger nearby protests, not distant ones\n")
} else {
  cat(sprintf("  ⚠️  Spatial bandwidth (%.0f km) is quite large\n", sigma_m1))
  cat("    May indicate broader regional coordination\n")
}

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  ESTIMATION COMPLETE                                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
