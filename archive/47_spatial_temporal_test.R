# ============================================================================
# 47: TEST SPATIAL-TEMPORAL HAWKES VS TEMPORAL-ONLY
# ============================================================================
#
# PURPOSE:
# Test whether adding a spatial component improves the Hawkes model.
# Hypothesis: Protests don't spread spatially, so spatial-temporal models
# won't improve log-likelihood over purely temporal models.
#
# COMPARISON:
# - Model A: Univariate temporal Hawkes (model1_basic_hawkes.rds, LL ≈ -79,441)
# - Model B: Univariate spatial-temporal Hawkes (to be fitted)
#
# SPATIAL SPECIFICATION:
# - Gaussian kernel: g(d) = exp(-d²/(2σ²))
# - Regional scale: σ ≈ 100 km (test 50-200 km range)
# - Events beyond 200 km contribute negligibly
#
# LR TEST:
# - H₀: Spatial parameter adds no explanatory power
# - df = 1 (adding σ parameter)
# - If LR not significant → spatial spread doesn't matter
#
# ============================================================================

library(dplyr)
library(Rcpp)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   SPATIAL-TEMPORAL VS TEMPORAL HAWKES: TEST FOR SPATIAL SPREAD          ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("HYPOTHESIS:\n")
cat("  H₀: Protests do not spread spatially (σ → ∞ or spatial term ≈ 0)\n")
cat("  H₁: Protests spread spatially with finite range σ\n\n")

cat("SPECIFICATION:\n")
cat("  Temporal:        λ(t) = μ(d,y) + Σ α·exp(-β·τ)\n")
cat("  Spatial-Temporal: λ(s,t) = μ(d,y) + Σ α·g(d)·exp(-β·τ)\n")
cat("  where g(d) = exp(-d²/(2σ²))  [Gaussian spatial kernel]\n\n")

# ====================================================================
# CONFIGURATION
# ====================================================================

TEMPORAL_CUTOFF <- 90   # Days (same as temporal model)
SPATIAL_CUTOFF <- 200   # km (regional scale)
SIGMA_INIT <- 100       # Initial spatial scale (km)

# ====================================================================
# HAVERSINE DISTANCE FUNCTION
# ====================================================================

haversine_km <- function(lon1, lat1, lon2, lat2) {
  # Convert degrees to radians
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
  d <- R * c

  return(d)
}

# ====================================================================
# PARAMETER TRANSFORMATION FUNCTIONS (same as model1)
# ====================================================================

transform_gamma <- function(gamma_raw) {
  0.01 + 1.99 * plogis(gamma_raw)
}
inverse_transform_gamma <- function(gamma) {
  qlogis((gamma - 0.01) / 1.99)
}

transform_delta <- function(delta_raw) {
  -20 + 40 * plogis(delta_raw)
}
inverse_transform_delta <- function(delta) {
  qlogis((delta + 20) / 40)
}

transform_beta_year <- function(beta_year_raw) {
  -3 + 6 * plogis(beta_year_raw)
}
inverse_transform_beta_year <- function(beta_year) {
  qlogis((beta_year + 3) / 6)
}

transform_beta_0_trig <- function(beta_0_trig_raw) {
  -10 + 20 * plogis(beta_0_trig_raw)
}
inverse_transform_beta_0_trig <- function(beta_0_trig) {
  qlogis((beta_0_trig + 10) / 20)
}

transform_decay <- function(decay_raw) {
  0.001 + 9.999 * plogis(decay_raw)
}
inverse_transform_decay <- function(decay) {
  qlogis((decay - 0.001) / 9.999)
}

# Transform sigma_raw → sigma ∈ [1, 500] km
transform_sigma <- function(sigma_raw) {
  1 + 499 * plogis(sigma_raw)
}
inverse_transform_sigma <- function(sigma) {
  qlogis((sigma - 1) / 499)
}

# ====================================================================
# 1. LOAD DATA AND TEMPORAL MODEL
# ====================================================================

cat("=== LOADING DATA AND TEMPORAL MODEL ===\n\n")

# Load data
protests <- readRDS("protests_with_poverty.rds")
cat(sprintf("Loaded %d events\n", nrow(protests)))

# Prepare data with coordinates
marks_data <- protests %>%
  select(
    event_id = event_id_cnty,
    time = days_since_start,
    district,
    year,
    log_pop,
    population,
    poverty_decimal,
    longitude,
    latitude
  ) %>%
  arrange(time) %>%
  filter(!is.na(poverty_decimal) & !is.na(longitude) & !is.na(latitude))

cat(sprintf("Events with complete data (poverty + coordinates): %d\n", nrow(marks_data)))

# Load temporal model
model_temporal <- readRDS("model1_basic_hawkes.rds")
cat(sprintf("\nTemporal model log-likelihood: %.2f\n", model_temporal$loglik))
cat(sprintf("Temporal model parameters: %d\n", model_temporal$n_params))

# Extract temporal model parameters
temporal_params <- model_temporal$params
gamma_temporal <- transform_gamma(temporal_params$gamma_raw)
delta_temporal <- transform_delta(temporal_params$delta_raw)
alpha_temporal <- exp(transform_beta_0_trig(temporal_params$beta_0_trig_raw))
decay_temporal <- transform_decay(temporal_params$decay_raw)

cat(sprintf("\nTemporal model estimates:\n"))
cat(sprintf("  γ (population): %.4f\n", gamma_temporal))
cat(sprintf("  δ (poverty): %.4f\n", delta_temporal))
cat(sprintf("  α (triggering): %.4f\n", alpha_temporal))
cat(sprintf("  β (decay): %.4f (half-life: %.1f days)\n\n",
            decay_temporal, log(2)/decay_temporal))

# Prepare district-years
T_min <- min(marks_data$time)
T_max <- max(marks_data$time)

district_years <- marks_data %>%
  select(district, year, population, log_pop, poverty_decimal) %>%
  distinct()
district_years$days_observed <- T_max - T_min

# ====================================================================
# 2. COMPUTE DISTANCE MATRIX
# ====================================================================

cat("=== COMPUTING DISTANCE MATRIX ===\n")

times <- marks_data$time
lons <- marks_data$longitude
lats <- marks_data$latitude
n <- length(times)

cat(sprintf("Computing pairwise distances for %d events...\n", n))
cat(sprintf("Using spatial cutoff: %d km\n", SPATIAL_CUTOFF))

# Pre-compute distance matrix (only store distances within cutoff)
# This is memory-efficient: store as sparse representation
start_time <- Sys.time()

# For efficiency, compute distances on-the-fly during likelihood
# but pre-identify potential spatial neighbors

cat("Pre-computing spatial neighbor indices...\n")

# Sample-based approach for efficiency: pre-compute for subset to estimate sparsity
set.seed(42)
sample_idx <- sample(min(500, n), min(500, n))
n_within_cutoff <- 0
n_total_pairs <- 0

for(i in sample_idx) {
  if(i <= 1) next
  for(j in 1:(i-1)) {
    if(is.na(lons[j]) || is.na(lats[j]) || is.na(lons[i]) || is.na(lats[i])) next
    d <- haversine_km(lons[j], lats[j], lons[i], lats[i])
    n_total_pairs <- n_total_pairs + 1
    if(!is.na(d) && d <= SPATIAL_CUTOFF) {
      n_within_cutoff <- n_within_cutoff + 1
    }
  }
}

sparsity <- n_within_cutoff / n_total_pairs
cat(sprintf("Estimated sparsity: %.2f%% pairs within %d km\n",
            100 * sparsity, SPATIAL_CUTOFF))

# Store distances in compressed format
cat("Computing full distance matrix (compressed)...\n")

# For each event i, store list of (j, distance) pairs where j < i and d <= cutoff
spatial_neighbors <- vector("list", n)
spatial_distances <- vector("list", n)

progress_interval <- max(1, n %/% 20)
for(i in 2:n) {
  neighbors_i <- c()
  distances_i <- c()

  # Skip if current event has NA coordinates
  if(is.na(lons[i]) || is.na(lats[i])) {
    spatial_neighbors[[i]] <- neighbors_i
    spatial_distances[[i]] <- distances_i
    next
  }

  for(j in 1:(i-1)) {
    # Skip if past event has NA coordinates
    if(is.na(lons[j]) || is.na(lats[j])) next

    # Only consider if temporal also within cutoff
    tau <- times[i] - times[j]
    if(tau > 0 && tau <= TEMPORAL_CUTOFF) {
      d <- haversine_km(lons[j], lats[j], lons[i], lats[i])
      if(!is.na(d) && d <= SPATIAL_CUTOFF) {
        neighbors_i <- c(neighbors_i, j)
        distances_i <- c(distances_i, d)
      }
    }
  }

  spatial_neighbors[[i]] <- neighbors_i
  spatial_distances[[i]] <- distances_i

  if(i %% progress_interval == 0) {
    cat(sprintf("  Progress: %d/%d (%.0f%%)  \r", i, n, 100*i/n))
  }
}
cat(sprintf("  Progress: %d/%d (100%%)  \n", n, n))

# Summary statistics
avg_neighbors <- mean(sapply(2:n, function(i) length(spatial_neighbors[[i]])))
cat(sprintf("\nAverage spatiotemporal neighbors per event: %.1f\n", avg_neighbors))

runtime_dist <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("Distance computation time: %.1f minutes\n\n", runtime_dist))

# ====================================================================
# 3. SPATIAL-TEMPORAL LIKELIHOOD FUNCTION
# ====================================================================

cat("=== IMPLEMENTING SPATIAL-TEMPORAL LIKELIHOOD ===\n\n")

# Global state for optimization
.optim_state <- new.env()
.optim_state$n_evals <- 0
.optim_state$best_ll <- -Inf

reset_optim_state <- function() {
  .optim_state$n_evals <- 0
  .optim_state$best_ll <- -Inf
}

loglik_spatial_temporal <- function(params, times, marks, district_years,
                                    spatial_neighbors, spatial_distances,
                                    verbose = TRUE) {

  .optim_state$n_evals <- .optim_state$n_evals + 1
  eval_num <- .optim_state$n_evals

  n <- length(times)

  # Extract and transform parameters
  beta_0_bg <- params[["beta_0_bg"]]
  gamma <- transform_gamma(params[["gamma_raw"]])
  delta <- transform_delta(params[["delta_raw"]])

  # Year effects
  beta_years <- c(
    `2015` = 0,
    `2016` = transform_beta_year(params[["beta_2016_raw"]]),
    `2017` = transform_beta_year(params[["beta_2017_raw"]]),
    `2018` = transform_beta_year(params[["beta_2018_raw"]]),
    `2019` = transform_beta_year(params[["beta_2019_raw"]]),
    `2020` = transform_beta_year(params[["beta_2020_raw"]]),
    `2021` = transform_beta_year(params[["beta_2021_raw"]]),
    `2022` = transform_beta_year(params[["beta_2022_raw"]]),
    `2023` = transform_beta_year(params[["beta_2023_raw"]]),
    `2024` = transform_beta_year(params[["beta_2024_raw"]])
  )

  # Triggering parameters
  alpha <- exp(transform_beta_0_trig(params[["beta_0_trig_raw"]]))
  decay <- transform_decay(params[["decay_raw"]])
  sigma <- transform_sigma(params[["sigma_raw"]])

  # Gaussian spatial kernel
  spatial_kernel <- function(d, sigma) {
    exp(-d^2 / (2 * sigma^2))
  }

  # ===== Part 1: Sum of log intensities at event times =====
  loglik <- 0

  for(i in 1:n) {
    # Background rate
    year_i <- as.character(marks$year[i])
    beta_year_i <- beta_years[year_i]
    if(is.na(beta_year_i)) beta_year_i <- 0

    mu <- exp(beta_0_bg + gamma * marks$log_pop[i] +
              delta * marks$poverty_decimal[i] + beta_year_i)

    # Triggering from past events (spatial-temporal)
    trigger_sum <- 0
    if(i > 1 && length(spatial_neighbors[[i]]) > 0) {
      neighbors_j <- spatial_neighbors[[i]]
      distances_j <- spatial_distances[[i]]

      for(k in seq_along(neighbors_j)) {
        j <- neighbors_j[k]
        d <- distances_j[k]
        tau <- times[i] - times[j]

        # Spatial kernel
        g_spatial <- spatial_kernel(d, sigma)

        # Temporal kernel
        h_temporal <- exp(-decay * tau)

        trigger_sum <- trigger_sum + alpha * g_spatial * h_temporal
      }
    }

    # Total intensity
    lambda <- mu + trigger_sum

    if(lambda > 0) {
      loglik <- loglik + log(lambda)
    } else {
      return(-Inf)
    }
  }

  # ===== Part 2: Compensator =====

  # Background compensator
  background_comp <- 0
  for(k in 1:nrow(district_years)) {
    year_k <- as.character(district_years$year[k])
    beta_year_k <- beta_years[year_k]
    if(is.na(beta_year_k)) beta_year_k <- 0

    mu_k <- exp(beta_0_bg + gamma * district_years$log_pop[k] +
                delta * district_years$poverty_decimal[k] + beta_year_k)

    background_comp <- background_comp + mu_k * district_years$days_observed[k]
  }

  # Triggering compensator
  # Key insight: the spatial kernel just weights WHERE triggered events occur,
  # not HOW MANY. The expected number of offspring per parent remains α/β.
  # This is the same as the temporal model.
  #
  # The spatial weighting redistributes triggered events spatially, but
  # the total expected count is the same. Therefore, we use the same
  # compensator structure as the temporal model.

  trigger_comp <- 0
  for(i in 1:n) {
    tau_remain <- min(T_max - times[i], TEMPORAL_CUTOFF)
    if(tau_remain > 0) {
      # Temporal integral only (same as temporal model)
      temporal_int <- (1 - exp(-decay * tau_remain)) / decay
      trigger_comp <- trigger_comp + alpha * temporal_int
    }
  }

  compensator <- background_comp + trigger_comp
  ll_final <- loglik - compensator

  # Track best
  if(ll_final > .optim_state$best_ll) {
    .optim_state$best_ll <- ll_final
  }

  if(verbose && eval_num %% 10 == 0) {
    cat(sprintf("  Eval %d: LL=%.2f | σ=%.1f km, α=%.4f, β=%.4f\n",
                eval_num, ll_final, sigma, alpha, decay))
  }

  return(ll_final)
}

# ====================================================================
# 4. FIT SPATIAL-TEMPORAL MODEL
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FITTING SPATIAL-TEMPORAL MODEL                             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Initialize from temporal model parameters + spatial parameter
init_params <- list(
  beta_0_bg = temporal_params$beta_0_bg,
  gamma_raw = temporal_params$gamma_raw,
  delta_raw = temporal_params$delta_raw,
  beta_2016_raw = temporal_params$beta_2016_raw,
  beta_2017_raw = temporal_params$beta_2017_raw,
  beta_2018_raw = temporal_params$beta_2018_raw,
  beta_2019_raw = temporal_params$beta_2019_raw,
  beta_2020_raw = temporal_params$beta_2020_raw,
  beta_2021_raw = temporal_params$beta_2021_raw,
  beta_2022_raw = temporal_params$beta_2022_raw,
  beta_2023_raw = temporal_params$beta_2023_raw,
  beta_2024_raw = temporal_params$beta_2024_raw,
  beta_0_trig_raw = temporal_params$beta_0_trig_raw,
  decay_raw = temporal_params$decay_raw,
  sigma_raw = inverse_transform_sigma(SIGMA_INIT)  # Start at 100 km
)

param_names <- names(init_params)
par0 <- unlist(init_params)

# Bounds
lower_bounds <- c(
  beta_0_bg = -15, gamma_raw = -10, delta_raw = -10,
  beta_2016_raw = -10, beta_2017_raw = -10, beta_2018_raw = -10,
  beta_2019_raw = -10, beta_2020_raw = -10, beta_2021_raw = -10,
  beta_2022_raw = -10, beta_2023_raw = -10, beta_2024_raw = -10,
  beta_0_trig_raw = -10, decay_raw = -10, sigma_raw = -10
)

upper_bounds <- c(
  beta_0_bg = 5, gamma_raw = 10, delta_raw = 10,
  beta_2016_raw = 10, beta_2017_raw = 10, beta_2018_raw = 10,
  beta_2019_raw = 10, beta_2020_raw = 10, beta_2021_raw = 10,
  beta_2022_raw = 10, beta_2023_raw = 10, beta_2024_raw = 10,
  beta_0_trig_raw = 10, decay_raw = 10, sigma_raw = 10
)

# Objective function
obj_fun <- function(par) {
  params_list <- as.list(par)
  names(params_list) <- param_names
  ll <- loglik_spatial_temporal(params_list, times, marks_data, district_years,
                                 spatial_neighbors, spatial_distances,
                                 verbose = TRUE)
  return(-ll)
}

cat("Starting optimization...\n")
cat(sprintf("Initial σ = %.1f km\n\n", transform_sigma(init_params$sigma_raw)))

reset_optim_state()
start_time <- Sys.time()

fit_result <- optim(
  par = par0,
  fn = obj_fun,
  method = "L-BFGS-B",
  lower = lower_bounds[param_names],
  upper = upper_bounds[param_names],
  control = list(
    maxit = 500,
    factr = 1e7,
    trace = 1,
    REPORT = 20
  )
)

runtime <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\n✓ Optimization complete in %.1f minutes\n", runtime))

# Extract results
params_hat <- as.list(fit_result$par)
names(params_hat) <- param_names
loglik_spatial <- -fit_result$value

sigma_hat <- transform_sigma(params_hat$sigma_raw)
gamma_hat <- transform_gamma(params_hat$gamma_raw)
delta_hat <- transform_delta(params_hat$delta_raw)
alpha_hat <- exp(transform_beta_0_trig(params_hat$beta_0_trig_raw))
decay_hat <- transform_decay(params_hat$decay_raw)

# ====================================================================
# 5. LIKELIHOOD RATIO TEST
# ====================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  LIKELIHOOD RATIO TEST: SPATIAL COMPONENT                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

ll_temporal <- model_temporal$loglik
ll_spatial <- loglik_spatial

# LR statistic
lr_statistic <- 2 * (ll_spatial - ll_temporal)
df <- 1  # Adding sigma parameter
p_value <- pchisq(lr_statistic, df, lower.tail = FALSE)

cat("Model Comparison:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Temporal Hawkes:        LL = %.2f   (%d params)\n",
            ll_temporal, model_temporal$n_params))
cat(sprintf("  Spatial-Temporal:       LL = %.2f   (%d params)\n",
            ll_spatial, model_temporal$n_params + 1))
cat(sprintf("  Improvement:            ΔLL = %.2f\n", ll_spatial - ll_temporal))
cat("\n")

cat("LR Test (spatial component):\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  LR statistic: %.2f\n", lr_statistic))
cat(sprintf("  df: %d\n", df))
cat(sprintf("  p-value: %.4f", p_value))
if(p_value < 0.001) cat(" ***")
else if(p_value < 0.01) cat(" **")
else if(p_value < 0.05) cat(" *")
cat("\n\n")

cat("Estimated spatial scale:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  σ = %.1f km\n", sigma_hat))
cat(sprintf("  At d = σ: kernel = %.3f (%.1f%% of max)\n",
            exp(-0.5), 100*exp(-0.5)))
cat(sprintf("  At d = 2σ: kernel = %.3f (%.1f%% of max)\n",
            exp(-2), 100*exp(-2)))
cat("\n")

# ====================================================================
# 6. CONCLUSION
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
if(p_value >= 0.05) {
  cat("║  CONCLUSION: SPATIAL SPREAD DOES NOT MATTER                 ║\n")
} else {
  cat("║  CONCLUSION: SPATIAL SPREAD IS SIGNIFICANT                  ║\n")
}
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

if(p_value >= 0.05) {
  cat("The LR test is NOT significant at α = 0.05.\n")
  cat("Adding a spatial component does not significantly improve model fit.\n")
  cat("This supports the hypothesis that Indonesian protests spread primarily\n")
  cat("through temporal contagion (national-level triggering) rather than\n")
  cat("spatial diffusion to nearby locations.\n")
} else {
  cat("The LR test IS significant at α = 0.05.\n")
  cat("Adding a spatial component significantly improves model fit.\n")
  cat("This suggests that protests do spread spatially to nearby areas.\n")
}

# ====================================================================
# 7. PARAMETER COMPARISON
# ====================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PARAMETER COMPARISON                                        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

comparison_params <- data.frame(
  Parameter = c("γ (population)", "δ (poverty)", "α (triggering)",
                "β (decay)", "σ (spatial km)"),
  Temporal = c(
    sprintf("%.4f", gamma_temporal),
    sprintf("%.4f", delta_temporal),
    sprintf("%.4f", alpha_temporal),
    sprintf("%.4f", decay_temporal),
    "—"
  ),
  `Spatial-Temporal` = c(
    sprintf("%.4f", gamma_hat),
    sprintf("%.4f", delta_hat),
    sprintf("%.4f", alpha_hat),
    sprintf("%.4f", decay_hat),
    sprintf("%.1f", sigma_hat)
  ),
  check.names = FALSE
)

print(comparison_params)

# ====================================================================
# 8. SAVE RESULTS
# ====================================================================

cat("\n=== SAVING RESULTS ===\n")

results <- list(
  # Model comparison
  ll_temporal = ll_temporal,
  ll_spatial = ll_spatial,
  lr_statistic = lr_statistic,
  df = df,
  p_value = p_value,

  # Spatial-temporal model
  params_spatial = params_hat,
  sigma_km = sigma_hat,
  convergence = fit_result$convergence,
  runtime_mins = runtime,

  # Original temporal model
  params_temporal = temporal_params,

  # Conclusion
  spatial_significant = p_value < 0.05
)

saveRDS(results, "spatial_temporal_test.rds")
cat("✓ Saved: spatial_temporal_test.rds\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SPATIAL-TEMPORAL TEST COMPLETE                             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
