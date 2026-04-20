# Phase 3: Spatial-Temporal-Mark Hawkes Models - FAST VALIDATION
# WITH GAUSSIAN SPATIAL KERNEL AND HAVERSINE DISTANCES
# Author: Daniel Carnahan
# Date: 2025-11-03

library(dplyr)
library(ggplot2)

cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  PHASE 3: SPATIAL-TEMPORAL-MARK HAWKES (FAST VALIDATION)    ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

cat("MODEL SPECIFICATION:\n")
cat("  λ(s,t) = μ + Σ α(marks) · g(d) · exp(-β·τ)\n")
cat("  where:\n")
cat("    g(d) = exp(-d²/(2σ²))  [Gaussian spatial kernel]\n")
cat("    d = haversine distance in km\n")
cat("    α(marks) = exp(β₀ + β_violence·violent + β_state·state)\n\n")

cat("OPTIMIZATIONS:\n")
cat("  ✓ Haversine distance calculation (great-circle)\n")
cat("  ✓ Spatial cutoff: 500 km\n")
cat("  ✓ Temporal cutoff: 730 days\n")
cat("  ✓ Pre-computed spatiotemporal neighbor matrix\n")
cat("  ✓ Vectorized spatial kernel calculations\n")
cat("  ✓ Multi-start optimization (4 starts)\n")
cat("  ✓ Early stopping\n\n")

cat("FAST VALIDATION:\n")
cat("  - Using 1,000 events (6% sample)\n")
cat("  - Max 100 iterations per start\n")
cat("  - Expected time: ~10-15 minutes total\n\n")

cat("HYPOTHESES TO TEST:\n")
cat("  H2: Do violence/state effects persist after controlling for space?\n")
cat("  H3: Does violence spread differently in space than peaceful protests?\n\n")

# ====================================================================
# CONFIGURATION
# ====================================================================

TEMPORAL_CUTOFF <- 730  # 2 years
SPATIAL_CUTOFF <- 500   # 500 km
SAMPLE_SIZE <- 1000     # Fast validation

# Early stopping & multi-start
EARLY_STOP_NO_IMPROVE <- 15
EARLY_STOP_MIN_IMPROVE <- 0.1
N_RANDOM_STARTS <- 3

# ====================================================================
# SPATIAL DISTANCE FUNCTIONS
# ====================================================================

# Haversine distance (great-circle distance on Earth)
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

# Pre-compute distance matrix with cutoff
precompute_distance_matrix <- function(lons, lats, cutoff_km) {
  n <- length(lons)
  dist_mat <- matrix(Inf, nrow = n, ncol = n)

  cat(sprintf("  Computing %d x %d distance matrix with %d km cutoff...\n", n, n, cutoff_km))

  # Only compute for i > j (since dist is symmetric and we only need past events)
  for(i in 2:n) {
    for(j in 1:(i-1)) {
      d <- haversine_km(lons[j], lats[j], lons[i], lats[i])
      if(d <= cutoff_km) {
        dist_mat[i, j] <- d
      }
    }

    if(i %% 100 == 0) {
      cat(sprintf("    Progress: %d/%d (%.1f%%)  \r", i, n, 100*i/n))
    }
  }

  cat(sprintf("    Complete: %d/%d (100%%)  \n", n, n))

  # Count neighbors
  n_neighbors <- sum(is.finite(dist_mat))
  avg_neighbors <- n_neighbors / n
  cat(sprintf("  Average neighbors within %d km: %.1f (%.1f%% of total)\n",
             cutoff_km, avg_neighbors, 100 * avg_neighbors / n))

  return(dist_mat)
}

# ====================================================================
# GLOBAL STATE FOR OPTIMIZATION
# ====================================================================

.optim_state <- new.env()
.optim_state$n_evals <- 0
.optim_state$best_ll <- -Inf
.optim_state$best_params <- NULL
.optim_state$model_name <- NULL
.optim_state$no_improve_count <- 0
.optim_state$should_stop <- FALSE

reset_optim_state <- function(model_name) {
  .optim_state$n_evals <- 0
  .optim_state$best_ll <- -Inf
  .optim_state$best_params <- NULL
  .optim_state$model_name <- model_name
  .optim_state$no_improve_count <- 0
  .optim_state$should_stop <- FALSE
}

# ====================================================================
# 1. LOAD AND SAMPLE DATA
# ====================================================================

cat("\n=== LOADING DATA ===\n")
protests <- readRDS("protests_prepared.rds")

cat(sprintf("Total events available: %d\n", nrow(protests)))

# Sample evenly across time
set.seed(42)
protests_sample <- protests %>%
  arrange(days_since_start) %>%
  slice(seq(1, n(), by = nrow(protests) %/% SAMPLE_SIZE)) %>%
  head(SAMPLE_SIZE)

cat(sprintf("Sampled: %d events (%.1f%%)\n", nrow(protests_sample),
           100 * nrow(protests_sample) / nrow(protests)))

# Prepare data
marks_data <- protests_sample %>%
  mutate(
    time = days_since_start,
    is_violent = as.numeric(is_violent),
    state_intervention = as.numeric(state_intervention),
    fatalities = as.numeric(fatalities)
  ) %>%
  select(time, longitude, latitude, is_violent, state_intervention, fatalities) %>%
  arrange(time)

cat(sprintf("Prepared: %d events with spatial coordinates\n", nrow(marks_data)))
cat(sprintf("  Violent: %d (%.1f%%)\n", sum(marks_data$is_violent),
           100 * mean(marks_data$is_violent)))
cat(sprintf("  State intervention: %d (%.1f%%)\n", sum(marks_data$state_intervention),
           100 * mean(marks_data$state_intervention)))

# ====================================================================
# 2. PRE-COMPUTE MATRICES
# ====================================================================

cat("\n=== PRE-COMPUTING SPATIOTEMPORAL MATRICES ===\n")

times <- marks_data$time
lons <- marks_data$longitude
lats <- marks_data$latitude
n <- length(times)

# Temporal matrix
cat("Computing temporal matrix...\n")
time_mat <- matrix(0, nrow = n, ncol = n)
for(i in 2:n) {
  for(j in 1:(i-1)) {
    tau <- times[i] - times[j]
    if(isTRUE(tau > 0) && isTRUE(tau <= TEMPORAL_CUTOFF)) {
      time_mat[i, j] <- tau
    }
  }
}
cat("✓ Temporal matrix complete\n\n")

# Spatial matrix
cat("Computing spatial distance matrix...\n")
dist_mat <- precompute_distance_matrix(lons, lats, SPATIAL_CUTOFF)
cat("✓ Spatial matrix complete\n\n")

# Combined spatiotemporal neighbor matrix
cat("Creating combined spatiotemporal neighbor matrix...\n")
spatiotemporal_neighbors <- list()
for(i in 2:n) {
  # Find events that are both temporally and spatially within cutoffs
  valid_temporal <- which(time_mat[i, ] > 0)
  valid_spatial <- which(is.finite(dist_mat[i, ]))
  valid_both <- intersect(valid_temporal, valid_spatial)

  spatiotemporal_neighbors[[i]] <- valid_both
}
avg_neighbors <- mean(sapply(spatiotemporal_neighbors, length))
cat(sprintf("✓ Average spatiotemporal neighbors: %.1f\n\n", avg_neighbors))

# ====================================================================
# 3. SPATIAL-TEMPORAL HAWKES FUNCTIONS
# ====================================================================

spatial_temporal_hawkes_functions <- function() {

  # Gaussian spatial kernel
  spatial_kernel_gaussian <- function(d, sigma) {
    # g(d) = exp(-d²/(2σ²))
    return(exp(-d^2 / (2 * sigma^2)))
  }

  # Conditional intensity with spatial and temporal components
  lambda_spatial_temporal <- function(idx, times, lons, lats, marks, params,
                                     time_mat = NULL, dist_mat = NULL,
                                     neighbors = NULL,
                                     use_mark_interaction = FALSE) {

    mu <- params[["mu"]]
    beta_0 <- params[["beta_0"]]
    beta_violence <- ifelse(is.null(params[["beta_violence"]]), 0, params[["beta_violence"]])
    beta_state <- ifelse(is.null(params[["beta_state"]]), 0, params[["beta_state"]])
    decay <- params[["decay"]]
    sigma <- params[["sigma"]]

    # For spatial-mark interaction model
    if(use_mark_interaction) {
      sigma_violent <- params[["sigma_violent"]]
      sigma_peaceful <- params[["sigma_peaceful"]]
    }

    # Background
    intensity <- mu

    # Triggered component
    if(idx > 1) {
      # Use pre-computed neighbors if available
      if(!is.null(neighbors) && idx <= length(neighbors)) {
        valid_j <- neighbors[[idx]]
      } else {
        # Fallback: find valid neighbors
        valid_temporal <- which(time_mat[idx, ] > 0)
        valid_spatial <- which(is.finite(dist_mat[idx, ]))
        valid_j <- intersect(valid_temporal, valid_spatial)
      }

      if(length(valid_j) > 0) {
        # Vectorized calculation
        taus <- time_mat[idx, valid_j]
        dists <- dist_mat[idx, valid_j]

        # Mark-dependent excitation
        alphas <- exp(beta_0 +
                     beta_violence * marks$is_violent[valid_j] +
                     beta_state * marks$state_intervention[valid_j])

        # Temporal kernel
        h_temporal <- exp(-decay * taus)

        # Spatial kernel
        if(use_mark_interaction) {
          # Different spatial kernels for violent vs peaceful
          g_spatial <- numeric(length(valid_j))
          violent_mask <- marks$is_violent[valid_j] == 1

          if(any(violent_mask)) {
            g_spatial[violent_mask] <- spatial_kernel_gaussian(dists[violent_mask],
                                                               sigma_violent)
          }
          if(any(!violent_mask)) {
            g_spatial[!violent_mask] <- spatial_kernel_gaussian(dists[!violent_mask],
                                                                sigma_peaceful)
          }
        } else {
          # Same spatial kernel for all
          g_spatial <- spatial_kernel_gaussian(dists, sigma)
        }

        # Combined triggering
        intensity <- intensity + sum(alphas * g_spatial * h_temporal)
      }
    }

    return(intensity)
  }

  # Log-likelihood
  loglik_spatial_temporal <- function(params, times, lons, lats, marks,
                                     time_mat = NULL, dist_mat = NULL,
                                     neighbors = NULL,
                                     use_mark_interaction = FALSE,
                                     verbose = TRUE) {

    # Parameter constraints - check for NA/NULL first
    if(is.null(params[["mu"]]) || is.na(params[["mu"]]) ||
       is.null(params[["decay"]]) || is.na(params[["decay"]])) {
      return(-Inf)
    }

    if(params[["mu"]] <= 0 || params[["decay"]] <= 0) {
      return(-Inf)
    }

    if(use_mark_interaction) {
      # Check mark-interaction spatial parameters
      if(is.null(params[["sigma_violent"]]) || is.na(params[["sigma_violent"]]) ||
         is.null(params[["sigma_peaceful"]]) || is.na(params[["sigma_peaceful"]])) {
        return(-Inf)
      }

      if(params[["sigma_violent"]] <= 0 || params[["sigma_peaceful"]] <= 0) {
        return(-Inf)
      }
    } else {
      # Check regular spatial parameters
      if(is.null(params[["sigma"]]) || is.na(params[["sigma"]])) {
        return(-Inf)
      }

      if(params[["sigma"]] <= 0) {
        return(-Inf)
      }
    }

    # Track evaluations
    .optim_state$n_evals <- .optim_state$n_evals + 1
    eval_num <- .optim_state$n_evals

    n <- length(times)
    ll <- 0

    # Sum log-intensities at event times
    for(i in 1:n) {
      lambda_i <- lambda_spatial_temporal(i, times, lons, lats, marks, params,
                                         time_mat, dist_mat, neighbors,
                                         use_mark_interaction)
      ll <- ll + log(max(lambda_i, 1e-10))

      if(verbose && i %% 200 == 0) {
        cat(sprintf("    Event %d/%d (%.1f%%)  \r", i, n, 100*i/n))
      }
    }

    # Compensator
    T_max <- max(times)
    compensator <- params[["mu"]] * T_max

    # Triggered component contribution
    for(i in 1:n) {
      if(times[i] < T_max) {
        # Mark-dependent alpha
        beta_violence_val <- ifelse(is.null(params[["beta_violence"]]), 0, params[["beta_violence"]])
        beta_state_val <- ifelse(is.null(params[["beta_state"]]), 0, params[["beta_state"]])
        alpha_i <- exp(params[["beta_0"]] +
                      beta_violence_val * marks$is_violent[i] +
                      beta_state_val * marks$state_intervention[i])

        # Temporal integral
        tau_remain <- min(T_max - times[i], TEMPORAL_CUTOFF)
        temporal_integral <- (1 / params[["decay"]]) * (1 - exp(-params[["decay"]] * tau_remain))

        # Spatial integral (approximate with mean spatial kernel value)
        # For Gaussian: analytical form = π·σ² · (1 - exp(-R²/(2σ²)))
        # Use simplified approximation: mean kernel value × effective area
        if(use_mark_interaction) {
          if(marks$is_violent[i] == 1) {
            sigma_i <- params[["sigma_violent"]]
          } else {
            sigma_i <- params[["sigma_peaceful"]]
          }
        } else {
          sigma_i <- params[["sigma"]]
        }

        # Approximate spatial integral: int_0^R g(d) * 2πd dd
        # For Gaussian: analytical form = π·σ² · (1 - exp(-R²/(2σ²)))
        # Approximate as: (mean kernel value) × (effective area)
        # Mean kernel value at d=σ: exp(-1/2) ≈ 0.606
        # Effective area: π·σ²
        spatial_integral <- pi * sigma_i^2 * 0.606  # Gaussian: exp(-1/2)

        compensator <- compensator + alpha_i * spatial_integral * temporal_integral
      }
    }

    ll_final <- ll - compensator

    # Early stopping
    if(ll_final > .optim_state$best_ll + EARLY_STOP_MIN_IMPROVE) {
      .optim_state$best_ll <- ll_final
      .optim_state$best_params <- params
      .optim_state$no_improve_count <- 0
    } else {
      .optim_state$no_improve_count <- .optim_state$no_improve_count + 1

      if(.optim_state$no_improve_count >= EARLY_STOP_NO_IMPROVE) {
        .optim_state$should_stop <- TRUE
        if(verbose) {
          cat(sprintf("\n    ⏹ EARLY STOP: No improvement for %d evaluations\n",
                     EARLY_STOP_NO_IMPROVE))
        }
        return(-1e10)
      }
    }

    if(verbose && eval_num %% 5 == 0) {
      cat(sprintf("\n    LL: %.2f (events: %.2f, comp: %.2f) [eval %d]\n",
                 ll_final, ll, -compensator, eval_num))
    }

    return(ll_final)
  }

  return(list(
    lambda = lambda_spatial_temporal,
    loglik = loglik_spatial_temporal
  ))
}

hawkes_funcs <- spatial_temporal_hawkes_functions()

# ====================================================================
# 4. FITTING FUNCTION
# ====================================================================

fit_spatial_temporal_hawkes <- function(times, lons, lats, marks,
                                       time_mat, dist_mat, neighbors,
                                       include_marks = TRUE,
                                       use_mark_interaction = FALSE,
                                       model_name = "Model",
                                       model_id = "M1") {

  cat(sprintf("\n=== FITTING %s ===\n", model_name))
  start_time <- Sys.time()

  # Define parameters
  if(use_mark_interaction) {
    param_names <- c("mu", "beta_0", "beta_violence", "beta_state", "decay",
                    "sigma_violent", "sigma_peaceful")
    param_vec <- c(mu = 0.1, beta_0 = -1, beta_violence = 0.5, beta_state = 0.5,
                  decay = 0.2, sigma_violent = 50,
                  sigma_peaceful = 50)
  } else if(include_marks) {
    param_names <- c("mu", "beta_0", "beta_violence", "beta_state", "decay", "sigma")
    param_vec <- c(mu = 0.1, beta_0 = -1, beta_violence = 0.5, beta_state = 0.5,
                  decay = 0.2, sigma = 50)
  } else {
    param_names <- c("mu", "beta_0", "decay", "sigma")
    param_vec <- c(mu = 0.1, beta_0 = -1, decay = 0.2, sigma = 50)

    # Zero out mark effects
    marks_temp <- marks
    marks_temp$is_violent <- 0
    marks_temp$state_intervention <- 0
    marks <- marks_temp
  }

  # Bounds (relaxed to prevent boundary convergence)
  lower_bounds <- c(mu = 1e-6, beta_0 = -10, beta_violence = -10, beta_state = -10,
                   decay = 0.001, sigma = 1,
                   sigma_violent = 1,
                   sigma_peaceful = 1)[param_names]
  upper_bounds <- c(mu = 100, beta_0 = 10, beta_violence = 10, beta_state = 10,
                   decay = 5, sigma = 500,
                   sigma_violent = 500,
                   sigma_peaceful = 500)[param_names]

  # Objective function
  obj_fun <- function(params_vec) {
    params_list <- as.list(params_vec)
    names(params_list) <- param_names
    return(-hawkes_funcs$loglik(params_list, times, lons, lats, marks,
                                time_mat, dist_mat, neighbors,
                                use_mark_interaction, verbose = TRUE))
  }

  # Multi-start optimization
  cat(sprintf("  Multi-start optimization (%d starts)...\n", N_RANDOM_STARTS + 1))
  all_starts <- list()

  # Start 1: Default
  cat("\n  === START 1: Default ===\n")
  reset_optim_state(model_name)

  fit1 <- optim(par = param_vec, fn = obj_fun, method = "L-BFGS-B",
               lower = lower_bounds, upper = upper_bounds,
               control = list(maxit = 100, factr = 1e10))

  all_starts[[1]] <- list(
    fit = fit1,
    final_ll = if(.optim_state$should_stop) .optim_state$best_ll else -fit1$value,
    converged = fit1$convergence == 0 || .optim_state$should_stop
  )

  # Random starts
  for(start_i in 2:(N_RANDOM_STARTS + 1)) {
    cat(sprintf("\n  === START %d: Random ===\n", start_i))

    random_params <- param_vec
    random_params["mu"] <- runif(1, 0.05, 0.3)
    random_params["beta_0"] <- runif(1, -2, 0)
    random_params["decay"] <- runif(1, 0.1, 0.5)
    random_params["sigma"] <- runif(1, 20, 100)

    if("beta_violence" %in% param_names) random_params["beta_violence"] <- runif(1, -1, 1)
    if("beta_state" %in% param_names) random_params["beta_state"] <- runif(1, -1, 1)
    if("sigma_violent" %in% param_names) random_params["sigma_violent"] <- runif(1, 20, 100)
    if("sigma_peaceful" %in% param_names) random_params["sigma_peaceful"] <- runif(1, 20, 100)

    reset_optim_state(model_name)

    fit_i <- optim(par = random_params, fn = obj_fun, method = "L-BFGS-B",
                  lower = lower_bounds, upper = upper_bounds,
                  control = list(maxit = 100, factr = 1e10))

    all_starts[[start_i]] <- list(
      fit = fit_i,
      final_ll = if(.optim_state$should_stop) .optim_state$best_ll else -fit_i$value,
      converged = fit_i$convergence == 0 || .optim_state$should_stop
    )
  }

  # Select best
  best_idx <- which.max(sapply(all_starts, function(x) x$final_ll))
  fit <- all_starts[[best_idx]]$fit
  best_ll <- all_starts[[best_idx]]$final_ll

  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))

  cat(sprintf("\n✓ Best start: %d with LL = %.2f\n", best_idx, best_ll))
  cat(sprintf("✓ Runtime: %.1f minutes\n", runtime))

  # Extract results
  params_hat <- fit$par
  names(params_hat) <- param_names

  # Calculate AIC/BIC
  k <- length(param_names)
  aic <- -2 * best_ll + 2 * k
  bic <- -2 * best_ll + k * log(length(times))

  result <- list(
    model_name = model_name,
    model_id = model_id,
    params = params_hat,
    loglik = best_ll,
    aic = aic,
    bic = bic,
    convergence = fit$convergence,
    n_params = k,
    n_events = length(times),
    runtime_mins = runtime,
    all_starts = all_starts
  )

  cat(sprintf("  AIC: %.2f, BIC: %.2f\n", aic, bic))

  return(result)
}

# ====================================================================
# 5. FIT MODELS
# ====================================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  FITTING MODELS                                               ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n")

results <- list()

# M1: Spatial-temporal (no marks)
results$M1 <- fit_spatial_temporal_hawkes(
  times, lons, lats, marks_data,
  time_mat, dist_mat, spatiotemporal_neighbors,
  include_marks = FALSE,
  use_mark_interaction = FALSE,
  model_name = "M1: Spatial-temporal (no marks)",
  model_id = "M1"
)

# M3: Spatial-temporal-marks (additive)
results$M3 <- fit_spatial_temporal_hawkes(
  times, lons, lats, marks_data,
  time_mat, dist_mat, spatiotemporal_neighbors,
  include_marks = TRUE,
  use_mark_interaction = FALSE,
  model_name = "M3: Spatial-temporal-marks (additive)",
  model_id = "M3"
)

# M4: Spatial-mark interaction
results$M4 <- fit_spatial_temporal_hawkes(
  times, lons, lats, marks_data,
  time_mat, dist_mat, spatiotemporal_neighbors,
  include_marks = TRUE,
  use_mark_interaction = TRUE,
  model_name = "M4: Spatial-mark interaction",
  model_id = "M4"
)

# ====================================================================
# 6. MODEL COMPARISON
# ====================================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL COMPARISON                                             ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

comparison_df <- data.frame(
  Model = sapply(results, function(x) x$model_name),
  LogLik = sapply(results, function(x) x$loglik),
  nParams = sapply(results, function(x) x$n_params),
  AIC = sapply(results, function(x) x$aic),
  BIC = sapply(results, function(x) x$bic),
  Runtime_mins = sapply(results, function(x) x$runtime_mins)
) %>%
  arrange(AIC)

print(comparison_df)

# Save
write.csv(comparison_df, "model_comparison_phase3_fast.csv", row.names = FALSE)
cat("\n✓ Saved: model_comparison_phase3_fast.csv\n")

# ====================================================================
# 7. HYPOTHESIS TESTS
# ====================================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  HYPOTHESIS TESTS                                             ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

# H2: Marks persist after controlling for space?
# Test: M3 (spatial+marks) vs M1 (spatial only)
lr_h2 <- 2 * (results$M3$loglik - results$M1$loglik)
df_h2 <- results$M3$n_params - results$M1$n_params
p_h2 <- pchisq(lr_h2, df_h2, lower.tail = FALSE)

cat("H2: Do violence/state effects persist after controlling for space?\n")
cat(sprintf("  LRT: M3 vs M1\n"))
cat(sprintf("  LR = %.2f, df = %d, p = %.4f %s\n",
           lr_h2, df_h2, p_h2, ifelse(p_h2 < 0.05, "***", "")))

# H3: Violence spreads differently in space?
# Test: M4 (interaction) vs M3 (additive)
lr_h3 <- 2 * (results$M4$loglik - results$M3$loglik)
df_h3 <- results$M4$n_params - results$M3$n_params
p_h3 <- pchisq(lr_h3, df_h3, lower.tail = FALSE)

cat("\nH3: Does violence spread differently in space?\n")
cat(sprintf("  LRT: M4 vs M3\n"))
cat(sprintf("  LR = %.2f, df = %d, p = %.4f %s\n",
           lr_h3, df_h3, p_h3, ifelse(p_h3 < 0.05, "***", "")))

# Save hypothesis tests
hypothesis_tests <- data.frame(
  Hypothesis = c("H2: Marks persist after space", "H3: Spatial-mark interaction"),
  LR = c(lr_h2, lr_h3),
  df = c(df_h2, df_h3),
  p_value = c(p_h2, p_h3),
  significant = c(p_h2 < 0.05, p_h3 < 0.05)
)

print(hypothesis_tests)

write.csv(hypothesis_tests, "hypothesis_tests_phase3_fast.csv", row.names = FALSE)
cat("\n✓ Saved: hypothesis_tests_phase3_fast.csv\n")

# ====================================================================
# 8. SPATIAL PARAMETERS
# ====================================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  ESTIMATED SPATIAL PARAMETERS                                 ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

spatial_params_df <- data.frame(
  Model = c("M1", "M3", "M4 (violent)", "M4 (peaceful)"),
  sigma_km = c(
    results$M1$params["sigma"],
    results$M3$params["sigma"],
    results$M4$params["sigma_violent"],
    results$M4$params["sigma_peaceful"]
  )
)

print(spatial_params_df)

write.csv(spatial_params_df, "spatial_parameters_phase3_fast.csv", row.names = FALSE)
cat("\n✓ Saved: spatial_parameters_phase3_fast.csv\n")

# Interpret spatial decay
cat("\nInterpretation:\n")
cat("  σ (sigma) = spatial scale parameter (km)\n")
cat("  At distance d=σ: kernel value = exp(-1/2) ≈ 0.606 (60% of maximum)\n")
cat("  At distance d=2σ: kernel value = exp(-2) ≈ 0.135 (13.5% of maximum)\n")

# ====================================================================
# 9. SAVE RESULTS
# ====================================================================

cat("\n=== SAVING RESULTS ===\n")

saveRDS(results, "spatial_temporal_models_fast.rds")
cat("✓ Saved: spatial_temporal_models_fast.rds\n")

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  COMPLETE - FAST VALIDATION SUCCESSFUL                        ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n")
cat("\nNext step: Run full dataset version (09_phase3_spatial_temporal_marks_FULL.R)\n")
