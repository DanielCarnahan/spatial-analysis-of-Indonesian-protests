# Phase 2: Heterogeneous Background Rates - FULL DATASET
# Author: Daniel Carnahan
# Date: 2025-11-04
#
# This script implements population-varying background rates:
# μ(district, year) = exp(β₀ + γ·log(population) + Σ β_year·I(year))
#
# Key changes from Phase 2 mark-dependent models:
# 1. Load protests_with_population.rds (includes log_pop variable)
# 2. Implement heterogeneous background rates (11 parameters)
# 3. Change compensator calculation to sum over district-years
# 4. Keep mark-dependent triggering (3 parameters) from Phase 2
#
# Total parameters: 14
# - Background: β₀_bg, γ, β_2016...β_2024 (11 params, 2015 is reference)
# - Triggering: β₀_trig, β_violence, β_state (3 params)
# - Decay: β (1 param, shared)

library(dplyr)
library(ggplot2)

cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║   PHASE 2: HETEROGENEOUS BACKGROUND RATES (FULL DATASET)    ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

cat("MODEL SPECIFICATION:\n")
cat("  Background: μ(district,year) = exp(β₀ + γ·log_pop + Σ β_year·I(year))\n")
cat("  Triggering: α(marks) = exp(β₀_trig + β_violence·violent + β_state·intervention)\n")
cat("  Kernel: g(t) = exp(-β·t)\n")
cat("  Intensity: λ(t) = μ(district,year) + Σ α(marks_i)·g(t-t_i)\n\n")

cat("PARAMETERS (14 total):\n")
cat("  Background (11): β₀, γ, β_2016, β_2017, ..., β_2024\n")
cat("  Triggering (3): β₀_trig, β_violence, β_state\n")
cat("  Decay (1): β (shared)\n\n")

cat("FULL DATASET:\n")
cat("  - Using all 15,914 events (100% of data)\n")
cat("  - Max 500 iterations\n")
cat("  - Expected time: 12-15 minutes\n")
cat("  - Full estimation with robust parameter estimates\n\n")

# ====================================================================
# CONFIGURATION
# ====================================================================

TEMPORAL_CUTOFF <- 730  # 2 years
CHECKPOINT_EVERY_N_EVALS <- 10

# ====================================================================
# CHECKPOINT FUNCTIONS
# ====================================================================

checkpoint_dir <- "checkpoints_phase2_heterogeneous_full"
if(!dir.exists(checkpoint_dir)) {
  dir.create(checkpoint_dir)
  cat(sprintf("✓ Created checkpoint directory: %s\n", checkpoint_dir))
}

.optim_state <- new.env()
.optim_state$n_evals <- 0
.optim_state$best_ll <- -Inf
.optim_state$best_params <- NULL
.optim_state$model_name <- NULL

reset_optim_state <- function(model_name) {
  .optim_state$n_evals <- 0
  .optim_state$best_ll <- -Inf
  .optim_state$best_params <- NULL
  .optim_state$model_name <- model_name
}

save_intermediate_checkpoint <- function(model_id, eval_num, params, ll) {
  checkpoint_file <- file.path(checkpoint_dir,
                                sprintf("intermediate-%s-eval%04d.rds", model_id, eval_num))

  saveRDS(list(
    model_id = model_id,
    eval_num = eval_num,
    params = params,
    ll = ll,
    best_ll = .optim_state$best_ll,
    best_params = .optim_state$best_params,
    timestamp = Sys.time(),
    status = "in_progress"
  ), checkpoint_file)

  invisible(checkpoint_file)
}

save_final_checkpoint <- function(model_fit, model_id) {
  checkpoint_file <- file.path(checkpoint_dir, sprintf("model-%s-FINAL.rds", model_id))

  model_fit$status <- "completed"
  model_fit$timestamp <- Sys.time()

  saveRDS(model_fit, checkpoint_file)
  cat(sprintf("✓ Final checkpoint saved: %s\n\n", checkpoint_file))

  # Clean up intermediate checkpoints
  intermediate_files <- list.files(checkpoint_dir,
                                    pattern = sprintf("intermediate-%s-", model_id),
                                    full.names = TRUE)
  if(length(intermediate_files) > 0) {
    file.remove(intermediate_files)
    cat(sprintf("  Cleaned up %d intermediate checkpoints\n\n", length(intermediate_files)))
  }
}

find_latest_checkpoint <- function(model_id) {
  final_file <- file.path(checkpoint_dir, sprintf("model-%s-FINAL.rds", model_id))
  if(file.exists(final_file)) {
    return(list(type = "final", file = final_file, data = readRDS(final_file)))
  }

  pattern <- sprintf("intermediate-%s-eval\\d{4}\\.rds", model_id)
  checkpoints <- list.files(checkpoint_dir, pattern = pattern, full.names = TRUE)

  if(length(checkpoints) > 0) {
    latest <- checkpoints[order(file.mtime(checkpoints), decreasing = TRUE)[1]]
    return(list(type = "intermediate", file = latest, data = readRDS(latest)))
  }

  return(NULL)
}

# ====================================================================
# 1. LOAD DATA WITH POPULATION
# ====================================================================

cat("\n=== LOADING DATA ===\n")
protests <- readRDS("protests_with_population.rds")

cat("Loaded", nrow(protests), "events with population data\n")
cat("Time range:", min(protests$year), "to", max(protests$year), "\n\n")

# Check data
cat("Data summary:\n")
cat(sprintf("  Unique districts: %d\n", length(unique(protests$district))))
cat(sprintf("  Year range: %d-%d\n", min(protests$year), max(protests$year)))
cat(sprintf("  Population range: %.0f - %.0f\n",
            min(protests$population), max(protests$population)))
cat(sprintf("  Log population range: %.2f - %.2f\n\n",
            min(protests$log_pop), max(protests$log_pop)))

# Prepare marks data
marks_data <- protests %>%
  select(
    event_id = event_id_cnty,
    time = days_since_start,
    is_violent,
    is_peaceful,
    state_intervention,
    longitude,
    latitude,
    district,
    year,
    log_pop,
    population
  ) %>%
  arrange(time)

# Use FULL dataset (no sampling)
cat(sprintf("Using full dataset: %d events (100.0%%)\n", nrow(marks_data)))
cat("NOTE: This is the complete dataset for robust parameter estimation\n\n")

# Summary of marks
cat("Mark Statistics:\n")
cat(sprintf("  Violent events: %d (%.1f%%)\n",
            sum(marks_data$is_violent),
            100*mean(marks_data$is_violent)))
cat(sprintf("  Events with state intervention: %d (%.1f%%)\n",
            sum(marks_data$state_intervention),
            100*mean(marks_data$state_intervention)))
cat(sprintf("  Mean log_pop: %.2f (SD: %.2f)\n",
            mean(marks_data$log_pop), sd(marks_data$log_pop)))
cat("\n")

# ====================================================================
# 2. CALCULATE DISTRICT-YEAR OBSERVATION TIMES
# ====================================================================

cat("=== CALCULATING DISTRICT-YEAR OBSERVATION TIMES ===\n")

# Get observation window for each district-year
# (days observed in the sample period)
T_min <- min(marks_data$time)
T_max <- max(marks_data$time)

# For this sample, we'll use the simplified assumption:
# Each district-year is observed for the full duration it overlaps with [T_min, T_max]

# Get unique district-year combinations in the sample
district_years <- marks_data %>%
  select(district, year, log_pop) %>%
  distinct()

cat(sprintf("Found %d unique district-year combinations in sample\n", nrow(district_years)))

# Calculate days observed for each district-year
# This is a simplified version - full version would account for actual observation periods
district_years$days_observed <- T_max - T_min

cat(sprintf("Observation period: %.1f days (%.2f years)\n\n",
            T_max - T_min, (T_max - T_min)/365.25))

# ====================================================================
# 3. HETEROGENEOUS HAWKES IMPLEMENTATION
# ====================================================================

cat("=== IMPLEMENTING HETEROGENEOUS HAWKES ===\n\n")

heterogeneous_hawkes_functions <- function() {

  # Pre-compute time differences with cutoff
  precompute_time_matrix <- function(times, cutoff = 730) {
    n <- length(times)
    time_mat <- matrix(0, nrow = n, ncol = n)

    for(i in 1:n) {
      if(i > 1) {
        for(j in 1:(i-1)) {
          tau <- times[i] - times[j]
          if(isTRUE(tau > 0) && isTRUE(tau <= cutoff)) {
            time_mat[i, j] <- tau
          }
        }
      }
    }

    return(time_mat)
  }

  # Calculate district-year-specific background rate
  calc_mu <- function(log_pop, year, params) {
    # μ(district,year) = exp(β₀ + γ·log_pop + Σ β_year·I(year))

    beta_0_bg <- params[["beta_0_bg"]]
    gamma <- params[["gamma"]]

    # Year effects (2015 is reference, omitted)
    year_effect <- 0
    if(year == 2016 && "beta_2016" %in% names(params)) year_effect <- params[["beta_2016"]]
    if(year == 2017 && "beta_2017" %in% names(params)) year_effect <- params[["beta_2017"]]
    if(year == 2018 && "beta_2018" %in% names(params)) year_effect <- params[["beta_2018"]]
    if(year == 2019 && "beta_2019" %in% names(params)) year_effect <- params[["beta_2019"]]
    if(year == 2020 && "beta_2020" %in% names(params)) year_effect <- params[["beta_2020"]]
    if(year == 2021 && "beta_2021" %in% names(params)) year_effect <- params[["beta_2021"]]
    if(year == 2022 && "beta_2022" %in% names(params)) year_effect <- params[["beta_2022"]]
    if(year == 2023 && "beta_2023" %in% names(params)) year_effect <- params[["beta_2023"]]
    if(year == 2024 && "beta_2024" %in% names(params)) year_effect <- params[["beta_2024"]]

    mu <- exp(beta_0_bg + gamma * log_pop + year_effect)
    return(mu)
  }

  # Conditional intensity with heterogeneous background
  lambda_heterogeneous <- function(idx, times, marks, params, time_mat = NULL) {

    # District-year-specific background rate
    mu_i <- calc_mu(marks$log_pop[idx], marks$year[idx], params)
    intensity <- mu_i

    # Mark-dependent triggering parameters
    beta_0_trig <- params[["beta_0_trig"]]
    beta_violence <- params[["beta_violence"]]
    beta_state <- params[["beta_state"]]
    decay <- params[["decay"]]

    # Triggered component with cutoff
    if(idx > 1) {

      if(!is.null(time_mat)) {
        valid_j <- which(time_mat[idx, ] > 0)

        if(length(valid_j) > 0) {
          # Vectorized alpha calculation (triggering strength)
          alphas <- exp(beta_0_trig +
                       beta_violence * marks$is_violent[valid_j] +
                       beta_state * marks$state_intervention[valid_j])

          # Vectorized kernel evaluation
          taus <- time_mat[idx, valid_j]
          intensity <- intensity + sum(alphas * exp(-decay * taus))
        }

      } else {
        # Fallback: loop with cutoff
        current_time <- times[idx]
        cutoff_time <- current_time - TEMPORAL_CUTOFF

        for(j in 1:(idx-1)) {
          if(times[j] < cutoff_time) next

          tau <- current_time - times[j]
          if(tau > 0) {
            alpha_j <- exp(beta_0_trig +
                          beta_violence * marks$is_violent[j] +
                          beta_state * marks$state_intervention[j])

            intensity <- intensity + alpha_j * exp(-decay * tau)
          }
        }
      }
    }

    return(intensity)
  }

  # Log-likelihood with heterogeneous background
  loglik_heterogeneous <- function(params, times, marks, district_years,
                                   time_mat = NULL, verbose = TRUE, model_id = NULL) {

    # Parameter constraints
    if(params[["decay"]] <= 0) return(-Inf)

    # Track function evaluations
    .optim_state$n_evals <- .optim_state$n_evals + 1
    eval_num <- .optim_state$n_evals

    n <- length(times)
    ll <- 0

    # Sum log-intensities at event times
    for(i in 1:n) {
      lambda_i <- lambda_heterogeneous(i, times, marks, params, time_mat)
      ll <- ll + log(max(lambda_i, 1e-10))

      if(verbose && i %% 200 == 0) {
        cat(sprintf("    Event %d/%d (%.1f%%)  \r", i, n, 100*i/n))
      }
    }

    # Compensator: sum over district-years
    compensator <- 0

    # Background component: Σ_{district,year} μ(district,year) × days_observed
    for(k in 1:nrow(district_years)) {
      mu_k <- calc_mu(district_years$log_pop[k], district_years$year[k], params)
      compensator <- compensator + mu_k * district_years$days_observed[k]
    }

    # Triggered component contribution
    for(i in 1:n) {
      if(times[i] < max(times)) {

        # Mark-dependent alpha for event i
        alpha_i <- exp(params[["beta_0_trig"]] +
                      params[["beta_violence"]] * marks$is_violent[i] +
                      params[["beta_state"]] * marks$state_intervention[i])

        # Integral of exponential kernel (with cutoff consideration)
        tau_remain <- min(max(times) - times[i], TEMPORAL_CUTOFF)
        compensator <- compensator +
          (alpha_i / params[["decay"]]) * (1 - exp(-params[["decay"]] * tau_remain))
      }
    }

    ll_final <- ll - compensator

    # Track best result
    if(ll_final > .optim_state$best_ll) {
      .optim_state$best_ll <- ll_final
      .optim_state$best_params <- params
    }

    # Intermediate checkpoint
    if(!is.null(model_id) && eval_num %% CHECKPOINT_EVERY_N_EVALS == 0) {
      save_intermediate_checkpoint(model_id, eval_num, params, ll_final)
      if(verbose) {
        cat(sprintf("\n    💾 Checkpoint %d: LL=%.2f (best: %.2f)\n",
                   eval_num, ll_final, .optim_state$best_ll))
      }
    }

    if(verbose) {
      cat(sprintf("\n    LL: %.2f (events: %.2f, comp: %.2f) [eval %d]\n",
                 ll_final, ll, -compensator, eval_num))
    }

    return(ll_final)
  }

  return(list(
    precompute_time_matrix = precompute_time_matrix,
    calc_mu = calc_mu,
    lambda = lambda_heterogeneous,
    loglik = loglik_heterogeneous
  ))
}

hawkes_funcs <- heterogeneous_hawkes_functions()

# ====================================================================
# 4. PREPARE DATA
# ====================================================================

cat("=== PREPARING DATA FOR ESTIMATION ===\n")

times_sample <- marks_data$time
marks_sample <- marks_data
district_years_sample <- district_years

cat(sprintf("Sample size: %d events\n", length(times_sample)))
cat(sprintf("District-years: %d\n", nrow(district_years_sample)))

# Pre-compute time difference matrix
cat("Pre-computing time difference matrix...\n")
time_mat_sample <- hawkes_funcs$precompute_time_matrix(times_sample, TEMPORAL_CUTOFF)
cat("✓ Time matrix computed\n\n")

# ====================================================================
# 5. FITTING FUNCTION
# ====================================================================

fit_heterogeneous_hawkes <- function(times, marks, district_years, time_mat = NULL,
                                     model_name = "Heterogeneous Background Model",
                                     model_id = "HET") {

  cat(sprintf("\n--- Fitting %s ---\n", model_name))

  # Check for existing checkpoint
  checkpoint <- find_latest_checkpoint(model_id)

  if(!is.null(checkpoint)) {
    if(checkpoint$type == "final") {
      cat("✓ Model already completed (loading from final checkpoint)\n")
      return(checkpoint$data)
    } else {
      cat(sprintf("⏩ Resuming from intermediate checkpoint (eval %d, LL=%.2f)\n",
                 checkpoint$data$eval_num, checkpoint$data$best_ll))

      reset_optim_state(model_name)
      .optim_state$n_evals <- checkpoint$data$eval_num
      .optim_state$best_ll <- checkpoint$data$best_ll
      .optim_state$best_params <- checkpoint$data$best_params

      params_init_list <- checkpoint$data$best_params
    }
  } else {
    cat("Starting fresh optimization...\n")
    reset_optim_state(model_name)

    # Initial parameters
    params_init_list <- list(
      # Background parameters (11)
      beta_0_bg = log(0.01),      # Base background rate
      gamma = 0.5,                 # Population effect (positive = more protests in populous districts)
      beta_2016 = 0,
      beta_2017 = 0,
      beta_2018 = 0,
      beta_2019 = 0,
      beta_2020 = 0,
      beta_2021 = 0,
      beta_2022 = 0,
      beta_2023 = 0,
      beta_2024 = 0,
      # Triggering parameters (3)
      beta_0_trig = log(0.1),
      beta_violence = 0,
      beta_state = 0,
      # Decay (1)
      decay = 0.2
    )
  }

  start_time <- Sys.time()

  # Parameter vector for optimization (all 14 parameters)
  param_names <- c("beta_0_bg", "gamma",
                   "beta_2016", "beta_2017", "beta_2018", "beta_2019", "beta_2020",
                   "beta_2021", "beta_2022", "beta_2023", "beta_2024",
                   "beta_0_trig", "beta_violence", "beta_state", "decay")

  param_vec <- unlist(params_init_list[param_names])
  names(param_vec) <- param_names

  # Objective function
  obj_fun <- function(par) {

    params_full <- as.list(par)

    ll <- hawkes_funcs$loglik(params_full, times, marks, district_years, time_mat,
                               verbose = TRUE, model_id = model_id)
    return(-ll)
  }

  # Optimize
  cat("  Starting optimization (max 500 iterations)...\n")
  cat(sprintf("  Expected time: 12-15 minutes\n"))

  fit <- optim(
    par = param_vec,
    fn = obj_fun,
    method = "L-BFGS-B",
    lower = c(beta_0_bg = -10, gamma = -2,
              beta_2016 = -5, beta_2017 = -5, beta_2018 = -5, beta_2019 = -5, beta_2020 = -5,
              beta_2021 = -5, beta_2022 = -5, beta_2023 = -5, beta_2024 = -5,
              beta_0_trig = -5, beta_violence = -3, beta_state = -3, decay = 0.01),
    upper = c(beta_0_bg = 5, gamma = 3,
              beta_2016 = 5, beta_2017 = 5, beta_2018 = 5, beta_2019 = 5, beta_2020 = 5,
              beta_2021 = 5, beta_2022 = 5, beta_2023 = 5, beta_2024 = 5,
              beta_0_trig = 5, beta_violence = 3, beta_state = 3, decay = 2),
    control = list(maxit = 500, factr = 1e10, trace = 1, REPORT = 10)
  )

  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")

  cat(sprintf("\n  ✓ Optimization complete in %.2f minutes\n", as.numeric(runtime)))
  cat(sprintf("  Log-likelihood: %.2f\n", -fit$value))
  cat(sprintf("  Best LL found: %.2f\n", .optim_state$best_ll))
  cat(sprintf("  Convergence: %d\n", fit$convergence))

  # Extract results
  final_params <- if(.optim_state$best_ll > -fit$value) {
    cat("  Using best parameters found during optimization\n")
    .optim_state$best_params
  } else {
    as.list(fit$par)
  }

  final_ll <- if(.optim_state$best_ll > -fit$value) {
    .optim_state$best_ll
  } else {
    -fit$value
  }

  results <- list(
    model_name = model_name,
    model_id = model_id,
    params = final_params,
    loglik = final_ll,
    convergence = fit$convergence,
    n_params = length(fit$par),
    n_events = length(times),
    AIC = 2*length(fit$par) - 2*final_ll,
    BIC = length(fit$par)*log(length(times)) - 2*final_ll,
    runtime_mins = as.numeric(runtime)
  )

  # Save final checkpoint
  save_final_checkpoint(results, model_id)

  return(results)
}

# ====================================================================
# 6. FIT MODEL
# ====================================================================

cat("\n╔══════════════════════════════════════════╗\n")
cat("║   FITTING HETEROGENEOUS MODEL           ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

het_fit <- fit_heterogeneous_hawkes(
  times_sample,
  marks_sample,
  district_years_sample,
  time_mat_sample,
  model_name = "Heterogeneous Background (Full)",
  model_id = "HET_FULL"
)

# ====================================================================
# 7. RESULTS AND INTERPRETATION
# ====================================================================

cat("\n╔══════════════════════════════════════════╗\n")
cat("║         RESULTS SUMMARY                 ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

cat("MODEL FIT:\n")
cat(sprintf("  Log-likelihood: %.2f\n", het_fit$loglik))
cat(sprintf("  AIC: %.2f\n", het_fit$AIC))
cat(sprintf("  BIC: %.2f\n", het_fit$BIC))
cat(sprintf("  Parameters: %d\n", het_fit$n_params))
cat(sprintf("  Runtime: %.2f minutes\n\n", het_fit$runtime_mins))

cat("BACKGROUND RATE PARAMETERS:\n")
cat("---------------------------\n")
cat(sprintf("  β₀ (baseline): %.4f\n", het_fit$params$beta_0_bg))
cat(sprintf("  γ (population): %.4f\n", het_fit$params$gamma))
cat(sprintf("    → exp(γ) = %.3f\n", exp(het_fit$params$gamma)))
cat(sprintf("    → A 1-unit increase in log(pop) multiplies μ by %.3f\n",
            exp(het_fit$params$gamma)))
cat("\n")

cat("YEAR EFFECTS (relative to 2015):\n")
for(yr in 2016:2024) {
  param_name <- paste0("beta_", yr)
  if(param_name %in% names(het_fit$params)) {
    beta_yr <- het_fit$params[[param_name]]
    cat(sprintf("  %d: β = %6.3f  (exp(β) = %.3f)\n", yr, beta_yr, exp(beta_yr)))
  }
}
cat("\n")

cat("TRIGGERING PARAMETERS:\n")
cat("----------------------\n")
cat(sprintf("  β₀_trig: %.4f\n", het_fit$params$beta_0_trig))
cat(sprintf("  β_violence: %.4f (exp = %.3f)\n",
            het_fit$params$beta_violence, exp(het_fit$params$beta_violence)))
cat(sprintf("  β_state: %.4f (exp = %.3f)\n",
            het_fit$params$beta_state, exp(het_fit$params$beta_state)))
cat(sprintf("  Decay: %.4f (half-life = %.1f days)\n",
            het_fit$params$decay, log(2)/het_fit$params$decay))
cat("\n")

# Save results
saveRDS(het_fit, "heterogeneous_hawkes_full.rds")
cat("✓ Saved: heterogeneous_hawkes_full.rds\n\n")

# ====================================================================
# 8. COMPARISON TO HOMOGENEOUS MODEL
# ====================================================================

cat("╔══════════════════════════════════════════╗\n")
cat("║   COMPARISON TO HOMOGENEOUS MODEL       ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

cat("Next steps:\n")
cat("  1. Fit homogeneous baseline on same sample for comparison\n")
cat("  2. Conduct likelihood ratio test\n")
cat("  3. If heterogeneous model improves fit, run on full dataset\n")
cat("  4. Generate predicted background rates for each district-year\n")
cat("  5. Visualize spatial-temporal variation in μ\n\n")

cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║   PHASE 2 HETEROGENEOUS BACKGROUND COMPLETE (FULL DATASET)   ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n")
