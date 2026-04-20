# ============================================================================
# MODEL A: POPULATION OFFSET (Constant Per-Capita Rates)
# ============================================================================
#
# Specification: μ(district, year) = population × exp(β₀ + Σ β_year·I(year))
#
# This model assumes:
# - Protest rates are PROPORTIONAL to population (constant per-capita rate)
# - Population enters as an OFFSET (exposure term), not a covariate
# - Year effects capture temporal variation in per-capita rates
#
# Parameters (10 total):
# - Background (10): β₀, β_2016...β_2024 (no γ parameter)
# - Triggering (3): β₀_trig, β_violence, β_state
# - Decay (1): β
#
# ============================================================================

library(dplyr)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   MODEL A: POPULATION OFFSET (Constant Per-Capita Rates)    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("SPECIFICATION:\n")
cat("  Background: μ(d,y) = population × exp(β₀ + Σ β_year)\n")
cat("  Triggering: α(marks) = exp(β₀_trig + β_violence·violent + β_state·intervention)\n")
cat("  Kernel: g(t) = exp(-β·t)\n")
cat("  Intensity: λ(t) = μ(d,y) + Σ α(marks_i)·g(t-t_i)\n\n")

cat("KEY ASSUMPTION:\n")
cat("  → Districts with 2× population have 2× the protest rate\n")
cat("  → Per-capita rates are constant across population sizes\n")
cat("  → Population is an EXPOSURE term, not a predictor\n\n")

cat("PARAMETERS (14 total):\n")
cat("  Background (10): β₀, β_2016...β_2024 (NO γ parameter)\n")
cat("  Triggering (3): β₀_trig, β_violence, β_state\n")
cat("  Decay (1): β\n\n")

# ====================================================================
# CONFIGURATION
# ====================================================================

TEMPORAL_CUTOFF <- 730  # 2 years
CHECKPOINT_EVERY_N_EVALS <- 10

# ====================================================================
# CHECKPOINT FUNCTIONS
# ====================================================================

checkpoint_dir <- "checkpoints_model_a_offset"
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
# 1. LOAD DATA
# ====================================================================

cat("\n=== LOADING DATA ===\n")
protests <- readRDS("protests_with_population.rds")

cat("Loaded", nrow(protests), "events with population data\n")
cat("Time range:", min(protests$year), "to", max(protests$year), "\n\n")

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

cat(sprintf("Using full dataset: %d events\n", nrow(marks_data)))
cat(sprintf("Violent events: %d (%.1f%%)\n",
            sum(marks_data$is_violent),
            100*mean(marks_data$is_violent)))
cat(sprintf("Events with state intervention: %d (%.1f%%)\n\n",
            sum(marks_data$state_intervention),
            100*mean(marks_data$state_intervention)))

# ====================================================================
# 2. CALCULATE DISTRICT-YEAR OBSERVATION TIMES
# ====================================================================

cat("=== CALCULATING DISTRICT-YEAR OBSERVATION TIMES ===\n")

T_min <- min(marks_data$time)
T_max <- max(marks_data$time)

# Get unique district-year combinations
district_years <- marks_data %>%
  select(district, year, population) %>%
  distinct()

cat(sprintf("Found %d unique district-year combinations\n", nrow(district_years)))

# Calculate days observed for each district-year
district_years$days_observed <- T_max - T_min

cat(sprintf("Observation period: %.1f days (%.2f years)\n\n",
            T_max - T_min, (T_max - T_min)/365.25))

# ====================================================================
# 3. MODEL A IMPLEMENTATION
# ====================================================================

cat("=== IMPLEMENTING MODEL A (OFFSET) ===\n\n")

model_a_functions <- function() {

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
  # MODEL A: μ(d,y) = population × exp(β₀ + Σ β_year)
  calc_mu <- function(population, year, params) {

    beta_0_bg <- params[["beta_0_bg"]]

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

    # OFFSET MODEL: population enters multiplicatively
    mu <- population * exp(beta_0_bg + year_effect)
    return(mu)
  }

  # Conditional intensity
  lambda_offset <- function(idx, times, marks, params, time_mat = NULL) {

    # District-year-specific background rate (OFFSET MODEL)
    mu_i <- calc_mu(marks$population[idx], marks$year[idx], params)
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
          # Vectorized alpha calculation
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

  # Log-likelihood
  loglik_offset <- function(params, times, marks, district_years,
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
      lambda_i <- lambda_offset(i, times, marks, params, time_mat)
      ll <- ll + log(max(lambda_i, 1e-10))

      if(verbose && i %% 200 == 0) {
        cat(sprintf("    Event %d/%d (%.1f%%)  \r", i, n, 100*i/n))
      }
    }

    # Compensator: sum over district-years
    compensator <- 0

    # Background component: Σ_{district,year} μ(district,year) × days_observed
    # MODEL A: μ = population × exp(β₀ + year_effect)
    for(k in 1:nrow(district_years)) {
      mu_k <- calc_mu(district_years$population[k], district_years$year[k], params)
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
    lambda = lambda_offset,
    loglik = loglik_offset
  ))
}

hawkes_funcs <- model_a_functions()

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

fit_model_a <- function(times, marks, district_years, time_mat = NULL,
                        model_name = "Model A: Population Offset",
                        model_id = "MODEL_A") {

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

    # Initial parameters (NO GAMMA - only 10 background params)
    params_init_list <- list(
      # Background parameters (10) - NO gamma
      beta_0_bg = log(0.01 / mean(marks$population)),  # Per-capita rate
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

  # Parameter vector for optimization (14 parameters - NO gamma)
  param_names <- c("beta_0_bg",
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
    lower = c(beta_0_bg = -15,
              beta_2016 = -5, beta_2017 = -5, beta_2018 = -5, beta_2019 = -5, beta_2020 = -5,
              beta_2021 = -5, beta_2022 = -5, beta_2023 = -5, beta_2024 = -5,
              beta_0_trig = -5, beta_violence = -3, beta_state = -3, decay = 0.01),
    upper = c(beta_0_bg = -5,
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
cat("║   FITTING MODEL A (POPULATION OFFSET)   ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

model_a_fit <- fit_model_a(
  times_sample,
  marks_sample,
  district_years_sample,
  time_mat_sample,
  model_name = "Model A: Population Offset",
  model_id = "MODEL_A"
)

# ====================================================================
# 7. RESULTS
# ====================================================================

cat("\n╔══════════════════════════════════════════╗\n")
cat("║         MODEL A RESULTS                 ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

cat("MODEL FIT:\n")
cat(sprintf("  Log-likelihood: %.2f\n", model_a_fit$loglik))
cat(sprintf("  AIC: %.2f\n", model_a_fit$AIC))
cat(sprintf("  BIC: %.2f\n", model_a_fit$BIC))
cat(sprintf("  Parameters: %d\n", model_a_fit$n_params))
cat(sprintf("  Runtime: %.2f minutes\n\n", model_a_fit$runtime_mins))

cat("BACKGROUND RATE PARAMETERS:\n")
cat("---------------------------\n")
cat(sprintf("  β₀ (baseline per-capita rate): %.4f\n", model_a_fit$params$beta_0_bg))
cat(sprintf("    → exp(β₀) = %.6f (events per person per day)\n",
            exp(model_a_fit$params$beta_0_bg)))
cat(sprintf("    → Per 100k population: %.2f events/day\n",
            exp(model_a_fit$params$beta_0_bg) * 100000))
cat("\n")

cat("YEAR EFFECTS (relative to 2015):\n")
for(yr in 2016:2024) {
  param_name <- paste0("beta_", yr)
  if(param_name %in% names(model_a_fit$params)) {
    beta_yr <- model_a_fit$params[[param_name]]
    cat(sprintf("  %d: β = %6.3f  (exp(β) = %.3f)\n", yr, beta_yr, exp(beta_yr)))
  }
}
cat("\n")

cat("TRIGGERING PARAMETERS:\n")
cat("----------------------\n")
cat(sprintf("  β₀_trig: %.4f\n", model_a_fit$params$beta_0_trig))
cat(sprintf("  β_violence: %.4f (exp = %.3f)\n",
            model_a_fit$params$beta_violence, exp(model_a_fit$params$beta_violence)))
cat(sprintf("  β_state: %.4f (exp = %.3f)\n",
            model_a_fit$params$beta_state, exp(model_a_fit$params$beta_state)))
cat(sprintf("  Decay: %.4f (half-life = %.1f days)\n",
            model_a_fit$params$decay, log(2)/model_a_fit$params$decay))
cat("\n")

# Save results
saveRDS(model_a_fit, "model_a_offset.rds")
cat("✓ Saved: model_a_offset.rds\n\n")

cat("╔══════════════════════════════════════════╗\n")
cat("║   MODEL A ESTIMATION COMPLETE           ║\n")
cat("╚══════════════════════════════════════════╝\n")
