#!/usr/bin/env Rscript
# ==============================================================================
# MODEL 5: EXTENDED MARKS HAWKES PROCESS (EXPONENTIAL KERNEL)
# ==============================================================================
#
# Specification: λ(t) = μ(d,y) + Σ α(marks) · exp(-β·(t-t_i))
#
# MARKS (5 total):
#   Violence indicators:
#     - has_fatalities: Fatality indicator (NON-EXCLUSIVE - can overlap)
#     - is_violent: Violence indicator (violent vs peaceful)
#
#   Actor-based categories (non-exclusive):
#     - is_student: Student-led indicator
#     - is_labor: Labor-led indicator
#     - is_papua: Papua-related indicator
#
# COMPARISON:
#   - Model 5 vs Model 2: LR test for mark effects
#   - 5 mark coefficients (violence-focused)
#
# ==============================================================================

library(dplyr)
library(Rcpp)

# Compile C++ likelihood function
cat("Compiling C++ likelihood function for extended marks Hawkes model...\n")
sourceCpp("hawkes_extended_marks_likelihood.cpp")
cat("✓ C++ compiled\n\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   MODEL 5: EXTENDED MARKS HAWKES (EXPONENTIAL)               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("SPECIFICATION:\n")
cat("  Background: μ(d,y) = exp(β₀ + γ·log(pop) + δ·poverty + Σ β_year)\n")
cat("  Triggering: α(marks) = exp(β₀_trig + β_fatal + β_violent\n")
cat("                              + β_student + β_labor + β_papua)\n")
cat("  Kernel: g(t) = exp(-β·t)  [EXPONENTIAL]\n")
cat("  Intensity: λ(t) = μ(d,y) + Σ α(marks_i) · g(t - t_i)\n\n")

cat("MARKS (violence-focused):\n")
cat("  Violence indicators (2 parameters):\n")
cat("    - has_fatalities (non-exclusive)\n")
cat("    - is_violent (violent vs peaceful)\n")
cat("  Actor-based (3 parameters, non-exclusive):\n")
cat("    - is_student\n")
cat("    - is_labor\n")
cat("    - is_papua\n\n")

# ====================================================================
# CONFIGURATION
# ====================================================================

TEST_MODE <- FALSE
TEMPORAL_CUTOFF <- 90  # 90 days
CHECKPOINT_EVERY_N_EVALS <- 1

# ====================================================================
# PARAMETER TRANSFORMATION FUNCTIONS
# ====================================================================

# Background parameters
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

# Triggering parameters
transform_beta_0_trig <- function(beta_0_trig_raw) {
  -10 + 20 * plogis(beta_0_trig_raw)
}
inverse_transform_beta_0_trig <- function(beta_0_trig) {
  qlogis((beta_0_trig + 10) / 20)
}

transform_beta_mark <- function(beta_mark_raw) {
  -5 + 10 * plogis(beta_mark_raw)
}
inverse_transform_beta_mark <- function(beta_mark) {
  qlogis((beta_mark + 5) / 10)
}

transform_decay <- function(decay_raw) {
  0.001 + 9.999 * plogis(decay_raw)
}
inverse_transform_decay <- function(decay) {
  qlogis((decay - 0.001) / 9.999)
}

# ====================================================================
# CHECKPOINT SYSTEM
# ====================================================================

checkpoint_dir <- "checkpoints_model5_extended"
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
protests <- readRDS("protests_with_poverty.rds")

cat("Loaded", nrow(protests), "events with poverty data\n")
cat("Time range:", min(protests$year), "to", max(protests$year), "\n")

# TEST MODE: Sample for quick testing
if(TEST_MODE) {
  set.seed(42)
  n_sample <- min(1000, nrow(protests))
  sample_idx <- sample(nrow(protests), n_sample)
  protests <- protests[sample_idx, ]
  cat(sprintf("⚠ TEST MODE: Using random sample of %d events\n", n_sample))
}

# Check poverty coverage
n_with_poverty <- sum(!is.na(protests$poverty_decimal))
pct_coverage <- 100 * n_with_poverty / nrow(protests)
cat(sprintf("Poverty data coverage: %d (%.1f%%)\n\n", n_with_poverty, pct_coverage))

# ====================================================================
# 2. CONSTRUCT EXTENDED MARK VARIABLES
# ====================================================================

cat("=== CONSTRUCTING MARK VARIABLES ===\n")

marks_data <- protests %>%
  select(
    event_id = event_id_cnty,
    time = days_since_start,
    event_type,
    sub_event_type,
    assoc_actor_1,
    actor2,
    admin1,
    fatalities,
    longitude,
    latitude,
    district,
    year,
    log_pop,
    population,
    poverty_decimal
  ) %>%
  arrange(time) %>%
  mutate(
    # Violence indicators
    has_fatalities = as.integer(fatalities > 0),

    # Violence indicator (violent vs peaceful)
    is_violent = as.integer(
      sub_event_type %in% c("Mob violence", "Violent demonstration") |
      sub_event_type == "Excessive force against protesters"
    ),

    # Actor-based categories (non-exclusive)
    is_student = as.integer(grepl("Student", assoc_actor_1, ignore.case = TRUE)),
    is_labor = as.integer(grepl("Labor|Worker", assoc_actor_1, ignore.case = TRUE)),
    is_papua = as.integer(
      grepl("Papua", admin1, ignore.case = TRUE) |
      grepl("Papua", assoc_actor_1, ignore.case = TRUE)
    )
  )

# Filter to complete cases
marks_data <- marks_data %>%
  filter(!is.na(poverty_decimal))

cat(sprintf("Using complete cases: %d events\n\n", nrow(marks_data)))

# Summary of marks
cat("MARK VARIABLE SUMMARY (violence-focused):\n")
cat("-----------------------------------------\n")
cat(sprintf("  Fatal events: %d (%.1f%%)\n",
            sum(marks_data$has_fatalities), 100*mean(marks_data$has_fatalities)))
cat(sprintf("  Violent events: %d (%.1f%%)\n",
            sum(marks_data$is_violent), 100*mean(marks_data$is_violent)))
cat(sprintf("  Student-led: %d (%.1f%%)\n",
            sum(marks_data$is_student), 100*mean(marks_data$is_student)))
cat(sprintf("  Labor-led: %d (%.1f%%)\n",
            sum(marks_data$is_labor), 100*mean(marks_data$is_labor)))
cat(sprintf("  Papua-related: %d (%.1f%%)\n\n",
            sum(marks_data$is_papua), 100*mean(marks_data$is_papua)))

# ====================================================================
# 3. CALCULATE DISTRICT-YEAR OBSERVATION TIMES
# ====================================================================

cat("=== CALCULATING DISTRICT-YEAR OBSERVATION TIMES ===\n")

T_min <- min(marks_data$time)
T_max <- max(marks_data$time)

district_years <- marks_data %>%
  select(district, year, population, log_pop, poverty_decimal) %>%
  distinct()

cat(sprintf("Found %d unique district-year combinations\n", nrow(district_years)))

district_years$days_observed <- T_max - T_min

cat(sprintf("Observation period: %.1f days (%.2f years)\n\n",
            T_max - T_min, (T_max - T_min)/365.25))

# ====================================================================
# 4. LIKELIHOOD FUNCTION WRAPPER
# ====================================================================

loglik_extended <- function(params, times, marks, district_years,
                            verbose = TRUE, model_id = NULL) {

  .optim_state$n_evals <- .optim_state$n_evals + 1
  eval_num <- .optim_state$n_evals

  # Prepare parameter vector for C++ (19 parameters)
  param_vec <- c(
    params[["beta_0_bg"]],
    params[["gamma_raw"]],
    params[["delta_raw"]],
    # Year effects (2016-2024)
    params[["beta_2016_raw"]],
    params[["beta_2017_raw"]],
    params[["beta_2018_raw"]],
    params[["beta_2019_raw"]],
    params[["beta_2020_raw"]],
    params[["beta_2021_raw"]],
    params[["beta_2022_raw"]],
    params[["beta_2023_raw"]],
    params[["beta_2024_raw"]],
    # Triggering parameters (7): baseline + 5 marks + decay
    params[["beta_0_trig_raw"]],
    params[["beta_fatal_raw"]],
    params[["beta_violent_raw"]],
    params[["beta_student_raw"]],
    params[["beta_labor_raw"]],
    params[["beta_papua_raw"]],
    params[["decay_raw"]]
  )

  # Call C++ function
  neg_ll <- hawkes_extended_marks_negloglik_cpp(
    times = times,
    log_pop = marks$log_pop,
    poverty_decimal = marks$poverty_decimal,
    year = as.integer(marks$year),
    # 5 mark vectors (violence-focused)
    has_fatalities = as.integer(marks$has_fatalities),
    is_violent = as.integer(marks$is_violent),
    is_student = as.integer(marks$is_student),
    is_labor = as.integer(marks$is_labor),
    is_papua = as.integer(marks$is_papua),
    # District-year data
    district_year_exposure = district_years$days_observed,
    district_year_log_pop = district_years$log_pop,
    district_year_poverty = district_years$poverty_decimal,
    district_year_year = as.integer(district_years$year),
    params = param_vec,
    temporal_cutoff = TEMPORAL_CUTOFF
  )

  ll_final <- -neg_ll

  # Track best result
  if(ll_final > .optim_state$best_ll) {
    .optim_state$best_ll <- ll_final
    .optim_state$best_params <- params
  }

  # Intermediate checkpoint
  if(!is.null(model_id) && eval_num %% CHECKPOINT_EVERY_N_EVALS == 0) {
    gamma_display <- transform_gamma(params[["gamma_raw"]])
    save_intermediate_checkpoint(model_id, eval_num, params, ll_final)
    if(verbose) {
      cat(sprintf("\n    💾 Checkpoint %d: LL=%.2f (best: %.2f) | γ=%.4f\n",
                 eval_num, ll_final, .optim_state$best_ll, gamma_display))
    }
  }

  if(verbose) {
    gamma_display <- transform_gamma(params[["gamma_raw"]])
    cat(sprintf("    LL: %.2f | γ=%.4f [eval %d]\n",
               ll_final, gamma_display, eval_num))
    flush.console()
  }

  return(ll_final)
}

# ====================================================================
# 5. FITTING FUNCTION
# ====================================================================

fit_model_extended <- function(times, marks, district_years,
                               model_name = "Model 5: Extended Marks",
                               model_id = "MODEL5_EXTENDED",
                               use_multistart = TRUE,
                               n_starts = 5) {

  cat(sprintf("\n--- Fitting %s ---\n", model_name))
  if(use_multistart) {
    cat(sprintf("Using MULTI-START optimization with %d starting points\n", n_starts))
  }

  # Check for existing checkpoint
  checkpoint <- find_latest_checkpoint(model_id)

  if(!is.null(checkpoint)) {
    if(checkpoint$type == "final") {
      cat("✓ Model already completed (loading from final checkpoint)\n")
      return(checkpoint$data)
    } else {
      gamma_resume <- transform_gamma(checkpoint$data$best_params$gamma_raw)
      cat(sprintf("⏩ Resuming from intermediate checkpoint (eval %d, LL=%.2f, γ=%.4f)\n",
                 checkpoint$data$eval_num, checkpoint$data$best_ll, gamma_resume))

      reset_optim_state(model_name)
      .optim_state$n_evals <- checkpoint$data$eval_num
      .optim_state$best_ll <- checkpoint$data$best_ll
      .optim_state$best_params <- checkpoint$data$best_params

      use_multistart <- FALSE
      starting_points <- list(list(name = "Checkpoint", params = checkpoint$data$best_params))
    }
  } else {
    cat("Starting fresh optimization...\n")

    # Calculate baseline values
    mean_log_pop <- mean(marks$log_pop)
    mean_poverty <- mean(marks$poverty_decimal, na.rm = TRUE)
    total_events <- nrow(marks)
    total_exposure <- nrow(district_years) * (max(times) - min(times))
    observed_rate <- total_events / total_exposure

    # Define starting points
    starting_points <- list()

    # 1. Data-driven (primary)
    gamma1 <- 0.35; delta1 <- 0
    gamma1_raw <- inverse_transform_gamma(gamma1)
    delta1_raw <- inverse_transform_delta(delta1)
    beta_0_trig1 <- log(0.1); decay1 <- 0.2
    beta_0_trig1_raw <- inverse_transform_beta_0_trig(beta_0_trig1)
    decay1_raw <- inverse_transform_decay(decay1)
    beta_0_1 <- log(observed_rate) - gamma1 * mean_log_pop - delta1 * mean_poverty

    starting_points[[1]] <- list(
      name = "Data-driven",
      params = list(
        beta_0_bg = beta_0_1, gamma_raw = gamma1_raw, delta_raw = delta1_raw,
        # Year effects
        beta_2016_raw = 0, beta_2017_raw = 0, beta_2018_raw = 0,
        beta_2019_raw = 0, beta_2020_raw = 0, beta_2021_raw = 0,
        beta_2022_raw = 0, beta_2023_raw = 0, beta_2024_raw = 0,
        # Triggering parameters (5 marks + decay)
        beta_0_trig_raw = beta_0_trig1_raw,
        beta_fatal_raw = 0, beta_violent_raw = 0,
        beta_student_raw = 0, beta_labor_raw = 0,
        beta_papua_raw = 0,
        decay_raw = decay1_raw
      )
    )

    if(use_multistart && n_starts >= 2) {
      # 2. Conservative
      gamma2 <- 0.20; delta2 <- 0.5
      gamma2_raw <- inverse_transform_gamma(gamma2)
      delta2_raw <- inverse_transform_delta(delta2)
      beta_0_trig2 <- log(0.15); decay2 <- 0.25
      beta_0_trig2_raw <- inverse_transform_beta_0_trig(beta_0_trig2)
      decay2_raw <- inverse_transform_decay(decay2)
      beta_0_2 <- log(observed_rate) - gamma2 * mean_log_pop - delta2 * mean_poverty

      starting_points[[2]] <- list(
        name = "Conservative",
        params = list(
          beta_0_bg = beta_0_2, gamma_raw = gamma2_raw, delta_raw = delta2_raw,
          beta_2016_raw = -0.1, beta_2017_raw = -0.05, beta_2018_raw = 0,
          beta_2019_raw = 0, beta_2020_raw = 0.05, beta_2021_raw = 0.1,
          beta_2022_raw = 0.1, beta_2023_raw = 0.05, beta_2024_raw = 0,
          beta_0_trig_raw = beta_0_trig2_raw,
          beta_fatal_raw = 0.2, beta_violent_raw = 0.1,
          beta_student_raw = 0.1, beta_labor_raw = 0.1,
          beta_papua_raw = 0.1,
          decay_raw = decay2_raw
        )
      )
    }

    if(use_multistart && n_starts >= 3) {
      # 3. Proportional
      gamma3 <- 0.80; delta3 <- 0.2
      gamma3_raw <- inverse_transform_gamma(gamma3)
      delta3_raw <- inverse_transform_delta(delta3)
      beta_0_trig3 <- log(0.08); decay3 <- 0.15
      beta_0_trig3_raw <- inverse_transform_beta_0_trig(beta_0_trig3)
      decay3_raw <- inverse_transform_decay(decay3)
      beta_0_3 <- log(observed_rate) - gamma3 * mean_log_pop - delta3 * mean_poverty

      starting_points[[3]] <- list(
        name = "Proportional",
        params = list(
          beta_0_bg = beta_0_3, gamma_raw = gamma3_raw, delta_raw = delta3_raw,
          beta_2016_raw = 0.2, beta_2017_raw = 0.15, beta_2018_raw = 0.1,
          beta_2019_raw = 0, beta_2020_raw = -0.1, beta_2021_raw = -0.15,
          beta_2022_raw = -0.1, beta_2023_raw = 0, beta_2024_raw = 0.1,
          beta_0_trig_raw = beta_0_trig3_raw,
          beta_fatal_raw = -0.2, beta_violent_raw = -0.1,
          beta_student_raw = -0.1, beta_labor_raw = -0.2,
          beta_papua_raw = 0,
          decay_raw = decay3_raw
        )
      )
    }

    if(use_multistart && n_starts >= 4) {
      # 4. High fatality effect
      gamma4 <- 0.50; delta4 <- 0
      gamma4_raw <- inverse_transform_gamma(gamma4)
      delta4_raw <- inverse_transform_delta(delta4)
      beta_0_trig4 <- log(0.12); decay4 <- 0.18
      beta_0_trig4_raw <- inverse_transform_beta_0_trig(beta_0_trig4)
      decay4_raw <- inverse_transform_decay(decay4)
      beta_0_4 <- log(observed_rate) - gamma4 * mean_log_pop - delta4 * mean_poverty

      starting_points[[4]] <- list(
        name = "High-fatality",
        params = list(
          beta_0_bg = beta_0_4, gamma_raw = gamma4_raw, delta_raw = delta4_raw,
          beta_2016_raw = 0.05, beta_2017_raw = -0.05, beta_2018_raw = 0.1,
          beta_2019_raw = -0.1, beta_2020_raw = 0, beta_2021_raw = 0.05,
          beta_2022_raw = -0.05, beta_2023_raw = 0.1, beta_2024_raw = -0.1,
          beta_0_trig_raw = beta_0_trig4_raw,
          beta_fatal_raw = 0.5, beta_violent_raw = 0.2,
          beta_student_raw = 0, beta_labor_raw = 0,
          beta_papua_raw = 0.2,
          decay_raw = decay4_raw
        )
      )
    }

    if(use_multistart && n_starts >= 5) {
      # 5. Papua-focused
      gamma5 <- 0.40; delta5 <- 0.3
      gamma5_raw <- inverse_transform_gamma(gamma5)
      delta5_raw <- inverse_transform_delta(delta5)
      beta_0_trig5 <- log(0.1); decay5 <- 0.2
      beta_0_trig5_raw <- inverse_transform_beta_0_trig(beta_0_trig5)
      decay5_raw <- inverse_transform_decay(decay5)
      beta_0_5 <- log(observed_rate) - gamma5 * mean_log_pop - delta5 * mean_poverty

      starting_points[[5]] <- list(
        name = "Papua-focused",
        params = list(
          beta_0_bg = beta_0_5, gamma_raw = gamma5_raw, delta_raw = delta5_raw,
          beta_2016_raw = 0, beta_2017_raw = 0, beta_2018_raw = 0,
          beta_2019_raw = 0.2, beta_2020_raw = 0.1, beta_2021_raw = 0,
          beta_2022_raw = 0, beta_2023_raw = 0, beta_2024_raw = 0,
          beta_0_trig_raw = beta_0_trig5_raw,
          beta_fatal_raw = 0.3, beta_violent_raw = 0.1,
          beta_student_raw = -0.2, beta_labor_raw = -0.3,
          beta_papua_raw = 0.5,
          decay_raw = decay5_raw
        )
      )
    }
  }

  # Parameter names (19 total: 12 background + 7 triggering)
  param_names <- c("beta_0_bg", "gamma_raw", "delta_raw",
                   "beta_2016_raw", "beta_2017_raw", "beta_2018_raw",
                   "beta_2019_raw", "beta_2020_raw", "beta_2021_raw",
                   "beta_2022_raw", "beta_2023_raw", "beta_2024_raw",
                   "beta_0_trig_raw", "beta_fatal_raw",
                   "beta_violent_raw", "beta_student_raw",
                   "beta_labor_raw", "beta_papua_raw", "decay_raw")

  # Bounds
  lower_bounds <- c(beta_0_bg = -15, gamma_raw = -10, delta_raw = -10,
                    beta_2016_raw = -10, beta_2017_raw = -10, beta_2018_raw = -10,
                    beta_2019_raw = -10, beta_2020_raw = -10, beta_2021_raw = -10,
                    beta_2022_raw = -10, beta_2023_raw = -10, beta_2024_raw = -10,
                    beta_0_trig_raw = -10, beta_fatal_raw = -10,
                    beta_violent_raw = -10, beta_student_raw = -10,
                    beta_labor_raw = -10, beta_papua_raw = -10, decay_raw = -10)

  upper_bounds <- c(beta_0_bg = 5, gamma_raw = 10, delta_raw = 10,
                    beta_2016_raw = 10, beta_2017_raw = 10, beta_2018_raw = 10,
                    beta_2019_raw = 10, beta_2020_raw = 10, beta_2021_raw = 10,
                    beta_2022_raw = 10, beta_2023_raw = 10, beta_2024_raw = 10,
                    beta_0_trig_raw = 10, beta_fatal_raw = 10,
                    beta_violent_raw = 10, beta_student_raw = 10,
                    beta_labor_raw = 10, beta_papua_raw = 10, decay_raw = 10)

  # Objective functions
  obj_fun_silent <- function(par) {
    params_full <- as.list(par)
    names(params_full) <- param_names
    ll <- loglik_extended(params_full, times, marks, district_years,
                          verbose = FALSE, model_id = NULL)
    return(-ll)
  }

  obj_fun_verbose <- function(par) {
    params_full <- as.list(par)
    names(params_full) <- param_names
    ll <- loglik_extended(params_full, times, marks, district_years,
                          verbose = TRUE, model_id = model_id)
    return(-ll)
  }

  # Run optimization from each starting point
  start_time <- Sys.time()
  all_results <- list()

  cat(sprintf("\n=== MULTI-START OPTIMIZATION (%d starts) ===\n\n", length(starting_points)))

  for(i in seq_along(starting_points)) {
    start_point <- starting_points[[i]]
    cat(sprintf("--- Start %d/%d: %s ---\n", i, length(starting_points), start_point$name))
    gamma_init <- transform_gamma(start_point$params$gamma_raw)
    cat(sprintf("  Initial: γ=%.3f\n", gamma_init))

    param_vec <- unlist(start_point$params[param_names])
    names(param_vec) <- param_names

    use_obj <- if(i == 1 && length(starting_points) == 1) obj_fun_verbose else obj_fun_silent

    if(i == 1) {
      reset_optim_state(model_name)
    }

    start_i_time <- Sys.time()

    fit_i <- optim(
      par = param_vec,
      fn = use_obj,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(maxit = if(TEST_MODE) 100 else 1000, factr = 1e7, trace = 0)
    )

    end_i_time <- Sys.time()
    runtime_i <- difftime(end_i_time, start_i_time, units = "mins")

    all_results[[i]] <- list(
      name = start_point$name,
      params = as.list(fit_i$par),
      loglik = -fit_i$value,
      convergence = fit_i$convergence,
      runtime = as.numeric(runtime_i)
    )

    gamma_result <- transform_gamma(all_results[[i]]$params$gamma_raw)
    cat(sprintf("  → LL: %.2f | γ=%.4f | Conv: %d | Time: %.1f min\n\n",
               all_results[[i]]$loglik, gamma_result,
               all_results[[i]]$convergence, all_results[[i]]$runtime))
  }

  # Find best result
  best_ll <- -Inf
  best_idx <- 1
  for(i in seq_along(all_results)) {
    if(all_results[[i]]$loglik > best_ll) {
      best_ll <- all_results[[i]]$loglik
      best_idx <- i
    }
  }

  cat(sprintf("╔══════════════════════════════════════════════════════════════╗\n"))
  cat(sprintf("║  BEST START: %s (LL=%.2f)\n", all_results[[best_idx]]$name, best_ll))
  cat(sprintf("╚══════════════════════════════════════════════════════════════╝\n\n"))

  # Re-run best with verbose if needed
  if(best_idx != 1 || length(starting_points) > 1) {
    cat("Re-running best initialization with full tracking...\n")
    reset_optim_state(model_name)

    param_vec <- unlist(all_results[[best_idx]]$params[param_names])
    names(param_vec) <- param_names

    fit <- optim(
      par = param_vec,
      fn = obj_fun_verbose,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(maxit = if(TEST_MODE) 100 else 1000, factr = 1e7, trace = 1, REPORT = 10)
    )

    final_ll <- -fit$value
    final_params <- as.list(fit$par)
    names(final_params) <- param_names
    final_convergence <- fit$convergence
  } else {
    final_ll <- all_results[[1]]$loglik
    final_params <- all_results[[1]]$params
    final_convergence <- all_results[[1]]$convergence
  }

  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")

  cat(sprintf("\n  ✓ Optimization complete in %.2f minutes\n", as.numeric(runtime)))
  cat(sprintf("  Log-likelihood: %.2f\n", final_ll))
  cat(sprintf("  Convergence: %d\n", final_convergence))

  # Assemble results
  n_params <- 19
  results <- list(
    model_name = model_name,
    model_id = model_id,
    params = final_params,
    loglik = final_ll,
    convergence = final_convergence,
    n_params = n_params,
    n_events = length(times),
    AIC = 2*n_params - 2*final_ll,
    BIC = n_params*log(length(times)) - 2*final_ll,
    runtime_mins = as.numeric(runtime),
    multistart_results = if(use_multistart) all_results else NULL,
    best_start = if(use_multistart) all_results[[best_idx]]$name else NULL
  )

  # Save final checkpoint
  save_final_checkpoint(results, model_id)

  return(results)
}

# ====================================================================
# 6. FIT MODEL
# ====================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║   FITTING MODEL 5: EXTENDED MARKS (EXPONENTIAL)              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

times_sample <- marks_data$time
marks_sample <- marks_data
district_years_sample <- district_years

model5_fit <- fit_model_extended(
  times_sample,
  marks_sample,
  district_years_sample,
  model_name = "Model 5: Extended Marks (Exponential)",
  model_id = "MODEL5_EXTENDED",
  use_multistart = !TEST_MODE,
  n_starts = if(TEST_MODE) 1 else 5
)

# ====================================================================
# 7. RESULTS
# ====================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║         MODEL 5 RESULTS                                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("MODEL FIT:\n")
cat(sprintf("  Log-likelihood: %.2f\n", model5_fit$loglik))
cat(sprintf("  AIC: %.2f\n", model5_fit$AIC))
cat(sprintf("  BIC: %.2f\n", model5_fit$BIC))
cat(sprintf("  Parameters: %d\n", model5_fit$n_params))
cat(sprintf("  Runtime: %.2f minutes\n\n", model5_fit$runtime_mins))

cat("BACKGROUND RATE PARAMETERS:\n")
cat("---------------------------\n")
gamma_final <- transform_gamma(model5_fit$params$gamma_raw)
delta_final <- transform_delta(model5_fit$params$delta_raw)

cat(sprintf("  β₀ (baseline): %.4f\n", model5_fit$params$beta_0_bg))
cat(sprintf("  γ (population elasticity): %.4f\n", gamma_final))
cat(sprintf("    → 10× population → %.1f× more protests\n", 10^gamma_final))
cat(sprintf("  δ (poverty effect): %.4f\n", delta_final))
cat(sprintf("    → +10%% poverty → %.3f× more protests\n\n", exp(delta_final * 0.10)))

cat("TRIGGERING PARAMETERS (violence-focused):\n")
cat("------------------------------------------\n")
beta_0_trig_final <- transform_beta_0_trig(model5_fit$params$beta_0_trig_raw)
beta_fatal_final <- transform_beta_mark(model5_fit$params$beta_fatal_raw)
beta_violent_final <- transform_beta_mark(model5_fit$params$beta_violent_raw)
beta_student_final <- transform_beta_mark(model5_fit$params$beta_student_raw)
beta_labor_final <- transform_beta_mark(model5_fit$params$beta_labor_raw)
beta_papua_final <- transform_beta_mark(model5_fit$params$beta_papua_raw)
decay_final <- transform_decay(model5_fit$params$decay_raw)

cat(sprintf("  β₀_trig (baseline): %.4f\n", beta_0_trig_final))
cat(sprintf("  β_fatal: %.4f (exp = %.3f)\n", beta_fatal_final, exp(beta_fatal_final)))
cat(sprintf("  β_violent: %.4f (exp = %.3f)\n", beta_violent_final, exp(beta_violent_final)))
cat(sprintf("  β_student: %.4f (exp = %.3f)\n", beta_student_final, exp(beta_student_final)))
cat(sprintf("  β_labor: %.4f (exp = %.3f)\n", beta_labor_final, exp(beta_labor_final)))
cat(sprintf("  β_papua: %.4f (exp = %.3f)\n", beta_papua_final, exp(beta_papua_final)))
cat(sprintf("  Decay: %.5f (half-life = %.1f days)\n\n", decay_final, log(2)/decay_final))

# Save results
saveRDS(model5_fit, "model5_extended_marks.rds")
cat("✓ Saved: model5_extended_marks.rds\n\n")

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   MODEL 5 ESTIMATION COMPLETE                                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
