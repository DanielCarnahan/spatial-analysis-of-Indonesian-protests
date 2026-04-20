#!/usr/bin/env Rscript
# ==============================================================================
# MODEL 1b: BASIC HAWKES WITH POWER-LAW KERNEL (NO MARKS)
# ==============================================================================
#
# Specification: λ(t) = μ(d,y) + α · (Δt + c)^(-β)
#
# PARAMETERS (15 total):
#   Background: β₀, γ, δ + 9 year effects
#   Triggering: α (CONSTANT - no mark dependence)
#   Kernel: power-law exponent β, offset c
#
# PURPOSE:
#   - Establish baseline for power-law kernel WITHOUT mark effects
#   - Compare to Model 4 to test if marks matter (both use power-law)
#   - Compare to Model 1 to test if power-law > exponential (both no marks)
#
# NESTED MODEL TESTS:
#   - Model 0 vs Model 1b: Does triggering exist? (LR test, df=3)
#   - Model 1b vs Model 4: Do marks matter? (LR test, df=4)
#
# ==============================================================================

library(dplyr)
library(Rcpp)

# Compile C++ likelihood function
cat("Compiling C++ likelihood function for basic power-law Hawkes model...\n")
sourceCpp("hawkes_powerlaw_basic_likelihood.cpp")
cat("✓ C++ compiled\n\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   MODEL 1b: BASIC HAWKES + POWER-LAW KERNEL                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("SPECIFICATION:\n")
cat("  Background: μ(d,y) = exp(β₀ + γ·log(pop) + δ·poverty + Σ β_year)\n")
cat("  Triggering: α = exp(β₀_trig)  [CONSTANT - NO MARKS]\n")
cat("  Kernel: g(Δt) = (Δt + c)^(-α)  [POWER-LAW]\n")
cat("  Intensity: λ(t) = μ(d,y) + α · Σ g(t - t_i)\n\n")

cat("HYPOTHESIS:\n")
cat("  → Tests if triggering exists with power-law kernel\n")
cat("  → Baseline for testing mark effects (vs Model 4)\n")
cat("  → Tests power-law vs exponential (vs Model 1)\n\n")

# ====================================================================
# CONFIGURATION
# ====================================================================

TEST_MODE <- FALSE
TEMPORAL_CUTOFF <- 90  # 90 days
CHECKPOINT_EVERY_N_EVALS <- 1

# ====================================================================
# PARAMETER TRANSFORMATION FUNCTIONS
# ====================================================================

# Background parameters (same as all models)
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

# Triggering base parameter
transform_beta_0_trig <- function(beta_0_trig_raw) {
  -10 + 20 * plogis(beta_0_trig_raw)
}
inverse_transform_beta_0_trig <- function(beta_0_trig) {
  qlogis((beta_0_trig + 10) / 20)
}

# Power-law kernel parameters
transform_alpha <- function(alpha_raw) {
  1.05 + 2.45 * plogis(alpha_raw)  # [1.05, 3.5]
}
inverse_transform_alpha <- function(alpha) {
  qlogis((alpha - 1.05) / 2.45)
}

transform_c <- function(c_raw) {
  0.01 + 0.99 * plogis(c_raw)  # [0.01, 1.0]
}
inverse_transform_c <- function(c) {
  qlogis((c - 0.01) / 0.99)
}

# ====================================================================
# CHECKPOINT SYSTEM
# ====================================================================

checkpoint_dir <- "checkpoints_model1b_powerlaw_basic"
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
    eval_num = eval_num,
    params = params,
    loglik = ll,
    timestamp = Sys.time()
  ), checkpoint_file)
}

# ====================================================================
# 1. LOAD AND PREPARE DATA
# ====================================================================

cat("=== LOADING DATA ===\n")

protests <- readRDS("protests_with_poverty.rds")
cat(sprintf("Loaded %d events with poverty data\n", nrow(protests)))
cat(sprintf("Time range: %d to %d \n", min(protests$year, na.rm=TRUE),
            max(protests$year, na.rm=TRUE)))

# Use test mode or full dataset
if (TEST_MODE) {
  protests <- protests %>% sample_n(min(500, nrow(protests)))
  cat(sprintf("✓ TEST MODE: Using %d events\n\n", nrow(protests)))
} else {
  cat(sprintf("✓ Using full dataset: %d events across %d years\n",
             nrow(protests), length(unique(protests$year))))
}

# Prepare data (NO mark indicators needed for Model 1b)
marks_data <- protests %>%
  select(event_id_cnty, time = days_since_start, district = event_id_cnty, year,
         log_pop, population, poverty_decimal) %>%
  filter(!is.na(time), !is.na(log_pop), !is.na(poverty_decimal))

cat(sprintf("\nUsing complete cases: %d events (%.1f%% of total)\n\n",
           nrow(marks_data), 100*nrow(marks_data)/nrow(protests)))

# ====================================================================
# 2. CALCULATE DISTRICT-YEAR OBSERVATION TIMES
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
# 3. MODEL IMPLEMENTATION
# ====================================================================

cat("=== IMPLEMENTING BASIC POWER-LAW HAWKES MODEL ===\n\n")

model_functions <- function() {

  loglik_powerlaw_basic <- function(params, times, marks, district_years,
                                  verbose = TRUE, model_id = NULL) {

    .optim_state$n_evals <- .optim_state$n_evals + 1
    eval_num <- .optim_state$n_evals

    # Prepare parameter vector for C++ (15 parameters)
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
      # Triggering base (constant, no marks)
      params[["beta_0_trig_raw"]],
      # Power-law kernel (2 total)
      params[["alpha_raw"]],
      params[["c_raw"]]
    )

    # Call C++ function (returns NEGATIVE log-likelihood)
    neg_ll <- hawkes_powerlaw_basic_negloglik_cpp(
      times = times,
      log_pop = marks$log_pop,
      poverty_decimal = marks$poverty_decimal,
      year = as.integer(marks$year),
      district_year_exposure = district_years$days_observed,
      district_year_log_pop = district_years$log_pop,
      district_year_poverty = district_years$poverty_decimal,
      district_year_year = as.integer(district_years$year),
      params = param_vec,
      temporal_cutoff = TEMPORAL_CUTOFF
    )

    ll_final <- -neg_ll

    # Update best result
    if (ll_final > .optim_state$best_ll) {
      .optim_state$best_ll <- ll_final
      .optim_state$best_params <- params
    }

    # Checkpoint
    if (eval_num %% CHECKPOINT_EVERY_N_EVALS == 0 && !is.null(model_id)) {
      save_intermediate_checkpoint(model_id, eval_num, params, ll_final)

      # Display progress
      alpha <- transform_alpha(params[["alpha_raw"]])
      c_val <- transform_c(params[["c_raw"]])

      cat(sprintf("\n    💾 Checkpoint %d: LL=%.2f (best: %.2f) | α=%.3f c=%.4f\n",
                 eval_num, ll_final, .optim_state$best_ll, alpha, c_val))

      if (verbose) {
        cat(sprintf("    LL: %.2f | α=%.3f c=%.4f [eval %d]\n",
                   ll_final, alpha, c_val, eval_num))
      }
    }

    return(neg_ll)
  }

  list(loglik = loglik_powerlaw_basic)
}

# ====================================================================
# 4. OPTIMIZATION WITH MULTI-START
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  STARTING OPTIMIZATION                                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("--- Fitting Model 1b: Basic Hawkes + Power-Law Kernel ---\n")

use_multistart <- TRUE
n_starts <- 3  # Fewer starts for basic model

model_funcs <- model_functions()
objective_function <- model_funcs$loglik

cat(sprintf("Using MULTI-START optimization with %d starting points\n", n_starts))
cat("Starting fresh optimization...\n\n")

# ================================================================
# STARTING POINTS
# ================================================================

starting_points <- list()

# Load Model 2 (or Model 0) for baseline background parameters
if (file.exists("model_poverty.rds")) {
  model2 <- readRDS("model_poverty.rds")
  baseline_params <- model2$params[1:12]  # background + year effects
} else {
  # If Model 2 doesn't exist, use neutral starting values
  baseline_params <- list(
    beta_0_bg = -5,
    gamma_raw = 0,
    delta_raw = 0,
    beta_2016_raw = 0, beta_2017_raw = 0, beta_2018_raw = 0,
    beta_2019_raw = 0, beta_2020_raw = 0, beta_2021_raw = 0,
    beta_2022_raw = 0, beta_2023_raw = 0, beta_2024_raw = 0
  )
}

# 1. Data-driven: α=2.0, c=0.1 (PRIMARY)
alpha1 <- 2.0; c1 <- 0.1
starting_points[[1]] <- list(
  name = "Data-driven (α=2.0, c=0.1)",
  params = c(
    baseline_params,
    list(
      beta_0_trig_raw = 0,  # Neutral starting value
      alpha_raw = inverse_transform_alpha(alpha1),
      c_raw = inverse_transform_c(c1)
    )
  )
)

if(use_multistart && n_starts >= 2) {
  # 2. Fast-decay: α=2.5, c=0.05
  alpha2 <- 2.5; c2 <- 0.05
  starting_points[[2]] <- list(
    name = "Fast-decay (α=2.5, c=0.05)",
    params = c(
      baseline_params,
      list(
        beta_0_trig_raw = 0,
        alpha_raw = inverse_transform_alpha(alpha2),
        c_raw = inverse_transform_c(c2)
      )
    )
  )
}

if(use_multistart && n_starts >= 3) {
  # 3. Slow-tail: α=1.5, c=0.3
  alpha3 <- 1.5; c3 <- 0.3
  starting_points[[3]] <- list(
    name = "Slow-tail (α=1.5, c=0.3)",
    params = c(
      baseline_params,
      list(
        beta_0_trig_raw = 0,
        alpha_raw = inverse_transform_alpha(alpha3),
        c_raw = inverse_transform_c(c3)
      )
    )
  )
}

# Parameter names (15 total)
param_names <- c("beta_0_bg", "gamma_raw", "delta_raw",
                "beta_2016_raw", "beta_2017_raw", "beta_2018_raw",
                "beta_2019_raw", "beta_2020_raw", "beta_2021_raw",
                "beta_2022_raw", "beta_2023_raw", "beta_2024_raw",
                "beta_0_trig_raw",
                "alpha_raw", "c_raw")

# Run multi-start optimization
all_results <- list()

for(s in seq_along(starting_points)) {
  sp <- starting_points[[s]]

  cat(sprintf("║  Starting Point %d/%d: %-42s ║\n", s, length(starting_points), sp$name))
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")

  # Extract starting values
  par0 <- unlist(sp$params[param_names])

  # Display starting values
  alpha_start <- transform_alpha(par0[["alpha_raw"]])
  c_start <- transform_c(par0[["c_raw"]])
  cat(sprintf("Starting values: α=%.3f, c=%.4f\n\n", alpha_start, c_start))

  # Reset optimization state
  reset_optim_state(sprintf("model1b_start%d", s))

  # Run optimization
  start_time <- Sys.time()

  result <- optim(
    par = par0,
    fn = function(params_vec) {
      params_list <- as.list(params_vec)
      names(params_list) <- param_names
      objective_function(params_list, marks_data$time, marks_data, district_years,
                        verbose = TRUE, model_id = sprintf("model1b_start%d", s))
    },
    method = "L-BFGS-B",
    control = list(maxit = 1000, factr = 1e7)
  )

  runtime <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

  # Transform parameters
  params_final <- as.list(result$par)
  names(params_final) <- param_names

  alpha_final <- transform_alpha(params_final$alpha_raw)
  c_final <- transform_c(params_final$c_raw)

  # Store result
  all_results[[s]] <- list(
    starting_point = sp$name,
    params = params_final,
    loglik = -result$value,
    convergence = result$convergence,
    runtime = runtime,
    alpha = alpha_final,
    c = c_final
  )

  cat(sprintf("\n✓ Optimization complete\n"))
  cat(sprintf("  Log-likelihood: %.2f\n", -result$value))
  cat(sprintf("  Convergence: %d\n", result$convergence))
  cat(sprintf("  Runtime: %.1f minutes\n", runtime))
  cat(sprintf("  Final: α=%.3f, c=%.4f\n\n", alpha_final, c_final))

  if (s < length(starting_points)) {
    cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  }
}

# ====================================================================
# 5. SELECT BEST RESULT AND SAVE
# ====================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MULTI-START RESULTS SUMMARY                                 ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

for(i in seq_along(all_results)) {
  res <- all_results[[i]]
  cat(sprintf("%d. %s\n", i, res$starting_point))
  cat(sprintf("   LL = %.2f | α=%.3f, c=%.4f | Conv=%d | Time=%.1f min\n\n",
             res$loglik, res$alpha, res$c, res$convergence, res$runtime))
}

# Find best by log-likelihood
logliks <- sapply(all_results, function(x) x$loglik)
best_idx <- which.max(logliks)
best_result <- all_results[[best_idx]]

cat(sprintf("✓ Best result: Starting point %d (%s)\n", best_idx, best_result$starting_point))
cat(sprintf("  Log-likelihood: %.2f\n", best_result$loglik))
cat(sprintf("  α = %.4f (power-law exponent)\n", best_result$alpha))
cat(sprintf("  c = %.4f days (~%.1f hours offset)\n\n", best_result$c, best_result$c * 24))

# Save best model
model1b_result <- list(
  model_name = "Model 1b: Basic Hawkes + Power-Law Kernel",
  loglik = best_result$loglik,
  convergence = best_result$convergence,
  n_params = length(param_names),
  n_events = nrow(marks_data),
  params = best_result$params,
  alpha = best_result$alpha,
  c = best_result$c,
  all_results = all_results,
  runtime = sum(sapply(all_results, function(x) x$runtime)),
  timestamp = Sys.time()
)

saveRDS(model1b_result, "model1b_powerlaw_basic.rds")
cat("✓ Saved: model1b_powerlaw_basic.rds\n\n")

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL 1b ESTIMATION COMPLETE                                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("FINAL PARAMETERS:\n")
cat(sprintf("  α (power-law exponent) = %.4f\n", best_result$alpha))
cat(sprintf("  c (offset) = %.4f days (%.1f hours)\n", best_result$c, best_result$c * 24))
cat(sprintf("  Log-likelihood = %.2f\n", best_result$loglik))
cat(sprintf("  Triggering: CONSTANT (no mark dependence)\n\n"))
