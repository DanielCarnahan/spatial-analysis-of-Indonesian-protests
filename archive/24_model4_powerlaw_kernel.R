#!/usr/bin/env Rscript
# ==============================================================================
# MODEL 4: POWER-LAW KERNEL HAWKES PROCESS
# ==============================================================================
#
# Specification: λ(t) = μ(d,y) + Σ α(marks) · (Δt + c)^(-α)
#
# NEW PARAMETERS (vs Model 2):
#   α ∈ [1.05, 3.5]: Power-law exponent (controls decay rate)
#   c ∈ [0.01, 1.0]: Offset parameter (prevents singularity at Δt=0)
#
# KERNEL: g(Δt) = (Δt + c)^(-α)
#
# INTEGRAL: ∫₀^Δt (s + c)^(-α) ds = [c^(1-α) - (Δt+c)^(1-α)] / (α - 1)
#
# PURPOSE:
#   - Address Model 2's unrealistic slow decay (β = 0.0014, half-life = 477 days)
#   - Model 3 (mixed exponential) FAILED - worse fit, boundary issues
#   - Power-law allows heavier tails than exponential (better long-term fit)
#   - Common in social contagion literature (Crane & Sornette 2008)
#
# COMPARISON:
#   - Model 2 vs Model 4: LR test (df=1, replaces β with α, c)
#   - Expected: LR > 0, p < 0.001, power-law improves over exponential
#
# ==============================================================================

library(dplyr)
library(Rcpp)

# Compile C++ likelihood function
cat("Compiling C++ likelihood function for power-law kernel Hawkes model...\n")
sourceCpp("hawkes_powerlaw_kernel_likelihood.cpp")
cat("✓ C++ compiled\n\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   MODEL 4: POWER-LAW KERNEL HAWKES                           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("SPECIFICATION:\n")
cat("  Background: μ(d,y) = exp(β₀ + γ·log(pop) + δ·poverty + Σ β_year)\n")
cat("  Triggering: α(marks) = exp(β₀_trig + β_riot + β_fatal + β_student + β_labor)\n")
cat("  Kernel: g(Δt) = (Δt + c)^(-α)  [POWER-LAW!]\n")
cat("  Intensity: λ(t) = μ(d,y) + Σ α(marks_i) · g(t - t_i)\n\n")

cat("HYPOTHESIS:\n")
cat("  → Power-law better captures empirical ACF than exponential\n")
cat("  → α ≈ 2.0 (moderate decay), c ≈ 0.15 days (~3.6 hours)\n")
cat("  → Resolves Model 2's 477-day half-life issue\n\n")

# ====================================================================
# CONFIGURATION
# ====================================================================

TEST_MODE <- FALSE
TEMPORAL_CUTOFF <- 90  # 90 days
CHECKPOINT_EVERY_N_EVALS <- 1

# ====================================================================
# PARAMETER TRANSFORMATION FUNCTIONS
# ====================================================================

# Background parameters (same as Model 2)
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

# Triggering mark parameters (same as Model 2)
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

# NEW: Power-law kernel parameters
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

checkpoint_dir <- "checkpoints_model4_powerlaw_kernel"
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

load_latest_checkpoint <- function(model_id) {
  checkpoint_files <- list.files(checkpoint_dir,
                                 pattern = sprintf("intermediate-%s-eval.*\\.rds", model_id),
                                 full.names = TRUE)

  if (length(checkpoint_files) == 0) {
    return(NULL)
  }

  latest_file <- checkpoint_files[order(file.mtime(checkpoint_files), decreasing = TRUE)[1]]

  checkpoint_data <- readRDS(latest_file)

  return(checkpoint_data)
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

# Check poverty data coverage
poverty_coverage <- sum(!is.na(protests$poverty_decimal))
cat(sprintf("Poverty data coverage: %d (%.1f%%)\n",
           poverty_coverage, 100*poverty_coverage/nrow(protests)))

# Create mark indicators
marks_data <- protests %>%
  mutate(
    is_riot = as.integer(event_type == "Riots"),
    is_fatal = as.integer(has_fatalities > 0),
    is_student = as.integer(grepl("Student", assoc_actor_1, ignore.case = TRUE)),
    is_labor = as.integer(grepl("Labor|Worker", assoc_actor_1, ignore.case = TRUE))
  ) %>%
  select(event_id_cnty, time = days_since_start, district = event_id_cnty, year,
         log_pop, population, poverty_decimal,
         is_riot, is_fatal, is_student, is_labor)

# Filter to complete cases
marks_data <- marks_data %>%
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

cat("=== IMPLEMENTING POWER-LAW KERNEL HAWKES MODEL ===\n\n")

model_powerlaw_hawkes_functions <- function() {

  loglik_powerlaw_hawkes <- function(params, times, marks, district_years,
                                  verbose = TRUE, model_id = NULL) {

    .optim_state$n_evals <- .optim_state$n_evals + 1
    eval_num <- .optim_state$n_evals

    # Prepare parameter vector for C++ (19 parameters: 12 background + 5 marks + 2 kernel)
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
      # Triggering parameters (5 total: base + 4 marks)
      params[["beta_0_trig_raw"]],
      params[["beta_riot_raw"]],
      params[["beta_fatal_raw"]],
      params[["beta_student_raw"]],
      params[["beta_labor_raw"]],
      # Power-law kernel (2 total)
      params[["alpha_raw"]],
      params[["c_raw"]]
    )

    # Call C++ function (returns NEGATIVE log-likelihood)
    neg_ll <- hawkes_powerlaw_negloglik_cpp(
      times = times,
      log_pop = marks$log_pop,
      poverty_decimal = marks$poverty_decimal,
      year = as.integer(marks$year),
      is_riot = as.integer(marks$is_riot),
      has_fatalities = as.integer(marks$is_fatal),
      is_student = as.integer(marks$is_student),
      is_labor = as.integer(marks$is_labor),
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

  list(loglik = loglik_powerlaw_hawkes)
}

# ====================================================================
# 4. PREPARE DATA FOR ESTIMATION
# ====================================================================

cat("=== PREPARING DATA FOR ESTIMATION ===\n")
cat(sprintf("Sample size: %d events\n", nrow(marks_data)))
cat(sprintf("District-years: %d\n", nrow(district_years)))

times_vec <- marks_data$time
cat("✓ Ready for optimization\n\n")

# ====================================================================
# 5. OPTIMIZATION WITH MULTI-START
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  STARTING OPTIMIZATION                                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("--- Fitting Model 4: Power-Law Kernel Hawkes ---\n")

use_multistart <- TRUE
n_starts <- 5

model_funcs <- model_powerlaw_hawkes_functions()
objective_function <- model_funcs$loglik

# Check for existing checkpoint
checkpoint <- load_latest_checkpoint("model4_powerlaw")

if (!is.null(checkpoint)) {
  cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  CHECKPOINT FOUND - RESUMING FROM PREVIOUS RUN               ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")

  cat(sprintf("  Checkpoint from: %s\n", checkpoint$timestamp))
  cat(sprintf("  Evaluation: %d\n", checkpoint$eval_num))
  cat(sprintf("  Log-likelihood: %.2f\n\n", checkpoint$loglik))

  # Resume from checkpoint
  starting_points <- list(list(name = "Checkpoint", params = checkpoint$params))

} else {
  cat(sprintf("Using MULTI-START optimization with %d starting points\n", n_starts))
  cat("Starting fresh optimization...\n\n")

  # ================================================================
  # STARTING POINTS
  # ================================================================

  starting_points <- list()

  # Load Model 2 estimates to use as baseline
  model2 <- readRDS("model_poverty.rds")

  # 1. Data-driven: α=2.0, c=0.1 (PRIMARY)
  alpha1 <- 2.0; c1 <- 0.1
  alpha1_raw <- inverse_transform_alpha(alpha1)
  c1_raw <- inverse_transform_c(c1)

  starting_points[[1]] <- list(
    name = "Data-driven (α=2.0, c=0.1)",
    params = c(
      model2$params[1:17],  # Use Model 2's background + mark estimates
      list(
        alpha_raw = alpha1_raw,
        c_raw = c1_raw
      )
    )
  )

  if(use_multistart && n_starts >= 2) {
    # 2. Fast-decay: α=2.5, c=0.05
    alpha2 <- 2.5; c2 <- 0.05
    starting_points[[2]] <- list(
      name = "Fast-decay (α=2.5, c=0.05)",
      params = c(
        model2$params[1:17],
        list(
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
        model2$params[1:17],
        list(
          alpha_raw = inverse_transform_alpha(alpha3),
          c_raw = inverse_transform_c(c3)
        )
      )
    )
  }

  if(use_multistart && n_starts >= 4) {
    # 4. Moderate: α=1.8, c=0.15
    alpha4 <- 1.8; c4 <- 0.15
    starting_points[[4]] <- list(
      name = "Moderate (α=1.8, c=0.15)",
      params = c(
        model2$params[1:17],
        list(
          alpha_raw = inverse_transform_alpha(alpha4),
          c_raw = inverse_transform_c(c4)
        )
      )
    )
  }

  if(use_multistart && n_starts >= 5) {
    # 5. High-c: α=2.2, c=0.5
    alpha5 <- 2.2; c5 <- 0.5
    starting_points[[5]] <- list(
      name = "High-c (α=2.2, c=0.5)",
      params = c(
        model2$params[1:17],
        list(
          alpha_raw = inverse_transform_alpha(alpha5),
          c_raw = inverse_transform_c(c5)
        )
      )
    )
  }
}

# Parameter names (19 total: 12 background + 5 triggering + 2 kernel)
param_names <- c("beta_0_bg", "gamma_raw", "delta_raw",
                "beta_2016_raw", "beta_2017_raw", "beta_2018_raw",
                "beta_2019_raw", "beta_2020_raw", "beta_2021_raw",
                "beta_2022_raw", "beta_2023_raw", "beta_2024_raw",
                "beta_0_trig_raw", "beta_riot_raw", "beta_fatal_raw",
                "beta_student_raw", "beta_labor_raw",
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
  reset_optim_state(sprintf("model4_start%d", s))

  # Run optimization
  start_time <- Sys.time()

  result <- optim(
    par = par0,
    fn = function(params_vec) {
      params_list <- as.list(params_vec)
      names(params_list) <- param_names
      objective_function(params_list, times_vec, marks_data, district_years,
                        verbose = TRUE, model_id = sprintf("model4_start%d", s))
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
  cat(sprintf("  Final: α=%.3f, c=%.4f\n", alpha_final, c_final))
  cat(sprintf("  Half-life equiv: %.1f days\n\n", log(2) / alpha_final))

  if (s < length(starting_points)) {
    cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  }
}

# ====================================================================
# 6. SELECT BEST RESULT
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
cat(sprintf("  c = %.4f days (~%.1f hours offset)\n", best_result$c, best_result$c * 24))
cat("\n")

# ====================================================================
# 7. SAVE RESULTS
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SAVING RESULTS                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Save best model
model4_result <- list(
  model_name = "Model 4: Power-Law Kernel Hawkes",
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

saveRDS(model4_result, "model4_powerlaw.rds")
cat("✓ Saved: model4_powerlaw.rds\n")

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL 4 ESTIMATION COMPLETE                                 ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("FINAL PARAMETERS:\n")
cat(sprintf("  α (power-law exponent) = %.4f\n", best_result$alpha))
cat(sprintf("  c (offset) = %.4f days (%.1f hours)\n", best_result$c, best_result$c * 24))
cat(sprintf("  Log-likelihood = %.2f\n", best_result$loglik))
cat("\n")
