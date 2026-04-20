#!/usr/bin/env Rscript
# ==============================================================================
# MODEL 3: MIXED KERNEL HAWKES PROCESS
# ==============================================================================
#
# Specification: λ(t) = μ(d,y) + Σ α(marks) · [p·exp(-β_fast·Δt) + (1-p)·exp(-β_slow·Δt)]
#
# NEW PARAMETERS (vs Model 2):
#   p ∈ [0,1]: Weight on fast component
#   β_fast ∈ [0.01, 5]: Fast decay rate (captures short-term contagion)
#   β_slow ∈ [0.001, 0.1]: Slow decay rate (captures long-term trends)
#
# CONSTRAINT: β_fast > β_slow (fast must decay faster than slow)
#
# PURPOSE:
#   - Resolve Model 2's slow decay issue (β = 0.0014, half-life = 477 days)
#   - Capture both short-term clustering (days) and long-term correlation (weeks/months)
#   - Test if two-timescale model fits data better than single exponential
#
# COMPARISON:
#   - Model 2 vs Model 3: LR test (df=2, adds p and second decay rate)
#   - Expected: LR >> 0, p < 0.001, strong rejection of single exponential
#
# ==============================================================================

library(dplyr)
library(Rcpp)

# Compile C++ likelihood function
cat("Compiling C++ likelihood function for mixed kernel Hawkes model...\n")
sourceCpp("hawkes_mixed_kernel_likelihood.cpp")
cat("✓ C++ compiled\n\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   MODEL 3: MIXED KERNEL HAWKES                               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("SPECIFICATION:\n")
cat("  Background: μ(d,y) = exp(β₀ + γ·log(pop) + δ·poverty + Σ β_year)\n")
cat("  Triggering: α(marks) = exp(β₀_trig + β_riot + β_fatal + β_student + β_labor)\n")
cat("  Kernel: g(Δt) = p·exp(-β_fast·Δt) + (1-p)·exp(-β_slow·Δt)  [MIXED!]\n")
cat("  Intensity: λ(t) = μ(d,y) + Σ α(marks_i) · g(t - t_i)\n\n")

cat("HYPOTHESIS:\n")
cat("  → Two distinct timescales: fast (days) + slow (weeks/months)\n")
cat("  → Most contagion is short-term (p ≈ 0.7-0.9)\n")
cat("  → Resolves single exponential's misspecification\n\n")

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

# NEW: Mixed kernel parameters
transform_p <- function(p_raw) {
  plogis(p_raw)  # [0, 1]
}
inverse_transform_p <- function(p) {
  qlogis(p)
}

transform_beta_fast <- function(beta_fast_raw) {
  0.1 + 9.9 * plogis(beta_fast_raw)  # [0.1, 10]
}
inverse_transform_beta_fast <- function(beta_fast) {
  qlogis((beta_fast - 0.1) / 9.9)
}

transform_beta_slow <- function(beta_slow_raw) {
  0.01 + 0.99 * plogis(beta_slow_raw)  # [0.01, 1]
}
inverse_transform_beta_slow <- function(beta_slow) {
  qlogis((beta_slow - 0.01) / 0.99)
}

# ====================================================================
# CHECKPOINT SYSTEM
# ====================================================================

checkpoint_dir <- "checkpoints_model3_mixed_kernel"
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

cat(sprintf("✓ Using full dataset: %d events across %d years\n",
            nrow(protests), length(unique(protests$year))))

if(TEST_MODE) {
  set.seed(42)
  n_sample <- min(1000, nrow(protests))
  sample_idx <- sample(nrow(protests), n_sample)
  protests <- protests[sample_idx, ]
  cat(sprintf("⚠ TEST MODE: Using random sample of %d events\n", n_sample))
}

n_with_poverty <- sum(!is.na(protests$poverty_decimal))
pct_coverage <- 100 * n_with_poverty / nrow(protests)
cat(sprintf("Poverty data coverage: %d (%.1f%%)\n\n", n_with_poverty, pct_coverage))

# Prepare data - create mark indicators first
marks_data <- protests %>%
  mutate(
    is_riot = as.integer(event_type == "Riots"),
    is_fatal = as.integer(has_fatalities > 0),
    is_student = as.integer(grepl("Student", assoc_actor_1, ignore.case = TRUE)),
    is_labor = as.integer(grepl("Labor|Worker", assoc_actor_1, ignore.case = TRUE))
  ) %>%
  select(
    event_id = event_id_cnty,
    time = days_since_start,
    district,
    year,
    log_pop,
    population,
    poverty_decimal,
    is_riot,
    is_fatal,
    is_student,
    is_labor
  ) %>%
  arrange(time) %>%
  filter(!is.na(poverty_decimal))

cat(sprintf("Using complete cases: %d events (%.1f%% of total)\n\n",
            nrow(marks_data),
            100 * nrow(marks_data) / nrow(protests)))

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

cat("=== IMPLEMENTING MIXED KERNEL HAWKES MODEL ===\n\n")

model_mixed_hawkes_functions <- function() {

  loglik_mixed_hawkes <- function(params, times, marks, district_years,
                                  verbose = TRUE, model_id = NULL) {

    .optim_state$n_evals <- .optim_state$n_evals + 1
    eval_num <- .optim_state$n_evals

    # Prepare parameter vector for C++ (20 parameters: 12 background + 8 triggering)
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
      # Triggering parameters (8 total: base + 4 marks + p + beta_fast + beta_slow)
      params[["beta_0_trig_raw"]],
      params[["beta_riot_raw"]],
      params[["beta_fatal_raw"]],
      params[["beta_student_raw"]],
      params[["beta_labor_raw"]],
      params[["p_raw"]],
      params[["beta_fast_raw"]],
      params[["beta_slow_raw"]]
    )

    # Call C++ function (returns NEGATIVE log-likelihood)
    neg_ll <- hawkes_mixed_negloglik_cpp(
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

    if(ll_final > .optim_state$best_ll) {
      .optim_state$best_ll <- ll_final
      .optim_state$best_params <- params
    }

    if(!is.null(model_id) && eval_num %% CHECKPOINT_EVERY_N_EVALS == 0) {
      p_display <- transform_p(params[["p_raw"]])
      beta_fast_display <- transform_beta_fast(params[["beta_fast_raw"]])
      beta_slow_display <- transform_beta_slow(params[["beta_slow_raw"]])

      save_intermediate_checkpoint(model_id, eval_num, params, ll_final)
      if(verbose) {
        cat(sprintf("\n    💾 Checkpoint %d: LL=%.2f (best: %.2f) | p=%.3f β_fast=%.4f β_slow=%.4f\n",
                   eval_num, ll_final, .optim_state$best_ll,
                   p_display, beta_fast_display, beta_slow_display))
      }
    }

    if(verbose) {
      p_display <- transform_p(params[["p_raw"]])
      beta_fast_display <- transform_beta_fast(params[["beta_fast_raw"]])
      beta_slow_display <- transform_beta_slow(params[["beta_slow_raw"]])

      cat(sprintf("    LL: %.2f | p=%.3f β_fast=%.4f β_slow=%.4f [eval %d]\n",
                 ll_final, p_display, beta_fast_display, beta_slow_display, eval_num))
      flush.console()
    }

    return(ll_final)
  }

  return(list(loglik = loglik_mixed_hawkes))
}

hawkes_mixed_funcs <- model_mixed_hawkes_functions()

# ====================================================================
# 4. PREPARE DATA
# ====================================================================

cat("=== PREPARING DATA FOR ESTIMATION ===\n")

times_sample <- marks_data$time
marks_sample <- marks_data
district_years_sample <- district_years

cat(sprintf("Sample size: %d events\n", length(times_sample)))
cat(sprintf("District-years: %d\n", nrow(district_years_sample)))
cat("✓ Ready for optimization\n\n")

# ====================================================================
# 5. FITTING FUNCTION
# ====================================================================

fit_model_mixed_hawkes <- function(times, marks, district_years,
                                   model_name = "Model 3: Mixed Kernel Hawkes",
                                   model_id = "MODEL3_MIXED_KERNEL",
                                   use_multistart = TRUE,
                                   n_starts = 5) {

  cat(sprintf("\n--- Fitting %s ---\n", model_name))
  if(use_multistart) {
    cat(sprintf("Using MULTI-START optimization with %d starting points\n", n_starts))
  }

  checkpoint <- find_latest_checkpoint(model_id)

  if(!is.null(checkpoint)) {
    if(checkpoint$type == "final") {
      cat("✓ Model already completed (loading from final checkpoint)\n")
      return(checkpoint$data)
    } else {
      p_resume <- transform_p(checkpoint$data$best_params$p_raw)
      beta_fast_resume <- transform_beta_fast(checkpoint$data$best_params$beta_fast_raw)
      cat(sprintf("⏩ Resuming from intermediate checkpoint (eval %d, LL=%.2f, p=%.3f, β_fast=%.4f)\n",
                 checkpoint$data$eval_num, checkpoint$data$best_ll,
                 p_resume, beta_fast_resume))

      reset_optim_state(model_name)
      .optim_state$n_evals <- checkpoint$data$eval_num
      .optim_state$best_ll <- checkpoint$data$best_ll
      .optim_state$best_params <- checkpoint$data$best_params

      use_multistart <- FALSE
      starting_points <- list(list(name = "Checkpoint", params = checkpoint$data$best_params))
    }
  } else {
    cat("Starting fresh optimization...\n")

    mean_log_pop <- mean(marks$log_pop)
    mean_poverty <- mean(marks$poverty_decimal, na.rm = TRUE)
    total_events <- nrow(marks)
    total_exposure <- nrow(district_years) * (max(times) - min(times))
    observed_rate <- total_events / total_exposure

    starting_points <- list()

    # Load Model 2 estimates to use as baseline
    model2 <- readRDS("model_poverty.rds")

    # 1. Data-driven: p=0.7, β_fast=1.0 (half-life=0.69 days), β_slow=0.05 (half-life=14 days)
    p1 <- 0.7; beta_fast1 <- 1.0; beta_slow1 <- 0.05
    p1_raw <- inverse_transform_p(p1)
    beta_fast1_raw <- inverse_transform_beta_fast(beta_fast1)
    beta_slow1_raw <- inverse_transform_beta_slow(beta_slow1)

    starting_points[[1]] <- list(
      name = "Data-driven (p=0.7)",
      params = c(
        model2$params[1:17],  # Use Model 2's background + mark estimates
        list(
          p_raw = p1_raw,
          beta_fast_raw = beta_fast1_raw,
          beta_slow_raw = beta_slow1_raw
        )
      )
    )

    if(use_multistart && n_starts >= 2) {
      # 2. Fast-dominant: p=0.9, β_fast=2.0 (half-life=0.35 days), β_slow=0.03 (half-life=23 days)
      p2 <- 0.9; beta_fast2 <- 2.0; beta_slow2 <- 0.03
      starting_points[[2]] <- list(
        name = "Fast-dominant (p=0.9)",
        params = c(
          model2$params[1:17],
          list(
            p_raw = inverse_transform_p(p2),
            beta_fast_raw = inverse_transform_beta_fast(beta_fast2),
            beta_slow_raw = inverse_transform_beta_slow(beta_slow2)
          )
        )
      )
    }

    if(use_multistart && n_starts >= 3) {
      # 3. Balanced: p=0.5, β_fast=1.5 (half-life=0.46 days), β_slow=0.1 (half-life=7 days)
      p3 <- 0.5; beta_fast3 <- 1.5; beta_slow3 <- 0.1
      starting_points[[3]] <- list(
        name = "Balanced (p=0.5)",
        params = c(
          model2$params[1:17],
          list(
            p_raw = inverse_transform_p(p3),
            beta_fast_raw = inverse_transform_beta_fast(beta_fast3),
            beta_slow_raw = inverse_transform_beta_slow(beta_slow3)
          )
        )
      )
    }

    if(use_multistart && n_starts >= 4) {
      # 4. Slow-dominant: p=0.3, β_fast=3.0 (half-life=0.23 days), β_slow=0.2 (half-life=3.5 days)
      p4 <- 0.3; beta_fast4 <- 3.0; beta_slow4 <- 0.2
      starting_points[[4]] <- list(
        name = "Slow-dominant (p=0.3)",
        params = c(
          model2$params[1:17],
          list(
            p_raw = inverse_transform_p(p4),
            beta_fast_raw = inverse_transform_beta_fast(beta_fast4),
            beta_slow_raw = inverse_transform_beta_slow(beta_slow4)
          )
        )
      )
    }

    if(use_multistart && n_starts >= 5) {
      # 5. Conservative: p=0.6, β_fast=0.7 (half-life=1 day), β_slow=0.07 (half-life=10 days)
      p5 <- 0.6; beta_fast5 <- 0.7; beta_slow5 <- 0.07
      starting_points[[5]] <- list(
        name = "Conservative (p=0.6)",
        params = c(
          model2$params[1:17],
          list(
            p_raw = inverse_transform_p(p5),
            beta_fast_raw = inverse_transform_beta_fast(beta_fast5),
            beta_slow_raw = inverse_transform_beta_slow(beta_slow5)
          )
        )
      )
    }
  }

  # Parameter names (20 total: 12 background + 8 triggering)
  param_names <- c("beta_0_bg", "gamma_raw", "delta_raw",
                   "beta_2016_raw", "beta_2017_raw", "beta_2018_raw",
                   "beta_2019_raw", "beta_2020_raw", "beta_2021_raw",
                   "beta_2022_raw", "beta_2023_raw", "beta_2024_raw",
                   "beta_0_trig_raw", "beta_riot_raw", "beta_fatal_raw",
                   "beta_student_raw", "beta_labor_raw",
                   "p_raw", "beta_fast_raw", "beta_slow_raw")

  lower_bounds <- c(beta_0_bg = -15, gamma_raw = -10, delta_raw = -10,
                    beta_2016_raw = -10, beta_2017_raw = -10, beta_2018_raw = -10,
                    beta_2019_raw = -10, beta_2020_raw = -10, beta_2021_raw = -10,
                    beta_2022_raw = -10, beta_2023_raw = -10, beta_2024_raw = -10,
                    beta_0_trig_raw = -10, beta_riot_raw = -10, beta_fatal_raw = -10,
                    beta_student_raw = -10, beta_labor_raw = -10,
                    p_raw = -10, beta_fast_raw = -10, beta_slow_raw = -10)

  upper_bounds <- c(beta_0_bg = 5, gamma_raw = 10, delta_raw = 10,
                    beta_2016_raw = 10, beta_2017_raw = 10, beta_2018_raw = 10,
                    beta_2019_raw = 10, beta_2020_raw = 10, beta_2021_raw = 10,
                    beta_2022_raw = 10, beta_2023_raw = 10, beta_2024_raw = 10,
                    beta_0_trig_raw = 10, beta_riot_raw = 10, beta_fatal_raw = 10,
                    beta_student_raw = 10, beta_labor_raw = 10,
                    p_raw = 10, beta_fast_raw = 10, beta_slow_raw = 10)

  obj_fun_silent <- function(par) {
    params_full <- as.list(par)
    names(params_full) <- param_names
    ll <- hawkes_mixed_funcs$loglik(params_full, times, marks, district_years,
                                    verbose = FALSE, model_id = NULL)
    return(-ll)
  }

  obj_fun_verbose <- function(par) {
    params_full <- as.list(par)
    names(params_full) <- param_names
    ll <- hawkes_mixed_funcs$loglik(params_full, times, marks, district_years,
                                    verbose = TRUE, model_id = model_id)
    return(-ll)
  }

  results_list <- list()
  start_time_total <- Sys.time()

  for(s in seq_along(starting_points)) {
    sp <- starting_points[[s]]
    cat(sprintf("\n╔══════════════════════════════════════════════════════════════╗\n"))
    cat(sprintf("║  Starting Point %d/%d: %-42s ║\n", s, length(starting_points), sp$name))
    cat(sprintf("╚══════════════════════════════════════════════════════════════╝\n\n"))

    par0 <- unlist(sp$params)[param_names]

    p_start <- transform_p(sp$params$p_raw)
    beta_fast_start <- transform_beta_fast(sp$params$beta_fast_raw)
    beta_slow_start <- transform_beta_slow(sp$params$beta_slow_raw)
    cat(sprintf("Starting values: p=%.3f, β_fast=%.4f, β_slow=%.4f\n\n",
               p_start, beta_fast_start, beta_slow_start))

    reset_optim_state(sprintf("%s (start %d)", model_name, s))
    start_time <- Sys.time()

    if(s == 1 || !use_multistart) {
      optim_result <- optim(
        par = par0,
        fn = obj_fun_verbose,
        method = "L-BFGS-B",
        lower = lower_bounds[param_names],
        upper = upper_bounds[param_names],
        control = list(
          maxit = if(TEST_MODE) 100 else 1000,
          factr = 1e7,
          trace = 1,
          REPORT = 10
        )
      )
    } else {
      optim_result <- optim(
        par = par0,
        fn = obj_fun_silent,
        method = "L-BFGS-B",
        lower = lower_bounds[param_names],
        upper = upper_bounds[param_names],
        control = list(
          maxit = if(TEST_MODE) 100 else 1000,
          factr = 1e7
        )
      )
    }

    runtime <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

    params_final <- as.list(optim_result$par)
    names(params_final) <- param_names
    loglik_final <- -optim_result$value

    p_final <- transform_p(params_final$p_raw)
    beta_fast_final <- transform_beta_fast(params_final$beta_fast_raw)
    beta_slow_final <- transform_beta_slow(params_final$beta_slow_raw)

    cat(sprintf("\n✓ Optimization complete\n"))
    cat(sprintf("  Log-likelihood: %.2f\n", loglik_final))
    cat(sprintf("  Convergence: %d\n", optim_result$convergence))
    cat(sprintf("  Runtime: %.1f minutes\n", runtime))
    cat(sprintf("  Final: p=%.3f, β_fast=%.4f, β_slow=%.4f\n",
               p_final, beta_fast_final, beta_slow_final))
    cat(sprintf("  Half-lives: fast=%.1f days, slow=%.1f days\n\n",
               log(2)/beta_fast_final, log(2)/beta_slow_final))

    results_list[[s]] <- list(
      starting_point = sp$name,
      params = params_final,
      loglik = loglik_final,
      convergence = optim_result$convergence,
      runtime = runtime,
      iterations = optim_result$counts[1]
    )
  }

  logliks <- sapply(results_list, function(x) x$loglik)
  best_idx <- which.max(logliks)
  best_result <- results_list[[best_idx]]

  total_runtime <- as.numeric(difftime(Sys.time(), start_time_total, units = "mins"))

  cat("\n╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  OPTIMIZATION SUMMARY                                        ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")

  for(s in seq_along(results_list)) {
    r <- results_list[[s]]
    marker <- if(s == best_idx) " ★ BEST" else ""
    cat(sprintf("  %d. %-30s  LL: %10.2f  Conv: %d%s\n",
               s, r$starting_point, r$loglik, r$convergence, marker))
  }

  cat(sprintf("\nTotal runtime: %.1f minutes\n\n", total_runtime))

  model_fit <- list(
    model_name = model_name,
    model_id = model_id,
    params = best_result$params,
    loglik = best_result$loglik,
    convergence = best_result$convergence,
    starting_point = best_result$starting_point,
    runtime = total_runtime,
    all_results = results_list,
    n_params = 20,
    n_events = nrow(marks),
    n_district_years = nrow(district_years),
    aic = -2 * best_result$loglik + 2 * 20,
    bic = -2 * best_result$loglik + 20 * log(nrow(marks))
  )

  save_final_checkpoint(model_fit, model_id)

  return(model_fit)
}

# ====================================================================
# 6. FIT MODEL
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  STARTING OPTIMIZATION                                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")

fit_model3 <- fit_model_mixed_hawkes(
  times = times_sample,
  marks = marks_sample,
  district_years = district_years_sample,
  model_name = "Model 3: Mixed Kernel Hawkes",
  model_id = "MODEL3_MIXED_KERNEL",
  use_multistart = TRUE,
  n_starts = 5
)

# ====================================================================
# 7. SAVE RESULTS
# ====================================================================

cat("\n=== SAVING RESULTS ===\n")

saveRDS(fit_model3, "model3_mixed_kernel.rds")
cat("✓ Saved: model3_mixed_kernel.rds\n")

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FINAL RESULTS: MODEL 3 (MIXED KERNEL)                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("Model: %s\n", fit_model3$model_name))
cat(sprintf("Log-likelihood: %.2f\n", fit_model3$loglik))
cat(sprintf("AIC: %.2f\n", fit_model3$aic))
cat(sprintf("BIC: %.2f\n", fit_model3$bic))
cat(sprintf("Parameters: %d\n", fit_model3$n_params))
cat(sprintf("Convergence: %d\n", fit_model3$convergence))
cat(sprintf("Runtime: %.1f minutes\n\n", fit_model3$runtime))

p_est <- transform_p(fit_model3$params$p_raw)
beta_fast_est <- transform_beta_fast(fit_model3$params$beta_fast_raw)
beta_slow_est <- transform_beta_slow(fit_model3$params$beta_slow_raw)

cat("Mixed Kernel Parameters:\n")
cat(sprintf("  p (weight on fast): %.3f\n", p_est))
cat(sprintf("  β_fast: %.4f (half-life: %.1f days)\n", beta_fast_est, log(2)/beta_fast_est))
cat(sprintf("  β_slow: %.4f (half-life: %.1f days)\n\n", beta_slow_est, log(2)/beta_slow_est))

cat("Interpretation:\n")
cat(sprintf("  %.1f%% of triggering follows fast decay (%.1f day half-life)\n",
           100*p_est, log(2)/beta_fast_est))
cat(sprintf("  %.1f%% follows slow decay (%.1f day half-life)\n\n",
           100*(1-p_est), log(2)/beta_slow_est))

cat("✓ Model 3 fitting complete!\n")
cat("  Next: Compare to Model 2, create visualizations\n\n")
