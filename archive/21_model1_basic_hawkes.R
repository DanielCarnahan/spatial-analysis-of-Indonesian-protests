# ============================================================================
# MODEL 1: BASIC HAWKES PROCESS (CONSTANT TRIGGERING, NO MARKS)
# ============================================================================
#
# Specification: λ(t) = μ(d,y) + Σ α·exp(-β·(t-t_i))
#               where μ(d,y) = exp(β₀ + γ·log(pop) + δ·poverty + Σ β_year)
#                     α = exp(β₀_trig) is CONSTANT (no mark dependence)
#
# PURPOSE:
# - Tests for self-excitation while assuming constant triggering
# - Intermediate model between Poisson (Model 0) and full Hawkes (Model 2)
# - Will be compared to both Model 0 and Model 2 via likelihood ratio tests
#
# COMPARISON STRATEGY:
# - Model 0 vs Model 1: Tests for self-excitation (LR test, df=2: β₀_trig, β_decay)
# - Model 1 vs Model 2: Tests for mark-dependent triggering (LR test, df=4: mark effects)
#
# ============================================================================

library(dplyr)
library(Rcpp)

# Compile C++ likelihood function for speed
cat("Compiling C++ likelihood function for basic Hawkes model...\n")
sourceCpp("hawkes_basic_likelihood.cpp")
cat("✓ C++ compiled\n\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   MODEL 1: BASIC HAWKES (CONSTANT TRIGGERING)               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("SPECIFICATION:\n")
cat("  Background: μ(d,y) = exp(β₀ + γ·log(pop) + δ·poverty_decimal + Σ β_year)\n")
cat("  Triggering: α = exp(β₀_trig) [CONSTANT, no mark dependence]\n")
cat("  Kernel: g(t) = exp(-β·t)\n")
cat("  Intensity: λ(t) = μ(d,y) + Σ α·g(t-t_i)\n\n")

cat("HYPOTHESIS:\n")
cat("  → Self-exciting process (protests trigger more protests)\n")
cat("  → Constant triggering effect (no heterogeneity by event characteristics)\n")
cat("  → Baseline for testing mark-dependent effects in Model 2\n\n")

# ====================================================================
# CONFIGURATION
# ====================================================================

TEST_MODE <- FALSE
TEMPORAL_CUTOFF <- 90  # 90 days for protest contagion
CHECKPOINT_EVERY_N_EVALS <- 1

# ====================================================================
# PARAMETER TRANSFORMATION FUNCTIONS
# ====================================================================

# Transform gamma_raw → gamma ∈ [0.001, 5]  (RELAXED from [0.01, 2])
transform_gamma <- function(gamma_raw) {
  0.001 + 4.999 * plogis(gamma_raw)
}
inverse_transform_gamma <- function(gamma) {
  qlogis((gamma - 0.001) / 4.999)
}

# Transform delta_raw → delta ∈ [-20, 20]
transform_delta <- function(delta_raw) {
  -20 + 40 * plogis(delta_raw)
}
inverse_transform_delta <- function(delta) {
  qlogis((delta + 20) / 40)
}

# Transform beta_year_raw → beta_year ∈ [-5, 5]  (RELAXED from [-3, 3])
transform_beta_year <- function(beta_year_raw) {
  -5 + 10 * plogis(beta_year_raw)
}
inverse_transform_beta_year <- function(beta_year) {
  qlogis((beta_year + 5) / 10)
}

# Transform beta_0_trig_raw → beta_0_trig ∈ [-15, 5]  (RELAXED from [-10, 10])
transform_beta_0_trig <- function(beta_0_trig_raw) {
  -15 + 20 * plogis(beta_0_trig_raw)
}
inverse_transform_beta_0_trig <- function(beta_0_trig) {
  qlogis((beta_0_trig + 15) / 20)
}

# Transform decay_raw → decay ∈ [0.01, 50]  (RELAXED from [0.001, 10])
transform_decay <- function(decay_raw) {
  0.01 + 49.99 * plogis(decay_raw)
}
inverse_transform_decay <- function(decay) {
  qlogis((decay - 0.01) / 49.99)
}

# ====================================================================
# CHECKPOINT SYSTEM
# ====================================================================

checkpoint_dir <- "checkpoints_model_basic_hawkes_v2"  # v2 with relaxed bounds
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
cat(sprintf("Poverty data coverage: %d (%.1f%%)\n", n_with_poverty, pct_coverage))

cat("\nPoverty rate summary:\n")
cat(sprintf("  Mean: %.4f (%.2f%%)\n", mean(protests$poverty_decimal, na.rm = TRUE),
            mean(protests$poverty_decimal, na.rm = TRUE)*100))
cat(sprintf("  Range: %.4f to %.4f\n\n",
            min(protests$poverty_decimal, na.rm = TRUE),
            max(protests$poverty_decimal, na.rm = TRUE)))

# Prepare data
marks_data <- protests %>%
  select(
    event_id = event_id_cnty,
    time = days_since_start,
    district,
    year,
    log_pop,
    population,
    poverty_decimal
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

cat("=== IMPLEMENTING BASIC HAWKES MODEL ===\n\n")

model_basic_hawkes_functions <- function() {

  loglik_basic_hawkes <- function(params, times, marks, district_years,
                                  verbose = TRUE, model_id = NULL) {

    .optim_state$n_evals <- .optim_state$n_evals + 1
    eval_num <- .optim_state$n_evals

    # Prepare parameter vector for C++ (14 parameters: 12 background + 2 triggering)
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
      # Triggering parameters (NO MARKS - constant α)
      params[["beta_0_trig_raw"]],
      params[["decay_raw"]]
    )

    # Call C++ function (returns NEGATIVE log-likelihood)
    neg_ll <- hawkes_basic_negloglik_cpp(
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

    if(ll_final > .optim_state$best_ll) {
      .optim_state$best_ll <- ll_final
      .optim_state$best_params <- params
    }

    if(!is.null(model_id) && eval_num %% CHECKPOINT_EVERY_N_EVALS == 0) {
      gamma_display <- transform_gamma(params[["gamma_raw"]])
      alpha_display <- exp(transform_beta_0_trig(params[["beta_0_trig_raw"]]))
      decay_display <- transform_decay(params[["decay_raw"]])

      save_intermediate_checkpoint(model_id, eval_num, params, ll_final)
      if(verbose) {
        cat(sprintf("\n    💾 Checkpoint %d: LL=%.2f (best: %.2f) | γ=%.4f α=%.4f β=%.4f\n",
                   eval_num, ll_final, .optim_state$best_ll,
                   gamma_display, alpha_display, decay_display))
      }
    }

    if(verbose) {
      gamma_display <- transform_gamma(params[["gamma_raw"]])
      alpha_display <- exp(transform_beta_0_trig(params[["beta_0_trig_raw"]]))
      decay_display <- transform_decay(params[["decay_raw"]])

      cat(sprintf("    LL: %.2f | γ=%.4f α=%.4f β=%.4f [eval %d]\n",
                 ll_final, gamma_display, alpha_display, decay_display, eval_num))
      flush.console()
    }

    return(ll_final)
  }

  return(list(loglik = loglik_basic_hawkes))
}

hawkes_basic_funcs <- model_basic_hawkes_functions()

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

fit_model_basic_hawkes <- function(times, marks, district_years,
                                   model_name = "Model 1: Basic Hawkes (Constant Triggering)",
                                   model_id = "MODEL1_BASIC_HAWKES",
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
      gamma_resume <- transform_gamma(checkpoint$data$best_params$gamma_raw)
      alpha_resume <- exp(transform_beta_0_trig(checkpoint$data$best_params$beta_0_trig_raw))
      cat(sprintf("⏩ Resuming from intermediate checkpoint (eval %d, LL=%.2f, γ=%.4f, α=%.4f)\n",
                 checkpoint$data$eval_num, checkpoint$data$best_ll,
                 gamma_resume, alpha_resume))

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

    # 1. Data-driven
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
        beta_2016_raw = 0, beta_2017_raw = 0, beta_2018_raw = 0,
        beta_2019_raw = 0, beta_2020_raw = 0, beta_2021_raw = 0,
        beta_2022_raw = 0, beta_2023_raw = 0, beta_2024_raw = 0,
        beta_0_trig_raw = beta_0_trig1_raw,
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
          decay_raw = decay3_raw
        )
      )
    }

    if(use_multistart && n_starts >= 4) {
      # 4. Fast decay
      gamma4 <- 0.50; delta4 <- -0.2
      gamma4_raw <- inverse_transform_gamma(gamma4)
      delta4_raw <- inverse_transform_delta(delta4)
      beta_0_trig4 <- log(0.2); decay4 <- 0.3
      beta_0_trig4_raw <- inverse_transform_beta_0_trig(beta_0_trig4)
      decay4_raw <- inverse_transform_decay(decay4)
      beta_0_4 <- log(observed_rate) - gamma4 * mean_log_pop - delta4 * mean_poverty
      starting_points[[4]] <- list(
        name = "Fast-decay",
        params = list(
          beta_0_bg = beta_0_4, gamma_raw = gamma4_raw, delta_raw = delta4_raw,
          beta_2016_raw = 0.05, beta_2017_raw = -0.05, beta_2018_raw = 0.1,
          beta_2019_raw = -0.1, beta_2020_raw = 0, beta_2021_raw = 0.05,
          beta_2022_raw = -0.05, beta_2023_raw = 0.1, beta_2024_raw = -0.1,
          beta_0_trig_raw = beta_0_trig4_raw,
          decay_raw = decay4_raw
        )
      )
    }

    if(use_multistart && n_starts >= 5) {
      # 5. Slow decay
      gamma5 <- 0.40; delta5 <- 0.5
      gamma5_raw <- inverse_transform_gamma(gamma5)
      delta5_raw <- inverse_transform_delta(delta5)
      beta_0_trig5 <- log(0.05); decay5 <- 0.1
      beta_0_trig5_raw <- inverse_transform_beta_0_trig(beta_0_trig5)
      decay5_raw <- inverse_transform_decay(decay5)
      beta_0_5 <- log(observed_rate) - gamma5 * mean_log_pop - delta5 * mean_poverty
      starting_points[[5]] <- list(
        name = "Slow-decay",
        params = list(
          beta_0_bg = beta_0_5, gamma_raw = gamma5_raw, delta_raw = delta5_raw,
          beta_2016_raw = -0.2, beta_2017_raw = -0.1, beta_2018_raw = 0,
          beta_2019_raw = 0.1, beta_2020_raw = 0.2, beta_2021_raw = 0.1,
          beta_2022_raw = 0, beta_2023_raw = -0.1, beta_2024_raw = -0.2,
          beta_0_trig_raw = beta_0_trig5_raw,
          decay_raw = decay5_raw
        )
      )
    }
  }

  # Parameter names (14 total: 12 background + 2 triggering)
  param_names <- c("beta_0_bg", "gamma_raw", "delta_raw",
                   "beta_2016_raw", "beta_2017_raw", "beta_2018_raw",
                   "beta_2019_raw", "beta_2020_raw", "beta_2021_raw",
                   "beta_2022_raw", "beta_2023_raw", "beta_2024_raw",
                   "beta_0_trig_raw", "decay_raw")

  lower_bounds <- c(beta_0_bg = -15, gamma_raw = -10, delta_raw = -10,
                    beta_2016_raw = -10, beta_2017_raw = -10, beta_2018_raw = -10,
                    beta_2019_raw = -10, beta_2020_raw = -10, beta_2021_raw = -10,
                    beta_2022_raw = -10, beta_2023_raw = -10, beta_2024_raw = -10,
                    beta_0_trig_raw = -10, decay_raw = -10)

  upper_bounds <- c(beta_0_bg = 5, gamma_raw = 10, delta_raw = 10,
                    beta_2016_raw = 10, beta_2017_raw = 10, beta_2018_raw = 10,
                    beta_2019_raw = 10, beta_2020_raw = 10, beta_2021_raw = 10,
                    beta_2022_raw = 10, beta_2023_raw = 10, beta_2024_raw = 10,
                    beta_0_trig_raw = 10, decay_raw = 10)

  obj_fun_silent <- function(par) {
    params_full <- as.list(par)
    ll <- hawkes_basic_funcs$loglik(params_full, times, marks, district_years,
                                    verbose = FALSE, model_id = NULL)
    return(-ll)
  }

  obj_fun_verbose <- function(par) {
    params_full <- as.list(par)
    ll <- hawkes_basic_funcs$loglik(params_full, times, marks, district_years,
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

    gamma_start <- transform_gamma(sp$params$gamma_raw)
    alpha_start <- exp(transform_beta_0_trig(sp$params$beta_0_trig_raw))
    decay_start <- transform_decay(sp$params$decay_raw)
    cat(sprintf("Starting values: γ=%.4f, α=%.4f, β=%.4f\n\n",
               gamma_start, alpha_start, decay_start))

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

    gamma_final <- transform_gamma(params_final$gamma_raw)
    alpha_final <- exp(transform_beta_0_trig(params_final$beta_0_trig_raw))
    decay_final <- transform_decay(params_final$decay_raw)

    cat(sprintf("\n✓ Optimization complete\n"))
    cat(sprintf("  Log-likelihood: %.2f\n", loglik_final))
    cat(sprintf("  Convergence: %d\n", optim_result$convergence))
    cat(sprintf("  Runtime: %.1f minutes\n", runtime))
    cat(sprintf("  Final: γ=%.4f, α=%.4f, β=%.4f\n\n",
               gamma_final, alpha_final, decay_final))

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
    cat(sprintf("  %d. %-20s  LL: %10.2f  Conv: %d%s\n",
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
    n_params = 14,
    n_events = nrow(marks),
    n_district_years = nrow(district_years),
    aic = -2 * best_result$loglik + 2 * 14,
    bic = -2 * best_result$loglik + 14 * log(nrow(marks))
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

fit_model1 <- fit_model_basic_hawkes(
  times = times_sample,
  marks = marks_sample,
  district_years = district_years_sample,
  model_name = "Model 1: Basic Hawkes (Constant Triggering)",
  model_id = "MODEL1_BASIC_HAWKES",
  use_multistart = TRUE,
  n_starts = 5
)

# ====================================================================
# 7. SAVE RESULTS
# ====================================================================

cat("\n=== SAVING RESULTS ===\n")

saveRDS(fit_model1, "model1_basic_hawkes_v2.rds")  # v2 with relaxed bounds
cat("✓ Saved: model1_basic_hawkes_v2.rds\n")

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FINAL RESULTS: MODEL 1 (BASIC HAWKES)                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("Model: %s\n", fit_model1$model_name))
cat(sprintf("Log-likelihood: %.2f\n", fit_model1$loglik))
cat(sprintf("AIC: %.2f\n", fit_model1$aic))
cat(sprintf("BIC: %.2f\n", fit_model1$bic))
cat(sprintf("Parameters: %d\n", fit_model1$n_params))
cat(sprintf("Convergence: %d\n", fit_model1$convergence))
cat(sprintf("Runtime: %.1f minutes\n\n", fit_model1$runtime))

gamma_est <- transform_gamma(fit_model1$params$gamma_raw)
delta_est <- transform_delta(fit_model1$params$delta_raw)
beta_0_trig_est <- transform_beta_0_trig(fit_model1$params$beta_0_trig_raw)
alpha_est <- exp(beta_0_trig_est)
decay_est <- transform_decay(fit_model1$params$decay_raw)

cat("Background parameters:\n")
cat(sprintf("  β₀ (baseline): %.4f\n", fit_model1$params$beta_0_bg))
cat(sprintf("  γ (population): %.4f\n", gamma_est))
cat(sprintf("  δ (poverty): %.4f\n\n", delta_est))

cat("Triggering parameters:\n")
cat(sprintf("  β₀_trig: %.4f → α = %.4f\n", beta_0_trig_est, alpha_est))
cat(sprintf("  β_decay: %.4f (half-life: %.1f days)\n\n",
           decay_est, log(2)/decay_est))

cat("Year effects:\n")
for(year in 2016:2024) {
  param_name <- sprintf("beta_%d_raw", year)
  beta_year <- transform_beta_year(fit_model1$params[[param_name]])
  cat(sprintf("  β_%d: %.4f\n", year, beta_year))
}

cat("\n✓ Model 1 fitting complete!\n")
cat("  Next: Run model comparison to test hypotheses\n\n")
