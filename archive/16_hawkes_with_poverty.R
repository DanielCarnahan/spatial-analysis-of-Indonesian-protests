# ============================================================================
# MODEL WITH POVERTY RATE: TESTING ECONOMIC DETERMINANTS
# ============================================================================
#
# Specification: μ(district, year) = exp(β₀ + γ·log(pop) + δ·poverty_rate + Σ β_year)
#
# RESEARCH QUESTION:
# - Does poverty rate explain the sub-proportional population effect?
# - Expected: γ should increase toward 1.0 if poverty explains the pattern
# - Larger districts have LOWER poverty (cor = -0.480)
#
# COMPARISON TO BASELINE:
# - Baseline Model B: γ ≈ 0.35 (10× pop → 2.2× more protests)
# - With poverty: expect γ ≈ 0.6-1.0 if poverty explains the discrepancy
#
# ============================================================================

library(dplyr)
library(Rcpp)

# Compile C++ likelihood function for speed
cat("Compiling C++ likelihood function...\n")
sourceCpp("hawkes_likelihood.cpp")
cat("✓ C++ compiled\n\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   HAWKES MODEL WITH POVERTY RATE COVARIATE                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("SPECIFICATION:\n")
cat("  Background: μ(d,y) = exp(β₀ + γ·log(pop) + δ·poverty_decimal + Σ β_year)\n")
cat("             (poverty_decimal = poverty_rate/100, range 0.017-0.437)\n")
cat("  Triggering: α(marks) = exp(β₀_trig + β_violence·violent + β_state·intervention)\n")
cat("  Kernel: g(t) = exp(-β·t)\n")
cat("  Intensity: λ(t) = μ(d,y) + Σ α(marks_i)·g(t-t_i)\n\n")

cat("HYPOTHESIS:\n")
cat("  → Larger districts have lower poverty (cor = -0.480)\n")
cat("  → If poverty explains sub-proportional effect, γ should increase\n")
cat("  → Baseline γ ≈ 0.35, with poverty expect γ ≈ 0.6-1.0\n\n")

# ====================================================================
# CONFIGURATION
# ====================================================================

TEST_MODE <- FALSE  # If TRUE, run on 1000-event sample for testing
TEMPORAL_CUTOFF <- 90  # 90 days (more realistic for protest contagion)
CHECKPOINT_EVERY_N_EVALS <- 1  # Checkpoint every evaluation for better tracking

# ====================================================================
# PARAMETER TRANSFORMATION FUNCTIONS
# ====================================================================
#
# To avoid boundary issues with L-BFGS-B, we reparameterize bounded
# parameters using sigmoid transformations:
#
# gamma ∈ [0.01, 2]  → gamma_raw ∈ (-∞, +∞)
# delta ∈ [-20, 20]  → delta_raw ∈ (-∞, +∞)
#
# This allows standard gradient-based optimization without boundaries.
# ====================================================================

# Transform gamma_raw (unconstrained) → gamma ∈ [0.01, 2]
transform_gamma <- function(gamma_raw) {
  0.01 + 1.99 * plogis(gamma_raw)
}

# Inverse: gamma ∈ [0.01, 2] → gamma_raw
inverse_transform_gamma <- function(gamma) {
  qlogis((gamma - 0.01) / 1.99)
}

# Transform delta_raw (unconstrained) → delta ∈ [-20, 20]
transform_delta <- function(delta_raw) {
  -20 + 40 * plogis(delta_raw)
}

# Inverse: delta ∈ [-20, 20] → delta_raw
inverse_transform_delta <- function(delta) {
  qlogis((delta + 20) / 40)
}

# Transform beta_0_trig_raw (unconstrained) → beta_0_trig ∈ [-10, 10]
transform_beta_0_trig <- function(beta_0_trig_raw) {
  -10 + 20 * plogis(beta_0_trig_raw)
}

# Inverse: beta_0_trig ∈ [-10, 10] → beta_0_trig_raw
inverse_transform_beta_0_trig <- function(beta_0_trig) {
  qlogis((beta_0_trig + 10) / 20)
}

# Transform beta_violence_raw (unconstrained) → beta_violence ∈ [-5, 5]
transform_beta_violence <- function(beta_violence_raw) {
  -5 + 10 * plogis(beta_violence_raw)
}

# Inverse: beta_violence ∈ [-5, 5] → beta_violence_raw
inverse_transform_beta_violence <- function(beta_violence) {
  qlogis((beta_violence + 5) / 10)
}

# Transform beta_state_raw (unconstrained) → beta_state ∈ [-5, 5]
transform_beta_state <- function(beta_state_raw) {
  -5 + 10 * plogis(beta_state_raw)
}

# Inverse: beta_state ∈ [-5, 5] → beta_state_raw
inverse_transform_beta_state <- function(beta_state) {
  qlogis((beta_state + 5) / 10)
}

# Transform decay_raw (unconstrained) → decay ∈ [0.001, 10]
transform_decay <- function(decay_raw) {
  0.001 + 9.999 * plogis(decay_raw)
}

# Inverse: decay ∈ [0.001, 10] → decay_raw
inverse_transform_decay <- function(decay) {
  qlogis((decay - 0.001) / 9.999)
}

# Transform beta_year_raw (unconstrained) → beta_year ∈ [-3, 3]
transform_beta_year <- function(beta_year_raw) {
  -3 + 6 * plogis(beta_year_raw)
}

# Inverse: beta_year ∈ [-3, 3] → beta_year_raw
inverse_transform_beta_year <- function(beta_year) {
  qlogis((beta_year + 3) / 6)
}

# Transform beta_riot_raw (unconstrained) → beta_riot ∈ [-5, 5]
transform_beta_riot <- function(beta_riot_raw) {
  -5 + 10 * plogis(beta_riot_raw)
}

# Inverse: beta_riot ∈ [-5, 5] → beta_riot_raw
inverse_transform_beta_riot <- function(beta_riot) {
  qlogis((beta_riot + 5) / 10)
}

# Transform beta_fatal_raw (unconstrained) → beta_fatal ∈ [-5, 5]
transform_beta_fatal <- function(beta_fatal_raw) {
  -5 + 10 * plogis(beta_fatal_raw)
}

# Inverse: beta_fatal ∈ [-5, 5] → beta_fatal_raw
inverse_transform_beta_fatal <- function(beta_fatal) {
  qlogis((beta_fatal + 5) / 10)
}

# Transform beta_student_raw (unconstrained) → beta_student ∈ [-5, 5]
transform_beta_student <- function(beta_student_raw) {
  -5 + 10 * plogis(beta_student_raw)
}

# Inverse: beta_student ∈ [-5, 5] → beta_student_raw
inverse_transform_beta_student <- function(beta_student) {
  qlogis((beta_student + 5) / 10)
}

# Transform beta_labor_raw (unconstrained) → beta_labor ∈ [-5, 5]
transform_beta_labor <- function(beta_labor_raw) {
  -5 + 10 * plogis(beta_labor_raw)
}

# Inverse: beta_labor ∈ [-5, 5] → beta_labor_raw
inverse_transform_beta_labor <- function(beta_labor) {
  qlogis((beta_labor + 5) / 10)
}

# Test transformations (run once to verify)
test_transforms <- function() {
  # Test gamma
  gamma_test <- c(0.01, 0.5, 1.0, 1.5, 2.0)
  for(g in gamma_test) {
    g_raw <- inverse_transform_gamma(g)
    g_back <- transform_gamma(g_raw)
    if(abs(g - g_back) > 1e-10) {
      stop(sprintf("Gamma transform failed: %.4f → %.4f → %.4f", g, g_raw, g_back))
    }
  }

  # Test delta
  delta_test <- c(-20, -10, 0, 10, 20)
  for(d in delta_test) {
    d_raw <- inverse_transform_delta(d)
    d_back <- transform_delta(d_raw)
    if(abs(d - d_back) > 1e-10) {
      stop(sprintf("Delta transform failed: %.4f → %.4f → %.4f", d, d_raw, d_back))
    }
  }

  # Test beta_0_trig
  beta_0_trig_test <- c(-10, -5, 0, 5, 10)
  for(b in beta_0_trig_test) {
    b_raw <- inverse_transform_beta_0_trig(b)
    b_back <- transform_beta_0_trig(b_raw)
    if(abs(b - b_back) > 1e-10) {
      stop(sprintf("Beta_0_trig transform failed: %.4f → %.4f → %.4f", b, b_raw, b_back))
    }
  }

  # Test beta_violence
  beta_violence_test <- c(-5, -2, 0, 2, 5)
  for(b in beta_violence_test) {
    b_raw <- inverse_transform_beta_violence(b)
    b_back <- transform_beta_violence(b_raw)
    if(abs(b - b_back) > 1e-10) {
      stop(sprintf("Beta_violence transform failed: %.4f → %.4f → %.4f", b, b_raw, b_back))
    }
  }

  # Test beta_state
  beta_state_test <- c(-5, -2, 0, 2, 5)
  for(b in beta_state_test) {
    b_raw <- inverse_transform_beta_state(b)
    b_back <- transform_beta_state(b_raw)
    if(abs(b - b_back) > 1e-10) {
      stop(sprintf("Beta_state transform failed: %.4f → %.4f → %.4f", b, b_raw, b_back))
    }
  }

  # Test decay
  decay_test <- c(0.001, 0.1, 1, 5, 10)
  for(dec in decay_test) {
    dec_raw <- inverse_transform_decay(dec)
    dec_back <- transform_decay(dec_raw)
    if(abs(dec - dec_back) > 1e-10) {
      stop(sprintf("Decay transform failed: %.4f → %.4f → %.4f", dec, dec_raw, dec_back))
    }
  }

  cat("✓ Parameter transformations verified\n")
}

# Run test
test_transforms()

# ====================================================================
# CHECKPOINT FUNCTIONS
# ====================================================================

checkpoint_dir <- "checkpoints_model_poverty"
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
# 1. LOAD DATA WITH POVERTY
# ====================================================================

cat("\n=== LOADING DATA WITH POVERTY ===\n")
protests <- readRDS("protests_with_poverty.rds")

cat("Loaded", nrow(protests), "events with poverty data\n")
cat("Time range:", min(protests$year), "to", max(protests$year), "\n")

# FULL DATASET MODE: Using all years 2015-2024
cat(sprintf("✓ Using full dataset: %d events across %d years\n",
            nrow(protests), length(unique(protests$year))))

# TEST MODE: Sample 1000 events for quick testing
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
cat(sprintf("Poverty data coverage: %d (%.1f%%)\n", n_with_poverty, pct_coverage))

if(pct_coverage < 95) {
  cat("⚠ WARNING: Low poverty coverage may affect results\n")
}

# Summary statistics
cat("\nPoverty rate summary (rescaled to decimal):\n")
cat(sprintf("  Mean: %.4f (%.2f%%)\n", mean(protests$poverty_decimal, na.rm = TRUE),
            mean(protests$poverty_decimal, na.rm = TRUE)*100))
cat(sprintf("  Median: %.4f (%.2f%%)\n", median(protests$poverty_decimal, na.rm = TRUE),
            median(protests$poverty_decimal, na.rm = TRUE)*100))
cat(sprintf("  Range: %.4f to %.4f (%.2f%% to %.2f%%)\n",
            min(protests$poverty_decimal, na.rm = TRUE),
            max(protests$poverty_decimal, na.rm = TRUE),
            min(protests$poverty_decimal, na.rm = TRUE)*100,
            max(protests$poverty_decimal, na.rm = TRUE)*100))

# Correlation check
cor_pop_poverty <- cor(protests$log_pop, protests$poverty_decimal, use = "complete.obs")
cat(sprintf("\nCorrelation (log(pop) vs poverty): %.3f\n", cor_pop_poverty))
if(abs(cor_pop_poverty) > 0.5) {
  cat("  ⚠ Moderate correlation - interpret γ carefully\n")
} else {
  cat("  ✓ Low to moderate correlation\n")
}
cat("\n")

# Prepare marks data
marks_data <- protests %>%
  select(
    event_id = event_id_cnty,
    time = days_since_start,
    is_violent,
    is_peaceful,
    state_intervention,
    event_type,
    assoc_actor_1,
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
    # New triggering marks
    is_riot = as.integer(event_type == "Riots"),  # Riots vs organized protests
    has_fatalities = as.integer(fatalities > 0),  # Any deaths
    is_student = as.integer(grepl("Student", assoc_actor_1, ignore.case = TRUE)),  # Student-led
    is_labor = as.integer(grepl("Labor|Worker", assoc_actor_1, ignore.case = TRUE))  # Labor-led
  )

# Filter to complete cases only
marks_data <- marks_data %>%
  filter(!is.na(poverty_decimal))

cat(sprintf("Using complete cases: %d events (%.1f%% of total)\n",
            nrow(marks_data),
            100 * nrow(marks_data) / nrow(protests)))
cat(sprintf("Violent events: %d (%.1f%%)\n",
            sum(marks_data$is_violent),
            100*mean(marks_data$is_violent)))
cat(sprintf("Events with state intervention: %d (%.1f%%)\n",
            sum(marks_data$state_intervention),
            100*mean(marks_data$state_intervention)))
cat(sprintf("Riot events: %d (%.1f%%)\n",
            sum(marks_data$is_riot),
            100*mean(marks_data$is_riot)))
cat(sprintf("Events with fatalities: %d (%.1f%%)\n",
            sum(marks_data$has_fatalities),
            100*mean(marks_data$has_fatalities)))
cat(sprintf("Student-led events: %d (%.1f%%)\n",
            sum(marks_data$is_student),
            100*mean(marks_data$is_student)))
cat(sprintf("Labor-led events: %d (%.1f%%)\n\n",
            sum(marks_data$is_labor),
            100*mean(marks_data$is_labor)))

# ====================================================================
# 2. CALCULATE DISTRICT-YEAR OBSERVATION TIMES
# ====================================================================

cat("=== CALCULATING DISTRICT-YEAR OBSERVATION TIMES ===\n")

T_min <- min(marks_data$time)
T_max <- max(marks_data$time)

# Get unique district-year combinations with poverty
district_years <- marks_data %>%
  select(district, year, population, log_pop, poverty_decimal) %>%
  distinct()

cat(sprintf("Found %d unique district-year combinations\n", nrow(district_years)))

# Calculate days observed for each district-year
district_years$days_observed <- T_max - T_min

cat(sprintf("Observation period: %.1f days (%.2f years)\n\n",
            T_max - T_min, (T_max - T_min)/365.25))

# ====================================================================
# 3. MODEL WITH POVERTY IMPLEMENTATION
# ====================================================================

cat("=== IMPLEMENTING MODEL WITH POVERTY ===\n\n")

model_poverty_functions <- function() {

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
  # MULTI-YEAR MODEL: μ(d,y) = exp(β₀ + γ·log_pop + δ·poverty_decimal + β_year)
  calc_mu <- function(log_pop, poverty_decimal, year, params) {

    beta_0_bg <- params[["beta_0_bg"]]

    # Transform raw (unconstrained) parameters to bounded parameters
    gamma <- transform_gamma(params[["gamma_raw"]])
    delta <- transform_delta(params[["delta_raw"]])

    # Lookup year effect based on discrete year (2015 is baseline)
    beta_year <- ifelse(year == 2015, 0.0,
                 ifelse(year == 2016, transform_beta_year(params[["beta_2016_raw"]]),
                 ifelse(year == 2017, transform_beta_year(params[["beta_2017_raw"]]),
                 ifelse(year == 2018, transform_beta_year(params[["beta_2018_raw"]]),
                 ifelse(year == 2019, transform_beta_year(params[["beta_2019_raw"]]),
                 ifelse(year == 2020, transform_beta_year(params[["beta_2020_raw"]]),
                 ifelse(year == 2021, transform_beta_year(params[["beta_2021_raw"]]),
                 ifelse(year == 2022, transform_beta_year(params[["beta_2022_raw"]]),
                 ifelse(year == 2023, transform_beta_year(params[["beta_2023_raw"]]),
                 ifelse(year == 2024, transform_beta_year(params[["beta_2024_raw"]]), 0.0))))))))))

    mu <- exp(beta_0_bg + gamma * log_pop + delta * poverty_decimal + beta_year)
    return(mu)
  }

  # Conditional intensity
  lambda_poverty <- function(idx, times, marks, params, time_mat = NULL) {

    # District-year-specific background rate (WITH POVERTY)
    mu_i <- calc_mu(marks$log_pop[idx], marks$poverty_decimal[idx],
                   marks$year_centered[idx], marks$election_year[idx], params)
    intensity <- mu_i

    # Mark-dependent triggering parameters (transform from raw/unconstrained)
    beta_0_trig <- transform_beta_0_trig(params[["beta_0_trig_raw"]])
    beta_violence <- transform_beta_violence(params[["beta_violence_raw"]])
    beta_state <- transform_beta_state(params[["beta_state_raw"]])
    decay <- transform_decay(params[["decay_raw"]])

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

  # Log-likelihood (C++ implementation for speed)
  loglik_poverty <- function(params, times, marks, district_years,
                            time_mat = NULL, verbose = TRUE, model_id = NULL) {

    # No parameter constraints needed - sigmoid transformations ensure valid ranges

    # Track function evaluations
    .optim_state$n_evals <- .optim_state$n_evals + 1
    eval_num <- .optim_state$n_evals

    # Prepare parameter vector for C++ (18 parameters: 3 background + 9 year effects + 6 triggering)
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
      # Triggering parameters
      params[["beta_0_trig_raw"]],
      params[["beta_riot_raw"]],
      params[["beta_fatal_raw"]],
      params[["beta_student_raw"]],
      params[["beta_labor_raw"]],
      params[["decay_raw"]]
    )

    # Call C++ function (returns NEGATIVE log-likelihood)
    neg_ll <- hawkes_negloglik_cpp(
      times = times,
      log_pop = marks$log_pop,
      poverty_decimal = marks$poverty_decimal,
      year = as.integer(marks$year),
      is_riot = as.integer(marks$is_riot),
      has_fatalities = as.integer(marks$has_fatalities),
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

    # Track best result
    if(ll_final > .optim_state$best_ll) {
      .optim_state$best_ll <- ll_final
      .optim_state$best_params <- params
    }

    # Intermediate checkpoint
    if(!is.null(model_id) && eval_num %% CHECKPOINT_EVERY_N_EVALS == 0) {
      # Transform parameters for display
      gamma_display <- transform_gamma(params[["gamma_raw"]])
      delta_display <- transform_delta(params[["delta_raw"]])

      save_intermediate_checkpoint(model_id, eval_num, params, ll_final)
      if(verbose) {
        cat(sprintf("\n    💾 Checkpoint %d: LL=%.2f (best: %.2f) | γ=%.4f δ=%.4f\n",
                   eval_num, ll_final, .optim_state$best_ll,
                   gamma_display, delta_display))
      }
    }

    if(verbose) {
      # Transform parameters for display
      gamma_display <- transform_gamma(params[["gamma_raw"]])
      delta_display <- transform_delta(params[["delta_raw"]])

      cat(sprintf("    LL: %.2f | γ=%.4f δ=%.4f [eval %d]\n",
                 ll_final, gamma_display, delta_display, eval_num))
      flush.console()  # Force output to display immediately
    }

    return(ll_final)
  }

  return(list(
    precompute_time_matrix = precompute_time_matrix,
    calc_mu = calc_mu,
    lambda = lambda_poverty,
    loglik = loglik_poverty
  ))
}

hawkes_funcs <- model_poverty_functions()

# ====================================================================
# 4. PREPARE DATA
# ====================================================================

cat("=== PREPARING DATA FOR ESTIMATION ===\n")

times_sample <- marks_data$time
marks_sample <- marks_data
district_years_sample <- district_years

cat(sprintf("Sample size: %d events\n", length(times_sample)))
cat(sprintf("District-years: %d\n", nrow(district_years_sample)))

# Pre-compute time difference matrix (SKIPPED - C++ computes on-the-fly)
cat("Skipping time matrix pre-computation (C++ handles this)...\n")
time_mat_sample <- NULL  # Not needed for C++ implementation
cat("✓ Ready for optimization\n\n")

# ====================================================================
# 5. FITTING FUNCTION
# ====================================================================

fit_model_poverty <- function(times, marks, district_years, time_mat = NULL,
                        model_name = "Model with Poverty Rate",
                        model_id = "MODEL_POVERTY",
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
      # Transform raw parameters back to original scale for display
      gamma_resume <- transform_gamma(checkpoint$data$best_params$gamma_raw)
      delta_resume <- transform_delta(checkpoint$data$best_params$delta_raw)
      cat(sprintf("⏩ Resuming from intermediate checkpoint (eval %d, LL=%.2f, γ=%.4f, δ=%.4f)\n",
                 checkpoint$data$eval_num, checkpoint$data$best_ll,
                 gamma_resume, delta_resume))

      reset_optim_state(model_name)
      .optim_state$n_evals <- checkpoint$data$eval_num
      .optim_state$best_ll <- checkpoint$data$best_ll
      .optim_state$best_params <- checkpoint$data$best_params

      # Use checkpoint params, disable multistart
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

    # Define multiple starting points
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
        # Year effects (2016-2024) - start at 0 (transforms to 0)
        beta_2016_raw = 0, beta_2017_raw = 0, beta_2018_raw = 0,
        beta_2019_raw = 0, beta_2020_raw = 0, beta_2021_raw = 0,
        beta_2022_raw = 0, beta_2023_raw = 0, beta_2024_raw = 0,
        # Triggering parameters (new marks, start at 0)
        beta_0_trig_raw = beta_0_trig1_raw,
        beta_riot_raw = 0, beta_fatal_raw = 0,
        beta_student_raw = 0, beta_labor_raw = 0,
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
          # Year effects - slight variation from baseline
          beta_2016_raw = -0.1, beta_2017_raw = -0.05, beta_2018_raw = 0,
          beta_2019_raw = 0, beta_2020_raw = 0.05, beta_2021_raw = 0.1,
          beta_2022_raw = 0.1, beta_2023_raw = 0.05, beta_2024_raw = 0,
          # Triggering parameters - small positive effects
          beta_0_trig_raw = beta_0_trig2_raw,
          beta_riot_raw = 0.2, beta_fatal_raw = 0.2,
          beta_student_raw = 0.1, beta_labor_raw = 0.1,
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
          # Year effects - larger variation
          beta_2016_raw = 0.2, beta_2017_raw = 0.15, beta_2018_raw = 0.1,
          beta_2019_raw = 0, beta_2020_raw = -0.1, beta_2021_raw = -0.15,
          beta_2022_raw = -0.1, beta_2023_raw = 0, beta_2024_raw = 0.1,
          # Triggering parameters - negative effects (test alternative)
          beta_0_trig_raw = beta_0_trig3_raw,
          beta_riot_raw = -0.2, beta_fatal_raw = -0.1,
          beta_student_raw = -0.1, beta_labor_raw = -0.2,
          decay_raw = decay3_raw
        )
      )
    }

    if(use_multistart && n_starts >= 4) {
      # 4. Negative poverty (wealth drives protests?)
      gamma4 <- 0.50; delta4 <- -0.5
      gamma4_raw <- inverse_transform_gamma(gamma4)
      delta4_raw <- inverse_transform_delta(delta4)
      beta_0_trig4 <- log(0.12); decay4 <- 0.18
      beta_0_trig4_raw <- inverse_transform_beta_0_trig(beta_0_trig4)
      decay4_raw <- inverse_transform_decay(decay4)
      beta_0_4 <- log(observed_rate) - gamma4 * mean_log_pop - delta4 * mean_poverty
      starting_points[[4]] <- list(
        name = "Negative-poverty",
        params = list(
          beta_0_bg = beta_0_4, gamma_raw = gamma4_raw, delta_raw = delta4_raw,
          # Year effects - random variation
          beta_2016_raw = 0.05, beta_2017_raw = -0.05, beta_2018_raw = 0.1,
          beta_2019_raw = -0.1, beta_2020_raw = 0, beta_2021_raw = 0.05,
          beta_2022_raw = -0.05, beta_2023_raw = 0.1, beta_2024_raw = -0.1,
          # Triggering parameters - mixed effects
          beta_0_trig_raw = beta_0_trig4_raw,
          beta_riot_raw = 0.1, beta_fatal_raw = -0.1,
          beta_student_raw = 0, beta_labor_raw = 0,
          decay_raw = decay4_raw
        )
      )
    }

    if(use_multistart && n_starts >= 5) {
      # 5. High poverty effect
      gamma5 <- 0.25; delta5 <- 1.0
      gamma5_raw <- inverse_transform_gamma(gamma5)
      delta5_raw <- inverse_transform_delta(delta5)
      beta_0_trig5 <- log(0.2); decay5 <- 0.3
      beta_0_trig5_raw <- inverse_transform_beta_0_trig(beta_0_trig5)
      decay5_raw <- inverse_transform_decay(decay5)
      beta_0_5 <- log(observed_rate) - gamma5 * mean_log_pop - delta5 * mean_poverty
      starting_points[[5]] <- list(
        name = "High-poverty",
        params = list(
          beta_0_bg = beta_0_5, gamma_raw = gamma5_raw, delta_raw = delta5_raw,
          # Year effects - strong positive trend
          beta_2016_raw = -0.3, beta_2017_raw = -0.2, beta_2018_raw = -0.1,
          beta_2019_raw = 0, beta_2020_raw = 0.1, beta_2021_raw = 0.2,
          beta_2022_raw = 0.3, beta_2023_raw = 0.4, beta_2024_raw = 0.5,
          # Triggering parameters - strong positive effects
          beta_0_trig_raw = beta_0_trig5_raw,
          beta_riot_raw = 0.3, beta_fatal_raw = 0.3,
          beta_student_raw = 0.2, beta_labor_raw = 0.2,
          decay_raw = decay5_raw
        )
      )
    }
  }

  # Parameter names (18 total: 3 background + 9 year effects + 6 triggering)
  param_names <- c("beta_0_bg", "gamma_raw", "delta_raw",
                   # Year effects (2016-2024)
                   "beta_2016_raw", "beta_2017_raw", "beta_2018_raw",
                   "beta_2019_raw", "beta_2020_raw", "beta_2021_raw",
                   "beta_2022_raw", "beta_2023_raw", "beta_2024_raw",
                   # Triggering parameters
                   "beta_0_trig_raw", "beta_riot_raw", "beta_fatal_raw",
                   "beta_student_raw", "beta_labor_raw", "decay_raw")

  # Bounds: All parameters use sigmoid transforms, so set wide bounds
  # to allow L-BFGS-B to explore freely without hitting constraints
  lower_bounds <- c(beta_0_bg = -15, gamma_raw = -10, delta_raw = -10,
                    # Year effects (transformed to [-3, 3])
                    beta_2016_raw = -10, beta_2017_raw = -10, beta_2018_raw = -10,
                    beta_2019_raw = -10, beta_2020_raw = -10, beta_2021_raw = -10,
                    beta_2022_raw = -10, beta_2023_raw = -10, beta_2024_raw = -10,
                    # Triggering (transformed to appropriate ranges)
                    beta_0_trig_raw = -10, beta_riot_raw = -10, beta_fatal_raw = -10,
                    beta_student_raw = -10, beta_labor_raw = -10, decay_raw = -10)

  upper_bounds <- c(beta_0_bg = 5, gamma_raw = 10, delta_raw = 10,
                    # Year effects
                    beta_2016_raw = 10, beta_2017_raw = 10, beta_2018_raw = 10,
                    beta_2019_raw = 10, beta_2020_raw = 10, beta_2021_raw = 10,
                    beta_2022_raw = 10, beta_2023_raw = 10, beta_2024_raw = 10,
                    # Triggering
                    beta_0_trig_raw = 10, beta_riot_raw = 10, beta_fatal_raw = 10,
                    beta_student_raw = 10, beta_labor_raw = 10, decay_raw = 10)

  # Objective function (without verbose for multi-start)
  obj_fun_silent <- function(par) {
    params_full <- as.list(par)
    ll <- hawkes_funcs$loglik(params_full, times, marks, district_years, time_mat,
                               verbose = FALSE, model_id = NULL)
    return(-ll)
  }

  # Objective function (with verbose for final run)
  obj_fun_verbose <- function(par) {
    params_full <- as.list(par)
    ll <- hawkes_funcs$loglik(params_full, times, marks, district_years, time_mat,
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
    # Transform raw parameters back to original scale for display
    gamma_init <- transform_gamma(start_point$params$gamma_raw)
    delta_init <- transform_delta(start_point$params$delta_raw)
    cat(sprintf("  Initial: γ=%.3f, δ=%.3f\n", gamma_init, delta_init))

    param_vec <- unlist(start_point$params[param_names])
    names(param_vec) <- param_names

    # Use silent objective for multi-start, verbose only for final
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

    # Transform raw parameters back to original scale for display
    gamma_result <- transform_gamma(all_results[[i]]$params$gamma_raw)
    delta_result <- transform_delta(all_results[[i]]$params$delta_raw)
    cat(sprintf("  → LL: %.2f | γ=%.4f, δ=%.4f | Conv: %d | Time: %.1f min\n\n",
               all_results[[i]]$loglik,
               gamma_result, delta_result,
               all_results[[i]]$convergence,
               all_results[[i]]$runtime))
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

  # If best start wasn't the first, re-run with verbose and checkpointing
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
  # Transform raw parameters back to original scale for display
  cat(sprintf("  Final γ: %.4f\n", transform_gamma(final_params$gamma_raw)))
  cat(sprintf("  Final δ: %.4f\n", transform_delta(final_params$delta_raw)))

  # Assemble results
  results <- list(
    model_name = model_name,
    model_id = model_id,
    params = final_params,
    loglik = final_ll,
    convergence = final_convergence,
    n_params = 7,
    n_events = length(times),
    AIC = 2*7 - 2*final_ll,
    BIC = 7*log(length(times)) - 2*final_ll,
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
cat("║   FITTING MODEL WITH POVERTY RATE                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

model_poverty_fit <- fit_model_poverty(
  times_sample,
  marks_sample,
  district_years_sample,
  time_mat_sample,
  model_name = "Model with Poverty Rate (Multi-Start)",
  model_id = "MODEL_POVERTY",
  use_multistart = !TEST_MODE,  # Disable multi-start in test mode
  n_starts = if(TEST_MODE) 1 else 5
)

# ====================================================================
# 7. RESULTS
# ====================================================================

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║         MODEL WITH POVERTY RESULTS                           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("MODEL FIT:\n")
cat(sprintf("  Log-likelihood: %.2f\n", model_poverty_fit$loglik))
cat(sprintf("  AIC: %.2f\n", model_poverty_fit$AIC))
cat(sprintf("  BIC: %.2f\n", model_poverty_fit$BIC))
cat(sprintf("  Parameters: %d\n", model_poverty_fit$n_params))
cat(sprintf("  Runtime: %.2f minutes\n", model_poverty_fit$runtime_mins))

# Print multi-start summary if available
if(!is.null(model_poverty_fit$multistart_results)) {
  cat(sprintf("  Best starting point: %s\n\n", model_poverty_fit$best_start))

  cat("MULTI-START COMPARISON:\n")
  cat("-----------------------\n")
  for(i in seq_along(model_poverty_fit$multistart_results)) {
    res <- model_poverty_fit$multistart_results[[i]]
    marker <- if(res$name == model_poverty_fit$best_start) "★" else " "
    # Transform raw parameters back to original scale for display
    gamma_display <- transform_gamma(res$params$gamma_raw)
    delta_display <- transform_delta(res$params$delta_raw)
    cat(sprintf("%s %d. %-18s LL=%.2f | γ=%.4f, δ=%.4f, trend=%.4f, elec=%.4f | Conv=%d | %.1f min\n",
               marker, i, res$name, res$loglik, gamma_display,
               delta_display, res$params$beta_trend, res$params$beta_election,
               res$convergence, res$runtime))
  }
  cat("\n")
} else {
  cat("\n")
}

cat("BACKGROUND RATE PARAMETERS:\n")
cat("---------------------------\n")
# Transform raw parameters back to original scale for interpretation
gamma_final <- transform_gamma(model_poverty_fit$params$gamma_raw)
delta_final <- transform_delta(model_poverty_fit$params$delta_raw)

cat(sprintf("  β₀ (baseline): %.4f\n", model_poverty_fit$params$beta_0_bg))
cat(sprintf("  γ (population elasticity): %.4f ★\n", gamma_final))
cat(sprintf("    → Interpretation: 10× population → %.1f× more protests\n",
            10^gamma_final))
cat(sprintf("  δ (poverty effect): %.4f ★\n", delta_final))
cat(sprintf("    → Interpretation: +0.01 poverty (1%%) → %.3f× more protests\n",
            exp(delta_final * 0.01)))
cat(sprintf("    → Or equivalently: +0.10 poverty (10%%) → %.3f× more protests\n",
            exp(delta_final * 0.10)))

if(delta_final > 0) {
  cat("    → ✓ Higher poverty → MORE protests\n")
} else {
  cat("    → Lower poverty → MORE protests (unexpected!)\n")
}

# Compare to baseline expectation
if(gamma_final > 0.35 && gamma_final < 1) {
  cat("\n  INTERPRETATION:\n")
  cat("  → γ increased from baseline (0.35), suggesting poverty partially\n")
  cat("    explains the sub-proportional population effect\n")
} else if(gamma_final >= 1) {
  cat("\n  INTERPRETATION:\n")
  cat("  → γ ≥ 1: Poverty fully explains sub-proportional effect!\n")
  cat("    Once poverty controlled, protests scale proportionally with population\n")
}
cat("\n")

cat("YEAR EFFECTS (2015 is baseline):\n")
cat("--------------------------------\n")
# Transform year effects back to original scale
beta_2016_final <- transform_beta_year(model_poverty_fit$params$beta_2016_raw)
beta_2017_final <- transform_beta_year(model_poverty_fit$params$beta_2017_raw)
beta_2018_final <- transform_beta_year(model_poverty_fit$params$beta_2018_raw)
beta_2019_final <- transform_beta_year(model_poverty_fit$params$beta_2019_raw)
beta_2020_final <- transform_beta_year(model_poverty_fit$params$beta_2020_raw)
beta_2021_final <- transform_beta_year(model_poverty_fit$params$beta_2021_raw)
beta_2022_final <- transform_beta_year(model_poverty_fit$params$beta_2022_raw)
beta_2023_final <- transform_beta_year(model_poverty_fit$params$beta_2023_raw)
beta_2024_final <- transform_beta_year(model_poverty_fit$params$beta_2024_raw)

cat(sprintf("  β_2016: %+.4f (exp = %.3f)\n", beta_2016_final, exp(beta_2016_final)))
cat(sprintf("  β_2017: %+.4f (exp = %.3f)\n", beta_2017_final, exp(beta_2017_final)))
cat(sprintf("  β_2018: %+.4f (exp = %.3f)\n", beta_2018_final, exp(beta_2018_final)))
cat(sprintf("  β_2019: %+.4f (exp = %.3f)\n", beta_2019_final, exp(beta_2019_final)))
cat(sprintf("  β_2020: %+.4f (exp = %.3f)\n", beta_2020_final, exp(beta_2020_final)))
cat(sprintf("  β_2021: %+.4f (exp = %.3f)\n", beta_2021_final, exp(beta_2021_final)))
cat(sprintf("  β_2022: %+.4f (exp = %.3f)\n", beta_2022_final, exp(beta_2022_final)))
cat(sprintf("  β_2023: %+.4f (exp = %.3f)\n", beta_2023_final, exp(beta_2023_final)))
cat(sprintf("  β_2024: %+.4f (exp = %.3f)\n", beta_2024_final, exp(beta_2024_final)))
cat("\n")

cat("TRIGGERING PARAMETERS:\n")
cat("----------------------\n")
# Transform raw parameters back to original scale for display
beta_0_trig_final <- transform_beta_0_trig(model_poverty_fit$params$beta_0_trig_raw)
beta_riot_final <- transform_beta_riot(model_poverty_fit$params$beta_riot_raw)
beta_fatal_final <- transform_beta_fatal(model_poverty_fit$params$beta_fatal_raw)
beta_student_final <- transform_beta_student(model_poverty_fit$params$beta_student_raw)
beta_labor_final <- transform_beta_labor(model_poverty_fit$params$beta_labor_raw)
decay_final <- transform_decay(model_poverty_fit$params$decay_raw)

cat(sprintf("  β₀_trig: %.4f\n", beta_0_trig_final))
cat(sprintf("  β_riot: %.4f (exp = %.3f)\n",
            beta_riot_final, exp(beta_riot_final)))
cat(sprintf("  β_fatal: %.4f (exp = %.3f)\n",
            beta_fatal_final, exp(beta_fatal_final)))
cat(sprintf("  β_student: %.4f (exp = %.3f)\n",
            beta_student_final, exp(beta_student_final)))
cat(sprintf("  β_labor: %.4f (exp = %.3f)\n",
            beta_labor_final, exp(beta_labor_final)))
cat(sprintf("  Decay: %.4f (half-life = %.1f days)\n",
            decay_final, log(2)/decay_final))
cat("\n")

# Save results
saveRDS(model_poverty_fit, "model_poverty.rds")
cat("✓ Saved: model_poverty.rds\n\n")

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   MODEL WITH POVERTY ESTIMATION COMPLETE                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("NEXT STEPS:\n")
cat("-----------\n")
cat("1. Compare AIC/BIC with baseline Model B (17_model_comparison.R)\n")
cat("2. Check if poverty explains population effect (γ increased?)\n")
cat("3. Residual analysis to validate improved fit\n")
cat("4. Consider interaction terms if needed\n")
cat("\n")
