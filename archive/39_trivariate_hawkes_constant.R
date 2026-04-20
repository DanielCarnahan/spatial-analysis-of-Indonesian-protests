################################################################################
#           TRIVARIATE HAWKES MODEL - CONSTANT BASELINE
#           Multi-starting point optimization with hypothesis testing
#
#           Model: Three interacting point processes (Fatal, Violent, Peaceful)
#           λ_F(t) = μ_F + Σ α_FF·g(t-t_j) + Σ α_VF·g(t-t_k) + Σ α_PF·g(t-t_l)
#           λ_V(t) = μ_V + Σ α_FV·g(t-t_j) + Σ α_VV·g(t-t_k) + Σ α_PV·g(t-t_l)
#           λ_P(t) = μ_P + Σ α_FP·g(t-t_j) + Σ α_VP·g(t-t_k) + Σ α_PP·g(t-t_l)
#           where g(s) = exp(-β·s)
################################################################################

library(tidyverse)
library(numDeriv)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   TRIVARIATE HAWKES MODEL - CONSTANT BASELINE                           ║\n")
cat("║   Event types: Fatal (F), Violent non-fatal (V), Peaceful (P)           ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

N_STARTS <- 5           # Number of starting points
TEMPORAL_CUTOFF <- 90   # Days - for computational tractability
MAX_ITER <- 2000        # Maximum iterations per optimization

cat("Configuration:\n")
cat(sprintf("  Starting points: %d\n", N_STARTS))
cat(sprintf("  Temporal cutoff: %d days\n", TEMPORAL_CUTOFF))
cat(sprintf("  Max iterations: %d\n\n", MAX_ITER))

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

cat("=== LOADING DATA ===\n")

protests <- readRDS("protests_trivariate.rds")

# Convert to numeric time (days since start)
start_date <- min(protests$event_date)
protests$time <- as.numeric(protests$event_date - start_date)
T_max <- max(protests$time) + 1

# Extract event counts by type
n_f <- sum(protests$event_type_3 == "F")
n_v <- sum(protests$event_type_3 == "V")
n_p <- sum(protests$event_type_3 == "P")
n_total <- n_f + n_v + n_p

cat(sprintf("Loaded %d events\n", nrow(protests)))
cat(sprintf("  Fatal events (F):    %d (%.1f%%)\n", n_f, 100*n_f/n_total))
cat(sprintf("  Violent events (V):  %d (%.1f%%)\n", n_v, 100*n_v/n_total))
cat(sprintf("  Peaceful events (P): %d (%.1f%%)\n", n_p, 100*n_p/n_total))
cat(sprintf("  Observation period: %.0f days (%.1f years)\n\n", T_max, T_max/365.25))

# =============================================================================
# 2. TRIVARIATE HAWKES LOG-LIKELIHOOD (EFFICIENT RECURSIVE COMPUTATION)
# =============================================================================

cat("=== SETTING UP LIKELIHOOD FUNCTION ===\n")

# Combine and sort all events for recursive computation
# type_numeric: 1=Fatal, 2=Violent, 3=Peaceful
all_times <- protests$time
all_types <- protests$type_numeric
ord <- order(all_times)
all_times <- all_times[ord]
all_types <- all_types[ord]
n <- length(all_times)

cat(sprintf("Combined %d events for recursive computation\n\n", n))

# Pre-compute time differences for efficiency
time_diffs <- diff(all_times)

trivariate_loglik <- function(params) {
  # Parameters (log scale for unconstrained optimization)
  # 13 total: 3 baselines + 9 alpha + 1 beta

  mu_f <- exp(params[1])
  mu_v <- exp(params[2])
  mu_p <- exp(params[3])

  # Triggering matrix (row = parent type, col = child type)
  # α_XY = effect of X parent on Y child intensity
  alpha_ff <- exp(params[4])   # F → F
  alpha_fv <- exp(params[5])   # F → V
  alpha_fp <- exp(params[6])   # F → P
  alpha_vf <- exp(params[7])   # V → F
  alpha_vv <- exp(params[8])   # V → V
  alpha_vp <- exp(params[9])   # V → P
  alpha_pf <- exp(params[10])  # P → F
  alpha_pv <- exp(params[11])  # P → V
  alpha_pp <- exp(params[12])  # P → P

  beta <- exp(params[13])

  # Check for valid parameters
  all_params <- c(mu_f, mu_v, mu_p, alpha_ff, alpha_fv, alpha_fp,
                  alpha_vf, alpha_vv, alpha_vp, alpha_pf, alpha_pv, alpha_pp, beta)
  if (any(!is.finite(all_params))) {
    return(-1e10)
  }
  if (beta < 1e-6 || beta > 100) return(-1e10)

  # Recursive computation of triggering sums
  # R_f[i] = sum of exp(-β*(t_i - t_j)) for fatal j < i
  # R_v[i] = sum of exp(-β*(t_i - t_j)) for violent j < i
  # R_p[i] = sum of exp(-β*(t_i - t_j)) for peaceful j < i

  R_f <- numeric(n)
  R_v <- numeric(n)
  R_p <- numeric(n)

  for (i in 2:n) {
    decay <- exp(-beta * time_diffs[i-1])

    # Limit decay to avoid numerical issues
    if (decay < 1e-100) decay <- 0

    if (all_types[i-1] == 1) {
      # Previous event was Fatal
      R_f[i] <- decay * (1 + R_f[i-1])
      R_v[i] <- decay * R_v[i-1]
      R_p[i] <- decay * R_p[i-1]
    } else if (all_types[i-1] == 2) {
      # Previous event was Violent
      R_f[i] <- decay * R_f[i-1]
      R_v[i] <- decay * (1 + R_v[i-1])
      R_p[i] <- decay * R_p[i-1]
    } else {
      # Previous event was Peaceful
      R_f[i] <- decay * R_f[i-1]
      R_v[i] <- decay * R_v[i-1]
      R_p[i] <- decay * (1 + R_p[i-1])
    }
  }

  # Log-likelihood: sum of log(λ) at event times
  loglik <- 0

  for (i in 1:n) {
    if (all_types[i] == 1) {
      # Fatal event: λ_F(t) = μ_F + α_FF*R_f + α_VF*R_v + α_PF*R_p
      lambda_i <- mu_f + alpha_ff * R_f[i] + alpha_vf * R_v[i] + alpha_pf * R_p[i]
    } else if (all_types[i] == 2) {
      # Violent event: λ_V(t) = μ_V + α_FV*R_f + α_VV*R_v + α_PV*R_p
      lambda_i <- mu_v + alpha_fv * R_f[i] + alpha_vv * R_v[i] + alpha_pv * R_p[i]
    } else {
      # Peaceful event: λ_P(t) = μ_P + α_FP*R_f + α_VP*R_v + α_PP*R_p
      lambda_i <- mu_p + alpha_fp * R_f[i] + alpha_vp * R_v[i] + alpha_pp * R_p[i]
    }

    if (lambda_i <= 1e-10) return(-1e10)
    loglik <- loglik + log(lambda_i)
  }

  # Integral terms (compensator)
  # For exponential kernel: ∫_{t_j}^{T} α·exp(-β·(t-t_j)) dt = (α/β)·(1 - exp(-β·(T-t_j)))

  integral_f <- mu_f * T_max
  integral_v <- mu_v * T_max
  integral_p <- mu_p * T_max

  for (j in 1:n) {
    contrib <- (1 - exp(-beta * (T_max - all_times[j]))) / beta

    if (all_types[j] == 1) {
      # Fatal parent contributes to all three child types
      integral_f <- integral_f + alpha_ff * contrib
      integral_v <- integral_v + alpha_fv * contrib
      integral_p <- integral_p + alpha_fp * contrib
    } else if (all_types[j] == 2) {
      # Violent parent
      integral_f <- integral_f + alpha_vf * contrib
      integral_v <- integral_v + alpha_vv * contrib
      integral_p <- integral_p + alpha_vp * contrib
    } else {
      # Peaceful parent
      integral_f <- integral_f + alpha_pf * contrib
      integral_v <- integral_v + alpha_pv * contrib
      integral_p <- integral_p + alpha_pp * contrib
    }
  }

  loglik <- loglik - integral_f - integral_v - integral_p

  if (!is.finite(loglik)) return(-1e10)

  return(loglik)
}

# Wrapper for optimization (minimization)
neg_loglik <- function(params) {
  -trivariate_loglik(params)
}

cat("✓ Likelihood function ready\n\n")

# =============================================================================
# 3. MULTI-STARTING POINT OPTIMIZATION
# =============================================================================

cat("=== MULTI-STARTING POINT OPTIMIZATION ===\n\n")

# Empirical estimates for starting values
emp_rate_f <- n_f / T_max
emp_rate_v <- n_v / T_max
emp_rate_p <- n_p / T_max

cat(sprintf("Empirical rates: μ_F ≈ %.5f, μ_V ≈ %.4f, μ_P ≈ %.4f events/day\n\n",
            emp_rate_f, emp_rate_v, emp_rate_p))

# Define starting points
starting_points <- list(

  # 1. Data-driven
  list(
    name = "Data-driven",
    params = c(
      log(emp_rate_f * 0.5),   # log_mu_f
      log(emp_rate_v * 0.5),   # log_mu_v
      log(emp_rate_p * 0.5),   # log_mu_p
      rep(log(0.02), 9),       # log_alpha (9 triggering params)
      log(0.1)                  # log_beta
    )
  ),

  # 2. Conservative (lower triggering, faster decay)
  list(
    name = "Conservative",
    params = c(
      log(emp_rate_f * 0.7),
      log(emp_rate_v * 0.7),
      log(emp_rate_p * 0.7),
      rep(log(0.01), 9),
      log(0.2)
    )
  ),

  # 3. Strong triggering (higher α, slower decay)
  list(
    name = "Strong-triggering",
    params = c(
      log(emp_rate_f * 0.3),
      log(emp_rate_v * 0.3),
      log(emp_rate_p * 0.3),
      rep(log(0.05), 9),
      log(0.05)
    )
  ),

  # 4. Asymmetric - Fatal triggering stronger
  list(
    name = "Asymmetric-Fatal",
    params = c(
      log(emp_rate_f * 0.5),
      log(emp_rate_v * 0.5),
      log(emp_rate_p * 0.5),
      log(0.10),  # α_FF
      log(0.10),  # α_FV
      log(0.10),  # α_FP
      log(0.02),  # α_VF
      log(0.02),  # α_VV
      log(0.02),  # α_VP
      log(0.01),  # α_PF
      log(0.01),  # α_PV
      log(0.01),  # α_PP
      log(0.15)
    )
  ),

  # 5. Asymmetric - Peaceful dominant
  list(
    name = "Asymmetric-Peaceful",
    params = c(
      log(emp_rate_f * 0.5),
      log(emp_rate_v * 0.5),
      log(emp_rate_p * 0.5),
      log(0.005), # α_FF
      log(0.005), # α_FV
      log(0.005), # α_FP
      log(0.01),  # α_VF
      log(0.01),  # α_VV
      log(0.01),  # α_VP
      log(0.05),  # α_PF
      log(0.05),  # α_PV
      log(0.05),  # α_PP
      log(0.1)
    )
  )
)

# Parameter names
param_names <- c("log_mu_f", "log_mu_v", "log_mu_p",
                 "log_alpha_ff", "log_alpha_fv", "log_alpha_fp",
                 "log_alpha_vf", "log_alpha_vv", "log_alpha_vp",
                 "log_alpha_pf", "log_alpha_pv", "log_alpha_pp",
                 "log_beta")

# Bounds for L-BFGS-B
lower_bounds <- c(rep(-15, 3), rep(-15, 9), -5)
upper_bounds <- c(rep(5, 3), rep(5, 9), 3)

# Track results from all starting points
all_results <- list()
best_ll <- -Inf
best_result <- NULL
best_start <- NULL

start_time_total <- Sys.time()

for (s in seq_along(starting_points)) {
  sp <- starting_points[[s]]

  cat(sprintf("╔══════════════════════════════════════════════════════════════╗\n"))
  cat(sprintf("║  Starting Point %d/%d: %-42s║\n", s, length(starting_points), sp$name))
  cat(sprintf("╚══════════════════════════════════════════════════════════════╝\n\n"))

  start_time <- Sys.time()

  # Run optimization
  result <- tryCatch({
    optim(
      par = sp$params,
      fn = neg_loglik,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(
        maxit = MAX_ITER,
        factr = 1e7,
        trace = 1,
        REPORT = 50
      )
    )
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", e$message))
    return(NULL)
  })

  runtime <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

  if (!is.null(result)) {
    ll <- -result$value

    cat(sprintf("\n✓ Completed in %.1f minutes\n", runtime))
    cat(sprintf("  Log-likelihood: %.2f\n", ll))
    cat(sprintf("  Convergence: %d\n\n", result$convergence))

    all_results[[s]] <- list(
      starting_point = sp$name,
      params = result$par,
      loglik = ll,
      convergence = result$convergence,
      runtime = runtime
    )

    if (ll > best_ll) {
      best_ll <- ll
      best_result <- result
      best_start <- sp$name
    }
  } else {
    cat(sprintf("\n✗ Optimization failed for this starting point\n\n"))
    all_results[[s]] <- list(
      starting_point = sp$name,
      params = NULL,
      loglik = -Inf,
      convergence = -1,
      runtime = runtime
    )
  }
}

total_runtime <- as.numeric(difftime(Sys.time(), start_time_total, units = "mins"))

# =============================================================================
# 4. OPTIMIZATION SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  OPTIMIZATION SUMMARY                                        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

for (s in seq_along(all_results)) {
  r <- all_results[[s]]
  marker <- if (r$starting_point == best_start) " ★ BEST" else ""
  cat(sprintf("  %d. %-22s  LL: %10.2f  Conv: %d  Time: %.1f min%s\n",
              s, r$starting_point, r$loglik, r$convergence, r$runtime, marker))
}

cat(sprintf("\nTotal runtime: %.1f minutes\n", total_runtime))

if (is.null(best_result)) {
  stop("All optimizations failed!")
}

# =============================================================================
# 5. EXTRACT BEST ESTIMATES
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PARAMETER ESTIMATES                                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

params <- best_result$par

# Transform to natural scale
mu_f <- exp(params[1])
mu_v <- exp(params[2])
mu_p <- exp(params[3])

alpha_ff <- exp(params[4])
alpha_fv <- exp(params[5])
alpha_fp <- exp(params[6])
alpha_vf <- exp(params[7])
alpha_vv <- exp(params[8])
alpha_vp <- exp(params[9])
alpha_pf <- exp(params[10])
alpha_pv <- exp(params[11])
alpha_pp <- exp(params[12])

beta <- exp(params[13])

cat("BASELINE RATES (events/day):\n")
cat(sprintf("  μ_F (fatal):    %.6f\n", mu_f))
cat(sprintf("  μ_V (violent):  %.6f\n", mu_v))
cat(sprintf("  μ_P (peaceful): %.6f\n\n", mu_p))

cat("TRIGGERING MATRIX (α coefficients):\n")
cat("                       → Fatal       → Violent     → Peaceful\n")
cat(sprintf("  Fatal parent:       %.6f     %.6f     %.6f\n", alpha_ff, alpha_fv, alpha_fp))
cat(sprintf("  Violent parent:     %.6f     %.6f     %.6f\n", alpha_vf, alpha_vv, alpha_vp))
cat(sprintf("  Peaceful parent:    %.6f     %.6f     %.6f\n\n", alpha_pf, alpha_pv, alpha_pp))

cat(sprintf("DECAY RATE: β = %.6f (half-life = %.1f days)\n\n", beta, log(2)/beta))

# =============================================================================
# 6. COMPUTE STANDARD ERRORS VIA HESSIAN
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  COMPUTING STANDARD ERRORS (Hessian)                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Computing Hessian matrix at MLE...\n")

hessian_mat <- tryCatch({
  hessian(neg_loglik, best_result$par)
}, error = function(e) {
  cat(sprintf("  Warning: Hessian computation error: %s\n", e$message))
  return(NULL)
})

if (!is.null(hessian_mat)) {
  # Variance-covariance matrix
  vcov_mat <- tryCatch({
    solve(hessian_mat)
  }, error = function(e) {
    cat("  Warning: Hessian not invertible, using pseudo-inverse\n")
    MASS::ginv(hessian_mat)
  })

  # Standard errors on log scale
  se_log <- sqrt(pmax(diag(vcov_mat), 0))

  cat("\nParameter estimates with standard errors (log scale):\n")
  cat("─────────────────────────────────────────────────────────────────\n")

  for (i in 1:13) {
    if (se_log[i] > 0 && !is.na(se_log[i]) && is.finite(se_log[i])) {
      z <- params[i] / se_log[i]
      p <- 2 * pnorm(-abs(z))
      sig <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
      cat(sprintf("  %-15s: %8.4f (SE: %6.4f, z: %6.2f, p: %.4f) %s\n",
                  param_names[i], params[i], se_log[i], z, p, sig))
    } else {
      cat(sprintf("  %-15s: %8.4f (SE: NA)\n", param_names[i], params[i]))
    }
  }
  cat("\n")

} else {
  se_log <- rep(NA, 13)
  vcov_mat <- NULL
  cat("  Could not compute standard errors\n\n")
}

# =============================================================================
# 7. MOBILIZATION COMPARISON WITH INFERENCE
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MOBILIZATION POTENTIAL COMPARISON                           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Expected offspring per parent
offspring_ff <- alpha_ff / beta
offspring_fv <- alpha_fv / beta
offspring_fp <- alpha_fp / beta
offspring_vf <- alpha_vf / beta
offspring_vv <- alpha_vv / beta
offspring_vp <- alpha_vp / beta
offspring_pf <- alpha_pf / beta
offspring_pv <- alpha_pv / beta
offspring_pp <- alpha_pp / beta

# Total mobilization by parent type
mobil_fatal <- offspring_ff + offspring_fv + offspring_fp
mobil_violent <- offspring_vf + offspring_vv + offspring_vp
mobil_peaceful <- offspring_pf + offspring_pv + offspring_pp

# Key ratios
ratio_fv <- mobil_fatal / mobil_violent      # Fatal vs Violent
ratio_vp <- mobil_violent / mobil_peaceful   # Violent vs Peaceful
ratio_fp <- mobil_fatal / mobil_peaceful     # Fatal vs Peaceful

cat("EXPECTED OFFSPRING PER PARENT EVENT:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("                       → Fatal    → Violent  → Peaceful  TOTAL\n")
cat(sprintf("  Fatal parent:       %.5f    %.5f    %.5f   %.5f\n",
            offspring_ff, offspring_fv, offspring_fp, mobil_fatal))
cat(sprintf("  Violent parent:     %.5f    %.5f    %.5f   %.5f\n",
            offspring_vf, offspring_vv, offspring_vp, mobil_violent))
cat(sprintf("  Peaceful parent:    %.5f    %.5f    %.5f   %.5f\n\n",
            offspring_pf, offspring_pv, offspring_pp, mobil_peaceful))

cat("MOBILIZATION RATIOS:\n")
cat(sprintf("  Fatal / Violent:    %.4f\n", ratio_fv))
cat(sprintf("  Violent / Peaceful: %.4f\n", ratio_vp))
cat(sprintf("  Fatal / Peaceful:   %.4f\n\n", ratio_fp))

# =============================================================================
# 8. DELTA METHOD FOR MOBILIZATION RATIO CIs
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  HYPOTHESIS TESTS                                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Store inference results
ci_fv_lower <- NA; ci_fv_upper <- NA; p_value_fv <- NA
ci_vp_lower <- NA; ci_vp_upper <- NA; p_value_vp <- NA
ci_fp_lower <- NA; ci_fp_upper <- NA; p_value_fp <- NA

if (!is.null(vcov_mat)) {

  # Function to compute log(mobilization ratio F/V) from params
  log_ratio_fv_func <- function(p) {
    alpha_ff <- exp(p[4]); alpha_fv <- exp(p[5]); alpha_fp <- exp(p[6])
    alpha_vf <- exp(p[7]); alpha_vv <- exp(p[8]); alpha_vp <- exp(p[9])
    beta <- exp(p[13])

    mobil_f <- (alpha_ff + alpha_fv + alpha_fp) / beta
    mobil_v <- (alpha_vf + alpha_vv + alpha_vp) / beta

    log(mobil_f / mobil_v)
  }

  # Function for V/P ratio
  log_ratio_vp_func <- function(p) {
    alpha_vf <- exp(p[7]); alpha_vv <- exp(p[8]); alpha_vp <- exp(p[9])
    alpha_pf <- exp(p[10]); alpha_pv <- exp(p[11]); alpha_pp <- exp(p[12])
    beta <- exp(p[13])

    mobil_v <- (alpha_vf + alpha_vv + alpha_vp) / beta
    mobil_p <- (alpha_pf + alpha_pv + alpha_pp) / beta

    log(mobil_v / mobil_p)
  }

  # TEST 1: Fatal vs Violent (H2)
  cat("TEST H2: Is mobilization ratio (Fatal/Violent) different from 1?\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat("  H0: Fatal and violent events have equal mobilization potential\n")
  cat("  H1: They differ in mobilization potential\n\n")

  gradient_fv <- grad(log_ratio_fv_func, params)
  var_log_fv <- as.numeric(t(gradient_fv) %*% vcov_mat %*% gradient_fv)
  se_log_fv <- sqrt(max(var_log_fv, 0))
  log_fv <- log(ratio_fv)

  if (se_log_fv > 0 && is.finite(se_log_fv)) {
    ci_fv_lower <- exp(log_fv - 1.96 * se_log_fv)
    ci_fv_upper <- exp(log_fv + 1.96 * se_log_fv)
    z_stat_fv <- log_fv / se_log_fv
    p_value_fv <- 2 * pnorm(-abs(z_stat_fv))

    cat(sprintf("  Mobilization ratio (F/V): %.4f\n", ratio_fv))
    cat(sprintf("  95%% CI: [%.4f, %.4f]\n", ci_fv_lower, ci_fv_upper))
    cat(sprintf("  z-statistic: %.4f, p-value: %.6f\n", z_stat_fv, p_value_fv))

    if (p_value_fv < 0.05 && ci_fv_lower > 1) {
      cat("  *** RESULT: Fatal events have HIGHER mobilization than violent ***\n\n")
    } else if (p_value_fv < 0.05 && ci_fv_upper < 1) {
      cat("  *** RESULT: Fatal events have LOWER mobilization than violent ***\n\n")
    } else {
      cat("  RESULT: No significant difference\n\n")
    }
  } else {
    cat("  Unable to compute test (SE = 0 or NA)\n\n")
  }

  # TEST 2: Violent vs Peaceful (H3)
  cat("TEST H3: Is mobilization ratio (Violent/Peaceful) different from 1?\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat("  H0: Violent and peaceful events have equal mobilization potential\n")
  cat("  H1: They differ in mobilization potential\n\n")

  gradient_vp <- grad(log_ratio_vp_func, params)
  var_log_vp <- as.numeric(t(gradient_vp) %*% vcov_mat %*% gradient_vp)
  se_log_vp <- sqrt(max(var_log_vp, 0))
  log_vp <- log(ratio_vp)

  if (se_log_vp > 0 && is.finite(se_log_vp)) {
    ci_vp_lower <- exp(log_vp - 1.96 * se_log_vp)
    ci_vp_upper <- exp(log_vp + 1.96 * se_log_vp)
    z_stat_vp <- log_vp / se_log_vp
    p_value_vp <- 2 * pnorm(-abs(z_stat_vp))

    cat(sprintf("  Mobilization ratio (V/P): %.4f\n", ratio_vp))
    cat(sprintf("  95%% CI: [%.4f, %.4f]\n", ci_vp_lower, ci_vp_upper))
    cat(sprintf("  z-statistic: %.4f, p-value: %.6f\n", z_stat_vp, p_value_vp))

    if (p_value_vp < 0.05 && ci_vp_lower > 1) {
      cat("  *** RESULT: Violent events have HIGHER mobilization than peaceful ***\n\n")
    } else if (p_value_vp < 0.05 && ci_vp_upper < 1) {
      cat("  *** RESULT: Violent events have LOWER mobilization than peaceful ***\n\n")
    } else {
      cat("  RESULT: No significant difference\n\n")
    }
  } else {
    cat("  Unable to compute test (SE = 0 or NA)\n\n")
  }

  # TEST 3: Check if escalation dominates (α_VV > α_VP) - H4
  cat("TEST H4: Does violence trigger violence more than peace? (α_VV vs α_VP)\n")
  cat("─────────────────────────────────────────────────────────────────\n")

  diff_log <- params[8] - params[9]  # log(α_VV) - log(α_VP)
  var_diff <- vcov_mat[8,8] + vcov_mat[9,9] - 2*vcov_mat[8,9]
  se_diff <- sqrt(max(var_diff, 0))

  if (se_diff > 0 && is.finite(se_diff)) {
    z_diff <- diff_log / se_diff
    p_diff <- 2 * pnorm(-abs(z_diff))

    cat(sprintf("  α_VV = %.6f, α_VP = %.6f\n", alpha_vv, alpha_vp))
    cat(sprintf("  Ratio α_VV/α_VP: %.4f\n", alpha_vv/alpha_vp))
    cat(sprintf("  z-statistic: %.4f, p-value: %.6f\n", z_diff, p_diff))

    if (p_diff < 0.05) {
      if (alpha_vv > alpha_vp) {
        cat("  *** RESULT: Violence ESCALATES (V→V > V→P) ***\n\n")
      } else {
        cat("  *** RESULT: Violence MOBILIZES (V→P > V→V) ***\n\n")
      }
    } else {
      cat("  RESULT: No significant difference between escalation and mobilization\n\n")
    }
  } else {
    cat("  Unable to compute test\n\n")
  }

} else {
  cat("Cannot perform hypothesis tests without standard errors\n\n")
}

# =============================================================================
# 9. FINAL SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FINAL SUMMARY                                               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("Model: Trivariate Hawkes Process (F/V/P)\n"))
cat(sprintf("Data: %d events (%d fatal, %d violent, %d peaceful) over %.0f days\n",
            n_total, n_f, n_v, n_p, T_max))
cat(sprintf("Log-likelihood: %.2f\n", best_ll))
cat(sprintf("Parameters: 13 (3 baselines + 9 triggering + 1 decay)\n"))
cat(sprintf("Best starting point: %s\n\n", best_start))

cat("TRIGGERING MATRIX (expected offspring):\n")
cat("┌─────────────────┬──────────┬──────────┬──────────┬──────────┐\n")
cat("│                 │ → Fatal  │ → Violent│ → Peaceful│  TOTAL  │\n")
cat("├─────────────────┼──────────┼──────────┼──────────┼──────────┤\n")
cat(sprintf("│ Fatal parent    │ %.5f  │ %.5f  │ %.5f  │ %.5f  │\n",
            offspring_ff, offspring_fv, offspring_fp, mobil_fatal))
cat(sprintf("│ Violent parent  │ %.5f  │ %.5f  │ %.5f  │ %.5f  │\n",
            offspring_vf, offspring_vv, offspring_vp, mobil_violent))
cat(sprintf("│ Peaceful parent │ %.5f  │ %.5f  │ %.5f  │ %.5f  │\n",
            offspring_pf, offspring_pv, offspring_pp, mobil_peaceful))
cat("└─────────────────┴──────────┴──────────┴──────────┴──────────┘\n\n")

cat("KEY MOBILIZATION RATIOS:\n")
cat(sprintf("  Fatal/Violent:    %.4f", ratio_fv))
if (!is.na(ci_fv_lower)) cat(sprintf(" [95%% CI: %.2f-%.2f]", ci_fv_lower, ci_fv_upper))
if (!is.na(p_value_fv)) cat(sprintf(" p=%.4f", p_value_fv))
cat("\n")

cat(sprintf("  Violent/Peaceful: %.4f", ratio_vp))
if (!is.na(ci_vp_lower)) cat(sprintf(" [95%% CI: %.2f-%.2f]", ci_vp_lower, ci_vp_upper))
if (!is.na(p_value_vp)) cat(sprintf(" p=%.4f", p_value_vp))
cat("\n\n")

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================

cat("=== SAVING RESULTS ===\n")

results <- list(
  # Model info
  model_name = "Trivariate Hawkes Process (F/V/P) - Constant Baseline",
  best_starting_point = best_start,
  all_results = all_results,

  # Data info
  n_fatal = n_f,
  n_violent = n_v,
  n_peaceful = n_p,
  n_total = n_total,
  T_max = T_max,

  # Parameter estimates (log scale)
  params = params,
  param_names = param_names,

  # Parameter estimates (natural scale)
  mu_f = mu_f, mu_v = mu_v, mu_p = mu_p,
  alpha_ff = alpha_ff, alpha_fv = alpha_fv, alpha_fp = alpha_fp,
  alpha_vf = alpha_vf, alpha_vv = alpha_vv, alpha_vp = alpha_vp,
  alpha_pf = alpha_pf, alpha_pv = alpha_pv, alpha_pp = alpha_pp,
  beta = beta,

  # Offspring (expected children per parent)
  offspring_ff = offspring_ff, offspring_fv = offspring_fv, offspring_fp = offspring_fp,
  offspring_vf = offspring_vf, offspring_vv = offspring_vv, offspring_vp = offspring_vp,
  offspring_pf = offspring_pf, offspring_pv = offspring_pv, offspring_pp = offspring_pp,

  # Mobilization
  mobil_fatal = mobil_fatal,
  mobil_violent = mobil_violent,
  mobil_peaceful = mobil_peaceful,
  ratio_fv = ratio_fv,
  ratio_vp = ratio_vp,
  ratio_fp = ratio_fp,

  # Inference - ratio F/V
  ci_fv_lower = ci_fv_lower, ci_fv_upper = ci_fv_upper, p_value_fv = p_value_fv,
  # Inference - ratio V/P
  ci_vp_lower = ci_vp_lower, ci_vp_upper = ci_vp_upper, p_value_vp = p_value_vp,

  # Variance-covariance
  se_log = se_log,
  vcov = vcov_mat,

  # Fit statistics
  loglik = best_ll,
  n_params = 13,
  aic = -2 * best_ll + 2 * 13,
  bic = -2 * best_ll + 13 * log(n_total),

  # Runtime
  total_runtime = total_runtime
)

saveRDS(results, "trivariate_hawkes_constant.rds")
cat("✓ Saved: trivariate_hawkes_constant.rds\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  TRIVARIATE HAWKES MODEL COMPLETE                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
