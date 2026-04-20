################################################################################
#           BIVARIATE HAWKES MODEL - FULL IMPLEMENTATION
#           Multi-starting point optimization with hypothesis testing
#
#           Model: Two interacting point processes (violent & peaceful)
#           λ_V(t) = μ_V + Σ α_VV·g(t-t_j) + Σ α_PV·g(t-t_k)
#           λ_P(t) = μ_P + Σ α_VP·g(t-t_j) + Σ α_PP·g(t-t_k)
#           where g(s) = exp(-β·s)
################################################################################

library(tidyverse)
library(numDeriv)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   BIVARIATE HAWKES MODEL - FULL IMPLEMENTATION                          ║\n")
cat("║   Multi-starting point optimization with hypothesis testing              ║\n")
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

protests <- readRDS("protests_prepared.rds")

protests <- protests %>%
  mutate(
    is_violent = sub_event_type %in% c("Mob violence", "Violent demonstration") |
                 sub_event_type == "Excessive force against protesters"
  ) %>%
  arrange(event_date)

# Convert to numeric time (days since start)
start_date <- min(protests$event_date)
protests$time <- as.numeric(protests$event_date - start_date)
T_max <- max(protests$time) + 1

# Separate streams
v_times <- protests$time[protests$is_violent]
p_times <- protests$time[!protests$is_violent]

n_v <- length(v_times)
n_p <- length(p_times)
n_total <- n_v + n_p

cat(sprintf("Loaded %d events\n", nrow(protests)))
cat(sprintf("  Violent events: %d (%.1f%%)\n", n_v, 100*n_v/n_total))
cat(sprintf("  Peaceful events: %d (%.1f%%)\n", n_p, 100*n_p/n_total))
cat(sprintf("  Observation period: %.0f days (%.1f years)\n\n", T_max, T_max/365.25))

# =============================================================================
# 2. BIVARIATE HAWKES LOG-LIKELIHOOD (EFFICIENT RECURSIVE COMPUTATION)
# =============================================================================

cat("=== SETTING UP LIKELIHOOD FUNCTION ===\n")

# Combine and sort all events for recursive computation
all_times <- c(v_times, p_times)
all_types <- c(rep(1, n_v), rep(0, n_p))  # 1 = violent, 0 = peaceful
ord <- order(all_times)
all_times <- all_times[ord]
all_types <- all_types[ord]
n <- length(all_times)

cat(sprintf("Combined %d events for recursive computation\n\n", n))

# Pre-compute time differences for efficiency
time_diffs <- diff(all_times)

bivariate_loglik <- function(params) {
  # Parameters (log scale for unconstrained optimization)
  mu_v <- exp(params[1])
  mu_p <- exp(params[2])
  alpha_vv <- exp(params[3])
  alpha_vp <- exp(params[4])
  alpha_pv <- exp(params[5])
  alpha_pp <- exp(params[6])
  beta <- exp(params[7])

  # Check for valid parameters
  if (any(!is.finite(c(mu_v, mu_p, alpha_vv, alpha_vp, alpha_pv, alpha_pp, beta)))) {
    return(-1e10)
  }
  if (beta < 1e-6 || beta > 100) return(-1e10)

  # Recursive computation of triggering sums
  # R_v[i] = sum of exp(-β*(t_i - t_j)) for violent j < i
  # R_p[i] = sum of exp(-β*(t_i - t_j)) for peaceful j < i

  R_v <- numeric(n)
  R_p <- numeric(n)

  for (i in 2:n) {
    decay <- exp(-beta * time_diffs[i-1])

    # Limit decay to avoid numerical issues
    if (decay < 1e-100) decay <- 0

    if (all_types[i-1] == 1) {
      # Previous event was violent
      R_v[i] <- decay * (1 + R_v[i-1])
      R_p[i] <- decay * R_p[i-1]
    } else {
      # Previous event was peaceful
      R_v[i] <- decay * R_v[i-1]
      R_p[i] <- decay * (1 + R_p[i-1])
    }
  }

  # Log-likelihood: sum of log(λ) at event times
  loglik <- 0

  for (i in 1:n) {
    if (all_types[i] == 1) {
      # Violent event
      lambda_i <- mu_v + alpha_vv * R_v[i] + alpha_pv * R_p[i]
    } else {
      # Peaceful event
      lambda_i <- mu_p + alpha_vp * R_v[i] + alpha_pp * R_p[i]
    }

    if (lambda_i <= 1e-10) return(-1e10)
    loglik <- loglik + log(lambda_i)
  }

  # Integral terms (compensator)
  # For exponential kernel: ∫_{t_j}^{T} α·exp(-β·(t-t_j)) dt = (α/β)·(1 - exp(-β·(T-t_j)))

  integral_v <- mu_v * T_max
  integral_p <- mu_p * T_max

  for (j in 1:n) {
    contrib <- (1 - exp(-beta * (T_max - all_times[j]))) / beta
    if (all_types[j] == 1) {
      # Violent parent
      integral_v <- integral_v + alpha_vv * contrib
      integral_p <- integral_p + alpha_vp * contrib
    } else {
      # Peaceful parent
      integral_v <- integral_v + alpha_pv * contrib
      integral_p <- integral_p + alpha_pp * contrib
    }
  }

  loglik <- loglik - integral_v - integral_p

  if (!is.finite(loglik)) return(-1e10)

  return(loglik)
}

# Wrapper for optimization (minimization)
neg_loglik <- function(params) {
  -bivariate_loglik(params)
}

cat("✓ Likelihood function ready\n\n")

# =============================================================================
# 3. MULTI-STARTING POINT OPTIMIZATION
# =============================================================================

cat("=== MULTI-STARTING POINT OPTIMIZATION ===\n\n")

# Empirical estimates for starting values
emp_rate_v <- n_v / T_max
emp_rate_p <- n_p / T_max

cat(sprintf("Empirical rates: μ_V ≈ %.4f, μ_P ≈ %.4f events/day\n\n", emp_rate_v, emp_rate_p))

# Define starting points
starting_points <- list(

  # 1. Data-driven
  list(
    name = "Data-driven",
    params = c(
      log(emp_rate_v * 0.5),   # log_mu_v
      log(emp_rate_p * 0.5),   # log_mu_p
      log(0.05),               # log_alpha_vv
      log(0.05),               # log_alpha_vp
      log(0.05),               # log_alpha_pv
      log(0.05),               # log_alpha_pp
      log(0.1)                 # log_beta
    )
  ),

  # 2. Conservative (lower triggering, faster decay)
  list(
    name = "Conservative",
    params = c(
      log(emp_rate_v * 0.7),
      log(emp_rate_p * 0.7),
      log(0.02),
      log(0.02),
      log(0.02),
      log(0.02),
      log(0.2)
    )
  ),

  # 3. Strong triggering (higher α, slower decay)
  list(
    name = "Strong-triggering",
    params = c(
      log(emp_rate_v * 0.3),
      log(emp_rate_p * 0.3),
      log(0.10),
      log(0.10),
      log(0.10),
      log(0.10),
      log(0.05)
    )
  ),

  # 4. Asymmetric (V→V stronger than P→P)
  list(
    name = "Asymmetric-VV",
    params = c(
      log(emp_rate_v * 0.5),
      log(emp_rate_p * 0.5),
      log(0.10),               # α_VV higher
      log(0.03),
      log(0.03),
      log(0.05),               # α_PP lower
      log(0.15)
    )
  ),

  # 5. Asymmetric (P→P stronger than V→V)
  list(
    name = "Asymmetric-PP",
    params = c(
      log(emp_rate_v * 0.5),
      log(emp_rate_p * 0.5),
      log(0.03),               # α_VV lower
      log(0.05),
      log(0.05),
      log(0.10),               # α_PP higher
      log(0.15)
    )
  )
)

# Parameter names
param_names <- c("log_mu_v", "log_mu_p", "log_alpha_vv", "log_alpha_vp",
                 "log_alpha_pv", "log_alpha_pp", "log_beta")

# Bounds for L-BFGS-B
lower_bounds <- c(-15, -15, -15, -15, -15, -15, -5)
upper_bounds <- c(5, 5, 5, 5, 5, 5, 3)

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

  # Display starting values
  mu_v_start <- exp(sp$params[1])
  mu_p_start <- exp(sp$params[2])
  alpha_vv_start <- exp(sp$params[3])
  beta_start <- exp(sp$params[7])

  cat(sprintf("Starting values: μ_V=%.4f, μ_P=%.4f, α_VV=%.4f, β=%.4f\n",
              mu_v_start, mu_p_start, alpha_vv_start, beta_start))

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

    # Extract final parameter values
    mu_v_final <- exp(result$par[1])
    mu_p_final <- exp(result$par[2])
    alpha_vv_final <- exp(result$par[3])
    alpha_vp_final <- exp(result$par[4])
    alpha_pv_final <- exp(result$par[5])
    alpha_pp_final <- exp(result$par[6])
    beta_final <- exp(result$par[7])

    cat(sprintf("\n✓ Completed in %.1f minutes\n", runtime))
    cat(sprintf("  Log-likelihood: %.2f\n", ll))
    cat(sprintf("  Convergence: %d\n", result$convergence))
    cat(sprintf("  Final: μ_V=%.4f, μ_P=%.4f, β=%.4f\n",
                mu_v_final, mu_p_final, beta_final))
    cat(sprintf("  Triggering: α_VV=%.5f, α_VP=%.5f, α_PV=%.5f, α_PP=%.5f\n\n",
                alpha_vv_final, alpha_vp_final, alpha_pv_final, alpha_pp_final))

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
  cat(sprintf("  %d. %-20s  LL: %10.2f  Conv: %d  Time: %.1f min%s\n",
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
mu_v <- exp(params[1])
mu_p <- exp(params[2])
alpha_vv <- exp(params[3])
alpha_vp <- exp(params[4])
alpha_pv <- exp(params[5])
alpha_pp <- exp(params[6])
beta <- exp(params[7])

cat("BASELINE RATES (events/day):\n")
cat(sprintf("  μ_V (violent):  %.6f\n", mu_v))
cat(sprintf("  μ_P (peaceful): %.6f\n\n", mu_p))

cat("TRIGGERING MATRIX (α coefficients):\n")
cat("                      → Violent     → Peaceful\n")
cat(sprintf("  Violent parent:     %.6f     %.6f\n", alpha_vv, alpha_vp))
cat(sprintf("  Peaceful parent:    %.6f     %.6f\n\n", alpha_pv, alpha_pp))

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

  for (i in 1:7) {
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
  se_log <- rep(NA, 7)
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
offspring_vv <- alpha_vv / beta
offspring_vp <- alpha_vp / beta
offspring_pv <- alpha_pv / beta
offspring_pp <- alpha_pp / beta

# Total mobilization
mobil_violent <- offspring_vv + offspring_vp
mobil_peaceful <- offspring_pv + offspring_pp
mobil_ratio <- mobil_violent / mobil_peaceful

cat("EXPECTED OFFSPRING PER PARENT EVENT:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("                        Violent     Peaceful    TOTAL\n")
cat("                        offspring   offspring   MOBILIZATION\n")
cat(sprintf("  From VIOLENT event:   %.5f     %.5f     %.5f\n",
            offspring_vv, offspring_vp, mobil_violent))
cat(sprintf("  From PEACEFUL event:  %.5f     %.5f     %.5f\n\n",
            offspring_pv, offspring_pp, mobil_peaceful))

cat(sprintf("MOBILIZATION RATIO (Violent / Peaceful): %.4f\n\n", mobil_ratio))

# =============================================================================
# 8. DELTA METHOD FOR MOBILIZATION RATIO CI
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  HYPOTHESIS TESTS                                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

if (!is.null(vcov_mat)) {

  # Function to compute log(mobilization ratio) from params
  log_mobil_ratio_func <- function(p) {
    alpha_vv <- exp(p[3])
    alpha_vp <- exp(p[4])
    alpha_pv <- exp(p[5])
    alpha_pp <- exp(p[6])
    beta <- exp(p[7])

    mobil_v <- (alpha_vv + alpha_vp) / beta
    mobil_p <- (alpha_pv + alpha_pp) / beta

    log(mobil_v / mobil_p)
  }

  # Compute gradient numerically
  gradient <- grad(log_mobil_ratio_func, params)

  # SE of log(ratio) using delta method
  var_log_ratio <- as.numeric(t(gradient) %*% vcov_mat %*% gradient)
  se_log_ratio <- sqrt(max(var_log_ratio, 0))

  log_ratio <- log(mobil_ratio)

  # 95% CI for log(ratio), then transform
  ci_log_lower <- log_ratio - 1.96 * se_log_ratio
  ci_log_upper <- log_ratio + 1.96 * se_log_ratio

  ci_lower <- exp(ci_log_lower)
  ci_upper <- exp(ci_log_upper)

  cat("TEST 1: Is mobilization ratio different from 1?\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat("  H0: Violent and peaceful events have equal mobilization potential\n")
  cat("  H1: They differ in mobilization potential\n\n")

  cat(sprintf("  Mobilization ratio: %.4f\n", mobil_ratio))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", ci_lower, ci_upper))

  if (se_log_ratio > 0 && is.finite(se_log_ratio)) {
    z_stat <- log_ratio / se_log_ratio
    p_value <- 2 * pnorm(-abs(z_stat))

    cat(sprintf("  z-statistic: %.4f\n", z_stat))
    cat(sprintf("  p-value: %.6f\n", p_value))

    if (p_value < 0.05) {
      if (ci_lower > 1) {
        cat("\n  *** RESULT: REJECT H0 (p < 0.05) ***\n")
        cat("  *** Violent events have SIGNIFICANTLY HIGHER mobilization potential ***\n")
      } else if (ci_upper < 1) {
        cat("\n  *** RESULT: REJECT H0 (p < 0.05) ***\n")
        cat("  *** Violent events have SIGNIFICANTLY LOWER mobilization potential ***\n")
      } else {
        cat("\n  Note: p < 0.05 but CI includes 1 (borderline significance)\n")
      }
    } else {
      cat("\n  RESULT: Cannot reject H0 (p >= 0.05)\n")
      cat("  No significant difference in mobilization potential\n")
    }
  } else {
    p_value <- NA
    cat("  Unable to compute z-test (SE = 0 or NA)\n")
  }

  # TEST 2: Individual triggering coefficients
  cat("\n\nTEST 2: Are individual triggering coefficients significant?\n")
  cat("─────────────────────────────────────────────────────────────────\n")

  coef_names <- c("α_VV (V→V)", "α_VP (V→P)", "α_PV (P→V)", "α_PP (P→P)")
  coef_values <- c(alpha_vv, alpha_vp, alpha_pv, alpha_pp)

  for (i in 3:6) {
    if (!is.na(se_log[i]) && se_log[i] > 0 && is.finite(se_log[i])) {
      z <- params[i] / se_log[i]
      p <- 2 * pnorm(-abs(z))
      sig <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
      cat(sprintf("  %s: %.6f (SE: %.4f, z: %.2f, p: %.4f) %s\n",
                  coef_names[i-2], coef_values[i-2], se_log[i], z, p, sig))
    } else {
      cat(sprintf("  %s: %.6f (SE: NA)\n", coef_names[i-2], coef_values[i-2]))
    }
  }

  # TEST 3: Is α_VV different from α_PP?
  cat("\n\nTEST 3: Is α_VV different from α_PP?\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat("  H0: Violence-to-violence = Peaceful-to-peaceful triggering\n\n")

  diff_log <- params[3] - params[6]  # log(α_VV) - log(α_PP)
  var_diff <- vcov_mat[3,3] + vcov_mat[6,6] - 2*vcov_mat[3,6]
  se_diff <- sqrt(max(var_diff, 0))

  if (se_diff > 0 && is.finite(se_diff)) {
    z_diff <- diff_log / se_diff
    p_diff <- 2 * pnorm(-abs(z_diff))

    cat(sprintf("  α_VV = %.6f, α_PP = %.6f\n", alpha_vv, alpha_pp))
    cat(sprintf("  Ratio α_VV/α_PP: %.4f\n", alpha_vv/alpha_pp))
    cat(sprintf("  z-statistic: %.4f\n", z_diff))
    cat(sprintf("  p-value: %.6f\n", p_diff))

    if (p_diff < 0.05) {
      cat("\n  *** α_VV IS significantly different from α_PP (p < 0.05) ***\n")
    } else {
      cat("\n  Cannot reject H0: α_VV and α_PP are not significantly different\n")
    }
  } else {
    p_diff <- NA
    cat("  Unable to compute test (SE = 0 or NA)\n")
  }

} else {
  p_value <- NA
  p_diff <- NA
  ci_lower <- NA
  ci_upper <- NA
  cat("Cannot perform hypothesis tests without standard errors\n")
}

# =============================================================================
# 9. SUMMARY TABLE
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FINAL SUMMARY                                               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("Model: Bivariate Hawkes Process\n"))
cat(sprintf("Data: %d events (%d violent, %d peaceful) over %.0f days\n",
            n_total, n_v, n_p, T_max))
cat(sprintf("Log-likelihood: %.2f\n", best_ll))
cat(sprintf("Best starting point: %s\n\n", best_start))

cat("TRIGGERING MATRIX (expected offspring):\n")
cat("┌─────────────────┬────────────────┬────────────────┬──────────────┐\n")
cat("│                 │  → Violent     │  → Peaceful    │    TOTAL     │\n")
cat("├─────────────────┼────────────────┼────────────────┼──────────────┤\n")
cat(sprintf("│ Violent parent  │    %.5f     │    %.5f     │    %.5f   │\n",
            offspring_vv, offspring_vp, mobil_violent))
cat(sprintf("│ Peaceful parent │    %.5f     │    %.5f     │    %.5f   │\n",
            offspring_pv, offspring_pp, mobil_peaceful))
cat("└─────────────────┴────────────────┴────────────────┴──────────────┘\n\n")

cat("KEY RESULT:\n")
cat(sprintf("  Mobilization ratio: %.4f", mobil_ratio))
if (!is.na(ci_lower) && !is.na(ci_upper)) {
  cat(sprintf(" (95%% CI: [%.4f, %.4f])", ci_lower, ci_upper))
}
if (!is.na(p_value)) {
  cat(sprintf(", p = %.4f", p_value))
}
cat("\n\n")

if (!is.na(p_value) && p_value < 0.05) {
  if (!is.na(ci_lower) && ci_lower > 1) {
    cat("══════════════════════════════════════════════════════════════════\n")
    cat("  CONCLUSION: Violent events trigger SIGNIFICANTLY MORE\n")
    cat("              total follow-on events than peaceful events\n")
    cat("══════════════════════════════════════════════════════════════════\n")
  } else if (!is.na(ci_upper) && ci_upper < 1) {
    cat("══════════════════════════════════════════════════════════════════\n")
    cat("  CONCLUSION: Violent events trigger SIGNIFICANTLY FEWER\n")
    cat("              total follow-on events than peaceful events\n")
    cat("══════════════════════════════════════════════════════════════════\n")
  }
} else {
  cat("══════════════════════════════════════════════════════════════════\n")
  cat("  CONCLUSION: No significant difference in mobilization potential\n")
  cat("              between violent and peaceful events\n")
  cat("══════════════════════════════════════════════════════════════════\n")
}

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================

cat("\n=== SAVING RESULTS ===\n")

results <- list(
  # Model info
  model_name = "Bivariate Hawkes Process",
  best_starting_point = best_start,
  all_results = all_results,

  # Data info
  n_violent = n_v,
  n_peaceful = n_p,
  n_total = n_total,
  T_max = T_max,

  # Parameter estimates (log scale)
  params = params,
  param_names = param_names,

  # Parameter estimates (natural scale)
  mu_v = mu_v, mu_p = mu_p,
  alpha_vv = alpha_vv, alpha_vp = alpha_vp,
  alpha_pv = alpha_pv, alpha_pp = alpha_pp,
  beta = beta,

  # Offspring
  offspring_vv = offspring_vv, offspring_vp = offspring_vp,
  offspring_pv = offspring_pv, offspring_pp = offspring_pp,
  mobil_violent = mobil_violent,
  mobil_peaceful = mobil_peaceful,
  mobil_ratio = mobil_ratio,

  # Inference
  se_log = se_log,
  vcov = vcov_mat,
  ci_lower = ci_lower,
  ci_upper = ci_upper,
  p_value = p_value,

  # Fit statistics
  loglik = best_ll,
  n_params = 7,
  aic = -2 * best_ll + 2 * 7,
  bic = -2 * best_ll + 7 * log(n_total),

  # Runtime
  total_runtime = total_runtime
)

saveRDS(results, "bivariate_hawkes_full.rds")
cat("✓ Saved: bivariate_hawkes_full.rds\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  BIVARIATE HAWKES MODEL COMPLETE                             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
