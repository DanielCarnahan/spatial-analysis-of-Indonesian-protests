################################################################################
#           BIVARIATE HAWKES MODEL WITH HYPOTHESIS TESTING
#           Standard errors via Hessian, tests for mobilization difference
################################################################################

library(tidyverse)
library(numDeriv)  # For Hessian computation

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   BIVARIATE HAWKES WITH INFERENCE                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Load and prepare data
protests <- readRDS("protests_prepared.rds")
protests <- protests %>%
  mutate(
    is_violent = sub_event_type %in% c("Mob violence", "Violent demonstration") |
                 sub_event_type == "Excessive force against protesters",
    year = lubridate::year(event_date)
  ) %>%
  arrange(event_date)

# Use 2019-2021 (same period as fast model that worked)
subset_data <- protests %>%
  filter(year >= 2019 & year <= 2021)

subset_data$time <- as.numeric(subset_data$event_date - min(subset_data$event_date))
T_max <- max(subset_data$time) + 1

v_times <- subset_data$time[subset_data$is_violent]
p_times <- subset_data$time[!subset_data$is_violent]
n_v <- length(v_times)
n_p <- length(p_times)

cat(sprintf("Data: %d violent, %d peaceful events over %.0f days\n\n", n_v, n_p, T_max))

# =============================================================================
# BIVARIATE HAWKES LOG-LIKELIHOOD
# =============================================================================

bivariate_loglik <- function(params, v_times, p_times, T_max) {

  mu_v <- exp(params[1])
  mu_p <- exp(params[2])
  alpha_vv <- exp(params[3])
  alpha_vp <- exp(params[4])
  alpha_pv <- exp(params[5])
  alpha_pp <- exp(params[6])
  beta <- exp(params[7])

  n_v <- length(v_times)
  n_p <- length(p_times)

  # Combine and sort
  all_times <- c(v_times, p_times)
  all_types <- c(rep(1, n_v), rep(0, n_p))
  ord <- order(all_times)
  all_times <- all_times[ord]
  all_types <- all_types[ord]
  n <- length(all_times)

  # Recursive computation of triggering sums
  R_v <- R_p <- numeric(n)
  for (i in 2:n) {
    decay <- exp(-beta * (all_times[i] - all_times[i-1]))
    if (all_types[i-1] == 1) {
      R_v[i] <- decay * (1 + R_v[i-1])
      R_p[i] <- decay * R_p[i-1]
    } else {
      R_v[i] <- decay * R_v[i-1]
      R_p[i] <- decay * (1 + R_p[i-1])
    }
  }

  # Log-likelihood: sum of log(lambda) at event times
  loglik <- 0
  for (i in 1:n) {
    if (all_types[i] == 1) {
      lambda_i <- mu_v + alpha_vv * R_v[i] + alpha_pv * R_p[i]
    } else {
      lambda_i <- mu_p + alpha_vp * R_v[i] + alpha_pp * R_p[i]
    }
    if (lambda_i <= 1e-10) return(-1e10)
    loglik <- loglik + log(lambda_i)
  }

  # Integral terms
  integral_v <- mu_v * T_max
  integral_p <- mu_p * T_max

  for (j in 1:n) {
    contrib <- (1 - exp(-beta * (T_max - all_times[j]))) / beta
    if (all_types[j] == 1) {
      integral_v <- integral_v + alpha_vv * contrib
      integral_p <- integral_p + alpha_vp * contrib
    } else {
      integral_v <- integral_v + alpha_pv * contrib
      integral_p <- integral_p + alpha_pp * contrib
    }
  }

  loglik <- loglik - integral_v - integral_p
  return(loglik)
}

# Wrapper for optimization
neg_loglik <- function(params) {
  ll <- bivariate_loglik(params, v_times, p_times, T_max)
  if (is.na(ll) || is.infinite(ll)) return(1e10)
  return(-ll)
}

# =============================================================================
# OPTIMIZATION
# =============================================================================

cat("Optimizing bivariate Hawkes model...\n")

# Starting values (log scale) - based on successful estimates from fast version
start_params <- c(
  log(n_v / T_max * 0.5),   # log_mu_v
  log(n_p / T_max * 0.5),   # log_mu_p
  log(0.05),                 # log_alpha_vv
  log(0.05),                 # log_alpha_vp
  log(0.05),                 # log_alpha_pv
  log(0.05),                 # log_alpha_pp
  log(0.1)                   # log_beta
)

# Multiple starting points
best_result <- NULL
best_ll <- Inf

for (trial in 1:10) {
  init <- start_params + rnorm(7, 0, 0.3)

  # Try BFGS first
  result <- tryCatch({
    optim(init, neg_loglik, method = "BFGS",
          control = list(maxit = 1000), hessian = FALSE)
  }, error = function(e) NULL)

  if (!is.null(result) && result$value < best_ll) {
    best_ll <- result$value
    best_result <- result
    cat(sprintf("  Trial %d (BFGS): LL = %.1f\n", trial, -result$value))
  }

  # Try Nelder-Mead as fallback
  result_nm <- tryCatch({
    optim(init, neg_loglik, method = "Nelder-Mead",
          control = list(maxit = 2000))
  }, error = function(e) NULL)

  if (!is.null(result_nm) && result_nm$value < best_ll) {
    best_ll <- result_nm$value
    best_result <- result_nm
    cat(sprintf("  Trial %d (NM): LL = %.1f\n", trial, -result_nm$value))
  }
}

if (is.null(best_result)) {
  stop("Optimization failed!")
}

cat(sprintf("Log-likelihood: %.2f\n", -best_result$value))

# =============================================================================
# COMPUTE HESSIAN FOR STANDARD ERRORS
# =============================================================================

cat("Computing Hessian for standard errors...\n")

# Compute Hessian at MLE
hessian_mat <- hessian(neg_loglik, best_result$par)

# Variance-covariance matrix (inverse of Hessian)
vcov_mat <- tryCatch({
  solve(hessian_mat)
}, error = function(e) {
  cat("Warning: Hessian not invertible, using pseudo-inverse\n")
  MASS::ginv(hessian_mat)
})

# Standard errors on log scale
se_log <- sqrt(diag(vcov_mat))

# =============================================================================
# PARAMETER ESTIMATES WITH STANDARD ERRORS
# =============================================================================

params <- best_result$par
param_names <- c("log_mu_v", "log_mu_p", "log_alpha_vv", "log_alpha_vp",
                 "log_alpha_pv", "log_alpha_pp", "log_beta")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   PARAMETER ESTIMATES (Log Scale)                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

for (i in 1:7) {
  z <- params[i] / se_log[i]
  p <- 2 * pnorm(-abs(z))
  cat(sprintf("  %s: %.4f (SE: %.4f, z: %.2f, p: %.4f)\n",
              param_names[i], params[i], se_log[i], z, p))
}

# Transform to natural scale
mu_v <- exp(params[1])
mu_p <- exp(params[2])
alpha_vv <- exp(params[3])
alpha_vp <- exp(params[4])
alpha_pv <- exp(params[5])
alpha_pp <- exp(params[6])
beta <- exp(params[7])

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   TRIGGERING MATRIX (Natural Scale)                          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("BASELINE RATES:\n")
cat(sprintf("  μ_V = %.4f events/day\n", mu_v))
cat(sprintf("  μ_P = %.4f events/day\n", mu_p))

cat("\nTRIGGERING COEFFICIENTS (α):\n")
cat("                    → Violent    → Peaceful\n")
cat(sprintf("  Violent parent:     %.5f      %.5f\n", alpha_vv, alpha_vp))
cat(sprintf("  Peaceful parent:    %.5f      %.5f\n", alpha_pv, alpha_pp))

cat(sprintf("\nDecay rate: β = %.4f (half-life = %.1f days)\n", beta, log(2)/beta))

# =============================================================================
# MOBILIZATION ESTIMATES WITH CONFIDENCE INTERVALS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   MOBILIZATION POTENTIAL WITH INFERENCE                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Expected offspring
off_vv <- alpha_vv / beta
off_vp <- alpha_vp / beta
off_pv <- alpha_pv / beta
off_pp <- alpha_pp / beta

mobil_v <- off_vv + off_vp
mobil_p <- off_pv + off_pp
mobil_ratio <- mobil_v / mobil_p

cat("EXPECTED OFFSPRING (point estimates):\n")
cat(sprintf("  From VIOLENT:  %.4f violent + %.4f peaceful = %.4f total\n",
            off_vv, off_vp, mobil_v))
cat(sprintf("  From PEACEFUL: %.4f violent + %.4f peaceful = %.4f total\n",
            off_pv, off_pp, mobil_p))
cat(sprintf("\n  MOBILIZATION RATIO: %.3f\n", mobil_ratio))

# =============================================================================
# DELTA METHOD FOR CONFIDENCE INTERVALS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   HYPOTHESIS TESTS                                           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Use delta method to get SE of mobilization ratio
# mobil_v = (exp(log_alpha_vv) + exp(log_alpha_vp)) / exp(log_beta)
# mobil_p = (exp(log_alpha_pv) + exp(log_alpha_pp)) / exp(log_beta)
# ratio = mobil_v / mobil_p

# Gradient of log(ratio) with respect to log-scale parameters
# log(ratio) = log(mobil_v) - log(mobil_p)
#            = log(exp(log_alpha_vv) + exp(log_alpha_vp)) - log_beta
#              - log(exp(log_alpha_pv) + exp(log_alpha_pp)) + log_beta
#            = log(exp(log_alpha_vv) + exp(log_alpha_vp))
#              - log(exp(log_alpha_pv) + exp(log_alpha_pp))

log_ratio <- log(mobil_ratio)

# Numerical gradient
grad_log_ratio <- function(params) {
  alpha_vv <- exp(params[3])
  alpha_vp <- exp(params[4])
  alpha_pv <- exp(params[5])
  alpha_pp <- exp(params[6])
  beta <- exp(params[7])

  mobil_v <- (alpha_vv + alpha_vp) / beta
  mobil_p <- (alpha_pv + alpha_pp) / beta

  log(mobil_v / mobil_p)
}

gradient <- grad(grad_log_ratio, params)

# SE of log(ratio) using delta method
var_log_ratio <- t(gradient) %*% vcov_mat %*% gradient
se_log_ratio <- sqrt(var_log_ratio)

# 95% CI for log(ratio)
ci_log_lower <- log_ratio - 1.96 * se_log_ratio
ci_log_upper <- log_ratio + 1.96 * se_log_ratio

# Transform to ratio scale
ci_lower <- exp(ci_log_lower)
ci_upper <- exp(ci_log_upper)

cat("TEST 1: Is mobilization ratio different from 1?\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Mobilization ratio: %.3f\n", mobil_ratio))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n", ci_lower, ci_upper))

# Z-test for log(ratio) = 0 (i.e., ratio = 1)
z_stat <- log_ratio / se_log_ratio
p_value <- 2 * pnorm(-abs(z_stat))

cat(sprintf("  z-statistic: %.3f\n", z_stat))
cat(sprintf("  p-value: %.6f\n", p_value))

if (p_value < 0.05) {
  cat("\n  *** REJECT H0: Mobilization ratio IS significantly different from 1 ***\n")
} else {
  cat("\n  Cannot reject H0: No significant difference in mobilization\n")
}

# =============================================================================
# TEST 2: IS VIOLENCE → VIOLENCE DIFFERENT FROM PEACEFUL → PEACEFUL?
# =============================================================================

cat("\n")
cat("TEST 2: Is α_VV different from α_PP?\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Test log(alpha_vv) - log(alpha_pp) = 0
diff_log <- params[3] - params[6]  # log_alpha_vv - log_alpha_pp
var_diff <- vcov_mat[3,3] + vcov_mat[6,6] - 2*vcov_mat[3,6]
se_diff <- sqrt(var_diff)

z_diff <- diff_log / se_diff
p_diff <- 2 * pnorm(-abs(z_diff))

cat(sprintf("  α_VV = %.5f, α_PP = %.5f\n", alpha_vv, alpha_pp))
cat(sprintf("  Ratio α_VV/α_PP: %.3f\n", alpha_vv/alpha_pp))
cat(sprintf("  z-statistic: %.3f\n", z_diff))
cat(sprintf("  p-value: %.6f\n", p_diff))

if (p_diff < 0.05) {
  cat("\n  *** α_VV IS significantly different from α_PP ***\n")
}

# =============================================================================
# TEST 3: IS VIOLENCE → PEACEFUL (SPILLOVER) SIGNIFICANT?
# =============================================================================

cat("\n")
cat("TEST 3: Is α_VP (violence → peaceful spillover) > 0?\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Already have z-test from individual parameter
z_vp <- params[4] / se_log[4]
p_vp <- 2 * pnorm(-abs(z_vp))

cat(sprintf("  α_VP = %.5f\n", alpha_vp))
cat(sprintf("  z-statistic: %.3f\n", z_vp))
cat(sprintf("  p-value: %.6f\n", p_vp))

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   SUMMARY: BIVARIATE HAWKES RESULTS                          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("TRIGGERING MATRIX:\n")
cat("┌─────────────────┬────────────────┬────────────────┐\n")
cat("│                 │  → Violent     │  → Peaceful    │\n")
cat("├─────────────────┼────────────────┼────────────────┤\n")
cat(sprintf("│ Violent parent  │    %.5f     │    %.5f     │\n", alpha_vv, alpha_vp))
cat(sprintf("│ Peaceful parent │    %.5f     │    %.5f     │\n", alpha_pv, alpha_pp))
cat("└─────────────────┴────────────────┴────────────────┘\n")

cat("\nMOBILIZATION COMPARISON:\n")
cat(sprintf("  Violent events:  %.4f total offspring\n", mobil_v))
cat(sprintf("  Peaceful events: %.4f total offspring\n", mobil_p))
cat(sprintf("  Ratio: %.3f (95%% CI: [%.3f, %.3f], p = %.4f)\n",
            mobil_ratio, ci_lower, ci_upper, p_value))

if (ci_lower > 1) {
  cat("\n  *** CONCLUSION: Violent events have SIGNIFICANTLY HIGHER ***\n")
  cat("  *** mobilization potential than peaceful events ***\n")
} else if (ci_upper < 1) {
  cat("\n  *** CONCLUSION: Violent events have SIGNIFICANTLY LOWER ***\n")
  cat("  *** mobilization potential than peaceful events ***\n")
} else {
  cat("\n  CONCLUSION: No significant difference in mobilization potential\n")
}

# Save results
results <- list(
  params = params,
  se = se_log,
  vcov = vcov_mat,
  mu_v = mu_v, mu_p = mu_p,
  alpha_vv = alpha_vv, alpha_vp = alpha_vp,
  alpha_pv = alpha_pv, alpha_pp = alpha_pp,
  beta = beta,
  mobil_v = mobil_v, mobil_p = mobil_p,
  mobil_ratio = mobil_ratio,
  ci_lower = ci_lower, ci_upper = ci_upper,
  p_value = p_value,
  loglik = -best_result$value,
  n_v = n_v, n_p = n_p, T_max = T_max
)

saveRDS(results, "bivariate_hawkes_inference.rds")
cat("\nResults saved to bivariate_hawkes_inference.rds\n")
