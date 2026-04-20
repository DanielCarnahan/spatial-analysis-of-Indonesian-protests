################################################################################
#           BIVARIATE HAWKES MODEL WITH HYPOTHESIS TESTING
#           Based on fast implementation that successfully optimizes
################################################################################

library(tidyverse)
library(numDeriv)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   BIVARIATE HAWKES WITH INFERENCE                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Load data
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

# Use 2019-2021 subset (known to work)
cat("Using 2019-2021 subset...\n")
subset_data <- protests %>%
  filter(event_date >= as.Date("2019-01-01") & event_date <= as.Date("2021-12-31"))

subset_data$time <- as.numeric(subset_data$event_date - min(subset_data$event_date))
T_max <- max(subset_data$time) + 1

# Separate streams
violent_times <- subset_data$time[subset_data$is_violent]
peaceful_times <- subset_data$time[!subset_data$is_violent]

n_v <- length(violent_times)
n_p <- length(peaceful_times)

cat(sprintf("  Violent events: %d\n", n_v))
cat(sprintf("  Peaceful events: %d\n", n_p))
cat(sprintf("  Total: %d over %.0f days\n\n", n_v + n_p, T_max))

# =============================================================================
# EFFICIENT BIVARIATE HAWKES LOG-LIKELIHOOD
# =============================================================================

bivariate_loglik <- function(params, v_times, p_times, T_max) {

  # Parameters (all on log scale for unconstrained optimization)
  mu_v <- exp(params[1])
  mu_p <- exp(params[2])
  alpha_vv <- exp(params[3])
  alpha_vp <- exp(params[4])
  alpha_pv <- exp(params[5])
  alpha_pp <- exp(params[6])
  beta <- exp(params[7])

  # Combine all events with type indicator
  all_times <- c(v_times, p_times)
  all_types <- c(rep(1, length(v_times)), rep(0, length(p_times)))  # 1 = violent
  ord <- order(all_times)
  all_times <- all_times[ord]
  all_types <- all_types[ord]
  n <- length(all_times)

  if (n == 0) return(-1e10)

  # Recursive sums for triggering intensity
  R_v <- numeric(n)
  R_p <- numeric(n)

  for (i in 2:n) {
    dt <- all_times[i] - all_times[i-1]
    decay <- exp(-beta * dt)

    if (all_types[i-1] == 1) {
      R_v[i] <- decay * (1 + R_v[i-1])
      R_p[i] <- decay * R_p[i-1]
    } else {
      R_v[i] <- decay * R_v[i-1]
      R_p[i] <- decay * (1 + R_p[i-1])
    }
  }

  # Log-likelihood
  loglik <- 0

  # Term 1: sum of log(lambda) at event times
  for (i in 1:n) {
    if (all_types[i] == 1) {
      lambda_i <- mu_v + alpha_vv * R_v[i] + alpha_pv * R_p[i]
    } else {
      lambda_i <- mu_p + alpha_vp * R_v[i] + alpha_pp * R_p[i]
    }

    if (lambda_i <= 1e-10) return(-1e10)
    loglik <- loglik + log(lambda_i)
  }

  # Term 2: integral of intensity
  integral_v <- mu_v * T_max
  for (j in which(all_types == 1)) {
    integral_v <- integral_v + (alpha_vv / beta) * (1 - exp(-beta * (T_max - all_times[j])))
  }
  for (j in which(all_types == 0)) {
    integral_v <- integral_v + (alpha_pv / beta) * (1 - exp(-beta * (T_max - all_times[j])))
  }

  integral_p <- mu_p * T_max
  for (j in which(all_types == 1)) {
    integral_p <- integral_p + (alpha_vp / beta) * (1 - exp(-beta * (T_max - all_times[j])))
  }
  for (j in which(all_types == 0)) {
    integral_p <- integral_p + (alpha_pp / beta) * (1 - exp(-beta * (T_max - all_times[j])))
  }

  loglik <- loglik - integral_v - integral_p

  return(loglik)
}

# =============================================================================
# OPTIMIZATION
# =============================================================================

cat("Optimizing bivariate Hawkes model...\n")

# Starting values
start_params <- c(
  log(n_v / T_max * 0.5),    # log_mu_v
  log(n_p / T_max * 0.5),    # log_mu_p
  log(0.05),                  # log_alpha_vv
  log(0.05),                  # log_alpha_vp
  log(0.05),                  # log_alpha_pv
  log(0.05),                  # log_alpha_pp
  log(0.1)                    # log_beta
)

neg_loglik <- function(params) {
  ll <- bivariate_loglik(params, violent_times, peaceful_times, T_max)
  if (is.na(ll) || is.infinite(ll)) return(1e10)
  return(-ll)
}

# Multiple starting points
best_result <- NULL
best_ll <- Inf

cat("  Running optimization with multiple starts...\n")

for (trial in 1:5) {
  init <- start_params + rnorm(7, 0, 0.5)

  result <- tryCatch({
    optim(par = init, fn = neg_loglik, method = "BFGS",
          control = list(maxit = 500))
  }, error = function(e) NULL)

  if (!is.null(result) && result$value < best_ll) {
    best_ll <- result$value
    best_result <- result
    cat(sprintf("    Trial %d: LL = %.1f\n", trial, -result$value))
  }
}

if (is.null(best_result)) {
  cat("Optimization failed!\n")
  quit(status = 1)
}

cat(sprintf("\nBest log-likelihood: %.1f\n", -best_result$value))

# =============================================================================
# COMPUTE HESSIAN FOR STANDARD ERRORS
# =============================================================================

cat("\nComputing Hessian for standard errors...\n")

hessian_mat <- hessian(neg_loglik, best_result$par)

# Variance-covariance matrix
vcov_mat <- tryCatch({
  solve(hessian_mat)
}, error = function(e) {
  cat("Warning: Hessian not invertible, using pseudo-inverse\n")
  MASS::ginv(hessian_mat)
})

# Standard errors on log scale
se_log <- sqrt(pmax(diag(vcov_mat), 0))  # Ensure non-negative before sqrt

# =============================================================================
# EXTRACT AND INTERPRET RESULTS
# =============================================================================

params <- best_result$par
param_names <- c("log_mu_v", "log_mu_p", "log_alpha_vv", "log_alpha_vp",
                 "log_alpha_pv", "log_alpha_pp", "log_beta")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   PARAMETER ESTIMATES (Log Scale)                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

for (i in 1:7) {
  if (se_log[i] > 0 && !is.na(se_log[i])) {
    z <- params[i] / se_log[i]
    p <- 2 * pnorm(-abs(z))
    cat(sprintf("  %s: %.4f (SE: %.4f, z: %.2f, p: %.4f)\n",
                param_names[i], params[i], se_log[i], z, p))
  } else {
    cat(sprintf("  %s: %.4f (SE: NA)\n", param_names[i], params[i]))
  }
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

cat("BASELINE RATES (events/day):\n")
cat(sprintf("  μ_V (violent baseline): %.4f\n", mu_v))
cat(sprintf("  μ_P (peaceful baseline): %.4f\n", mu_p))

cat("\nTRIGGERING MATRIX (α coefficients):\n")
cat("                     → Violent    → Peaceful\n")
cat(sprintf("  Violent parent:      %.5f      %.5f\n", alpha_vv, alpha_vp))
cat(sprintf("  Peaceful parent:     %.5f      %.5f\n", alpha_pv, alpha_pp))

cat(sprintf("\nDECAY RATE: β = %.4f (half-life = %.1f days)\n", beta, log(2)/beta))

# =============================================================================
# MOBILIZATION COMPARISON WITH INFERENCE
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   MOBILIZATION POTENTIAL WITH INFERENCE                      ║\n")
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
cat(sprintf("  From VIOLENT event:   %.4f      %.4f      %.4f\n",
            offspring_vv, offspring_vp, mobil_violent))
cat(sprintf("  From PEACEFUL event:  %.4f      %.4f      %.4f\n",
            offspring_pv, offspring_pp, mobil_peaceful))

cat(sprintf("\n  MOBILIZATION RATIO (Violent / Peaceful): %.3f\n\n", mobil_ratio))

# =============================================================================
# DELTA METHOD FOR CONFIDENCE INTERVAL ON MOBILIZATION RATIO
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   HYPOTHESIS TESTS                                           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Function to compute log(mobilization ratio) from params
# mobil_v = (exp(log_alpha_vv) + exp(log_alpha_vp)) / exp(log_beta)
# mobil_p = (exp(log_alpha_pv) + exp(log_alpha_pp)) / exp(log_beta)
# ratio = mobil_v / mobil_p
# log(ratio) = log(mobil_v) - log(mobil_p)

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

# Gradient using numerical differentiation
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

cat(sprintf("  Mobilization ratio: %.3f\n", mobil_ratio))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n", ci_lower, ci_upper))

# Z-test for log(ratio) = 0 (i.e., ratio = 1)
if (se_log_ratio > 0) {
  z_stat <- log_ratio / se_log_ratio
  p_value <- 2 * pnorm(-abs(z_stat))

  cat(sprintf("  z-statistic: %.3f\n", z_stat))
  cat(sprintf("  p-value: %.6f\n", p_value))

  if (p_value < 0.05) {
    if (mobil_ratio > 1) {
      cat("\n  *** REJECT H0: Violent events have SIGNIFICANTLY HIGHER ***\n")
      cat("  *** mobilization potential than peaceful events (p < 0.05) ***\n")
    } else {
      cat("\n  *** REJECT H0: Violent events have SIGNIFICANTLY LOWER ***\n")
      cat("  *** mobilization potential than peaceful events (p < 0.05) ***\n")
    }
  } else {
    cat("\n  Cannot reject H0: No significant difference in mobilization (p >= 0.05)\n")
  }
} else {
  cat("  Unable to compute z-test (SE = 0 or NA)\n")
  p_value <- NA
}

# =============================================================================
# TEST 2: Individual coefficient tests
# =============================================================================

cat("\n\nTEST 2: Are individual triggering coefficients significant?\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Test each alpha coefficient (indices 3-6)
coef_names <- c("α_VV (V→V)", "α_VP (V→P)", "α_PV (P→V)", "α_PP (P→P)")
for (i in 3:6) {
  if (se_log[i] > 0 && !is.na(se_log[i])) {
    z <- params[i] / se_log[i]
    p <- 2 * pnorm(-abs(z))
    sig <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
    cat(sprintf("  %s: %.5f (SE: %.4f, z: %.2f, p: %.4f) %s\n",
                coef_names[i-2], exp(params[i]), se_log[i], z, p, sig))
  }
}

# =============================================================================
# TEST 3: Is α_VV different from α_PP?
# =============================================================================

cat("\n\nTEST 3: Is α_VV different from α_PP?\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("  H0: Violence-to-violence and peaceful-to-peaceful triggering are equal\n\n")

# Difference in log scale: log(alpha_vv) - log(alpha_pp)
diff_log <- params[3] - params[6]
var_diff <- vcov_mat[3,3] + vcov_mat[6,6] - 2*vcov_mat[3,6]
se_diff <- sqrt(max(var_diff, 0))

if (se_diff > 0) {
  z_diff <- diff_log / se_diff
  p_diff <- 2 * pnorm(-abs(z_diff))

  cat(sprintf("  α_VV = %.5f, α_PP = %.5f\n", alpha_vv, alpha_pp))
  cat(sprintf("  Ratio α_VV/α_PP: %.3f\n", alpha_vv/alpha_pp))
  cat(sprintf("  z-statistic: %.3f\n", z_diff))
  cat(sprintf("  p-value: %.6f\n", p_diff))

  if (p_diff < 0.05) {
    cat("\n  *** α_VV IS significantly different from α_PP (p < 0.05) ***\n")
  } else {
    cat("\n  Cannot reject H0: α_VV and α_PP are not significantly different\n")
  }
} else {
  cat("  Unable to compute test (SE = 0 or NA)\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   SUMMARY: BIVARIATE HAWKES RESULTS                          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("TRIGGERING MATRIX (expected offspring):\n")
cat("┌─────────────────┬────────────────┬────────────────┬──────────────┐\n")
cat("│                 │  → Violent     │  → Peaceful    │    TOTAL     │\n")
cat("├─────────────────┼────────────────┼────────────────┼──────────────┤\n")
cat(sprintf("│ Violent parent  │    %.4f      │    %.4f      │    %.4f    │\n",
            offspring_vv, offspring_vp, mobil_violent))
cat(sprintf("│ Peaceful parent │    %.4f      │    %.4f      │    %.4f    │\n",
            offspring_pv, offspring_pp, mobil_peaceful))
cat("└─────────────────┴────────────────┴────────────────┴──────────────┘\n")

cat("\nKEY RESULT:\n")
cat(sprintf("  Mobilization ratio: %.3f (95%% CI: [%.3f, %.3f])\n",
            mobil_ratio, ci_lower, ci_upper))

if (!is.na(p_value) && p_value < 0.05) {
  if (mobil_ratio > 1 && ci_lower > 1) {
    cat("\n  *** CONCLUSION: Violent events trigger SIGNIFICANTLY MORE ***\n")
    cat("  *** total follow-on events than peaceful events ***\n")
  } else if (mobil_ratio < 1 && ci_upper < 1) {
    cat("\n  *** CONCLUSION: Violent events trigger SIGNIFICANTLY FEWER ***\n")
    cat("  *** total follow-on events than peaceful events ***\n")
  } else {
    cat("\n  Note: Point estimate differs from 1 but CI includes 1\n")
  }
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
  offspring_vv = offspring_vv, offspring_vp = offspring_vp,
  offspring_pv = offspring_pv, offspring_pp = offspring_pp,
  mobil_violent = mobil_violent,
  mobil_peaceful = mobil_peaceful,
  mobil_ratio = mobil_ratio,
  ci_lower = ci_lower, ci_upper = ci_upper,
  p_value = p_value,
  loglik = -best_result$value,
  n_violent = n_v, n_peaceful = n_p, T_max = T_max
)

saveRDS(results, "bivariate_hawkes_inference.rds")
cat("\nResults saved to bivariate_hawkes_inference.rds\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   BIVARIATE MODEL WITH INFERENCE COMPLETE                    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
