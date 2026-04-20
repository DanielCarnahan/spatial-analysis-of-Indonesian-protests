################################################################################
#           BIVARIATE HAWKES MODEL WITH COVARIATE-DEPENDENT BACKGROUND RATE
#
#           Background: μ_V(d,y) = exp(β₀_V + γ·log(pop) + δ·poverty + Σβ_year)
#                       μ_P(d,y) = exp(β₀_P + γ·log(pop) + δ·poverty + Σβ_year)
#           Triggering: 4 cross-coefficients (VV, VP, PV, PP) with shared decay
################################################################################

library(tidyverse)
library(numDeriv)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   BIVARIATE HAWKES WITH COVARIATE-DEPENDENT BACKGROUND                  ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

N_STARTS <- 5
MAX_ITER <- 2000

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

cat("=== LOADING DATA ===\n")

protests <- readRDS("protests_with_poverty.rds")

protests <- protests %>%
  mutate(
    is_violent = sub_event_type %in% c("Mob violence", "Violent demonstration") |
                 sub_event_type == "Excessive force against protesters"
  ) %>%
  arrange(event_date) %>%
  filter(!is.na(poverty_decimal) & !is.na(log_pop))

# Convert to numeric time
start_date <- min(protests$event_date)
protests$time <- as.numeric(protests$event_date - start_date)
T_max <- max(protests$time) + 1

n_total <- nrow(protests)
n_v <- sum(protests$is_violent)
n_p <- sum(!protests$is_violent)

cat(sprintf("Loaded %d events with complete covariate data\n", n_total))
cat(sprintf("  Violent events: %d (%.1f%%)\n", n_v, 100*n_v/n_total))
cat(sprintf("  Peaceful events: %d (%.1f%%)\n", n_p, 100*n_p/n_total))
cat(sprintf("  Observation period: %.0f days\n", T_max))
cat(sprintf("  Years: %s\n\n", paste(sort(unique(protests$year)), collapse=", ")))

# Prepare data vectors
times <- protests$time
types <- as.integer(protests$is_violent)  # 1 = violent, 0 = peaceful
log_pop <- protests$log_pop
poverty <- protests$poverty_decimal
years <- protests$year

# Pre-compute time differences
time_diffs <- diff(times)

# District-year exposure for integral
district_years <- protests %>%
  select(district, year, log_pop, poverty_decimal) %>%
  distinct()
district_years$exposure <- T_max  # Simplified: same exposure for all

cat(sprintf("District-years: %d unique combinations\n\n", nrow(district_years)))

# =============================================================================
# 2. PARAMETER STRUCTURE
# =============================================================================

# Parameters (18 total):
# 1: β₀_V (violent intercept)
# 2: β₀_P (peaceful intercept)
# 3: γ (log population effect, shared)
# 4: δ (poverty effect, shared)
# 5-13: β_year for 2016-2024 (2015 is reference)
# 14: log_α_VV
# 15: log_α_VP
# 16: log_α_PV
# 17: log_α_PP
# 18: log_β (decay)

param_names <- c(
  "beta_0_V", "beta_0_P", "gamma", "delta",
  paste0("beta_", 2016:2024),
  "log_alpha_VV", "log_alpha_VP", "log_alpha_PV", "log_alpha_PP",
  "log_beta"
)

n_params <- length(param_names)
cat(sprintf("Model has %d parameters\n\n", n_params))

# =============================================================================
# 3. LOG-LIKELIHOOD FUNCTION
# =============================================================================

bivariate_loglik_cov <- function(params) {

  # Extract parameters
  beta_0_V <- params[1]
  beta_0_P <- params[2]
  gamma <- params[3]
  delta <- params[4]
  beta_years <- params[5:13]  # 2016-2024

  alpha_vv <- exp(params[14])
  alpha_vp <- exp(params[15])
  alpha_pv <- exp(params[16])
  alpha_pp <- exp(params[17])
  beta_decay <- exp(params[18])

  # Check bounds
  if (beta_decay < 1e-6 || beta_decay > 100) return(-1e10)
  if (any(!is.finite(c(alpha_vv, alpha_vp, alpha_pv, alpha_pp)))) return(-1e10)

  n <- length(times)

  # Compute background rates for each event
  year_effects <- numeric(n)
  for (i in 1:n) {
    yr <- years[i]
    if (yr >= 2016 && yr <= 2024) {
      year_effects[i] <- beta_years[yr - 2015]
    }
  }

  # Background rate at each event location/time
  mu_v <- exp(beta_0_V + gamma * log_pop + delta * poverty + year_effects)
  mu_p <- exp(beta_0_P + gamma * log_pop + delta * poverty + year_effects)

  # Recursive computation of triggering sums
  R_v <- numeric(n)
  R_p <- numeric(n)

  for (i in 2:n) {
    decay <- exp(-beta_decay * time_diffs[i-1])
    if (decay < 1e-100) decay <- 0

    if (types[i-1] == 1) {
      R_v[i] <- decay * (1 + R_v[i-1])
      R_p[i] <- decay * R_p[i-1]
    } else {
      R_v[i] <- decay * R_v[i-1]
      R_p[i] <- decay * (1 + R_p[i-1])
    }
  }

  # Log-likelihood: sum of log(λ) at event times
  loglik <- 0

  for (i in 1:n) {
    if (types[i] == 1) {
      lambda_i <- mu_v[i] + alpha_vv * R_v[i] + alpha_pv * R_p[i]
    } else {
      lambda_i <- mu_p[i] + alpha_vp * R_v[i] + alpha_pp * R_p[i]
    }

    if (lambda_i <= 1e-10) return(-1e10)
    loglik <- loglik + log(lambda_i)
  }

  # Integral terms (compensator)
  # For background: integrate over all district-years
  integral_bg_v <- 0
  integral_bg_p <- 0

  for (j in 1:nrow(district_years)) {
    yr <- district_years$year[j]
    yr_effect <- if (yr >= 2016 && yr <= 2024) beta_years[yr - 2015] else 0

    mu_v_j <- exp(beta_0_V + gamma * district_years$log_pop[j] +
                   delta * district_years$poverty_decimal[j] + yr_effect)
    mu_p_j <- exp(beta_0_P + gamma * district_years$log_pop[j] +
                   delta * district_years$poverty_decimal[j] + yr_effect)

    # Exposure is days observed (simplified to T_max/n_years per district-year)
    exposure_j <- 365.25  # ~1 year per district-year
    integral_bg_v <- integral_bg_v + mu_v_j * exposure_j
    integral_bg_p <- integral_bg_p + mu_p_j * exposure_j
  }

  # For triggering: integrate kernel contributions
  integral_trig_v <- 0
  integral_trig_p <- 0

  for (j in 1:n) {
    contrib <- (1 - exp(-beta_decay * (T_max - times[j]))) / beta_decay
    if (types[j] == 1) {
      integral_trig_v <- integral_trig_v + alpha_vv * contrib
      integral_trig_p <- integral_trig_p + alpha_vp * contrib
    } else {
      integral_trig_v <- integral_trig_v + alpha_pv * contrib
      integral_trig_p <- integral_trig_p + alpha_pp * contrib
    }
  }

  loglik <- loglik - integral_bg_v - integral_bg_p - integral_trig_v - integral_trig_p

  if (!is.finite(loglik)) return(-1e10)

  return(loglik)
}

neg_loglik <- function(params) -bivariate_loglik_cov(params)

cat("✓ Likelihood function ready\n\n")

# =============================================================================
# 4. MULTI-STARTING POINT OPTIMIZATION
# =============================================================================

cat("=== MULTI-STARTING POINT OPTIMIZATION ===\n\n")

# Empirical estimates
emp_rate_v <- n_v / T_max
emp_rate_p <- n_p / T_max
mean_log_pop <- mean(log_pop)
mean_poverty <- mean(poverty)

cat(sprintf("Empirical rates: violent=%.4f, peaceful=%.4f per day\n", emp_rate_v, emp_rate_p))
cat(sprintf("Mean log(pop)=%.2f, mean poverty=%.3f\n\n", mean_log_pop, mean_poverty))

# Starting points
starting_points <- list(
  list(
    name = "Data-driven",
    params = c(
      log(emp_rate_v) - 0.3 * mean_log_pop,  # β₀_V
      log(emp_rate_p) - 0.3 * mean_log_pop,  # β₀_P
      0.3,   # γ
      0.0,   # δ
      rep(0, 9),  # year effects
      log(0.05), log(0.05), log(0.05), log(0.05),  # α's
      log(0.1)   # β_decay
    )
  ),
  list(
    name = "Strong-pop",
    params = c(
      log(emp_rate_v) - 0.5 * mean_log_pop,
      log(emp_rate_p) - 0.5 * mean_log_pop,
      0.5,   # stronger pop effect
      0.5,   # positive poverty
      rep(0, 9),
      log(0.03), log(0.03), log(0.03), log(0.03),
      log(0.15)
    )
  ),
  list(
    name = "Fast-decay",
    params = c(
      log(emp_rate_v) - 0.3 * mean_log_pop,
      log(emp_rate_p) - 0.3 * mean_log_pop,
      0.3,
      -0.5,  # negative poverty
      rep(0, 9),
      log(0.10), log(0.10), log(0.10), log(0.10),
      log(0.3)   # faster decay
    )
  ),
  list(
    name = "Asymmetric-VV",
    params = c(
      log(emp_rate_v) - 0.4 * mean_log_pop,
      log(emp_rate_p) - 0.4 * mean_log_pop,
      0.4,
      0.2,
      rep(0, 9),
      log(0.15), log(0.05), log(0.02), log(0.05),  # Higher VV
      log(0.1)
    )
  ),
  list(
    name = "Asymmetric-VP",
    params = c(
      log(emp_rate_v) - 0.4 * mean_log_pop,
      log(emp_rate_p) - 0.4 * mean_log_pop,
      0.4,
      0.2,
      rep(0, 9),
      log(0.05), log(0.15), log(0.02), log(0.05),  # Higher VP
      log(0.1)
    )
  )
)

# Bounds
lower_bounds <- c(
  -20, -20,      # intercepts
  0.01, -10,    # γ (positive), δ
  rep(-5, 9),   # year effects
  -15, -15, -15, -15,  # log α's
  -5            # log β
)
upper_bounds <- c(
  10, 10,       # intercepts
  2, 10,        # γ, δ
  rep(5, 9),    # year effects
  5, 5, 5, 5,   # log α's
  5             # log β
)

# Run optimization
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
# 5. RESULTS SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  OPTIMIZATION SUMMARY                                        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

for (s in seq_along(all_results)) {
  r <- all_results[[s]]
  marker <- if (r$starting_point == best_start) " ★ BEST" else ""
  cat(sprintf("  %d. %-20s  LL: %10.2f  Conv: %d%s\n",
              s, r$starting_point, r$loglik, r$convergence, marker))
}
cat(sprintf("\nTotal runtime: %.1f minutes\n", total_runtime))

if (is.null(best_result)) {
  stop("All optimizations failed!")
}

# =============================================================================
# 6. EXTRACT ESTIMATES
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PARAMETER ESTIMATES                                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

params <- best_result$par

# Background parameters
beta_0_V <- params[1]
beta_0_P <- params[2]
gamma <- params[3]
delta <- params[4]
beta_years <- params[5:13]

# Triggering parameters
alpha_vv <- exp(params[14])
alpha_vp <- exp(params[15])
alpha_pv <- exp(params[16])
alpha_pp <- exp(params[17])
beta_decay <- exp(params[18])

cat("BACKGROUND RATE PARAMETERS:\n")
cat(sprintf("  β₀_V (violent intercept): %.4f\n", beta_0_V))
cat(sprintf("  β₀_P (peaceful intercept): %.4f\n", beta_0_P))
cat(sprintf("  γ (log population): %.4f\n", gamma))
cat(sprintf("  δ (poverty): %.4f\n\n", delta))

cat("YEAR EFFECTS (relative to 2015):\n")
for (i in 1:9) {
  cat(sprintf("  β_%d: %.4f\n", 2015 + i, beta_years[i]))
}

cat("\nTRIGGERING PARAMETERS:\n")
cat("                      → Violent     → Peaceful\n")
cat(sprintf("  Violent parent:     %.6f     %.6f\n", alpha_vv, alpha_vp))
cat(sprintf("  Peaceful parent:    %.6f     %.6f\n", alpha_pv, alpha_pp))
cat(sprintf("\nDecay rate: β = %.6f (half-life = %.1f days)\n", beta_decay, log(2)/beta_decay))

# =============================================================================
# 7. MOBILIZATION COMPARISON
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MOBILIZATION POTENTIAL                                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

offspring_vv <- alpha_vv / beta_decay
offspring_vp <- alpha_vp / beta_decay
offspring_pv <- alpha_pv / beta_decay
offspring_pp <- alpha_pp / beta_decay

mobil_violent <- offspring_vv + offspring_vp
mobil_peaceful <- offspring_pv + offspring_pp
mobil_ratio <- mobil_violent / mobil_peaceful

cat("EXPECTED OFFSPRING PER PARENT:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("                        Violent     Peaceful    TOTAL\n")
cat(sprintf("  From VIOLENT:         %.5f     %.5f     %.5f\n",
            offspring_vv, offspring_vp, mobil_violent))
cat(sprintf("  From PEACEFUL:        %.5f     %.5f     %.5f\n",
            offspring_pv, offspring_pp, mobil_peaceful))
cat(sprintf("\nMOBILIZATION RATIO: %.4f\n", mobil_ratio))

# =============================================================================
# 8. COMPUTE STANDARD ERRORS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  COMPUTING STANDARD ERRORS                                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Computing Hessian...\n")

hessian_mat <- tryCatch({
  hessian(neg_loglik, best_result$par)
}, error = function(e) {
  cat(sprintf("  Warning: %s\n", e$message))
  return(NULL)
})

if (!is.null(hessian_mat)) {
  vcov_mat <- tryCatch({
    solve(hessian_mat)
  }, error = function(e) {
    cat("  Using pseudo-inverse\n")
    MASS::ginv(hessian_mat)
  })

  se <- sqrt(pmax(diag(vcov_mat), 0))

  # Delta method for mobilization ratio CI
  log_mobil_ratio_func <- function(p) {
    alpha_vv <- exp(p[14])
    alpha_vp <- exp(p[15])
    alpha_pv <- exp(p[16])
    alpha_pp <- exp(p[17])
    beta <- exp(p[18])

    mobil_v <- (alpha_vv + alpha_vp) / beta
    mobil_p <- (alpha_pv + alpha_pp) / beta

    log(mobil_v / mobil_p)
  }

  gradient <- grad(log_mobil_ratio_func, params)
  var_log_ratio <- as.numeric(t(gradient) %*% vcov_mat %*% gradient)
  se_log_ratio <- sqrt(max(var_log_ratio, 0))

  log_ratio <- log(mobil_ratio)
  ci_lower <- exp(log_ratio - 1.96 * se_log_ratio)
  ci_upper <- exp(log_ratio + 1.96 * se_log_ratio)

  if (se_log_ratio > 0) {
    z_stat <- log_ratio / se_log_ratio
    p_value <- 2 * pnorm(-abs(z_stat))
  } else {
    z_stat <- NA
    p_value <- NA
  }

  cat("\nHYPOTHESIS TEST: Mobilization ratio = 1?\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat(sprintf("  Mobilization ratio: %.4f\n", mobil_ratio))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", ci_lower, ci_upper))
  if (!is.na(z_stat)) {
    cat(sprintf("  z-statistic: %.4f\n", z_stat))
    cat(sprintf("  p-value: %.6f\n", p_value))
  }

  if (!is.na(p_value) && p_value < 0.05 && ci_lower > 1) {
    cat("\n  *** VIOLENT EVENTS HAVE SIGNIFICANTLY HIGHER MOBILIZATION ***\n")
  } else if (!is.na(p_value) && p_value < 0.05 && ci_upper < 1) {
    cat("\n  *** VIOLENT EVENTS HAVE SIGNIFICANTLY LOWER MOBILIZATION ***\n")
  }

} else {
  se <- rep(NA, n_params)
  vcov_mat <- NULL
  ci_lower <- ci_upper <- p_value <- NA
  cat("  Could not compute standard errors\n")
}

# =============================================================================
# 9. FINAL SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FINAL SUMMARY                                               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("Model: Bivariate Hawkes with Covariate-Dependent Background\n"))
cat(sprintf("Data: %d events (%d violent, %d peaceful)\n", n_total, n_v, n_p))
cat(sprintf("Log-likelihood: %.2f\n", best_ll))
cat(sprintf("Parameters: %d\n", n_params))
cat(sprintf("AIC: %.2f\n", -2 * best_ll + 2 * n_params))
cat(sprintf("BIC: %.2f\n\n", -2 * best_ll + n_params * log(n_total)))

cat("TRIGGERING MATRIX (expected offspring):\n")
cat("┌─────────────────┬────────────────┬────────────────┬──────────────┐\n")
cat("│                 │  → Violent     │  → Peaceful    │    TOTAL     │\n")
cat("├─────────────────┼────────────────┼────────────────┼──────────────┤\n")
cat(sprintf("│ Violent parent  │    %.5f     │    %.5f     │    %.5f   │\n",
            offspring_vv, offspring_vp, mobil_violent))
cat(sprintf("│ Peaceful parent │    %.5f     │    %.5f     │    %.5f   │\n",
            offspring_pv, offspring_pp, mobil_peaceful))
cat("└─────────────────┴────────────────┴────────────────┴──────────────┘\n\n")

cat(sprintf("KEY RESULT: Mobilization ratio = %.4f", mobil_ratio))
if (!is.na(ci_lower)) {
  cat(sprintf(" (95%% CI: [%.4f, %.4f])", ci_lower, ci_upper))
}
if (!is.na(p_value)) {
  cat(sprintf(", p = %.4f", p_value))
}
cat("\n")

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================

cat("\n=== SAVING RESULTS ===\n")

results <- list(
  model_name = "Bivariate Hawkes with Covariates",

  # Data info
  n_violent = n_v,
  n_peaceful = n_p,
  n_total = n_total,
  T_max = T_max,

  # Parameters
  params = params,
  param_names = param_names,

  # Background
  beta_0_V = beta_0_V, beta_0_P = beta_0_P,
  gamma = gamma, delta = delta,
  beta_years = beta_years,

  # Triggering
  alpha_vv = alpha_vv, alpha_vp = alpha_vp,
  alpha_pv = alpha_pv, alpha_pp = alpha_pp,
  beta_decay = beta_decay,

  # Offspring
  offspring_vv = offspring_vv, offspring_vp = offspring_vp,
  offspring_pv = offspring_pv, offspring_pp = offspring_pp,
  mobil_violent = mobil_violent,
  mobil_peaceful = mobil_peaceful,
  mobil_ratio = mobil_ratio,

  # Inference
  se = se,
  vcov = vcov_mat,
  ci_lower = ci_lower,
  ci_upper = ci_upper,
  p_value = p_value,

  # Fit
  loglik = best_ll,
  n_params = n_params,
  aic = -2 * best_ll + 2 * n_params,
  bic = -2 * best_ll + n_params * log(n_total),

  # Optimization
  all_results = all_results,
  best_starting_point = best_start,
  total_runtime = total_runtime
)

saveRDS(results, "bivariate_hawkes_covariates.rds")
cat("✓ Saved: bivariate_hawkes_covariates.rds\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  BIVARIATE MODEL WITH COVARIATES COMPLETE                    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
