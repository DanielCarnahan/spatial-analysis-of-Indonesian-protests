################################################################################
#           TRIVARIATE HAWKES MODEL WITH COVARIATE-DEPENDENT BACKGROUND RATE
#
#           Background: μ_k(d,y) = exp(β₀_k + γ·log(pop) + δ·poverty + Σβ_year)
#           where k ∈ {F, V, P}
#           Triggering: 9 cross-coefficients with shared decay
################################################################################

library(tidyverse)
library(numDeriv)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   TRIVARIATE HAWKES WITH COVARIATE-DEPENDENT BACKGROUND                 ║\n")
cat("║   Event types: Fatal (F), Violent non-fatal (V), Peaceful (P)           ║\n")
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

# Create trivariate classification
protests <- protests %>%
  mutate(
    is_violent = sub_event_type %in% c("Mob violence", "Violent demonstration") |
                 sub_event_type == "Excessive force against protesters",
    event_type_3 = case_when(
      fatalities > 0 ~ "F",
      is_violent & fatalities == 0 ~ "V",
      TRUE ~ "P"
    ),
    type_numeric = case_when(
      event_type_3 == "F" ~ 1L,
      event_type_3 == "V" ~ 2L,
      event_type_3 == "P" ~ 3L
    )
  ) %>%
  arrange(event_date) %>%
  filter(!is.na(poverty_decimal) & !is.na(log_pop))

# Convert to numeric time
start_date <- min(protests$event_date)
protests$time <- as.numeric(protests$event_date - start_date)
T_max <- max(protests$time) + 1

n_total <- nrow(protests)
n_f <- sum(protests$event_type_3 == "F")
n_v <- sum(protests$event_type_3 == "V")
n_p <- sum(protests$event_type_3 == "P")

cat(sprintf("Loaded %d events with complete covariate data\n", n_total))
cat(sprintf("  Fatal events (F):    %d (%.1f%%)\n", n_f, 100*n_f/n_total))
cat(sprintf("  Violent events (V):  %d (%.1f%%)\n", n_v, 100*n_v/n_total))
cat(sprintf("  Peaceful events (P): %d (%.1f%%)\n", n_p, 100*n_p/n_total))
cat(sprintf("  Observation period: %.0f days\n", T_max))
cat(sprintf("  Years: %s\n\n", paste(sort(unique(protests$year)), collapse=", ")))

# Prepare data vectors
times <- protests$time
types <- protests$type_numeric  # 1=F, 2=V, 3=P
log_pop <- protests$log_pop
poverty <- protests$poverty_decimal
years <- protests$year

# Pre-compute time differences
time_diffs <- diff(times)

# District-year exposure for integral
district_years <- protests %>%
  select(district, year, log_pop, poverty_decimal) %>%
  distinct()

cat(sprintf("District-years: %d unique combinations\n\n", nrow(district_years)))

# =============================================================================
# 2. PARAMETER STRUCTURE
# =============================================================================

# Parameters (24 total):
# 1: β₀_F (fatal intercept)
# 2: β₀_V (violent intercept)
# 3: β₀_P (peaceful intercept)
# 4: γ (log population effect, shared)
# 5: δ (poverty effect, shared)
# 6-14: β_year for 2016-2024 (2015 is reference)
# 15-23: log_α (9 triggering: FF, FV, FP, VF, VV, VP, PF, PV, PP)
# 24: log_β (decay)

param_names <- c(
  "beta_0_F", "beta_0_V", "beta_0_P", "gamma", "delta",
  paste0("beta_", 2016:2024),
  "log_alpha_FF", "log_alpha_FV", "log_alpha_FP",
  "log_alpha_VF", "log_alpha_VV", "log_alpha_VP",
  "log_alpha_PF", "log_alpha_PV", "log_alpha_PP",
  "log_beta"
)

n_params <- length(param_names)
cat(sprintf("Model has %d parameters\n\n", n_params))

# =============================================================================
# 3. LOG-LIKELIHOOD FUNCTION
# =============================================================================

trivariate_loglik_cov <- function(params) {

  # Extract parameters
  beta_0_F <- params[1]
  beta_0_V <- params[2]
  beta_0_P <- params[3]
  gamma <- params[4]
  delta <- params[5]
  beta_years <- params[6:14]  # 2016-2024

  # Triggering matrix
  alpha_ff <- exp(params[15])
  alpha_fv <- exp(params[16])
  alpha_fp <- exp(params[17])
  alpha_vf <- exp(params[18])
  alpha_vv <- exp(params[19])
  alpha_vp <- exp(params[20])
  alpha_pf <- exp(params[21])
  alpha_pv <- exp(params[22])
  alpha_pp <- exp(params[23])

  beta_decay <- exp(params[24])

  # Check bounds
  if (beta_decay < 1e-6 || beta_decay > 100) return(-1e10)
  all_alpha <- c(alpha_ff, alpha_fv, alpha_fp, alpha_vf, alpha_vv, alpha_vp, alpha_pf, alpha_pv, alpha_pp)
  if (any(!is.finite(all_alpha))) return(-1e10)

  n <- length(times)

  # Compute year effects for each event
  year_effects <- numeric(n)
  for (i in 1:n) {
    yr <- years[i]
    if (yr >= 2016 && yr <= 2024) {
      year_effects[i] <- beta_years[yr - 2015]
    }
  }

  # Background rate at each event location/time
  mu_f <- exp(beta_0_F + gamma * log_pop + delta * poverty + year_effects)
  mu_v <- exp(beta_0_V + gamma * log_pop + delta * poverty + year_effects)
  mu_p <- exp(beta_0_P + gamma * log_pop + delta * poverty + year_effects)

  # Recursive computation of triggering sums
  R_f <- numeric(n)
  R_v <- numeric(n)
  R_p <- numeric(n)

  for (i in 2:n) {
    decay <- exp(-beta_decay * time_diffs[i-1])
    if (decay < 1e-100) decay <- 0

    if (types[i-1] == 1) {
      # Previous was Fatal
      R_f[i] <- decay * (1 + R_f[i-1])
      R_v[i] <- decay * R_v[i-1]
      R_p[i] <- decay * R_p[i-1]
    } else if (types[i-1] == 2) {
      # Previous was Violent
      R_f[i] <- decay * R_f[i-1]
      R_v[i] <- decay * (1 + R_v[i-1])
      R_p[i] <- decay * R_p[i-1]
    } else {
      # Previous was Peaceful
      R_f[i] <- decay * R_f[i-1]
      R_v[i] <- decay * R_v[i-1]
      R_p[i] <- decay * (1 + R_p[i-1])
    }
  }

  # Log-likelihood: sum of log(λ) at event times
  loglik <- 0

  for (i in 1:n) {
    if (types[i] == 1) {
      # Fatal event
      lambda_i <- mu_f[i] + alpha_ff * R_f[i] + alpha_vf * R_v[i] + alpha_pf * R_p[i]
    } else if (types[i] == 2) {
      # Violent event
      lambda_i <- mu_v[i] + alpha_fv * R_f[i] + alpha_vv * R_v[i] + alpha_pv * R_p[i]
    } else {
      # Peaceful event
      lambda_i <- mu_p[i] + alpha_fp * R_f[i] + alpha_vp * R_v[i] + alpha_pp * R_p[i]
    }

    if (lambda_i <= 1e-10) return(-1e10)
    loglik <- loglik + log(lambda_i)
  }

  # Integral terms (compensator)
  # For background: integrate over all district-years
  integral_bg_f <- 0
  integral_bg_v <- 0
  integral_bg_p <- 0

  for (j in 1:nrow(district_years)) {
    yr <- district_years$year[j]
    yr_effect <- if (yr >= 2016 && yr <= 2024) beta_years[yr - 2015] else 0

    mu_f_j <- exp(beta_0_F + gamma * district_years$log_pop[j] +
                   delta * district_years$poverty_decimal[j] + yr_effect)
    mu_v_j <- exp(beta_0_V + gamma * district_years$log_pop[j] +
                   delta * district_years$poverty_decimal[j] + yr_effect)
    mu_p_j <- exp(beta_0_P + gamma * district_years$log_pop[j] +
                   delta * district_years$poverty_decimal[j] + yr_effect)

    exposure_j <- 365.25
    integral_bg_f <- integral_bg_f + mu_f_j * exposure_j
    integral_bg_v <- integral_bg_v + mu_v_j * exposure_j
    integral_bg_p <- integral_bg_p + mu_p_j * exposure_j
  }

  # For triggering: integrate kernel contributions
  integral_trig_f <- 0
  integral_trig_v <- 0
  integral_trig_p <- 0

  for (j in 1:n) {
    contrib <- (1 - exp(-beta_decay * (T_max - times[j]))) / beta_decay

    if (types[j] == 1) {
      # Fatal parent
      integral_trig_f <- integral_trig_f + alpha_ff * contrib
      integral_trig_v <- integral_trig_v + alpha_fv * contrib
      integral_trig_p <- integral_trig_p + alpha_fp * contrib
    } else if (types[j] == 2) {
      # Violent parent
      integral_trig_f <- integral_trig_f + alpha_vf * contrib
      integral_trig_v <- integral_trig_v + alpha_vv * contrib
      integral_trig_p <- integral_trig_p + alpha_vp * contrib
    } else {
      # Peaceful parent
      integral_trig_f <- integral_trig_f + alpha_pf * contrib
      integral_trig_v <- integral_trig_v + alpha_pv * contrib
      integral_trig_p <- integral_trig_p + alpha_pp * contrib
    }
  }

  loglik <- loglik - integral_bg_f - integral_bg_v - integral_bg_p
  loglik <- loglik - integral_trig_f - integral_trig_v - integral_trig_p

  if (!is.finite(loglik)) return(-1e10)

  return(loglik)
}

neg_loglik <- function(params) -trivariate_loglik_cov(params)

cat("✓ Likelihood function ready\n\n")

# =============================================================================
# 4. MULTI-STARTING POINT OPTIMIZATION
# =============================================================================

cat("=== MULTI-STARTING POINT OPTIMIZATION ===\n\n")

# Empirical estimates
emp_rate_f <- n_f / T_max
emp_rate_v <- n_v / T_max
emp_rate_p <- n_p / T_max
mean_log_pop <- mean(log_pop)
mean_poverty <- mean(poverty)

cat(sprintf("Empirical rates: F=%.5f, V=%.4f, P=%.4f per day\n",
            emp_rate_f, emp_rate_v, emp_rate_p))
cat(sprintf("Mean log(pop)=%.2f, mean poverty=%.3f\n\n", mean_log_pop, mean_poverty))

# Starting points
starting_points <- list(
  list(
    name = "Data-driven",
    params = c(
      log(emp_rate_f) - 0.3 * mean_log_pop,  # β₀_F
      log(emp_rate_v) - 0.3 * mean_log_pop,  # β₀_V
      log(emp_rate_p) - 0.3 * mean_log_pop,  # β₀_P
      0.3,   # γ
      0.0,   # δ
      rep(0, 9),  # year effects
      rep(log(0.02), 9),  # 9 α's
      log(0.1)   # β_decay
    )
  ),
  list(
    name = "Strong-pop",
    params = c(
      log(emp_rate_f) - 0.5 * mean_log_pop,
      log(emp_rate_v) - 0.5 * mean_log_pop,
      log(emp_rate_p) - 0.5 * mean_log_pop,
      0.5,   # stronger pop effect
      0.5,   # positive poverty
      rep(0, 9),
      rep(log(0.01), 9),
      log(0.15)
    )
  ),
  list(
    name = "Fast-decay",
    params = c(
      log(emp_rate_f) - 0.3 * mean_log_pop,
      log(emp_rate_v) - 0.3 * mean_log_pop,
      log(emp_rate_p) - 0.3 * mean_log_pop,
      0.3,
      -0.5,  # negative poverty
      rep(0, 9),
      rep(log(0.05), 9),
      log(0.3)
    )
  ),
  list(
    name = "Asymmetric-Fatal",
    params = c(
      log(emp_rate_f) - 0.4 * mean_log_pop,
      log(emp_rate_v) - 0.4 * mean_log_pop,
      log(emp_rate_p) - 0.4 * mean_log_pop,
      0.4,
      0.2,
      rep(0, 9),
      log(0.10), log(0.10), log(0.10),  # Fatal→ high
      log(0.02), log(0.02), log(0.02),  # Violent→ low
      log(0.01), log(0.01), log(0.01),  # Peaceful→ lowest
      log(0.1)
    )
  ),
  list(
    name = "Asymmetric-Peaceful",
    params = c(
      log(emp_rate_f) - 0.4 * mean_log_pop,
      log(emp_rate_v) - 0.4 * mean_log_pop,
      log(emp_rate_p) - 0.4 * mean_log_pop,
      0.4,
      0.2,
      rep(0, 9),
      log(0.005), log(0.005), log(0.005),  # Fatal→ low
      log(0.01), log(0.01), log(0.01),     # Violent→ low
      log(0.05), log(0.05), log(0.05),     # Peaceful→ highest
      log(0.1)
    )
  )
)

# Bounds
lower_bounds <- c(
  -20, -20, -20,  # 3 intercepts
  0.01, -10,      # γ (positive), δ
  rep(-5, 9),     # year effects
  rep(-15, 9),    # 9 log α's
  -5              # log β
)
upper_bounds <- c(
  10, 10, 10,     # 3 intercepts
  2, 10,          # γ, δ
  rep(5, 9),      # year effects
  rep(5, 9),      # 9 log α's
  5               # log β
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

    # Only consider converged results for best
    if (result$convergence == 0 && ll > best_ll) {
      best_ll <- ll
      best_result <- result
      best_start <- paste0(sp$name, " (converged)")
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
  cat(sprintf("  %d. %-22s  LL: %10.2f  Conv: %d%s\n",
              s, r$starting_point, r$loglik, r$convergence, marker))
}
cat(sprintf("\nTotal runtime: %.1f minutes\n", total_runtime))

if (is.null(best_result)) {
  cat("\n⚠️  WARNING: No starting points converged!\n")
  cat("Using best non-converged result (use with caution).\n\n")

  # Find best non-converged
  for (r in all_results) {
    if (!is.null(r$params) && r$loglik > best_ll) {
      best_ll <- r$loglik
      best_result <- list(par = r$params, value = -r$loglik, convergence = r$convergence)
      best_start <- paste0(r$starting_point, " (non-converged)")
    }
  }

  if (is.null(best_result)) {
    stop("All optimizations failed!")
  }
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
beta_0_F <- params[1]
beta_0_V <- params[2]
beta_0_P <- params[3]
gamma <- params[4]
delta <- params[5]
beta_years <- params[6:14]

# Triggering parameters
alpha_ff <- exp(params[15])
alpha_fv <- exp(params[16])
alpha_fp <- exp(params[17])
alpha_vf <- exp(params[18])
alpha_vv <- exp(params[19])
alpha_vp <- exp(params[20])
alpha_pf <- exp(params[21])
alpha_pv <- exp(params[22])
alpha_pp <- exp(params[23])

beta_decay <- exp(params[24])

cat("BACKGROUND RATE PARAMETERS:\n")
cat(sprintf("  β₀_F (fatal intercept): %.4f\n", beta_0_F))
cat(sprintf("  β₀_V (violent intercept): %.4f\n", beta_0_V))
cat(sprintf("  β₀_P (peaceful intercept): %.4f\n", beta_0_P))
cat(sprintf("  γ (log population): %.4f\n", gamma))
cat(sprintf("  δ (poverty): %.4f\n\n", delta))

cat("YEAR EFFECTS (relative to 2015):\n")
for (i in 1:9) {
  cat(sprintf("  β_%d: %.4f\n", 2015 + i, beta_years[i]))
}

cat("\nTRIGGERING PARAMETERS:\n")
cat("                       → Fatal       → Violent     → Peaceful\n")
cat(sprintf("  Fatal parent:       %.6f     %.6f     %.6f\n", alpha_ff, alpha_fv, alpha_fp))
cat(sprintf("  Violent parent:     %.6f     %.6f     %.6f\n", alpha_vf, alpha_vv, alpha_vp))
cat(sprintf("  Peaceful parent:    %.6f     %.6f     %.6f\n\n", alpha_pf, alpha_pv, alpha_pp))

cat(sprintf("Decay rate: β = %.6f (half-life = %.1f days)\n", beta_decay, log(2)/beta_decay))

# =============================================================================
# 7. MOBILIZATION COMPARISON
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MOBILIZATION POTENTIAL                                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

offspring_ff <- alpha_ff / beta_decay
offspring_fv <- alpha_fv / beta_decay
offspring_fp <- alpha_fp / beta_decay
offspring_vf <- alpha_vf / beta_decay
offspring_vv <- alpha_vv / beta_decay
offspring_vp <- alpha_vp / beta_decay
offspring_pf <- alpha_pf / beta_decay
offspring_pv <- alpha_pv / beta_decay
offspring_pp <- alpha_pp / beta_decay

mobil_fatal <- offspring_ff + offspring_fv + offspring_fp
mobil_violent <- offspring_vf + offspring_vv + offspring_vp
mobil_peaceful <- offspring_pf + offspring_pv + offspring_pp

ratio_fv <- mobil_fatal / mobil_violent
ratio_vp <- mobil_violent / mobil_peaceful
ratio_fp <- mobil_fatal / mobil_peaceful

cat("EXPECTED OFFSPRING PER PARENT:\n")
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
cat(sprintf("  Fatal / Peaceful:   %.4f\n", ratio_fp))

# =============================================================================
# 8. COMPUTE STANDARD ERRORS & HYPOTHESIS TESTS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  COMPUTING STANDARD ERRORS & HYPOTHESIS TESTS                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Computing Hessian...\n")

hessian_mat <- tryCatch({
  hessian(neg_loglik, best_result$par)
}, error = function(e) {
  cat(sprintf("  Warning: %s\n", e$message))
  return(NULL)
})

# Initialize CI/p-value storage
ci_fv_lower <- NA; ci_fv_upper <- NA; p_value_fv <- NA
ci_vp_lower <- NA; ci_vp_upper <- NA; p_value_vp <- NA

if (!is.null(hessian_mat)) {
  vcov_mat <- tryCatch({
    solve(hessian_mat)
  }, error = function(e) {
    cat("  Using pseudo-inverse\n")
    MASS::ginv(hessian_mat)
  })

  se <- sqrt(pmax(diag(vcov_mat), 0))

  # Delta method for F/V ratio
  log_ratio_fv_func <- function(p) {
    alpha_ff <- exp(p[15]); alpha_fv <- exp(p[16]); alpha_fp <- exp(p[17])
    alpha_vf <- exp(p[18]); alpha_vv <- exp(p[19]); alpha_vp <- exp(p[20])
    beta <- exp(p[24])
    mobil_f <- (alpha_ff + alpha_fv + alpha_fp) / beta
    mobil_v <- (alpha_vf + alpha_vv + alpha_vp) / beta
    log(mobil_f / mobil_v)
  }

  log_ratio_vp_func <- function(p) {
    alpha_vf <- exp(p[18]); alpha_vv <- exp(p[19]); alpha_vp <- exp(p[20])
    alpha_pf <- exp(p[21]); alpha_pv <- exp(p[22]); alpha_pp <- exp(p[23])
    beta <- exp(p[24])
    mobil_v <- (alpha_vf + alpha_vv + alpha_vp) / beta
    mobil_p <- (alpha_pf + alpha_pv + alpha_pp) / beta
    log(mobil_v / mobil_p)
  }

  # H2: Fatal vs Violent
  gradient_fv <- grad(log_ratio_fv_func, params)
  var_log_fv <- as.numeric(t(gradient_fv) %*% vcov_mat %*% gradient_fv)
  se_log_fv <- sqrt(max(var_log_fv, 0))
  log_fv <- log(ratio_fv)

  cat("\nTEST H2: Fatal vs Violent mobilization\n")
  cat("─────────────────────────────────────────────────────────────────\n")

  if (se_log_fv > 0 && is.finite(se_log_fv)) {
    ci_fv_lower <- exp(log_fv - 1.96 * se_log_fv)
    ci_fv_upper <- exp(log_fv + 1.96 * se_log_fv)
    z_stat_fv <- log_fv / se_log_fv
    p_value_fv <- 2 * pnorm(-abs(z_stat_fv))

    cat(sprintf("  Ratio (F/V): %.4f [95%% CI: %.4f, %.4f]\n", ratio_fv, ci_fv_lower, ci_fv_upper))
    cat(sprintf("  z = %.4f, p = %.6f\n", z_stat_fv, p_value_fv))

    if (p_value_fv < 0.05 && ci_fv_lower > 1) {
      cat("  *** Fatal has HIGHER mobilization than violent ***\n")
    } else if (p_value_fv < 0.05 && ci_fv_upper < 1) {
      cat("  *** Fatal has LOWER mobilization than violent ***\n")
    } else {
      cat("  No significant difference\n")
    }
  } else {
    cat("  Unable to compute CI\n")
  }

  # H3: Violent vs Peaceful
  gradient_vp <- grad(log_ratio_vp_func, params)
  var_log_vp <- as.numeric(t(gradient_vp) %*% vcov_mat %*% gradient_vp)
  se_log_vp <- sqrt(max(var_log_vp, 0))
  log_vp <- log(ratio_vp)

  cat("\nTEST H3: Violent vs Peaceful mobilization\n")
  cat("─────────────────────────────────────────────────────────────────\n")

  if (se_log_vp > 0 && is.finite(se_log_vp)) {
    ci_vp_lower <- exp(log_vp - 1.96 * se_log_vp)
    ci_vp_upper <- exp(log_vp + 1.96 * se_log_vp)
    z_stat_vp <- log_vp / se_log_vp
    p_value_vp <- 2 * pnorm(-abs(z_stat_vp))

    cat(sprintf("  Ratio (V/P): %.4f [95%% CI: %.4f, %.4f]\n", ratio_vp, ci_vp_lower, ci_vp_upper))
    cat(sprintf("  z = %.4f, p = %.6f\n", z_stat_vp, p_value_vp))

    if (p_value_vp < 0.05 && ci_vp_lower > 1) {
      cat("  *** Violent has HIGHER mobilization than peaceful ***\n")
    } else if (p_value_vp < 0.05 && ci_vp_upper < 1) {
      cat("  *** Violent has LOWER mobilization than peaceful ***\n")
    } else {
      cat("  No significant difference\n")
    }
  } else {
    cat("  Unable to compute CI\n")
  }

} else {
  se <- rep(NA, n_params)
  vcov_mat <- NULL
  cat("  Could not compute standard errors\n")
}

# =============================================================================
# 9. FINAL SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  FINAL SUMMARY                                               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("Model: Trivariate Hawkes with Covariate-Dependent Background\n"))
cat(sprintf("Data: %d events (%d fatal, %d violent, %d peaceful)\n", n_total, n_f, n_v, n_p))
cat(sprintf("Log-likelihood: %.2f\n", best_ll))
cat(sprintf("Parameters: %d\n", n_params))
cat(sprintf("AIC: %.2f\n", -2 * best_ll + 2 * n_params))
cat(sprintf("BIC: %.2f\n\n", -2 * best_ll + n_params * log(n_total)))

cat("COVARIATE EFFECTS:\n")
cat(sprintf("  γ (population): %.4f\n", gamma))
cat(sprintf("  δ (poverty): %.4f\n\n", delta))

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
cat("\n")

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================

cat("\n=== SAVING RESULTS ===\n")

results <- list(
  model_name = "Trivariate Hawkes with Covariates",

  # Data info
  n_fatal = n_f,
  n_violent = n_v,
  n_peaceful = n_p,
  n_total = n_total,
  T_max = T_max,

  # Parameters
  params = params,
  param_names = param_names,

  # Background
  beta_0_F = beta_0_F, beta_0_V = beta_0_V, beta_0_P = beta_0_P,
  gamma = gamma, delta = delta,
  beta_years = beta_years,

  # Triggering
  alpha_ff = alpha_ff, alpha_fv = alpha_fv, alpha_fp = alpha_fp,
  alpha_vf = alpha_vf, alpha_vv = alpha_vv, alpha_vp = alpha_vp,
  alpha_pf = alpha_pf, alpha_pv = alpha_pv, alpha_pp = alpha_pp,
  beta_decay = beta_decay,

  # Offspring
  offspring_ff = offspring_ff, offspring_fv = offspring_fv, offspring_fp = offspring_fp,
  offspring_vf = offspring_vf, offspring_vv = offspring_vv, offspring_vp = offspring_vp,
  offspring_pf = offspring_pf, offspring_pv = offspring_pv, offspring_pp = offspring_pp,
  mobil_fatal = mobil_fatal,
  mobil_violent = mobil_violent,
  mobil_peaceful = mobil_peaceful,

  # Ratios
  ratio_fv = ratio_fv, ratio_vp = ratio_vp, ratio_fp = ratio_fp,

  # Inference
  ci_fv_lower = ci_fv_lower, ci_fv_upper = ci_fv_upper, p_value_fv = p_value_fv,
  ci_vp_lower = ci_vp_lower, ci_vp_upper = ci_vp_upper, p_value_vp = p_value_vp,
  se = se,
  vcov = vcov_mat,

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

saveRDS(results, "trivariate_hawkes_covariates.rds")
cat("✓ Saved: trivariate_hawkes_covariates.rds\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  TRIVARIATE MODEL WITH COVARIATES COMPLETE                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
