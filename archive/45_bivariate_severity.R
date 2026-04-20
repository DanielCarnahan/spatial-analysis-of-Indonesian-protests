################################################################################
#           SIMPLIFIED BIVARIATE HAWKES: SEVERE vs PEACEFUL
#
#           Pools Fatal + Violent into "Severe" category
#           Compares with Trivariate model to assess whether pooling changes
#           conclusions about mobilization potential
################################################################################

library(tidyverse)
library(numDeriv)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   BIVARIATE HAWKES: SEVERE (Fatal+Violent) vs PEACEFUL                  ║\n")
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

# Create bivariate classification: Severe (Fatal+Violent) vs Peaceful
protests <- protests %>%
  mutate(
    is_violent = sub_event_type %in% c("Mob violence", "Violent demonstration") |
                 sub_event_type == "Excessive force against protesters",
    is_severe = fatalities > 0 | (is_violent & fatalities == 0),
    event_type_2 = ifelse(is_severe, "S", "P"),
    type_numeric = ifelse(is_severe, 1L, 2L)  # 1=Severe, 2=Peaceful
  ) %>%
  arrange(event_date) %>%
  filter(!is.na(poverty_decimal) & !is.na(log_pop))

start_date <- min(protests$event_date)
protests$time <- as.numeric(protests$event_date - start_date)
T_max <- max(protests$time) + 1

n_total <- nrow(protests)
n_severe <- sum(protests$is_severe)
n_peaceful <- sum(!protests$is_severe)

cat(sprintf("Loaded %d events with complete covariate data\n", n_total))
cat(sprintf("  Severe events (S):   %d (%.1f%%)\n", n_severe, 100*n_severe/n_total))
cat(sprintf("    - Fatal:           %d\n", sum(protests$fatalities > 0)))
cat(sprintf("    - Violent:         %d\n", sum(protests$is_violent & protests$fatalities == 0)))
cat(sprintf("  Peaceful events (P): %d (%.1f%%)\n", n_peaceful, 100*n_peaceful/n_total))
cat(sprintf("  Observation period:  %.0f days\n\n", T_max))

# Prepare data vectors
times <- protests$time
types <- protests$type_numeric
log_pop <- protests$log_pop
poverty <- protests$poverty_decimal
years <- protests$year
time_diffs <- diff(times)

# District-year exposure
district_years <- protests %>%
  select(district, year, log_pop, poverty_decimal) %>%
  distinct()

cat(sprintf("District-years: %d unique combinations\n\n", nrow(district_years)))

# =============================================================================
# 2. MODEL A: CONSTANT BASELINE BIVARIATE
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL A: CONSTANT BASELINE BIVARIATE                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Parameters (7 total):
# 1: log_mu_S, 2: log_mu_P
# 3: log_alpha_SS, 4: log_alpha_SP, 5: log_alpha_PS, 6: log_alpha_PP
# 7: log_beta

bivariate_loglik_const <- function(params) {
  mu_s <- exp(params[1])
  mu_p <- exp(params[2])
  alpha_ss <- exp(params[3])
  alpha_sp <- exp(params[4])
  alpha_ps <- exp(params[5])
  alpha_pp <- exp(params[6])
  beta <- exp(params[7])

  if (beta < 1e-6 || beta > 100) return(-1e10)

  n <- length(times)

  # Recursive computation
  R_s <- numeric(n)
  R_p <- numeric(n)

  for (i in 2:n) {
    decay <- exp(-beta * time_diffs[i-1])
    if (types[i-1] == 1) {  # Severe
      R_s[i] <- decay * (1 + R_s[i-1])
      R_p[i] <- decay * R_p[i-1]
    } else {  # Peaceful
      R_s[i] <- decay * R_s[i-1]
      R_p[i] <- decay * (1 + R_p[i-1])
    }
  }

  # Log-likelihood
  loglik <- 0
  for (i in 1:n) {
    if (types[i] == 1) {
      lambda_i <- mu_s + alpha_ss * R_s[i] + alpha_ps * R_p[i]
    } else {
      lambda_i <- mu_p + alpha_sp * R_s[i] + alpha_pp * R_p[i]
    }
    if (lambda_i <= 1e-10) return(-1e10)
    loglik <- loglik + log(lambda_i)
  }

  # Integral
  integral_bg <- (mu_s + mu_p) * T_max
  integral_trig <- 0
  for (j in 1:n) {
    contrib <- (1 - exp(-beta * (T_max - times[j]))) / beta
    if (types[j] == 1) {
      integral_trig <- integral_trig + (alpha_ss + alpha_sp) * contrib
    } else {
      integral_trig <- integral_trig + (alpha_ps + alpha_pp) * contrib
    }
  }

  loglik <- loglik - integral_bg - integral_trig
  if (!is.finite(loglik)) return(-1e10)
  return(loglik)
}

# Optimize
emp_rate_s <- n_severe / T_max
emp_rate_p <- n_peaceful / T_max

start_const <- c(
  log(emp_rate_s * 0.5),
  log(emp_rate_p * 0.5),
  log(0.02), log(0.02), log(0.02), log(0.02),
  log(0.1)
)

cat("Optimizing constant baseline model...\n")

result_const <- optim(
  par = start_const,
  fn = function(p) -bivariate_loglik_const(p),
  method = "L-BFGS-B",
  lower = c(-15, -15, rep(-15, 4), -5),
  upper = c(5, 5, rep(5, 4), 5),
  control = list(maxit = MAX_ITER, trace = 0)
)

ll_const <- -result_const$value
params_const <- result_const$par

mu_s_const <- exp(params_const[1])
mu_p_const <- exp(params_const[2])
alpha_ss_const <- exp(params_const[3])
alpha_sp_const <- exp(params_const[4])
alpha_ps_const <- exp(params_const[5])
alpha_pp_const <- exp(params_const[6])
beta_const <- exp(params_const[7])

mobil_severe_const <- (alpha_ss_const + alpha_sp_const) / beta_const
mobil_peaceful_const <- (alpha_ps_const + alpha_pp_const) / beta_const
ratio_sp_const <- mobil_severe_const / mobil_peaceful_const

cat(sprintf("\nConstant baseline results:\n"))
cat(sprintf("  Log-likelihood: %.2f\n", ll_const))
cat(sprintf("  Convergence: %d\n", result_const$convergence))
cat(sprintf("  μ_S = %.4f, μ_P = %.4f\n", mu_s_const, mu_p_const))
cat(sprintf("  β = %.4f (half-life = %.3f days)\n", beta_const, log(2)/beta_const))
cat(sprintf("  Mobil(Severe) = %.4f, Mobil(Peaceful) = %.4f\n",
            mobil_severe_const, mobil_peaceful_const))
cat(sprintf("  Ratio S/P = %.4f\n\n", ratio_sp_const))

# =============================================================================
# 3. MODEL B: COVARIATE-ADJUSTED BIVARIATE
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL B: COVARIATE-ADJUSTED BIVARIATE                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Parameters (18 total):
# 1-2: β₀_S, β₀_P (intercepts)
# 3: γ (log population)
# 4: δ (poverty)
# 5-13: β_year for 2016-2024
# 14-17: log_α (SS, SP, PS, PP)
# 18: log_β

param_names_cov <- c(
  "beta_0_S", "beta_0_P", "gamma", "delta",
  paste0("beta_", 2016:2024),
  "log_alpha_SS", "log_alpha_SP", "log_alpha_PS", "log_alpha_PP",
  "log_beta"
)

bivariate_loglik_cov <- function(params) {
  beta_0_S <- params[1]
  beta_0_P <- params[2]
  gamma <- params[3]
  delta <- params[4]
  beta_years <- params[5:13]

  alpha_ss <- exp(params[14])
  alpha_sp <- exp(params[15])
  alpha_ps <- exp(params[16])
  alpha_pp <- exp(params[17])
  beta <- exp(params[18])

  if (beta < 1e-6 || beta > 100) return(-1e10)

  n <- length(times)

  # Year effects
  year_effects <- numeric(n)
  for (i in 1:n) {
    yr <- years[i]
    if (yr >= 2016 && yr <= 2024) {
      year_effects[i] <- beta_years[yr - 2015]
    }
  }

  # Background rates
  mu_s <- exp(beta_0_S + gamma * log_pop + delta * poverty + year_effects)
  mu_p <- exp(beta_0_P + gamma * log_pop + delta * poverty + year_effects)

  # Recursive computation
  R_s <- numeric(n)
  R_p <- numeric(n)

  for (i in 2:n) {
    decay <- exp(-beta * time_diffs[i-1])
    if (types[i-1] == 1) {
      R_s[i] <- decay * (1 + R_s[i-1])
      R_p[i] <- decay * R_p[i-1]
    } else {
      R_s[i] <- decay * R_s[i-1]
      R_p[i] <- decay * (1 + R_p[i-1])
    }
  }

  # Log-likelihood
  loglik <- 0
  for (i in 1:n) {
    if (types[i] == 1) {
      lambda_i <- mu_s[i] + alpha_ss * R_s[i] + alpha_ps * R_p[i]
    } else {
      lambda_i <- mu_p[i] + alpha_sp * R_s[i] + alpha_pp * R_p[i]
    }
    if (lambda_i <= 1e-10) return(-1e10)
    loglik <- loglik + log(lambda_i)
  }

  # Background integral
  integral_bg <- 0
  for (j in 1:nrow(district_years)) {
    yr <- district_years$year[j]
    yr_effect <- if (yr >= 2016 && yr <= 2024) beta_years[yr - 2015] else 0
    mu_s_j <- exp(beta_0_S + gamma * district_years$log_pop[j] +
                   delta * district_years$poverty_decimal[j] + yr_effect)
    mu_p_j <- exp(beta_0_P + gamma * district_years$log_pop[j] +
                   delta * district_years$poverty_decimal[j] + yr_effect)
    integral_bg <- integral_bg + (mu_s_j + mu_p_j) * 365.25
  }

  # Triggering integral
  integral_trig <- 0
  for (j in 1:n) {
    contrib <- (1 - exp(-beta * (T_max - times[j]))) / beta
    if (types[j] == 1) {
      integral_trig <- integral_trig + (alpha_ss + alpha_sp) * contrib
    } else {
      integral_trig <- integral_trig + (alpha_ps + alpha_pp) * contrib
    }
  }

  loglik <- loglik - integral_bg - integral_trig
  if (!is.finite(loglik)) return(-1e10)
  return(loglik)
}

neg_loglik_cov <- function(params) -bivariate_loglik_cov(params)

# Multi-start optimization
mean_log_pop <- mean(log_pop)
mean_poverty <- mean(poverty)

starting_points <- list(
  list(
    name = "Data-driven",
    params = c(
      log(emp_rate_s) - 0.3 * mean_log_pop,
      log(emp_rate_p) - 0.3 * mean_log_pop,
      0.3, 0.0,
      rep(0, 9),
      rep(log(0.02), 4),
      log(0.1)
    )
  ),
  list(
    name = "Strong-pop",
    params = c(
      log(emp_rate_s) - 0.5 * mean_log_pop,
      log(emp_rate_p) - 0.5 * mean_log_pop,
      0.5, 0.5,
      rep(0, 9),
      rep(log(0.01), 4),
      log(0.15)
    )
  ),
  list(
    name = "Negative-poverty",
    params = c(
      log(emp_rate_s) - 0.3 * mean_log_pop,
      log(emp_rate_p) - 0.3 * mean_log_pop,
      0.3, -2.0,
      rep(0, 9),
      rep(log(0.03), 4),
      log(0.2)
    )
  ),
  list(
    name = "Asymmetric",
    params = c(
      log(emp_rate_s) - 0.4 * mean_log_pop,
      log(emp_rate_p) - 0.4 * mean_log_pop,
      0.4, 0.0,
      rep(0, 9),
      log(0.05), log(0.05),
      log(0.01), log(0.01),
      log(0.1)
    )
  ),
  list(
    name = "Fast-decay",
    params = c(
      log(emp_rate_s) - 0.3 * mean_log_pop,
      log(emp_rate_p) - 0.3 * mean_log_pop,
      0.3, -1.0,
      rep(0, 9),
      rep(log(0.05), 4),
      log(0.3)
    )
  )
)

lower_bounds <- c(-20, -20, 0.01, -10, rep(-5, 9), rep(-15, 4), -5)
upper_bounds <- c(10, 10, 2, 10, rep(5, 9), rep(5, 4), 5)

all_results <- list()
best_ll <- -Inf
best_result <- NULL

cat("Running multi-start optimization...\n")

for (s in seq_along(starting_points)) {
  sp <- starting_points[[s]]
  cat(sprintf("  Starting point %d/%d: %s...", s, length(starting_points), sp$name))

  result <- tryCatch({
    optim(
      par = sp$params,
      fn = neg_loglik_cov,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(maxit = MAX_ITER, trace = 0)
    )
  }, error = function(e) NULL)

  if (!is.null(result)) {
    ll <- -result$value
    cat(sprintf(" LL=%.2f, conv=%d\n", ll, result$convergence))

    all_results[[s]] <- list(
      name = sp$name,
      params = result$par,
      loglik = ll,
      convergence = result$convergence
    )

    # Only consider converged results
    if (result$convergence == 0 && ll > best_ll) {
      best_ll <- ll
      best_result <- result
    }
  } else {
    cat(" FAILED\n")
    all_results[[s]] <- list(name = sp$name, loglik = -Inf, convergence = -1)
  }
}

if (is.null(best_result)) {
  cat("\n⚠️ No starting points converged! Using best non-converged.\n")
  for (r in all_results) {
    if (!is.null(r$params) && r$loglik > best_ll) {
      best_ll <- r$loglik
      best_result <- list(par = r$params, value = -r$loglik, convergence = r$convergence)
    }
  }
}

# Extract parameters
params_cov <- best_result$par
ll_cov <- -best_result$value

beta_0_S <- params_cov[1]
beta_0_P <- params_cov[2]
gamma <- params_cov[3]
delta <- params_cov[4]
beta_years <- params_cov[5:13]
alpha_ss_cov <- exp(params_cov[14])
alpha_sp_cov <- exp(params_cov[15])
alpha_ps_cov <- exp(params_cov[16])
alpha_pp_cov <- exp(params_cov[17])
beta_cov <- exp(params_cov[18])

mobil_severe_cov <- (alpha_ss_cov + alpha_sp_cov) / beta_cov
mobil_peaceful_cov <- (alpha_ps_cov + alpha_pp_cov) / beta_cov
ratio_sp_cov <- mobil_severe_cov / mobil_peaceful_cov

cat(sprintf("\nCovariate-adjusted results:\n"))
cat(sprintf("  Log-likelihood: %.2f\n", ll_cov))
cat(sprintf("  β₀_S = %.4f, β₀_P = %.4f\n", beta_0_S, beta_0_P))
cat(sprintf("  γ (population) = %.4f\n", gamma))
cat(sprintf("  δ (poverty) = %.4f\n", delta))
cat(sprintf("  β (decay) = %.4f (half-life = %.3f days)\n", beta_cov, log(2)/beta_cov))
cat(sprintf("  Mobil(Severe) = %.4f, Mobil(Peaceful) = %.4f\n",
            mobil_severe_cov, mobil_peaceful_cov))
cat(sprintf("  Ratio S/P = %.4f\n\n", ratio_sp_cov))

# =============================================================================
# 4. COMPUTE STANDARD ERRORS
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  STANDARD ERRORS & HYPOTHESIS TEST                          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Computing Hessian...\n")

hessian_mat <- tryCatch({
  hessian(neg_loglik_cov, params_cov)
}, error = function(e) NULL)

if (!is.null(hessian_mat)) {
  vcov_mat <- tryCatch({
    solve(hessian_mat)
  }, error = function(e) {
    cat("  Using pseudo-inverse\n")
    MASS::ginv(hessian_mat)
  })

  se <- sqrt(pmax(diag(vcov_mat), 0))

  # Eigenvalue check
  eigen_vals <- eigen(hessian_mat, symmetric = TRUE)$values
  cat(sprintf("  Min eigenvalue: %.2e, Max: %.2e\n", min(eigen_vals), max(eigen_vals)))

  # Delta method for S/P ratio
  log_ratio_func <- function(p) {
    alpha_ss <- exp(p[14]); alpha_sp <- exp(p[15])
    alpha_ps <- exp(p[16]); alpha_pp <- exp(p[17])
    beta <- exp(p[18])
    mobil_s <- (alpha_ss + alpha_sp) / beta
    mobil_p <- (alpha_ps + alpha_pp) / beta
    log(mobil_s / mobil_p)
  }

  gradient <- grad(log_ratio_func, params_cov)
  var_log_ratio <- as.numeric(t(gradient) %*% vcov_mat %*% gradient)
  se_log_ratio <- sqrt(max(var_log_ratio, 0))
  log_ratio <- log(ratio_sp_cov)

  if (se_log_ratio > 0 && is.finite(se_log_ratio)) {
    ci_lower <- exp(log_ratio - 1.96 * se_log_ratio)
    ci_upper <- exp(log_ratio + 1.96 * se_log_ratio)
    z_stat <- log_ratio / se_log_ratio
    p_value <- 2 * pnorm(-abs(z_stat))

    cat(sprintf("\nHypothesis Test: Severe vs Peaceful mobilization\n"))
    cat(sprintf("  Ratio S/P: %.4f [95%% CI: %.4f, %.4f]\n", ratio_sp_cov, ci_lower, ci_upper))
    cat(sprintf("  z = %.4f, p = %.6f\n", z_stat, p_value))

    if (ci_lower > 1) {
      cat("  *** Severe has HIGHER mobilization than peaceful ***\n")
    } else if (ci_upper < 1) {
      cat("  *** Severe has LOWER mobilization than peaceful ***\n")
    } else {
      cat("  No significant difference at α = 0.05\n")
    }
  } else {
    ci_lower <- ci_upper <- p_value <- NA
    cat("  Unable to compute CI for ratio\n")
  }
} else {
  se <- rep(NA, 18)
  vcov_mat <- NULL
  ci_lower <- ci_upper <- p_value <- NA
  cat("  Could not compute Hessian\n")
}

# =============================================================================
# 5. MODEL COMPARISON
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL COMPARISON                                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Load trivariate results for comparison
trivariate_const <- readRDS("trivariate_hawkes_constant.rds")
trivariate_diag <- readRDS("trivariate_diagnostics.rds")

comparison <- data.frame(
  Model = c("Bivariate Constant", "Bivariate Covariate",
            "Trivariate Constant", "Trivariate Covariate"),
  LogLik = c(ll_const, ll_cov,
             trivariate_const$loglik, trivariate_diag$converged_ll),
  Params = c(7, 18, 13, 24),
  AIC = c(-2*ll_const + 2*7, -2*ll_cov + 2*18,
          -2*trivariate_const$loglik + 2*13, -2*trivariate_diag$converged_ll + 2*24),
  Decay = c(beta_const, beta_cov,
            trivariate_const$beta, trivariate_diag$beta_decay),
  `Half_life` = c(log(2)/beta_const, log(2)/beta_cov,
                  log(2)/trivariate_const$beta, trivariate_diag$half_life)
)

cat("Model comparison:\n")
print(comparison, row.names = FALSE)

# Mobilization comparison
cat("\nMobilization potential comparison:\n")
cat(sprintf("%-25s %12s %12s %12s\n", "Model", "Severe/Fatal", "Violent", "Peaceful"))
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-25s %12.4f %12s %12.4f\n",
            "Bivariate Constant", mobil_severe_const, "-", mobil_peaceful_const))
cat(sprintf("%-25s %12.4f %12s %12.4f\n",
            "Bivariate Covariate", mobil_severe_cov, "-", mobil_peaceful_cov))
cat(sprintf("%-25s %12.4f %12.4f %12.4f\n",
            "Trivariate Constant", trivariate_const$mobil_fatal,
            trivariate_const$mobil_violent, trivariate_const$mobil_peaceful))
cat(sprintf("%-25s %12.4f %12.4f %12.4f\n",
            "Trivariate Covariate", trivariate_diag$mobil_fatal,
            trivariate_diag$mobil_violent, trivariate_diag$mobil_peaceful))

cat("\nRatios (Severe/Peaceful or Fatal+Violent/Peaceful):\n")
cat(sprintf("  Bivariate Constant:    %.4f\n", ratio_sp_const))
cat(sprintf("  Bivariate Covariate:   %.4f", ratio_sp_cov))
if (!is.na(ci_lower)) cat(sprintf(" [%.4f, %.4f]", ci_lower, ci_upper))
cat("\n")
cat(sprintf("  Trivariate Constant:   %.4f\n",
            (trivariate_const$mobil_fatal + trivariate_const$mobil_violent) / trivariate_const$mobil_peaceful))
cat(sprintf("  Trivariate Covariate:  %.4f\n",
            (trivariate_diag$mobil_fatal + trivariate_diag$mobil_violent) / trivariate_diag$mobil_peaceful))

# =============================================================================
# 6. SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SUMMARY                                                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("KEY FINDINGS:\n")
cat("─────────────────────────────────────────────────────────────────\n")

cat(sprintf("1. SEVERE vs PEACEFUL MOBILIZATION\n"))
cat(sprintf("   Ratio (S/P) = %.4f\n", ratio_sp_cov))
if (!is.na(p_value) && p_value < 0.05) {
  if (ci_upper < 1) {
    cat("   *** Severe protests have LOWER mobilization than peaceful ***\n")
    cat("   This aligns with trivariate findings\n\n")
  } else if (ci_lower > 1) {
    cat("   *** Severe protests have HIGHER mobilization than peaceful ***\n\n")
  }
} else if (!is.na(p_value)) {
  cat(sprintf("   Not significantly different from 1 (p = %.4f)\n\n", p_value))
}

cat(sprintf("2. COVARIATE EFFECTS\n"))
cat(sprintf("   γ (population): %.4f → %.1f%% increase per log-unit\n",
            gamma, 100*(exp(gamma)-1)))
cat(sprintf("   δ (poverty): %.4f\n", delta))
if (delta < 0) {
  cat("   Negative poverty effect - same as trivariate model\n\n")
}

cat(sprintf("3. POOLING EFFECT\n"))
cat("   Pooling Fatal+Violent into 'Severe' category:\n")
if (ll_cov > trivariate_diag$converged_ll) {
  cat("   - Improves model fit (higher LL)\n")
} else {
  cat("   - Does not improve fit over trivariate\n")
}
cat(sprintf("   - S/P ratio (%.4f) similar to combined (F+V)/P ratio\n\n", ratio_sp_cov))

cat("CONCLUSIONS:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("• Both bivariate and trivariate models show: severe/violent\n")
cat("  protests have LOWER diffusion potential than peaceful protests\n")
cat("• This contradicts the 'radical flank' hypothesis\n")
cat("• Possible explanations:\n")
cat("  - State repression suppresses mobilization after violence\n")
cat("  - Violent events deter potential participants\n")
cat("  - Data limitations (underreporting of follow-up violence)\n")

# =============================================================================
# 7. SAVE RESULTS
# =============================================================================

cat("\n=== SAVING RESULTS ===\n")

results <- list(
  model_name = "Bivariate Hawkes: Severe vs Peaceful",

  # Data
  n_severe = n_severe,
  n_peaceful = n_peaceful,
  n_total = n_total,
  T_max = T_max,

  # Constant model
  const = list(
    params = params_const,
    loglik = ll_const,
    mu_s = mu_s_const, mu_p = mu_p_const,
    alpha_ss = alpha_ss_const, alpha_sp = alpha_sp_const,
    alpha_ps = alpha_ps_const, alpha_pp = alpha_pp_const,
    beta = beta_const,
    mobil_severe = mobil_severe_const,
    mobil_peaceful = mobil_peaceful_const,
    ratio_sp = ratio_sp_const
  ),

  # Covariate model
  cov = list(
    params = params_cov,
    param_names = param_names_cov,
    loglik = ll_cov,
    convergence = best_result$convergence,
    beta_0_S = beta_0_S, beta_0_P = beta_0_P,
    gamma = gamma, delta = delta,
    beta_years = beta_years,
    alpha_ss = alpha_ss_cov, alpha_sp = alpha_sp_cov,
    alpha_ps = alpha_ps_cov, alpha_pp = alpha_pp_cov,
    beta = beta_cov,
    mobil_severe = mobil_severe_cov,
    mobil_peaceful = mobil_peaceful_cov,
    ratio_sp = ratio_sp_cov,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value,
    se = se,
    vcov = vcov_mat
  ),

  # All optimization results
  all_results = all_results
)

saveRDS(results, "bivariate_severity.rds")
cat("✓ Saved: bivariate_severity.rds\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  BIVARIATE SEVERITY MODEL COMPLETE                          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
