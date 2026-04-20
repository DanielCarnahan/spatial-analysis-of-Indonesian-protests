################################################################################
#           TRIVARIATE HAWKES MODEL DIAGNOSTICS
#
#           Investigation of model behavior, convergence, and covariate effects
################################################################################

library(tidyverse)
library(numDeriv)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   TRIVARIATE HAWKES MODEL DIAGNOSTICS                                    ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. LOAD SAVED RESULTS
# =============================================================================

cat("=== LOADING SAVED RESULTS ===\n\n")

results_cov <- readRDS("trivariate_hawkes_covariates.rds")
results_const <- readRDS("trivariate_hawkes_constant.rds")

cat("Covariate model:\n")
cat(sprintf("  Best starting point: %s\n", results_cov$best_starting_point))
cat(sprintf("  Log-likelihood: %.2f\n", results_cov$loglik))
cat(sprintf("  γ (population): %.4f\n", results_cov$gamma))
cat(sprintf("  δ (poverty): %.4f\n\n", results_cov$delta))

cat("Constant model:\n")
cat(sprintf("  Log-likelihood: %.2f\n\n", results_const$loglik))

# =============================================================================
# 2. CHECK CONVERGENCE OF ALL STARTING POINTS
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  CONVERGENCE ANALYSIS                                        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Starting point results:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("%-25s %12s %8s %s\n", "Starting Point", "Log-lik", "Conv", "Status"))
cat("─────────────────────────────────────────────────────────────────\n")

converged_results <- list()
for (r in results_cov$all_results) {
  conv_status <- if (r$convergence == 0) "✓ Converged" else sprintf("✗ Code %d", r$convergence)
  marker <- if (r$starting_point == results_cov$best_starting_point) " ← SELECTED" else ""
  cat(sprintf("%-25s %12.2f %8d %s%s\n",
              r$starting_point, r$loglik, r$convergence, conv_status, marker))

  if (r$convergence == 0 && !is.null(r$params)) {
    converged_results[[length(converged_results) + 1]] <- r
  }
}

cat("\n")
cat(sprintf("Total starting points: %d\n", length(results_cov$all_results)))
cat(sprintf("Converged (code 0): %d\n", length(converged_results)))

# Find best CONVERGED result
if (length(converged_results) > 0) {
  best_converged_ll <- -Inf
  best_converged <- NULL
  for (r in converged_results) {
    if (r$loglik > best_converged_ll) {
      best_converged_ll <- r$loglik
      best_converged <- r
    }
  }

  cat(sprintf("\n*** BEST CONVERGED: %s (LL = %.2f) ***\n",
              best_converged$starting_point, best_converged$loglik))

  if (best_converged$starting_point != results_cov$best_starting_point) {
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════════╗\n")
    cat("║  ⚠️  WARNING: SELECTED RESULT DID NOT CONVERGE!              ║\n")
    cat("║  The saved model used a non-converged optimization.         ║\n")
    cat("║  Results below use the best CONVERGED estimate.             ║\n")
    cat("╚══════════════════════════════════════════════════════════════╝\n\n")

    params_to_use <- best_converged$params
    ll_to_use <- best_converged$loglik
    use_converged <- TRUE
  } else {
    params_to_use <- results_cov$params
    ll_to_use <- results_cov$loglik
    use_converged <- FALSE
  }
} else {
  cat("\n⚠️  No starting points converged!\n")
  params_to_use <- results_cov$params
  ll_to_use <- results_cov$loglik
  use_converged <- FALSE
}

# =============================================================================
# 3. EXTRACT PARAMETERS FROM CONVERGED MODEL
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  CONVERGED MODEL PARAMETERS                                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Background parameters
beta_0_F <- params_to_use[1]
beta_0_V <- params_to_use[2]
beta_0_P <- params_to_use[3]
gamma <- params_to_use[4]
delta <- params_to_use[5]
beta_years <- params_to_use[6:14]

# Triggering parameters
alpha_ff <- exp(params_to_use[15])
alpha_fv <- exp(params_to_use[16])
alpha_fp <- exp(params_to_use[17])
alpha_vf <- exp(params_to_use[18])
alpha_vv <- exp(params_to_use[19])
alpha_vp <- exp(params_to_use[20])
alpha_pf <- exp(params_to_use[21])
alpha_pv <- exp(params_to_use[22])
alpha_pp <- exp(params_to_use[23])
beta_decay <- exp(params_to_use[24])

cat("BACKGROUND RATE:\n")
cat(sprintf("  β₀_F (fatal intercept): %.4f\n", beta_0_F))
cat(sprintf("  β₀_V (violent intercept): %.4f\n", beta_0_V))
cat(sprintf("  β₀_P (peaceful intercept): %.4f\n", beta_0_P))
cat(sprintf("  γ (log population): %.4f\n", gamma))
cat(sprintf("  δ (poverty): %.4f\n\n", delta))

cat("YEAR EFFECTS (relative to 2015):\n")
for (i in 1:9) {
  cat(sprintf("  β_%d: %+.4f\n", 2015 + i, beta_years[i]))
}

cat("\nDECAY RATE:\n")
cat(sprintf("  β = %.4f\n", beta_decay))
cat(sprintf("  Half-life = %.3f days (%.1f hours)\n\n", log(2)/beta_decay, log(2)/beta_decay * 24))

# Offspring matrix
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

cat("OFFSPRING MATRIX (α/β):\n")
cat("                       → Fatal    → Violent  → Peaceful  TOTAL\n")
cat(sprintf("  Fatal parent:       %.5f    %.5f    %.5f   %.5f\n",
            offspring_ff, offspring_fv, offspring_fp, mobil_fatal))
cat(sprintf("  Violent parent:     %.5f    %.5f    %.5f   %.5f\n",
            offspring_vf, offspring_vv, offspring_vp, mobil_violent))
cat(sprintf("  Peaceful parent:    %.5f    %.5f    %.5f   %.5f\n\n",
            offspring_pf, offspring_pv, offspring_pp, mobil_peaceful))

cat("MOBILIZATION RATIOS:\n")
cat(sprintf("  Fatal / Violent:    %.4f\n", mobil_fatal / mobil_violent))
cat(sprintf("  Violent / Peaceful: %.4f\n", mobil_violent / mobil_peaceful))
cat(sprintf("  Fatal / Peaceful:   %.4f\n", mobil_fatal / mobil_peaceful))

# =============================================================================
# 4. COMPARE CONVERGED VS NON-CONVERGED
# =============================================================================

if (use_converged) {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  CONVERGED vs NON-CONVERGED COMPARISON                       ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")

  params_nonconv <- results_cov$params

  cat(sprintf("%-25s %15s %15s\n", "Parameter", "Converged", "Non-converged"))
  cat("─────────────────────────────────────────────────────────────────\n")
  cat(sprintf("%-25s %15.2f %15.2f\n", "Log-likelihood", ll_to_use, results_cov$loglik))
  cat(sprintf("%-25s %15.4f %15.4f\n", "β₀_F", beta_0_F, params_nonconv[1]))
  cat(sprintf("%-25s %15.4f %15.4f\n", "β₀_V", beta_0_V, params_nonconv[2]))
  cat(sprintf("%-25s %15.4f %15.4f\n", "β₀_P", beta_0_P, params_nonconv[3]))
  cat(sprintf("%-25s %15.4f %15.4f\n", "γ (population)", gamma, params_nonconv[4]))
  cat(sprintf("%-25s %15.4f %15.4f\n", "δ (poverty)", delta, params_nonconv[5]))
  cat(sprintf("%-25s %15.4f %15.4f\n", "β (decay)", beta_decay, exp(params_nonconv[24])))

  cat("\n")
  cat("The non-converged result has a higher log-likelihood because\n")
  cat("the optimizer got stuck or found an invalid region.\n")
  cat("Always use converged results for inference!\n")
}

# =============================================================================
# 5. HESSIAN ANALYSIS (AT CONVERGED PARAMETERS)
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  HESSIAN ANALYSIS                                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Need to reload data and recompute likelihood for Hessian
cat("Loading data and recomputing likelihood...\n")

protests <- readRDS("protests_with_poverty.rds")
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

start_date <- min(protests$event_date)
protests$time <- as.numeric(protests$event_date - start_date)
T_max <- max(protests$time) + 1

times <- protests$time
types <- protests$type_numeric
log_pop <- protests$log_pop
poverty <- protests$poverty_decimal
years <- protests$year
time_diffs <- diff(times)

district_years <- protests %>%
  select(district, year, log_pop, poverty_decimal) %>%
  distinct()

# Define likelihood function (same as in main script)
trivariate_loglik_cov <- function(params) {
  beta_0_F <- params[1]; beta_0_V <- params[2]; beta_0_P <- params[3]
  gamma <- params[4]; delta <- params[5]
  beta_years <- params[6:14]

  alpha_ff <- exp(params[15]); alpha_fv <- exp(params[16]); alpha_fp <- exp(params[17])
  alpha_vf <- exp(params[18]); alpha_vv <- exp(params[19]); alpha_vp <- exp(params[20])
  alpha_pf <- exp(params[21]); alpha_pv <- exp(params[22]); alpha_pp <- exp(params[23])
  beta_decay <- exp(params[24])

  if (beta_decay < 1e-6 || beta_decay > 100) return(-1e10)
  all_alpha <- c(alpha_ff, alpha_fv, alpha_fp, alpha_vf, alpha_vv, alpha_vp, alpha_pf, alpha_pv, alpha_pp)
  if (any(!is.finite(all_alpha))) return(-1e10)

  n <- length(times)
  year_effects <- numeric(n)
  for (i in 1:n) {
    yr <- years[i]
    if (yr >= 2016 && yr <= 2024) year_effects[i] <- beta_years[yr - 2015]
  }

  mu_f <- exp(beta_0_F + gamma * log_pop + delta * poverty + year_effects)
  mu_v <- exp(beta_0_V + gamma * log_pop + delta * poverty + year_effects)
  mu_p <- exp(beta_0_P + gamma * log_pop + delta * poverty + year_effects)

  R_f <- numeric(n); R_v <- numeric(n); R_p <- numeric(n)
  for (i in 2:n) {
    decay <- exp(-beta_decay * time_diffs[i-1])
    if (decay < 1e-100) decay <- 0
    if (types[i-1] == 1) {
      R_f[i] <- decay * (1 + R_f[i-1]); R_v[i] <- decay * R_v[i-1]; R_p[i] <- decay * R_p[i-1]
    } else if (types[i-1] == 2) {
      R_f[i] <- decay * R_f[i-1]; R_v[i] <- decay * (1 + R_v[i-1]); R_p[i] <- decay * R_p[i-1]
    } else {
      R_f[i] <- decay * R_f[i-1]; R_v[i] <- decay * R_v[i-1]; R_p[i] <- decay * (1 + R_p[i-1])
    }
  }

  loglik <- 0
  for (i in 1:n) {
    if (types[i] == 1) {
      lambda_i <- mu_f[i] + alpha_ff * R_f[i] + alpha_vf * R_v[i] + alpha_pf * R_p[i]
    } else if (types[i] == 2) {
      lambda_i <- mu_v[i] + alpha_fv * R_f[i] + alpha_vv * R_v[i] + alpha_pv * R_p[i]
    } else {
      lambda_i <- mu_p[i] + alpha_fp * R_f[i] + alpha_vp * R_v[i] + alpha_pp * R_p[i]
    }
    if (lambda_i <= 1e-10) return(-1e10)
    loglik <- loglik + log(lambda_i)
  }

  # Background integral
  integral_bg <- 0
  for (j in 1:nrow(district_years)) {
    yr <- district_years$year[j]
    yr_effect <- if (yr >= 2016 && yr <= 2024) beta_years[yr - 2015] else 0
    mu_f_j <- exp(beta_0_F + gamma * district_years$log_pop[j] + delta * district_years$poverty_decimal[j] + yr_effect)
    mu_v_j <- exp(beta_0_V + gamma * district_years$log_pop[j] + delta * district_years$poverty_decimal[j] + yr_effect)
    mu_p_j <- exp(beta_0_P + gamma * district_years$log_pop[j] + delta * district_years$poverty_decimal[j] + yr_effect)
    integral_bg <- integral_bg + (mu_f_j + mu_v_j + mu_p_j) * 365.25
  }

  # Triggering integral
  integral_trig <- 0
  for (j in 1:n) {
    contrib <- (1 - exp(-beta_decay * (T_max - times[j]))) / beta_decay
    if (types[j] == 1) {
      integral_trig <- integral_trig + (alpha_ff + alpha_fv + alpha_fp) * contrib
    } else if (types[j] == 2) {
      integral_trig <- integral_trig + (alpha_vf + alpha_vv + alpha_vp) * contrib
    } else {
      integral_trig <- integral_trig + (alpha_pf + alpha_pv + alpha_pp) * contrib
    }
  }

  loglik <- loglik - integral_bg - integral_trig
  if (!is.finite(loglik)) return(-1e10)
  return(loglik)
}

neg_loglik <- function(params) -trivariate_loglik_cov(params)

# Compute Hessian at converged parameters
cat("Computing Hessian at converged MLE...\n\n")

hessian_mat <- tryCatch({
  hessian(neg_loglik, params_to_use)
}, error = function(e) {
  cat(sprintf("  Error: %s\n", e$message))
  NULL
})

if (!is.null(hessian_mat)) {
  # Eigenvalue analysis
  eigen_result <- eigen(hessian_mat, symmetric = TRUE)
  eigenvalues <- eigen_result$values

  cat("HESSIAN EIGENVALUES:\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat(sprintf("  Min eigenvalue: %.4e\n", min(eigenvalues)))
  cat(sprintf("  Max eigenvalue: %.4e\n", max(eigenvalues)))
  cat(sprintf("  Condition number: %.4e\n\n", max(eigenvalues) / max(min(eigenvalues), 1e-10)))

  # Check for positive definiteness
  n_negative <- sum(eigenvalues < 0)
  n_small <- sum(abs(eigenvalues) < 1e-6)

  if (n_negative > 0) {
    cat(sprintf("  ⚠️  %d NEGATIVE eigenvalues - NOT at a local maximum!\n", n_negative))
  } else if (n_small > 0) {
    cat(sprintf("  ⚠️  %d near-zero eigenvalues - possible identification issues\n", n_small))
  } else {
    cat("  ✓ All eigenvalues positive - valid maximum\n")
  }

  cat("\nEigenvalue distribution:\n")
  eigen_breaks <- c(-Inf, -1e-6, 1e-6, 1, 100, 1000, Inf)
  eigen_labels <- c("Negative", "Near zero", "Small", "Medium", "Large", "Very large")
  eigen_table <- cut(eigenvalues, breaks = eigen_breaks, labels = eigen_labels)
  print(table(eigen_table))

  # Compute variance-covariance matrix
  vcov_mat <- tryCatch({
    solve(hessian_mat)
  }, error = function(e) {
    cat("\n  Using pseudo-inverse due to singular Hessian\n")
    MASS::ginv(hessian_mat)
  })

  # Standard errors
  se <- sqrt(pmax(diag(vcov_mat), 0))

  param_names <- c(
    "beta_0_F", "beta_0_V", "beta_0_P", "gamma", "delta",
    paste0("beta_", 2016:2024),
    "log_alpha_FF", "log_alpha_FV", "log_alpha_FP",
    "log_alpha_VF", "log_alpha_VV", "log_alpha_VP",
    "log_alpha_PF", "log_alpha_PV", "log_alpha_PP",
    "log_beta"
  )

  cat("\n")
  cat("PARAMETER STANDARD ERRORS:\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  cat(sprintf("%-18s %12s %12s %8s\n", "Parameter", "Estimate", "SE", "t-stat"))
  cat("─────────────────────────────────────────────────────────────────\n")

  for (i in 1:24) {
    t_stat <- if (se[i] > 0) params_to_use[i] / se[i] else NA
    flag <- if (is.na(t_stat) || abs(t_stat) < 1.96) "" else "*"
    cat(sprintf("%-18s %12.4f %12.4f %8.2f %s\n",
                param_names[i], params_to_use[i], se[i], ifelse(is.na(t_stat), NA, t_stat), flag))
  }
  cat("\n* indicates |t| > 1.96 (significant at 5%)\n")

  # Parameter correlation matrix
  cat("\n")
  cat("PARAMETER CORRELATIONS (|r| > 0.8):\n")
  cat("─────────────────────────────────────────────────────────────────\n")

  cor_mat <- cov2cor(vcov_mat)
  high_cor <- which(abs(cor_mat) > 0.8 & upper.tri(cor_mat), arr.ind = TRUE)

  if (nrow(high_cor) > 0) {
    for (k in 1:nrow(high_cor)) {
      i <- high_cor[k, 1]
      j <- high_cor[k, 2]
      cat(sprintf("  %s ↔ %s: r = %.3f\n",
                  param_names[i], param_names[j], cor_mat[i, j]))
    }
  } else {
    cat("  No correlations > 0.8\n")
  }

} else {
  cat("Could not compute Hessian\n")
  se <- rep(NA, 24)
}

# =============================================================================
# 6. POVERTY AND POPULATION ANALYSIS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  COVARIATE ANALYSIS: POVERTY AND POPULATION                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("COVARIATE SUMMARY:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Population (log): mean=%.2f, sd=%.2f, range=[%.2f, %.2f]\n",
            mean(log_pop), sd(log_pop), min(log_pop), max(log_pop)))
cat(sprintf("  Poverty (decimal): mean=%.3f, sd=%.3f, range=[%.3f, %.3f]\n",
            mean(poverty), sd(poverty), min(poverty), max(poverty)))
cat(sprintf("  Correlation(log_pop, poverty): r = %.3f\n\n", cor(log_pop, poverty)))

# Cross-tabulate poverty with event types
protests$poverty_quartile <- cut(protests$poverty_decimal,
                                  breaks = quantile(protests$poverty_decimal, c(0, 0.25, 0.5, 0.75, 1)),
                                  labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"),
                                  include.lowest = TRUE)

cat("EVENT COUNTS BY POVERTY QUARTILE:\n")
poverty_table <- protests %>%
  group_by(poverty_quartile, event_type_3) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = event_type_3, values_from = n, values_fill = 0) %>%
  mutate(Total = F + V + P,
         Pct_Fatal = 100 * F / Total,
         Pct_Violent = 100 * V / Total)

print(poverty_table)

# Event rates by poverty quartile (controlling for exposure)
cat("\nEVENT RATES PER DISTRICT-YEAR BY POVERTY QUARTILE:\n")

district_year_counts <- protests %>%
  group_by(district, year, poverty_quartile) %>%
  summarise(
    n_fatal = sum(event_type_3 == "F"),
    n_violent = sum(event_type_3 == "V"),
    n_peaceful = sum(event_type_3 == "P"),
    n_total = n(),
    log_pop = first(log_pop),
    .groups = "drop"
  )

rate_by_poverty <- district_year_counts %>%
  group_by(poverty_quartile) %>%
  summarise(
    n_district_years = n(),
    mean_fatal = mean(n_fatal),
    mean_violent = mean(n_violent),
    mean_peaceful = mean(n_peaceful),
    mean_total = mean(n_total),
    .groups = "drop"
  )

print(rate_by_poverty)

# Population quartile analysis
protests$pop_quartile <- cut(protests$log_pop,
                              breaks = quantile(protests$log_pop, c(0, 0.25, 0.5, 0.75, 1)),
                              labels = c("Q1 (Small)", "Q2", "Q3", "Q4 (Large)"),
                              include.lowest = TRUE)

cat("\nEVENT COUNTS BY POPULATION QUARTILE:\n")
pop_table <- protests %>%
  group_by(pop_quartile, event_type_3) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = event_type_3, values_from = n, values_fill = 0) %>%
  mutate(Total = F + V + P)

print(pop_table)

# =============================================================================
# 7. YEAR EFFECTS VISUALIZATION
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  YEAR EFFECTS                                                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

year_df <- data.frame(
  year = 2015:2024,
  effect = c(0, beta_years),  # 2015 is reference (effect = 0)
  se = c(0, if(exists("se") && !all(is.na(se))) se[6:14] else rep(NA, 9))
)

year_df$lower <- year_df$effect - 1.96 * year_df$se
year_df$upper <- year_df$effect + 1.96 * year_df$se

cat("Year effects (relative to 2015):\n")
print(year_df)

# =============================================================================
# 8. MODEL COMPARISON: CONSTANT VS COVARIATE
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MODEL COMPARISON: CONSTANT vs COVARIATE                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Model comparison table
comparison <- data.frame(
  Metric = c("Log-likelihood", "Parameters", "AIC", "BIC",
             "Mobil (Fatal)", "Mobil (Violent)", "Mobil (Peaceful)",
             "Ratio F/V", "Ratio V/P", "Decay (β)", "Half-life (days)"),
  Constant = c(
    results_const$loglik,
    results_const$n_params,
    -2 * results_const$loglik + 2 * results_const$n_params,
    -2 * results_const$loglik + results_const$n_params * log(results_const$n_total),
    results_const$mobil_fatal,
    results_const$mobil_violent,
    results_const$mobil_peaceful,
    results_const$ratio_fv,
    results_const$ratio_vp,
    results_const$beta,
    log(2) / results_const$beta
  ),
  Covariate = c(
    ll_to_use,
    24,
    -2 * ll_to_use + 2 * 24,
    -2 * ll_to_use + 24 * log(nrow(protests)),
    mobil_fatal,
    mobil_violent,
    mobil_peaceful,
    mobil_fatal / mobil_violent,
    mobil_violent / mobil_peaceful,
    beta_decay,
    log(2) / beta_decay
  )
)

print(comparison, row.names = FALSE)

# Likelihood ratio test
lr_stat <- 2 * (ll_to_use - results_const$loglik)
df_diff <- 24 - results_const$n_params
p_value <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)

cat(sprintf("\nLikelihood ratio test (covariate vs constant):\n"))
cat(sprintf("  LR statistic: %.2f\n", lr_stat))
cat(sprintf("  Degrees of freedom: %d\n", df_diff))
cat(sprintf("  p-value: %.4e\n", p_value))

if (lr_stat > 0 && p_value < 0.05) {
  cat("  *** Covariate model significantly better ***\n")
} else if (lr_stat < 0) {
  cat("  ⚠️ Constant model has HIGHER log-likelihood!\n")
  cat("  This suggests issues with the covariate model.\n")
}

# =============================================================================
# 9. INTERPRETATION SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  DIAGNOSTIC SUMMARY                                          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("KEY FINDINGS:\n")
cat("─────────────────────────────────────────────────────────────────\n")

if (use_converged) {
  cat("1. CONVERGENCE ISSUE DETECTED\n")
  cat("   - The saved model used non-converged parameters\n")
  cat("   - Diagnostics above use the best CONVERGED estimate\n\n")
}

cat(sprintf("2. POPULATION EFFECT (γ = %.4f)\n", gamma))
if (gamma > 0) {
  cat("   - Positive: larger populations have more protests\n")
  cat(sprintf("   - 1 unit increase in log(pop) → %.1f%% increase in rate\n\n", 100*(exp(gamma)-1)))
} else {
  cat("   - Negative or zero: unexpected for protest data\n\n")
}

cat(sprintf("3. POVERTY EFFECT (δ = %.4f)\n", delta))
if (delta < 0) {
  cat("   - NEGATIVE: Higher poverty → FEWER protests?\n")
  cat("   - This may reflect:\n")
  cat("     a) Underreporting in poor areas\n")
  cat("     b) Lower capacity for mobilization\n")
  cat("     c) Urban bias in data (poor rural areas have few events)\n\n")
} else {
  cat("   - Positive: Higher poverty → more protests\n\n")
}

cat(sprintf("4. DECAY RATE (β = %.2f, half-life = %.3f days)\n", beta_decay, log(2)/beta_decay))
if (beta_decay > 5) {
  cat("   - VERY FAST decay - triggering effects dissipate within hours\n")
  cat("   - May indicate weak temporal clustering\n")
  cat("   - Or: model compensating for misspecification\n\n")
}

cat(sprintf("5. MOBILIZATION RATIOS\n"))
cat(sprintf("   - Fatal/Violent: %.4f\n", mobil_fatal / mobil_violent))
cat(sprintf("   - Violent/Peaceful: %.4f\n", mobil_violent / mobil_peaceful))
if (mobil_violent / mobil_peaceful < 1) {
  cat("   - Violent protests have LOWER mobilization than peaceful\n")
  cat("   - This contradicts 'radical flank' hypothesis\n")
  cat("   - May reflect: selection effects, measurement issues, or true effect\n\n")
}

cat("6. MODEL FIT\n")
if (ll_to_use < results_const$loglik) {
  cat("   - ⚠️ Covariate model has LOWER likelihood than constant!\n")
  cat("   - This is unusual and suggests model issues\n")
  cat("   - Consider: simplifying the model\n")
} else {
  cat(sprintf("   - Covariate model improves fit (ΔLL = %.1f)\n", ll_to_use - results_const$loglik))
}

# =============================================================================
# 10. SAVE DIAGNOSTIC RESULTS
# =============================================================================

cat("\n=== SAVING DIAGNOSTIC RESULTS ===\n")

diagnostics <- list(
  # Convergence info
  use_converged = use_converged,
  converged_params = params_to_use,
  converged_ll = ll_to_use,
  best_converged_name = if(exists("best_converged")) best_converged$starting_point else NA,

  # Parameters
  beta_0_F = beta_0_F, beta_0_V = beta_0_V, beta_0_P = beta_0_P,
  gamma = gamma, delta = delta,
  beta_years = beta_years,
  beta_decay = beta_decay,
  half_life = log(2) / beta_decay,

  # Offspring
  mobil_fatal = mobil_fatal,
  mobil_violent = mobil_violent,
  mobil_peaceful = mobil_peaceful,
  ratio_fv = mobil_fatal / mobil_violent,
  ratio_vp = mobil_violent / mobil_peaceful,

  # Hessian analysis
  hessian = if(exists("hessian_mat")) hessian_mat else NULL,
  eigenvalues = if(exists("eigenvalues")) eigenvalues else NULL,
  se = se,
  vcov = if(exists("vcov_mat")) vcov_mat else NULL,

  # Model comparison
  constant_ll = results_const$loglik,
  covariate_ll = ll_to_use,
  lr_stat = lr_stat,
  lr_pvalue = p_value,

  # Data summaries
  poverty_table = poverty_table,
  rate_by_poverty = rate_by_poverty,
  pop_table = pop_table,
  year_effects = year_df,

  # Correlation
  cor_logpop_poverty = cor(log_pop, poverty)
)

saveRDS(diagnostics, "trivariate_diagnostics.rds")
cat("✓ Saved: trivariate_diagnostics.rds\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  DIAGNOSTICS COMPLETE                                        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
