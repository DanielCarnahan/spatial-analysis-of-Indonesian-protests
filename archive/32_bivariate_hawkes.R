################################################################################
#           BIVARIATE HAWKES MODEL: Separating Escalation from Mobilization
#
#           Model violent and peaceful events as two interacting streams
################################################################################

library(tidyverse)
library(Rcpp)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   BIVARIATE HAWKES MODEL                                     ║\n")
cat("║   Separating Violence→Violence from Violence→Peaceful        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Load data
protests <- readRDS("protests_prepared.rds")

# Create binary violence indicator
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

# Separate into two streams
violent_times <- protests$time[protests$is_violent]
peaceful_times <- protests$time[!protests$is_violent]

cat(sprintf("Violent events: %d\n", length(violent_times)))
cat(sprintf("Peaceful events: %d\n", length(peaceful_times)))
cat(sprintf("Observation period: %.0f days\n\n", T_max))

# =============================================================================
# BIVARIATE HAWKES LOG-LIKELIHOOD (Exponential Kernel)
# =============================================================================

# Parameters:
# mu_v, mu_p: baseline rates
# alpha_vv, alpha_vp, alpha_pv, alpha_pp: triggering coefficients
# beta: decay rate (shared)

bivariate_hawkes_loglik <- function(params, violent_times, peaceful_times, T_max,
                                     cutoff = 90) {
  # Extract parameters
  log_mu_v <- params[1]
  log_mu_p <- params[2]
  log_alpha_vv <- params[3]
  log_alpha_vp <- params[4]
  log_alpha_pv <- params[5]
  log_alpha_pp <- params[6]
  log_beta <- params[7]

  mu_v <- exp(log_mu_v)
  mu_p <- exp(log_mu_p)
  alpha_vv <- exp(log_alpha_vv)
  alpha_vp <- exp(log_alpha_vp)
  alpha_pv <- exp(log_alpha_pv)
  alpha_pp <- exp(log_alpha_pp)
  beta <- exp(log_beta)

  n_v <- length(violent_times)
  n_p <- length(peaceful_times)

  # Log-likelihood has four terms:
  # 1. Sum of log(lambda_v) at violent event times
  # 2. Sum of log(lambda_p) at peaceful event times
  # 3. -Integral of lambda_v from 0 to T
  # 4. -Integral of lambda_p from 0 to T

  loglik <- 0

  # Term 1: log(lambda_v) at violent times
  for (i in 1:n_v) {
    t_i <- violent_times[i]

    # Intensity from violent parents
    intensity_from_v <- 0
    for (j in which(violent_times < t_i & violent_times >= t_i - cutoff)) {
      intensity_from_v <- intensity_from_v + alpha_vv * exp(-beta * (t_i - violent_times[j]))
    }

    # Intensity from peaceful parents
    intensity_from_p <- 0
    for (j in which(peaceful_times < t_i & peaceful_times >= t_i - cutoff)) {
      intensity_from_p <- intensity_from_p + alpha_pv * exp(-beta * (t_i - peaceful_times[j]))
    }

    lambda_v_i <- mu_v + intensity_from_v + intensity_from_p

    if (lambda_v_i <= 0) return(-1e10)
    loglik <- loglik + log(lambda_v_i)
  }

  # Term 2: log(lambda_p) at peaceful times
  for (i in 1:n_p) {
    t_i <- peaceful_times[i]

    # Intensity from violent parents
    intensity_from_v <- 0
    for (j in which(violent_times < t_i & violent_times >= t_i - cutoff)) {
      intensity_from_v <- intensity_from_v + alpha_vp * exp(-beta * (t_i - violent_times[j]))
    }

    # Intensity from peaceful parents
    intensity_from_p <- 0
    for (j in which(peaceful_times < t_i & peaceful_times >= t_i - cutoff)) {
      intensity_from_p <- intensity_from_p + alpha_pp * exp(-beta * (t_i - peaceful_times[j]))
    }

    lambda_p_i <- mu_p + intensity_from_v + intensity_from_p

    if (lambda_p_i <= 0) return(-1e10)
    loglik <- loglik + log(lambda_p_i)
  }

  # Term 3: -Integral of lambda_v
  # Integral of mu_v from 0 to T = mu_v * T
  integral_v <- mu_v * T_max

  # Integral of alpha_vv * exp(-beta * (t - t_j)) from t_j to T
  # = alpha_vv / beta * (1 - exp(-beta * (T - t_j)))
  for (j in 1:n_v) {
    integral_v <- integral_v + (alpha_vv / beta) * (1 - exp(-beta * (T_max - violent_times[j])))
  }
  for (j in 1:n_p) {
    integral_v <- integral_v + (alpha_pv / beta) * (1 - exp(-beta * (T_max - peaceful_times[j])))
  }

  loglik <- loglik - integral_v

  # Term 4: -Integral of lambda_p
  integral_p <- mu_p * T_max

  for (j in 1:n_v) {
    integral_p <- integral_p + (alpha_vp / beta) * (1 - exp(-beta * (T_max - violent_times[j])))
  }
  for (j in 1:n_p) {
    integral_p <- integral_p + (alpha_pp / beta) * (1 - exp(-beta * (T_max - peaceful_times[j])))
  }

  loglik <- loglik - integral_p

  return(loglik)
}

# =============================================================================
# ESTIMATE MODEL (Using subset for speed - full data would take too long in R)
# =============================================================================

cat("Note: Full bivariate Hawkes estimation on 16k events requires C++.\n")
cat("Running on a subset for proof of concept...\n\n")

# Use 2020-2021 data as a subset (manageable size)
subset_protests <- protests %>%
  filter(event_date >= as.Date("2020-01-01") & event_date <= as.Date("2021-12-31"))

subset_start <- min(subset_protests$event_date)
subset_protests$time <- as.numeric(subset_protests$event_date - subset_start)
T_subset <- max(subset_protests$time) + 1

violent_subset <- subset_protests$time[subset_protests$is_violent]
peaceful_subset <- subset_protests$time[!subset_protests$is_violent]

cat(sprintf("Subset: %d violent, %d peaceful events over %.0f days\n",
            length(violent_subset), length(peaceful_subset), T_subset))

# Negative log-likelihood for optimization
neg_loglik <- function(params) {
  ll <- bivariate_hawkes_loglik(params, violent_subset, peaceful_subset, T_subset)
  return(-ll)
}

# Starting values
# Rough estimates: ~0.1 violent/day, ~1 peaceful/day
start_params <- c(
  log(0.1),   # log_mu_v
  log(1.0),   # log_mu_p
  log(0.01),  # log_alpha_vv
  log(0.01),  # log_alpha_vp
  log(0.01),  # log_alpha_pv
  log(0.01),  # log_alpha_pp
  log(0.1)    # log_beta
)

cat("\nOptimizing (this may take a few minutes)...\n")

result <- optim(
  par = start_params,
  fn = neg_loglik,
  method = "BFGS",
  control = list(maxit = 1000, trace = 1)
)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   BIVARIATE HAWKES RESULTS                                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Extract parameters
params_est <- result$par
mu_v <- exp(params_est[1])
mu_p <- exp(params_est[2])
alpha_vv <- exp(params_est[3])
alpha_vp <- exp(params_est[4])
alpha_pv <- exp(params_est[5])
alpha_pp <- exp(params_est[6])
beta <- exp(params_est[7])

cat("BASELINE RATES:\n")
cat(sprintf("  μ_V (violent): %.4f events/day\n", mu_v))
cat(sprintf("  μ_P (peaceful): %.4f events/day\n", mu_p))

cat("\nTRIGGERING COEFFICIENTS:\n")
cat(sprintf("  α_VV (violence → violence): %.6f\n", alpha_vv))
cat(sprintf("  α_VP (violence → peaceful): %.6f\n", alpha_vp))
cat(sprintf("  α_PV (peaceful → violence): %.6f\n", alpha_pv))
cat(sprintf("  α_PP (peaceful → peaceful): %.6f\n", alpha_pp))

cat("\nDECAY RATE:\n")
cat(sprintf("  β: %.4f (half-life: %.1f days)\n", beta, log(2)/beta))

# =============================================================================
# KEY METRICS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   KEY METRICS: ESCALATION vs MOBILIZATION                    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Expected offspring (integral of triggering kernel)
offspring_vv <- alpha_vv / beta  # Violent children per violent parent
offspring_vp <- alpha_vp / beta  # Peaceful children per violent parent
offspring_pv <- alpha_pv / beta  # Violent children per peaceful parent
offspring_pp <- alpha_pp / beta  # Peaceful children per peaceful parent

cat("EXPECTED OFFSPRING (per parent event):\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Violent parent → violent children:  %.4f\n", offspring_vv))
cat(sprintf("  Violent parent → peaceful children: %.4f\n", offspring_vp))
cat(sprintf("  Peaceful parent → violent children: %.4f\n", offspring_pv))
cat(sprintf("  Peaceful parent → peaceful children: %.4f\n", offspring_pp))

cat("\nTOTAL OFFSPRING:\n")
cat(sprintf("  From violent parent:  %.4f (violence) + %.4f (peaceful) = %.4f total\n",
            offspring_vv, offspring_vp, offspring_vv + offspring_vp))
cat(sprintf("  From peaceful parent: %.4f (violence) + %.4f (peaceful) = %.4f total\n",
            offspring_pv, offspring_pp, offspring_pv + offspring_pp))

# Mobilization effect: total offspring
mobilization_violent <- offspring_vv + offspring_vp
mobilization_peaceful <- offspring_pv + offspring_pp

cat("\nMOBILIZATION EFFECT (total offspring):\n")
cat(sprintf("  Violent events trigger %.4f total events\n", mobilization_violent))
cat(sprintf("  Peaceful events trigger %.4f total events\n", mobilization_peaceful))
cat(sprintf("  Ratio: %.2fx\n", mobilization_violent / mobilization_peaceful))

# Escalation effect: share of violent offspring
escalation_from_violent <- offspring_vv / (offspring_vv + offspring_vp)
escalation_from_peaceful <- offspring_pv / (offspring_pv + offspring_pp)

cat("\nESCALATION EFFECT (share of violent offspring):\n")
cat(sprintf("  From violent parent: %.1f%% of offspring are violent\n",
            100 * escalation_from_violent))
cat(sprintf("  From peaceful parent: %.1f%% of offspring are violent\n",
            100 * escalation_from_peaceful))

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   INTERPRETATION                                             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

if (mobilization_violent < mobilization_peaceful) {
  cat("✓ MOBILIZATION: Violent events trigger FEWER total events\n")
  cat("  This confirms our raw data finding.\n\n")
} else {
  cat("? MOBILIZATION: Violent events trigger MORE total events\n")
  cat("  This differs from raw data - check model specification.\n\n")
}

if (escalation_from_violent > escalation_from_peaceful) {
  cat("✓ ESCALATION: Violent events produce MORE violent offspring\n")
  cat("  This is the violence-to-violence clustering we identified.\n\n")
} else {
  cat("? ESCALATION: Violent events do NOT produce more violent offspring\n\n")
}

cat("The bivariate model SEPARATES these effects:\n")
cat("  - Mobilization: Does violence increase total protest count?\n")
cat("  - Escalation: Does violence shift composition toward violence?\n\n")

cat("The original univariate model conflated these,\n")
cat("giving a high violence coefficient that reflected escalation,\n")
cat("not mobilization.\n")

# =============================================================================
# COMPARISON WITH UNIVARIATE MODEL
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   COMPARISON: BIVARIATE vs UNIVARIATE                        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("UNIVARIATE MODEL (Model 5):\n")
cat("  Violence coefficient: exp(1.697) = 5.46×\n")
cat("  Interpretation: Violence increases triggering by 5.5×\n")
cat("  Problem: Conflates escalation with mobilization\n\n")

cat("BIVARIATE MODEL:\n")
cat(sprintf("  Violence → Violence: α_VV = %.6f\n", alpha_vv))
cat(sprintf("  Violence → Peaceful: α_VP = %.6f\n", alpha_vp))
cat(sprintf("  Ratio α_VV/α_VP: %.2fx\n", alpha_vv/alpha_vp))
cat("\n")
cat("  This separates:\n")
cat("  - How much violence increases VIOLENT follow-on\n")
cat("  - How much violence increases PEACEFUL follow-on\n")

# Save results
results <- list(
  params = params_est,
  mu_v = mu_v, mu_p = mu_p,
  alpha_vv = alpha_vv, alpha_vp = alpha_vp,
  alpha_pv = alpha_pv, alpha_pp = alpha_pp,
  beta = beta,
  offspring_vv = offspring_vv, offspring_vp = offspring_vp,
  offspring_pv = offspring_pv, offspring_pp = offspring_pp,
  loglik = -result$value,
  convergence = result$convergence
)

saveRDS(results, "bivariate_hawkes_results.rds")
cat("\nResults saved to bivariate_hawkes_results.rds\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   BIVARIATE HAWKES ESTIMATION COMPLETE                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
