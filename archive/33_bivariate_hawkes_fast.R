################################################################################
#           BIVARIATE HAWKES MODEL - Efficient Implementation
#           Comparing mobilization potential: Violent vs Peaceful
################################################################################

library(tidyverse)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   BIVARIATE HAWKES: MOBILIZATION COMPARISON                  ║\n")
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

# Use a subset for tractable computation (2019-2021: peak activity period)
cat("Using 2019-2021 subset for tractable computation...\n")
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

# Use recursive computation for exponential kernel
# R(t) = sum_{t_j < t} exp(-beta * (t - t_j))
# R(t_i) = exp(-beta * (t_i - t_{i-1})) * (1 + R(t_{i-1}))

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
  # R_v[i] = sum of exp(-beta * (t_i - t_j)) for violent j < i
  # R_p[i] = sum of exp(-beta * (t_i - t_j)) for peaceful j < i

  R_v <- numeric(n)
  R_p <- numeric(n)

  for (i in 2:n) {
    dt <- all_times[i] - all_times[i-1]
    decay <- exp(-beta * dt)

    # Update recursive sums
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
      # Violent event
      lambda_i <- mu_v + alpha_vv * R_v[i] + alpha_pv * R_p[i]
    } else {
      # Peaceful event
      lambda_i <- mu_p + alpha_vp * R_v[i] + alpha_pp * R_p[i]
    }

    if (lambda_i <= 1e-10) return(-1e10)
    loglik <- loglik + log(lambda_i)
  }

  # Term 2: integral of intensity
  # For exponential kernel: integral from t_j to T of alpha * exp(-beta*(t-t_j)) dt
  # = alpha/beta * (1 - exp(-beta*(T - t_j)))

  # Integral of violent intensity
  integral_v <- mu_v * T_max
  for (j in which(all_types == 1)) {
    integral_v <- integral_v + (alpha_vv / beta) * (1 - exp(-beta * (T_max - all_times[j])))
  }
  for (j in which(all_types == 0)) {
    integral_v <- integral_v + (alpha_pv / beta) * (1 - exp(-beta * (T_max - all_times[j])))
  }

  # Integral of peaceful intensity
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
  # Perturb starting values
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
# EXTRACT AND INTERPRET RESULTS
# =============================================================================

params <- best_result$par
mu_v <- exp(params[1])
mu_p <- exp(params[2])
alpha_vv <- exp(params[3])
alpha_vp <- exp(params[4])
alpha_pv <- exp(params[5])
alpha_pp <- exp(params[6])
beta <- exp(params[7])

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   BIVARIATE HAWKES RESULTS                                   ║\n")
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
# MOBILIZATION COMPARISON
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   MOBILIZATION POTENTIAL COMPARISON                          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Expected offspring per parent
offspring_vv <- alpha_vv / beta
offspring_vp <- alpha_vp / beta
offspring_pv <- alpha_pv / beta
offspring_pp <- alpha_pp / beta

cat("EXPECTED OFFSPRING PER PARENT EVENT:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("                        Violent     Peaceful    TOTAL\n")
cat("                        offspring   offspring   MOBILIZATION\n")
cat(sprintf("  From VIOLENT event:   %.4f      %.4f      %.4f\n",
            offspring_vv, offspring_vp, offspring_vv + offspring_vp))
cat(sprintf("  From PEACEFUL event:  %.4f      %.4f      %.4f\n",
            offspring_pv, offspring_pp, offspring_pv + offspring_pp))

# Total mobilization
mobil_violent <- offspring_vv + offspring_vp
mobil_peaceful <- offspring_pv + offspring_pp

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   KEY RESULT: MOBILIZATION RATIO                             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("Total mobilization from VIOLENT event:  %.4f offspring\n", mobil_violent))
cat(sprintf("Total mobilization from PEACEFUL event: %.4f offspring\n", mobil_peaceful))
cat(sprintf("\nMOBILIZATION RATIO (Violent / Peaceful): %.2fx\n",
            mobil_violent / mobil_peaceful))

if (mobil_violent > mobil_peaceful) {
  cat("\n→ Violent events trigger MORE total offspring\n")
} else if (mobil_violent < mobil_peaceful) {
  cat("\n→ Violent events trigger FEWER total offspring\n")
} else {
  cat("\n→ Similar mobilization potential\n")
}

# Decomposition
cat("\n")
cat("DECOMPOSITION:\n")
cat(sprintf("  Violence → Violence (escalation):  %.4f (%.0f%% of violent's offspring)\n",
            offspring_vv, 100 * offspring_vv / mobil_violent))
cat(sprintf("  Violence → Peaceful (spillover):   %.4f (%.0f%% of violent's offspring)\n",
            offspring_vp, 100 * offspring_vp / mobil_violent))

# Comparison with univariate model
cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   COMPARISON WITH UNIVARIATE MODEL                           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("UNIVARIATE MODEL (treats all events as one stream):\n")
cat("  Violence coefficient: 5.5×\n")
cat("  This CONFLATES escalation with mobilization\n\n")

cat("BIVARIATE MODEL (separates streams):\n")
cat(sprintf("  Violence → Violence: %.4f (escalation)\n", offspring_vv))
cat(sprintf("  Violence → Peaceful: %.4f (spillover/suppression)\n", offspring_vp))
cat(sprintf("  TOTAL mobilization ratio: %.2fx\n", mobil_violent / mobil_peaceful))
cat("\n  This SEPARATES escalation from overall mobilization!\n")

# Save results
results <- list(
  params = params,
  mu_v = mu_v, mu_p = mu_p,
  alpha_vv = alpha_vv, alpha_vp = alpha_vp,
  alpha_pv = alpha_pv, alpha_pp = alpha_pp,
  beta = beta,
  offspring_vv = offspring_vv, offspring_vp = offspring_vp,
  offspring_pv = offspring_pv, offspring_pp = offspring_pp,
  mobil_violent = mobil_violent,
  mobil_peaceful = mobil_peaceful,
  mobil_ratio = mobil_violent / mobil_peaceful,
  loglik = -best_result$value,
  n_violent = n_v, n_peaceful = n_p, T_max = T_max
)

saveRDS(results, "bivariate_hawkes_results.rds")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   BIVARIATE MODEL ESTIMATION COMPLETE                        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
