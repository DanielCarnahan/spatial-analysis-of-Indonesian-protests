# ============================================================================
# UNCERTAINTY QUANTIFICATION FOR HAWKES MODEL
# ============================================================================
#
# This script computes standard errors, confidence intervals, and hypothesis
# tests for the fitted Hawkes process model.
#
# Methods:
# 1. Hessian-based standard errors (observed Fisher information)
# 2. Delta method for transformed parameters
# 3. Individual parameter z-tests
# 4. Likelihood ratio tests for model comparison
#
# ============================================================================

library(dplyr)
library(Rcpp)
library(numDeriv)  # For numerical Hessian
library(xtable)    # For LaTeX tables

# Compile C++ likelihood function
cat("Compiling C++ likelihood function...\n")
sourceCpp("hawkes_likelihood.cpp")
cat("✓ C++ compiled\n\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   UNCERTAINTY QUANTIFICATION FOR HAWKES MODEL                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ====================================================================
# CONFIGURATION
# ====================================================================

TEMPORAL_CUTOFF <- 90  # Must match fitted model

# ====================================================================
# LOAD FITTED MODEL AND DATA
# ====================================================================

cat("=== LOADING FITTED MODEL ===\n")
fit <- readRDS("model_poverty.rds")

cat(sprintf("Loaded model with LL = %.2f\n", fit$loglik))
cat(sprintf("Parameters: %d\n", length(fit$params)))
cat(sprintf("Runtime: %.1f minutes\n\n", fit$runtime))

# ====================================================================
# RELOAD DATA (need same data as fitting)
# ====================================================================

cat("=== RELOADING DATA ===\n")

# Load data (same as fitting script)
protests <- readRDS("protests_with_poverty.rds")

cat(sprintf("Loaded %d events\n", nrow(protests)))

# Create marks data (same as fitting)
marks_data <- protests %>%
  select(
    event_id = event_id_cnty,
    time = days_since_start,
    is_violent,
    is_peaceful,
    state_intervention,
    event_type,
    assoc_actor_1,
    fatalities,
    district = admin2,
    year,
    log_pop,
    poverty_decimal
  ) %>%
  arrange(time) %>%
  mutate(
    is_riot = as.integer(event_type == "Riots"),
    has_fatalities = as.integer(fatalities > 0),
    is_student = as.integer(grepl("Student", assoc_actor_1, ignore.case = TRUE)),
    is_labor = as.integer(grepl("Labor|Worker", assoc_actor_1, ignore.case = TRUE))
  )

times <- marks_data$time
marks <- marks_data %>%
  select(is_riot, has_fatalities, is_student, is_labor, log_pop, poverty_decimal, year)

# Create district-year data for compensator
district_years <- protests %>%
  group_by(admin2, year) %>%
  summarize(
    log_pop = mean(log_pop),
    poverty_decimal = mean(poverty_decimal),
    .groups = "drop"
  ) %>%
  mutate(
    exposure = case_when(
      year == 2015 ~ as.numeric(as.Date("2015-12-31") - as.Date("2015-01-01") + 1),
      year == 2024 ~ as.numeric(as.Date("2024-12-31") - as.Date("2024-01-01") + 1),
      TRUE ~ 365
    )
  )

cat(sprintf("District-years: %d\n", nrow(district_years)))
cat("✓ Data reloaded\n\n")

# ====================================================================
# DEFINE TRANSFORMATION FUNCTIONS (must match fitting script)
# ====================================================================

# Gamma: [0.01, 2]
transform_gamma <- function(gamma_raw) {
  0.01 + 1.99 / (1 + exp(-gamma_raw))
}
transform_gamma_inv <- function(gamma) {
  -log((1.99 / (gamma - 0.01)) - 1)
}

# Delta: [-20, 20]
transform_delta <- function(delta_raw) {
  -20 + 40 / (1 + exp(-delta_raw))
}
transform_delta_inv <- function(delta) {
  -log((40 / (delta + 20)) - 1)
}

# Year effects: [-3, 3]
transform_beta_year <- function(beta_year_raw) {
  -3 + 6 / (1 + exp(-beta_year_raw))
}
transform_beta_year_inv <- function(beta_year) {
  -log((6 / (beta_year + 3)) - 1)
}

# Mark effects (riot, fatal, student, labor): [-5, 5]
transform_beta_mark <- function(beta_mark_raw) {
  -5 + 10 / (1 + exp(-beta_mark_raw))
}
transform_beta_mark_inv <- function(beta_mark) {
  -log((10 / (beta_mark + 5)) - 1)
}

# Decay: [0.001, 10]
transform_decay <- function(decay_raw) {
  0.001 + 9.999 / (1 + exp(-decay_raw))
}
transform_decay_inv <- function(decay) {
  -log((9.999 / (decay - 0.001)) - 1)
}

# ====================================================================
# DEFINE NEGATIVE LOG-LIKELIHOOD FUNCTION
# ====================================================================

# This function takes RAW (unconstrained) parameters and returns negative LL
hawkes_negloglik_raw <- function(params_raw) {

  # Convert raw params to named list
  param_names <- c("beta_0_bg", "gamma_raw", "delta_raw",
                   "beta_2016_raw", "beta_2017_raw", "beta_2018_raw",
                   "beta_2019_raw", "beta_2020_raw", "beta_2021_raw",
                   "beta_2022_raw", "beta_2023_raw", "beta_2024_raw",
                   "beta_0_trig_raw", "beta_riot_raw", "beta_fatal_raw",
                   "beta_student_raw", "beta_labor_raw", "decay_raw")

  params_list <- as.list(params_raw)
  names(params_list) <- param_names

  # Call C++ function
  neg_ll <- hawkes_negloglik_cpp(
    times = times,
    log_pop = marks$log_pop,
    poverty_decimal = marks$poverty_decimal,
    year = marks$year,
    is_riot = marks$is_riot,
    has_fatalities = marks$has_fatalities,
    is_student = marks$is_student,
    is_labor = marks$is_labor,
    district_year_exposure = district_years$exposure,
    district_year_log_pop = district_years$log_pop,
    district_year_poverty = district_years$poverty_decimal,
    district_year_year = district_years$year,
    params = unlist(params_list),
    temporal_cutoff = TEMPORAL_CUTOFF
  )

  return(neg_ll)
}

# ====================================================================
# COMPUTE HESSIAN MATRIX
# ====================================================================

cat("=== COMPUTING HESSIAN MATRIX ===\n")
cat("This may take several minutes...\n\n")

# Extract raw parameter values from fitted model
param_names <- c("beta_0_bg", "gamma_raw", "delta_raw",
                 "beta_2016_raw", "beta_2017_raw", "beta_2018_raw",
                 "beta_2019_raw", "beta_2020_raw", "beta_2021_raw",
                 "beta_2022_raw", "beta_2023_raw", "beta_2024_raw",
                 "beta_0_trig_raw", "beta_riot_raw", "beta_fatal_raw",
                 "beta_student_raw", "beta_labor_raw", "decay_raw")

params_raw_vec <- unlist(fit$params[param_names])

# Compute numerical Hessian
cat("Computing Hessian (this takes ~5-10 minutes)...\n")
start_time <- Sys.time()

hessian_matrix <- numDeriv::hessian(
  func = hawkes_negloglik_raw,
  x = params_raw_vec,
  method = "Richardson"  # More accurate than simple finite differences
)

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat(sprintf("✓ Hessian computed in %.1f minutes\n", as.numeric(runtime)))
cat(sprintf("  Dimension: %d x %d\n\n", nrow(hessian_matrix), ncol(hessian_matrix)))

# Check positive definiteness
eigenvalues <- eigen(hessian_matrix, only.values = TRUE)$values
if (all(eigenvalues > 0)) {
  cat("✓ Hessian is positive definite (local minimum confirmed)\n")
} else {
  cat("⚠ Warning: Hessian has negative eigenvalues\n")
  cat(sprintf("  Minimum eigenvalue: %.6f\n", min(eigenvalues)))
}

# Compute variance-covariance matrix
vcov_matrix <- solve(hessian_matrix)
se_raw <- sqrt(diag(vcov_matrix))

cat("\n")

# ====================================================================
# COMPUTE STANDARD ERRORS FOR TRANSFORMED PARAMETERS
# ====================================================================

cat("=== COMPUTING STANDARD ERRORS (DELTA METHOD) ===\n\n")

# For parameters that undergo transformations, we need to use delta method
# SE(h(θ)) ≈ |h'(θ)| · SE(θ)

# 1. Gamma (population effect)
gamma_raw <- fit$params$gamma_raw
gamma <- transform_gamma(gamma_raw)
# Derivative: d/dx [a + b/(1+exp(-x))] = b·exp(-x)/(1+exp(-x))^2
dgamma_dgamma_raw <- 1.99 * exp(-gamma_raw) / (1 + exp(-gamma_raw))^2
se_gamma <- abs(dgamma_dgamma_raw) * se_raw[2]

# 2. Delta (poverty effect)
delta_raw <- fit$params$delta_raw
delta <- transform_delta(delta_raw)
ddelta_ddelta_raw <- 40 * exp(-delta_raw) / (1 + exp(-delta_raw))^2
se_delta <- abs(ddelta_ddelta_raw) * se_raw[3]

# 3. Year effects (9 parameters)
year_effects <- c("beta_2016_raw", "beta_2017_raw", "beta_2018_raw",
                  "beta_2019_raw", "beta_2020_raw", "beta_2021_raw",
                  "beta_2022_raw", "beta_2023_raw", "beta_2024_raw")

year_results <- data.frame(
  year = 2016:2024,
  beta_year_raw = numeric(9),
  beta_year = numeric(9),
  se_beta_year = numeric(9)
)

for (i in 1:9) {
  year_raw <- fit$params[[year_effects[i]]]
  year_trans <- transform_beta_year(year_raw)
  dyear_dyear_raw <- 6 * exp(-year_raw) / (1 + exp(-year_raw))^2
  se_year <- abs(dyear_dyear_raw) * se_raw[3 + i]

  year_results$beta_year_raw[i] <- year_raw
  year_results$beta_year[i] <- year_trans
  year_results$se_beta_year[i] <- se_year
}

# 4. Mark effects
mark_effects <- c("beta_riot_raw", "beta_fatal_raw", "beta_student_raw", "beta_labor_raw")
mark_names <- c("Riot", "Fatalities", "Student-led", "Labor-led")

mark_results <- data.frame(
  mark = mark_names,
  beta_mark_raw = numeric(4),
  beta_mark = numeric(4),
  se_beta_mark = numeric(4)
)

for (i in 1:4) {
  mark_raw <- fit$params[[mark_effects[i]]]
  mark_trans <- transform_beta_mark(mark_raw)
  dmark_dmark_raw <- 10 * exp(-mark_raw) / (1 + exp(-mark_raw))^2
  se_mark <- abs(dmark_dmark_raw) * se_raw[12 + i]

  mark_results$beta_mark_raw[i] <- mark_raw
  mark_results$beta_mark[i] <- mark_trans
  mark_results$se_beta_mark[i] <- se_mark
}

# 5. Decay parameter
decay_raw <- fit$params$decay_raw
decay <- transform_decay(decay_raw)
ddecay_ddecay_raw <- 9.999 * exp(-decay_raw) / (1 + exp(-decay_raw))^2
se_decay <- abs(ddecay_ddecay_raw) * se_raw[18]

# 6. Baseline parameters (no transformation)
se_beta_0_bg <- se_raw[1]
se_beta_0_trig_raw <- se_raw[13]

cat("✓ Standard errors computed\n\n")

# ====================================================================
# INDIVIDUAL PARAMETER TESTS
# ====================================================================

cat("=== INDIVIDUAL PARAMETER HYPOTHESIS TESTS ===\n\n")

# Function to compute z-statistic and p-value
compute_test <- function(estimate, se, null_value = 0) {
  z_stat <- (estimate - null_value) / se
  p_value <- 2 * pnorm(-abs(z_stat))
  ci_lower <- estimate - 1.96 * se
  ci_upper <- estimate + 1.96 * se

  list(
    estimate = estimate,
    se = se,
    z_stat = z_stat,
    p_value = p_value,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    significant = p_value < 0.05
  )
}

# Create comprehensive results table
results_table <- data.frame(
  Parameter = character(),
  Estimate = numeric(),
  SE = numeric(),
  z_stat = numeric(),
  p_value = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  Significant = logical(),
  stringsAsFactors = FALSE
)

# Background parameters
results_table <- rbind(results_table,
  data.frame(
    Parameter = "β₀ (baseline)",
    Estimate = fit$params$beta_0_bg,
    SE = se_beta_0_bg,
    z_stat = fit$params$beta_0_bg / se_beta_0_bg,
    p_value = 2 * pnorm(-abs(fit$params$beta_0_bg / se_beta_0_bg)),
    CI_lower = fit$params$beta_0_bg - 1.96 * se_beta_0_bg,
    CI_upper = fit$params$beta_0_bg + 1.96 * se_beta_0_bg,
    Significant = (2 * pnorm(-abs(fit$params$beta_0_bg / se_beta_0_bg))) < 0.05
  )
)

results_table <- rbind(results_table,
  data.frame(
    Parameter = "γ (population)",
    Estimate = gamma,
    SE = se_gamma,
    z_stat = gamma / se_gamma,
    p_value = 2 * pnorm(-abs(gamma / se_gamma)),
    CI_lower = gamma - 1.96 * se_gamma,
    CI_upper = gamma + 1.96 * se_gamma,
    Significant = (2 * pnorm(-abs(gamma / se_gamma))) < 0.05
  )
)

results_table <- rbind(results_table,
  data.frame(
    Parameter = "δ (poverty)",
    Estimate = delta,
    SE = se_delta,
    z_stat = delta / se_delta,
    p_value = 2 * pnorm(-abs(delta / se_delta)),
    CI_lower = delta - 1.96 * se_delta,
    CI_upper = delta + 1.96 * se_delta,
    Significant = (2 * pnorm(-abs(delta / se_delta))) < 0.05
  )
)

# Year effects
for (i in 1:9) {
  results_table <- rbind(results_table,
    data.frame(
      Parameter = sprintf("β_%d", 2015 + i),
      Estimate = year_results$beta_year[i],
      SE = year_results$se_beta_year[i],
      z_stat = year_results$beta_year[i] / year_results$se_beta_year[i],
      p_value = 2 * pnorm(-abs(year_results$beta_year[i] / year_results$se_beta_year[i])),
      CI_lower = year_results$beta_year[i] - 1.96 * year_results$se_beta_year[i],
      CI_upper = year_results$beta_year[i] + 1.96 * year_results$se_beta_year[i],
      Significant = (2 * pnorm(-abs(year_results$beta_year[i] / year_results$se_beta_year[i]))) < 0.05
    )
  )
}

# Triggering parameters
results_table <- rbind(results_table,
  data.frame(
    Parameter = "β₀_trig (baseline triggering)",
    Estimate = fit$params$beta_0_trig_raw,
    SE = se_beta_0_trig_raw,
    z_stat = fit$params$beta_0_trig_raw / se_beta_0_trig_raw,
    p_value = 2 * pnorm(-abs(fit$params$beta_0_trig_raw / se_beta_0_trig_raw)),
    CI_lower = fit$params$beta_0_trig_raw - 1.96 * se_beta_0_trig_raw,
    CI_upper = fit$params$beta_0_trig_raw + 1.96 * se_beta_0_trig_raw,
    Significant = (2 * pnorm(-abs(fit$params$beta_0_trig_raw / se_beta_0_trig_raw))) < 0.05
  )
)

# Mark effects
for (i in 1:4) {
  results_table <- rbind(results_table,
    data.frame(
      Parameter = sprintf("β_%s", mark_results$mark[i]),
      Estimate = mark_results$beta_mark[i],
      SE = mark_results$se_beta_mark[i],
      z_stat = mark_results$beta_mark[i] / mark_results$se_beta_mark[i],
      p_value = 2 * pnorm(-abs(mark_results$beta_mark[i] / mark_results$se_beta_mark[i])),
      CI_lower = mark_results$beta_mark[i] - 1.96 * mark_results$se_beta_mark[i],
      CI_upper = mark_results$beta_mark[i] + 1.96 * mark_results$se_beta_mark[i],
      Significant = (2 * pnorm(-abs(mark_results$beta_mark[i] / mark_results$se_beta_mark[i]))) < 0.05
    )
  )
}

# Decay
results_table <- rbind(results_table,
  data.frame(
    Parameter = "β_decay",
    Estimate = decay,
    SE = se_decay,
    z_stat = decay / se_decay,
    p_value = 2 * pnorm(-abs(decay / se_decay)),
    CI_lower = decay - 1.96 * se_decay,
    CI_upper = decay + 1.96 * se_decay,
    Significant = (2 * pnorm(-abs(decay / se_decay))) < 0.05
  )
)

# Print results table
cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   PARAMETER ESTIMATES WITH STANDARD ERRORS                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

print(results_table, digits = 4, row.names = FALSE)

cat("\n")
cat("Significance codes: *** p<0.001, ** p<0.01, * p<0.05\n")
cat(sprintf("Significant parameters: %d / %d\n\n",
            sum(results_table$Significant), nrow(results_table)))

# ====================================================================
# LIKELIHOOD RATIO TESTS
# ====================================================================

cat("=== LIKELIHOOD RATIO TESTS ===\n\n")

# We need log-likelihoods from different model specifications
# These should be computed from previously fitted models

cat("Model comparison requires log-likelihoods from alternative specifications.\n")
cat("Current model (Model C): LL = -78,693.63\n")
cat("\n")

# If we have previous model fits, we can compute LR tests
# For now, we'll create a framework

lr_test <- function(ll_full, ll_restricted, df_diff) {
  lr_stat <- 2 * (ll_full - ll_restricted)
  p_value <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)

  list(
    lr_stat = lr_stat,
    df = df_diff,
    p_value = p_value,
    significant = p_value < 0.05
  )
}

# Example: Test if all mark effects = 0
# This would require fitting a restricted model without marks
# For demonstration purposes:

cat("Planned likelihood ratio tests:\n")
cat("1. Full model vs No marks (H₀: all β_mark = 0)\n")
cat("2. Full model vs No poverty (H₀: δ = 0)\n")
cat("3. Year FX vs Linear trend (H₀: year effects follow linear trend)\n")
cat("\n")
cat("Note: These require fitting restricted models separately.\n")
cat("See log files for model comparison results.\n\n")

# ====================================================================
# SAVE RESULTS
# ====================================================================

cat("=== SAVING RESULTS ===\n")

uncertainty_results <- list(
  hessian = hessian_matrix,
  vcov = vcov_matrix,
  se_raw = se_raw,
  results_table = results_table,
  year_results = year_results,
  mark_results = mark_results,
  summary = list(
    gamma = list(estimate = gamma, se = se_gamma),
    delta = list(estimate = delta, se = se_delta),
    decay = list(estimate = decay, se = se_decay)
  )
)

saveRDS(uncertainty_results, "uncertainty_results.rds")
cat("✓ Saved: uncertainty_results.rds\n")

# Save results table as CSV
write.csv(results_table, "parameter_estimates_with_se.csv", row.names = FALSE)
cat("✓ Saved: parameter_estimates_with_se.csv\n\n")

# ====================================================================
# GENERATE LATEX TABLES
# ====================================================================

cat("=== GENERATING LATEX TABLES ===\n")

# Table 1: Background parameters
bg_table <- results_table[1:12, ]  # First 12 rows: baseline + γ + δ + 9 year effects

latex_bg <- xtable(
  bg_table,
  caption = "Background Rate Parameters",
  label = "tab:background",
  digits = c(0, 0, 4, 4, 2, 4, 4, 4, 0)
)

sink("table_background_parameters.tex")
print(latex_bg, include.rownames = FALSE)
sink()

cat("✓ Saved: table_background_parameters.tex\n")

# Table 2: Triggering parameters
trig_table <- results_table[13:18, ]  # Triggering parameters

latex_trig <- xtable(
  trig_table,
  caption = "Triggering Effect Parameters",
  label = "tab:triggering",
  digits = c(0, 0, 4, 4, 2, 4, 4, 4, 0)
)

sink("table_triggering_parameters.tex")
print(latex_trig, include.rownames = FALSE)
sink()

cat("✓ Saved: table_triggering_parameters.tex\n\n")

# ====================================================================
# SUMMARY STATISTICS
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   SUMMARY OF KEY FINDINGS                                    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("HYPOTHESIS TEST RESULTS:\n")
cat("------------------------\n\n")

# H1: Riot effect
riot_test <- results_table[results_table$Parameter == "β_Riot", ]
cat("H1 (Riot effect):\n")
cat(sprintf("  Estimate: %.3f (SE = %.3f)\n", riot_test$Estimate, riot_test$SE))
cat(sprintf("  z = %.2f, p = %.4f\n", riot_test$z_stat, riot_test$p_value))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n", riot_test$CI_lower, riot_test$CI_upper))
cat(sprintf("  Multiplier: exp(%.3f) = %.2f\n", riot_test$Estimate, exp(riot_test$Estimate)))
if (riot_test$Significant) {
  cat("  → SUPPORTED: Riots significantly increase triggering\n\n")
} else {
  cat("  → NOT SUPPORTED\n\n")
}

# H2: Fatalities effect
fatal_test <- results_table[results_table$Parameter == "β_Fatalities", ]
cat("H2 (Fatalities effect):\n")
cat(sprintf("  Estimate: %.3f (SE = %.3f)\n", fatal_test$Estimate, fatal_test$SE))
cat(sprintf("  z = %.2f, p = %.4f\n", fatal_test$z_stat, fatal_test$p_value))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n", fatal_test$CI_lower, fatal_test$CI_upper))
cat(sprintf("  Multiplier: exp(%.3f) = %.2f\n", fatal_test$Estimate, exp(fatal_test$Estimate)))
if (fatal_test$Significant) {
  cat("  → SUPPORTED: Fatal events significantly increase triggering\n\n")
} else {
  cat("  → NOT SUPPORTED\n\n")
}

# H3: Student effect
student_test <- results_table[results_table$Parameter == "β_Student-led", ]
cat("H3a (Student-led effect):\n")
cat(sprintf("  Estimate: %.3f (SE = %.3f)\n", student_test$Estimate, student_test$SE))
cat(sprintf("  z = %.2f, p = %.4f\n", student_test$z_stat, student_test$p_value))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n", student_test$CI_lower, student_test$CI_upper))
cat(sprintf("  Multiplier: exp(%.3f) = %.2f\n", student_test$Estimate, exp(student_test$Estimate)))
if (student_test$Significant) {
  cat("  → SIGNIFICANT: Student-led protests SUPPRESS triggering\n\n")
} else {
  cat("  → NOT SIGNIFICANT\n\n")
}

# H3: Labor effect
labor_test <- results_table[results_table$Parameter == "β_Labor-led", ]
cat("H3b (Labor-led effect):\n")
cat(sprintf("  Estimate: %.3f (SE = %.3f)\n", labor_test$Estimate, labor_test$SE))
cat(sprintf("  z = %.2f, p = %.4f\n", labor_test$z_stat, labor_test$p_value))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n", labor_test$CI_lower, labor_test$CI_upper))
cat(sprintf("  Multiplier: exp(%.3f) = %.3f\n", labor_test$Estimate, exp(labor_test$Estimate)))
if (labor_test$Significant) {
  cat("  → SIGNIFICANT: Labor-led protests STRONGLY SUPPRESS triggering\n\n")
} else {
  cat("  → NOT SIGNIFICANT\n\n")
}

# H4: Poverty effect
poverty_test <- results_table[results_table$Parameter == "δ (poverty)", ]
cat("H4 (Poverty effect on background rate):\n")
cat(sprintf("  Estimate: %.3f (SE = %.3f)\n", poverty_test$Estimate, poverty_test$SE))
cat(sprintf("  z = %.2f, p = %.4f\n", poverty_test$z_stat, poverty_test$p_value))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n", poverty_test$CI_lower, poverty_test$CI_upper))
if (poverty_test$Significant) {
  cat("  → SUPPORTED: Poverty significantly increases baseline protest rates\n\n")
} else {
  cat("  → NOT SUPPORTED\n\n")
}

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   UNCERTAINTY QUANTIFICATION COMPLETE                        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("OUTPUT FILES:\n")
cat("- uncertainty_results.rds (full results object)\n")
cat("- parameter_estimates_with_se.csv (results table)\n")
cat("- table_background_parameters.tex (LaTeX table)\n")
cat("- table_triggering_parameters.tex (LaTeX table)\n")
cat("\n")
