# ============================================================================
# MODEL DIAGNOSTICS AND COLLINEARITY ANALYSIS
# ============================================================================
#
# This script investigates why the Hessian matrix has negative eigenvalues
# and identifies potential collinearity issues.
#
# ============================================================================

library(dplyr)
library(ggplot2)
library(Rcpp)

# Try to load corrplot, but don't fail if not available
has_corrplot <- requireNamespace("corrplot", quietly = TRUE)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   MODEL DIAGNOSTICS AND COLLINEARITY ANALYSIS                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ====================================================================
# LOAD RESULTS
# ====================================================================

cat("=== LOADING RESULTS ===\n")

# Load uncertainty quantification results
uq_results <- readRDS("uncertainty_results.rds")
fit <- readRDS("model_poverty.rds")

hessian <- uq_results$hessian
vcov <- uq_results$vcov
results_table <- uq_results$results_table

param_names <- c("beta_0_bg", "gamma_raw", "delta_raw",
                 "beta_2016_raw", "beta_2017_raw", "beta_2018_raw",
                 "beta_2019_raw", "beta_2020_raw", "beta_2021_raw",
                 "beta_2022_raw", "beta_2023_raw", "beta_2024_raw",
                 "beta_0_trig_raw", "beta_riot_raw", "beta_fatal_raw",
                 "beta_student_raw", "beta_labor_raw", "decay_raw")

cat(sprintf("Loaded Hessian: %d x %d\n", nrow(hessian), ncol(hessian)))
cat(sprintf("Parameters: %d\n\n", length(param_names)))

# ====================================================================
# 1. HESSIAN EIGENVALUE ANALYSIS
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   1. HESSIAN EIGENVALUE ANALYSIS                             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Compute eigenvalues
eigen_decomp <- eigen(hessian)
eigenvalues <- eigen_decomp$values
eigenvectors <- eigen_decomp$vectors

cat("Eigenvalue Summary:\n")
cat(sprintf("  Minimum: %.6f\n", min(eigenvalues)))
cat(sprintf("  Maximum: %.6f\n", max(eigenvalues)))
cat(sprintf("  Condition number: %.2e\n", max(eigenvalues) / abs(min(eigenvalues))))
cat(sprintf("  Negative eigenvalues: %d / %d\n\n", sum(eigenvalues < 0), length(eigenvalues)))

# List eigenvalues
cat("All eigenvalues:\n")
for (i in 1:length(eigenvalues)) {
  cat(sprintf("  λ_%d = %12.6f %s\n", i, eigenvalues[i],
              if (eigenvalues[i] < 0) "[NEGATIVE]" else ""))
}
cat("\n")

# Identify parameters associated with negative eigenvalues
if (sum(eigenvalues < 0) > 0) {
  cat("Parameters contributing to negative eigenvalues:\n")
  for (i in which(eigenvalues < 0)) {
    ev <- eigenvectors[, i]
    # Find parameters with largest absolute loadings
    top_idx <- order(abs(ev), decreasing = TRUE)[1:3]
    cat(sprintf("  Eigenvalue %d (λ = %.6f):\n", i, eigenvalues[i]))
    for (j in top_idx) {
      cat(sprintf("    %s: %.4f\n", param_names[j], ev[j]))
    }
  }
  cat("\n")
}

# ====================================================================
# 2. PARAMETER CORRELATION MATRIX
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   2. PARAMETER CORRELATION MATRIX                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Compute correlation matrix from covariance matrix
# Handle NaNs by replacing with 0
vcov_clean <- vcov
vcov_clean[is.nan(vcov_clean)] <- 0
vcov_clean[is.infinite(vcov_clean)] <- 0

# Get standard deviations
sds <- sqrt(diag(vcov_clean))
sds[sds == 0] <- 1  # Avoid division by zero

# Compute correlation
cor_matrix <- vcov_clean / outer(sds, sds)
rownames(cor_matrix) <- param_names
colnames(cor_matrix) <- param_names

# Find highly correlated pairs
high_cor_threshold <- 0.8
cat(sprintf("Highly correlated parameter pairs (|r| > %.2f):\n", high_cor_threshold))

found_high_cor <- FALSE
for (i in 1:(nrow(cor_matrix) - 1)) {
  for (j in (i + 1):ncol(cor_matrix)) {
    if (!is.nan(cor_matrix[i, j]) && abs(cor_matrix[i, j]) > high_cor_threshold) {
      cat(sprintf("  %s <-> %s: r = %.4f\n",
                  param_names[i], param_names[j], cor_matrix[i, j]))
      found_high_cor <- TRUE
    }
  }
}

if (!found_high_cor) {
  cat("  (None found - but many correlations are NaN)\n")
}
cat("\n")

# Save correlation matrix visualization
png("correlation_matrix.png", width = 1200, height = 1000, res = 120)
if (has_corrplot) {
  corrplot::corrplot(cor_matrix, method = "color", type = "lower",
                     tl.col = "black", tl.srt = 45, tl.cex = 0.7,
                     cl.cex = 0.7, mar = c(0, 0, 2, 0),
                     title = "Parameter Correlation Matrix",
                     na.label = "X", na.label.col = "gray")
} else {
  # Fallback to base R heatmap
  heatmap(cor_matrix, Rowv = NA, Colv = NA, scale = "none",
          main = "Parameter Correlation Matrix", margins = c(12, 12))
}
dev.off()
cat("✓ Saved: correlation_matrix.png\n\n")

# ====================================================================
# 3. PARAMETER IDENTIFICATION ISSUES
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   3. PARAMETER IDENTIFICATION ISSUES                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Check which parameters have NaN standard errors
nan_se_params <- results_table$Parameter[is.nan(results_table$SE)]

cat(sprintf("Parameters with NaN standard errors (%d / %d):\n",
            length(nan_se_params), nrow(results_table)))
for (param in nan_se_params) {
  estimate <- results_table$Estimate[results_table$Parameter == param]
  cat(sprintf("  %s: estimate = %.4f\n", param, estimate))
}
cat("\n")

# Check if parameters are hitting bounds
cat("Parameter bounds check:\n")
cat("  (Sigmoid transformations map unbounded → bounded)\n\n")

# Gamma: [0.01, 2]
gamma_raw <- fit$params$gamma_raw
gamma_trans <- 0.01 + 1.99 / (1 + exp(-gamma_raw))
cat(sprintf("  γ (population): %.4f (raw = %.4f, range [0.01, 2.00])\n",
            gamma_trans, gamma_raw))
if (abs(gamma_trans - 0.01) < 0.01) cat("    ⚠ Near lower bound\n")
if (abs(gamma_trans - 2.00) < 0.01) cat("    ⚠ Near upper bound\n")

# Delta: [-20, 20]
delta_raw <- fit$params$delta_raw
delta_trans <- -20 + 40 / (1 + exp(-delta_raw))
cat(sprintf("  δ (poverty): %.4f (raw = %.4f, range [-20, 20])\n",
            delta_trans, delta_raw))
if (abs(delta_trans - (-20)) < 1) cat("    ⚠ Near lower bound\n")
if (abs(delta_trans - 20) < 1) cat("    ⚠ Near upper bound\n")

# Year effects: [-3, 3]
year_params <- c("beta_2016_raw", "beta_2017_raw", "beta_2018_raw",
                 "beta_2019_raw", "beta_2020_raw", "beta_2021_raw",
                 "beta_2022_raw", "beta_2023_raw", "beta_2024_raw")

cat("\n  Year effects (range [-3, 3]):\n")
for (param in year_params) {
  raw_val <- fit$params[[param]]
  trans_val <- -3 + 6 / (1 + exp(-raw_val))
  year_num <- as.numeric(gsub("beta_|_raw", "", param))
  cat(sprintf("    %d: %.4f", year_num, trans_val))
  if (abs(trans_val - (-3)) < 0.1) cat(" ⚠ Near lower bound")
  if (abs(trans_val - 3) < 0.1) cat(" ⚠ Near upper bound")
  cat("\n")
}

# Mark effects: [-5, 5]
mark_params <- c("beta_riot_raw", "beta_fatal_raw", "beta_student_raw", "beta_labor_raw")
mark_names <- c("Riot", "Fatalities", "Student", "Labor")

cat("\n  Mark effects (range [-5, 5]):\n")
for (i in 1:length(mark_params)) {
  raw_val <- fit$params[[mark_params[i]]]
  trans_val <- -5 + 10 / (1 + exp(-raw_val))
  cat(sprintf("    %s: %.4f", mark_names[i], trans_val))
  if (abs(trans_val - (-5)) < 0.1) cat(" ⚠ Near lower bound")
  if (abs(trans_val - 5) < 0.1) cat(" ⚠ Near upper bound")
  cat("\n")
}

cat("\n")

# ====================================================================
# 4. DATA CORRELATION CHECK
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   4. DATA CORRELATION CHECK                                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Load data
protests <- readRDS("protests_with_poverty.rds")

# Correlation between log(pop) and poverty
cor_log_pop_poverty <- cor(protests$log_pop, protests$poverty_decimal, use = "complete.obs")
cat(sprintf("Correlation (log_pop, poverty_decimal): %.3f\n", cor_log_pop_poverty))

# Check correlation by year
cat("\nCorrelation by year:\n")
year_cors <- protests %>%
  group_by(year) %>%
  summarize(
    cor = cor(log_pop, poverty_decimal, use = "complete.obs"),
    n = n(),
    .groups = "drop"
  )

for (i in 1:nrow(year_cors)) {
  cat(sprintf("  %d: r = %.3f (n = %d)\n",
              year_cors$year[i], year_cors$cor[i], year_cors$n[i]))
}

cat("\n")

# ====================================================================
# 5. PARAMETER STABILITY ACROSS STARTING POINTS
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   5. PARAMETER STABILITY                                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("Note: Full multi-start comparison requires loading all checkpoints.\n")
cat("Currently showing best fit only.\n")
cat(sprintf("  Best LL: %.2f\n", fit$loglik))
cat(sprintf("  Convergence: %d\n", fit$convergence))
cat("\n")

# ====================================================================
# 6. MODEL SIMPLIFICATION RECOMMENDATIONS
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   6. RECOMMENDATIONS                                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("Based on diagnostic analysis:\n\n")

# Issue 1: Negative eigenvalues
cat("1. NEGATIVE EIGENVALUES:\n")
cat("   Problem: Hessian is not positive definite\n")
cat("   Likely cause: Model overparameterization\n")
cat("   → Some parameters are not well-identified\n")
cat("   → Flat likelihood surface in some directions\n\n")

# Issue 2: NaN standard errors
cat("2. NaN STANDARD ERRORS:\n")
cat(sprintf("   %d / %d parameters have undefined SEs\n", length(nan_se_params), nrow(results_table)))
cat("   Affected: ")
if (length(nan_se_params) > 0) {
  cat(paste(head(nan_se_params, 5), collapse = ", "))
  if (length(nan_se_params) > 5) cat(", ...")
  cat("\n")
}
cat("   → These parameters cannot be reliably estimated\n\n")

# Issue 3: Large SEs
large_se_threshold <- 10
large_se_params <- results_table %>%
  filter(!is.nan(SE), SE > large_se_threshold)

cat("3. PARAMETERS WITH LARGE STANDARD ERRORS (SE > 10):\n")
if (nrow(large_se_params) > 0) {
  for (i in 1:nrow(large_se_params)) {
    cat(sprintf("   %s: SE = %.2f\n",
                large_se_params$Parameter[i], large_se_params$SE[i]))
  }
  cat("   → High uncertainty, poorly identified\n\n")
} else {
  cat("   (None beyond those with NaN)\n\n")
}

# Recommendations
cat("RECOMMENDED ACTIONS:\n")
cat("-------------------\n\n")

cat("Option A: SIMPLIFY TEMPORAL STRUCTURE\n")
cat("  • Replace 9 year fixed effects with linear trend\n")
cat("  • Reduces parameters from 18 → 10\n")
cat("  • Previous test: LL = -113,113 (much worse)\n")
cat("  • But may improve identifiability\n\n")

cat("Option B: DROP POVERTY\n")
cat("  • Poverty effect has huge SE (33.15) and is non-significant\n")
cat("  • May be collinear with log(pop) despite r = -0.48\n")
cat("  • Test: Fit model without delta parameter\n")
cat("  • Reduces parameters from 18 → 17\n\n")

cat("Option C: COMBINE YEAR EFFECTS\n")
cat("  • Group years into periods (e.g., 2015-2017, 2018-2020, 2021-2024)\n")
cat("  • Reduces 9 year effects → 2-3 period effects\n")
cat("  • Reduces parameters from 18 → 11-12\n\n")

cat("Option D: SIMPLIFY MARKS\n")
cat("  • Student and labor both suppress triggering\n")
cat("  • Could combine into single \"organized\" indicator\n")
cat("  • Reduces parameters from 18 → 17\n\n")

cat("Option E: USE BOOTSTRAP (NO SIMPLIFICATION)\n")
cat("  • Accept that Hessian-based SEs don't work\n")
cat("  • Use bootstrap resampling for inference\n")
cat("  • Computationally expensive but robust\n")
cat("  • Recommended if you want to keep full model\n\n")

cat("PRIORITY RECOMMENDATION:\n")
cat("  1. Try Option B (drop poverty) - smallest change, may fix issue\n")
cat("  2. If still problematic, try Option C (combine year effects)\n")
cat("  3. If inference still needed, use Option E (bootstrap)\n\n")

# ====================================================================
# SAVE DIAGNOSTIC REPORT
# ====================================================================

cat("=== SAVING DIAGNOSTIC REPORT ===\n")

diagnostic_summary <- list(
  eigenvalues = eigenvalues,
  condition_number = max(eigenvalues) / abs(min(eigenvalues)),
  negative_eigenvalues = sum(eigenvalues < 0),
  correlation_matrix = cor_matrix,
  nan_se_parameters = nan_se_params,
  large_se_parameters = large_se_params,
  data_correlation = cor_log_pop_poverty,
  recommendations = c(
    "Option A: Linear trend instead of year FX",
    "Option B: Drop poverty parameter",
    "Option C: Combine year effects into periods",
    "Option D: Combine student/labor marks",
    "Option E: Use bootstrap for inference"
  )
)

saveRDS(diagnostic_summary, "model_diagnostics.rds")
cat("✓ Saved: model_diagnostics.rds\n")
cat("✓ Saved: correlation_matrix.png\n\n")

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   DIAGNOSTICS COMPLETE                                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("Next steps:\n")
cat("1. Review recommendations above\n")
cat("2. Decide on model simplification strategy\n")
cat("3. Re-fit simplified model\n")
cat("4. If needed, use bootstrap for final inference\n")
cat("\n")
