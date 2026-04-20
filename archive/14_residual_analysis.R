# ============================================================================
# RESIDUAL ANALYSIS FOR HETEROGENEOUS HAWKES MODEL
# ============================================================================
#
# Purpose: Diagnose model fit issues by examining residuals
#          - By district population
#          - By year
#          - By district characteristics
#          - Systematic over/under-prediction patterns
#
# ============================================================================

library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║     RESIDUAL ANALYSIS: HETEROGENEOUS HAWKES MODEL            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ====================================================================
# 1. LOAD DATA AND MODEL
# ====================================================================

cat("=== LOADING DATA AND MODEL RESULTS ===\n")

# Load model results
het_model <- readRDS("heterogeneous_hawkes_full.rds")
protests <- readRDS("protests_with_population.rds")

cat(sprintf("Loaded heterogeneous Hawkes model (full dataset)\n"))
cat(sprintf("  Log-likelihood: %.2f\n", het_model$loglik))
cat(sprintf("  AIC: %.2f\n", het_model$AIC))
cat(sprintf("  BIC: %.2f\n", het_model$BIC))
cat(sprintf("  Convergence: %d\n\n", het_model$convergence))

# ====================================================================
# 2. CALCULATE PREDICTED BACKGROUND RATES
# ====================================================================

cat("=== CALCULATING PREDICTED BACKGROUND RATES ===\n")

# Extract parameters
params <- het_model$params
beta_0 <- params$beta_0_bg
gamma <- params$gamma

cat(sprintf("Population elasticity (γ): %.4f\n", gamma))
cat(sprintf("Baseline intercept (β₀): %.4f\n\n", beta_0))

# Calculate predicted background rate for each district-year
district_year_data <- protests %>%
  group_by(district, year, population, log_pop) %>%
  summarise(
    n_observed = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Background rate: μ(d,y) = exp(β₀ + γ·log(pop) + β_year)
    year_effect = case_when(
      year == 2015 ~ 0,  # baseline
      year == 2016 ~ params$beta_2016,
      year == 2017 ~ params$beta_2017,
      year == 2018 ~ params$beta_2018,
      year == 2019 ~ params$beta_2019,
      year == 2020 ~ params$beta_2020,
      year == 2021 ~ params$beta_2021,
      year == 2022 ~ params$beta_2022,
      year == 2023 ~ params$beta_2023,
      year == 2024 ~ params$beta_2024,
      TRUE ~ 0
    ),
    mu_predicted = exp(beta_0 + gamma * log_pop + year_effect),
    # Approximate observation time (in days)
    # Simplified: assume full year observation
    days_obs = 365.25,
    # Expected counts = rate × time
    n_expected = mu_predicted * days_obs,
    # Residuals
    residual_raw = n_observed - n_expected,
    residual_pearson = (n_observed - n_expected) / sqrt(n_expected + 1),
    residual_deviance = sign(n_observed - n_expected) *
      sqrt(2 * (n_observed * log((n_observed + 0.5) / (n_expected + 0.5)) -
                (n_observed - n_expected)))
  )

cat(sprintf("Calculated predictions for %d district-year combinations\n\n",
            nrow(district_year_data)))

# ====================================================================
# 3. AGGREGATE BY DISTRICT
# ====================================================================

cat("=== AGGREGATING BY DISTRICT ===\n")

district_data <- district_year_data %>%
  group_by(district, population, log_pop) %>%
  summarise(
    n_years = n(),
    total_observed = sum(n_observed),
    total_expected = sum(n_expected),
    mean_residual = mean(residual_raw),
    mean_pearson = mean(residual_pearson),
    .groups = "drop"
  ) %>%
  mutate(
    residual_total = total_observed - total_expected,
    pct_error = 100 * (total_observed - total_expected) / total_expected,
    per_capita_observed = total_observed / (population / 100000),
    per_capita_expected = total_expected / (population / 100000)
  )

cat(sprintf("Aggregated to %d districts\n\n", nrow(district_data)))

# ====================================================================
# 4. SUMMARY STATISTICS
# ====================================================================

cat("RESIDUAL SUMMARY STATISTICS:\n")
cat("----------------------------\n")

cat(sprintf("Mean residual: %.2f\n", mean(district_year_data$residual_raw)))
cat(sprintf("Median residual: %.2f\n", median(district_year_data$residual_raw)))
cat(sprintf("SD residual: %.2f\n\n", sd(district_year_data$residual_raw)))

# Correlation between residuals and population
cor_resid_pop <- cor(district_data$residual_total, log(district_data$population))
cat(sprintf("Correlation (residual vs log(pop)): %.3f\n", cor_resid_pop))

if (abs(cor_resid_pop) > 0.3) {
  cat("  ⚠ STRONG systematic bias by population size!\n")
} else if (abs(cor_resid_pop) > 0.15) {
  cat("  ⚠ Moderate systematic bias by population size\n")
} else {
  cat("  ✓ Low bias by population size\n")
}
cat("\n")

# Top under-predicted districts
cat("TOP 10 UNDER-PREDICTED DISTRICTS (model predicts too few):\n")
under <- district_data %>%
  arrange(pct_error) %>%
  head(10) %>%
  select(district, population, total_observed, total_expected, pct_error)

print(as.data.frame(under), row.names = FALSE)
cat("\n")

# Top over-predicted districts
cat("TOP 10 OVER-PREDICTED DISTRICTS (model predicts too many):\n")
over <- district_data %>%
  arrange(desc(pct_error)) %>%
  head(10) %>%
  select(district, population, total_observed, total_expected, pct_error)

print(as.data.frame(over), row.names = FALSE)
cat("\n")

# ====================================================================
# 5. CREATE DIAGNOSTIC PLOTS
# ====================================================================

cat("=== CREATING RESIDUAL DIAGNOSTIC PLOTS ===\n")

# Plot 1: Residuals vs Population
p1 <- ggplot(district_data, aes(x = population, y = residual_total)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.6, size = 2, color = "#2C3E50") +
  geom_smooth(method = "loess", se = TRUE, color = "#E74C3C", linewidth = 1.2) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "A. Residuals vs District Population",
    subtitle = sprintf("Correlation = %.3f (systematic bias %s)",
                      cor_resid_pop,
                      ifelse(abs(cor_resid_pop) > 0.3, "PRESENT", "minimal")),
    x = "District Population (log scale)",
    y = "Total Residual (Observed - Expected)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray30")
  )

# Plot 2: Percent error vs Population
p2 <- ggplot(district_data, aes(x = population, y = pct_error)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.6, size = 2, color = "#2C3E50") +
  geom_smooth(method = "loess", se = TRUE, color = "#E74C3C", linewidth = 1.2) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "B. Percent Prediction Error vs Population",
    subtitle = "Relative bias across population sizes",
    x = "District Population (log scale)",
    y = "Prediction Error (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray30")
  )

# Plot 3: Observed vs Expected
p3 <- ggplot(district_data, aes(x = total_expected, y = total_observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = log10(population)), alpha = 0.6, size = 2) +
  scale_color_viridis_c(name = "log10(Pop)", option = "plasma") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "C. Observed vs Expected Protest Counts",
    subtitle = "Perfect prediction = diagonal line",
    x = "Expected Protests (model prediction)",
    y = "Observed Protests (actual data)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    legend.position = "right"
  )

# Plot 4: Per-capita rates - observed vs expected
p4 <- ggplot(district_data, aes(x = population)) +
  geom_point(aes(y = per_capita_observed, color = "Observed"),
             alpha = 0.5, size = 2) +
  geom_point(aes(y = per_capita_expected, color = "Expected (Model)"),
             alpha = 0.5, size = 2) +
  geom_smooth(aes(y = per_capita_observed, color = "Observed"),
              method = "loess", se = FALSE, linewidth = 1.2) +
  geom_smooth(aes(y = per_capita_expected, color = "Expected (Model)"),
              method = "loess", se = FALSE, linewidth = 1.2) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = c("Observed" = "#E74C3C",
                                 "Expected (Model)" = "#3498DB")) +
  labs(
    title = "D. Per-Capita Rates: Observed vs Model Predictions",
    subtitle = "Does the model capture population effect correctly?",
    x = "District Population (log scale)",
    y = "Protests per 100,000 Population",
    color = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    legend.position = "bottom"
  )

# Combine and save
pdf("heterogeneous_model_residuals.pdf", width = 12, height = 14)
grid.arrange(p1, p2, p3, p4, ncol = 1,
             top = "Residual Diagnostics: Heterogeneous Hawkes Model\nIndonesian Protests 2015-2024")
dev.off()

cat("✓ Saved: heterogeneous_model_residuals.pdf\n\n")

# ====================================================================
# 6. YEAR EFFECTS ANALYSIS
# ====================================================================

cat("=== TEMPORAL RESIDUAL PATTERNS ===\n")

year_residuals <- district_year_data %>%
  group_by(year) %>%
  summarise(
    n_districts = n(),
    total_observed = sum(n_observed),
    total_expected = sum(n_expected),
    mean_residual = mean(residual_raw),
    median_residual = median(residual_raw),
    mean_pct_error = mean((n_observed - n_expected) / (n_expected + 1) * 100),
    .groups = "drop"
  )

cat("\nRESIDUALS BY YEAR:\n")
print(as.data.frame(year_residuals), row.names = FALSE)
cat("\n")

# ====================================================================
# 7. CONCLUSIONS
# ====================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║                  DIAGNOSTIC CONCLUSIONS                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("POPULATION EFFECT DIAGNOSIS:\n")
cat("----------------------------\n")
if (abs(cor_resid_pop) > 0.3) {
  cat("⚠ STRONG systematic bias detected!\n")
  cat("  - Residuals significantly correlated with population\n")
  cat("  - Current power-law specification (γ·log(pop)) is inadequate\n")
  cat("  - Recommend: Model C (offset + interaction) to allow flexible\n")
  cat("    per-capita rates that vary with population\n\n")
} else {
  cat("✓ Population effect reasonably well-captured\n")
  cat("  - Residuals show weak correlation with population\n")
  cat("  - Current specification may be adequate\n\n")
}

cat("NEXT STEPS:\n")
cat("-----------\n")
cat("1. Fit alternative Model C (offset + interaction)\n")
cat("2. Compare AIC/BIC across models\n")
cat("3. Re-examine residuals from best-fitting model\n")
cat("4. Consider additional covariates if systematic patterns remain\n\n")

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║           RESIDUAL ANALYSIS COMPLETE                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
