# ============================================================================
# POPULATION DIAGNOSTICS FOR HETEROGENEOUS BACKGROUND MODELS
# ============================================================================
#
# Purpose: Investigate relationship between population and protest counts
#          to inform background rate specification
#
# Author: Analysis Script
# Date: 2025-11-06
# ============================================================================

library(dplyr)
library(ggplot2)
library(gridExtra)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║        POPULATION-PROTEST DIAGNOSTIC ANALYSIS                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Load data
cat("Loading protest data...\n")
protests <- readRDS("protests_with_population.rds")

cat(sprintf("  Total protests: %s\n", format(nrow(protests), big.mark = ",")))
cat(sprintf("  Unique districts: %d\n", length(unique(protests$district))))
cat(sprintf("  Time period: %d-%d\n", min(protests$year), max(protests$year)))
cat("\n")

# ============================================================================
# Aggregate by District
# ============================================================================

cat("Aggregating protests by district...\n")
district_summary <- protests %>%
  group_by(district, population) %>%
  summarise(
    n_protests = n(),
    log_pop = first(log_pop),
    .groups = "drop"
  ) %>%
  mutate(
    per_capita_rate = n_protests / (population / 100000),  # Per 100k population
    log_protests = log(n_protests)
  )

cat(sprintf("  Districts with protest data: %d\n", nrow(district_summary)))
cat(sprintf("  Population range: %s to %s\n",
            format(min(district_summary$population), big.mark = ","),
            format(max(district_summary$population), big.mark = ",")))
cat(sprintf("  Protest count range: %d to %d\n",
            min(district_summary$n_protests),
            max(district_summary$n_protests)))
cat("\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("SUMMARY STATISTICS:\n")
cat("-------------------\n")

# Correlation between log(population) and log(protests)
cor_log <- cor(district_summary$log_pop, district_summary$log_protests)
cat(sprintf("  Correlation (log-log): %.3f\n", cor_log))

# Simple linear regression in log-log space
lm_log <- lm(log_protests ~ log_pop, data = district_summary)
beta_log <- coef(lm_log)[2]
cat(sprintf("  Log-log regression slope: %.3f\n", beta_log))

if (abs(beta_log - 1) < 0.1) {
  cat("  → Close to 1: Suggests proportional relationship (offset model)\n")
} else if (beta_log > 1) {
  cat("  → Greater than 1: Larger districts have disproportionately MORE protests\n")
} else {
  cat("  → Less than 1: Larger districts have disproportionately FEWER protests per capita\n")
}

# Correlation between population and per-capita rate
cor_rate <- cor(log(district_summary$population), district_summary$per_capita_rate)
cat(sprintf("\n  Correlation (log(pop) vs per-capita rate): %.3f\n", cor_rate))

if (abs(cor_rate) < 0.2) {
  cat("  → Weak: Per-capita rates roughly constant (offset model appropriate)\n")
} else if (cor_rate > 0) {
  cat("  → Positive: Larger districts have HIGHER per-capita rates\n")
} else {
  cat("  → Negative: Larger districts have LOWER per-capita rates\n")
}
cat("\n")

# ============================================================================
# Create Diagnostic Plots
# ============================================================================

cat("Creating diagnostic plots...\n")

# Plot 1: Log(protests) vs Log(population)
# -----------------------------------------
p1 <- ggplot(district_summary, aes(x = population, y = n_protests)) +
  geom_point(alpha = 0.6, size = 2, color = "#2C3E50") +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C", linewidth = 1.2) +
  geom_abline(slope = 1, intercept = log(mean(district_summary$n_protests/district_summary$population)),
              linetype = "dashed", color = "#27AE60", linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  labs(
    title = "A. Protest Count vs Population (log-log scale)",
    subtitle = sprintf("OLS slope = %.3f (red line) | Proportional line (slope=1, green dashed)", beta_log),
    x = "District Population",
    y = "Total Protest Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    panel.grid.minor = element_blank()
  )

# Plot 2: Per-capita rate vs Population
# --------------------------------------
p2 <- ggplot(district_summary, aes(x = population, y = per_capita_rate)) +
  geom_point(alpha = 0.6, size = 2, color = "#2C3E50") +
  geom_smooth(method = "loess", se = TRUE, color = "#E74C3C", linewidth = 1.2) +
  geom_hline(yintercept = median(district_summary$per_capita_rate),
             linetype = "dashed", color = "#27AE60", linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "B. Per-Capita Protest Rate vs Population",
    subtitle = sprintf("Median rate = %.2f per 100k (green dashed line)",
                      median(district_summary$per_capita_rate)),
    x = "District Population (log scale)",
    y = "Protests per 100,000 Population"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    panel.grid.minor = element_blank()
  )

# Plot 3: Residuals from proportional model
# ------------------------------------------
# Expected counts under proportional model
district_summary <- district_summary %>%
  mutate(
    expected_proportional = mean(n_protests / population) * population,
    residual_proportional = n_protests - expected_proportional
  )

p3 <- ggplot(district_summary, aes(x = population, y = residual_proportional)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.6, size = 2, color = "#2C3E50") +
  geom_smooth(method = "loess", se = TRUE, color = "#E74C3C", linewidth = 1.2) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "C. Residuals from Proportional Model",
    subtitle = "Deviation from constant per-capita rate assumption",
    x = "District Population (log scale)",
    y = "Observed - Expected (Proportional Model)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    panel.grid.minor = element_blank()
  )

# Combine plots
combined <- grid.arrange(p1, p2, p3, ncol = 1,
                        top = "Population-Protest Diagnostic Plots\nIndonesian Protests 2015-2024")

# Save to PDF
pdf("population_diagnostics.pdf", width = 10, height = 12)
grid.arrange(p1, p2, p3, ncol = 1,
             top = "Population-Protest Diagnostic Plots\nIndonesian Protests 2015-2024")
dev.off()

cat("✓ Saved: population_diagnostics.pdf\n")
cat("\n")

# ============================================================================
# Model Recommendation
# ============================================================================

cat("MODEL RECOMMENDATIONS:\n")
cat("----------------------\n")

if (abs(beta_log - 1) < 0.15 && abs(cor_rate) < 0.2) {
  cat("✓ OFFSET MODEL (Option A) strongly recommended\n")
  cat("  - Log-log slope ≈ 1\n")
  cat("  - Per-capita rates roughly constant\n")
  cat("  - Use: μ(d,y) = population × exp(β₀ + year effects)\n")
} else if (beta_log > 1.15) {
  cat("→ OFFSET + INTERACTION MODEL (Option C) recommended\n")
  cat("  - Larger cities have disproportionately MORE protests\n")
  cat("  - Expect positive δ in: μ(d,y) = pop × exp(β₀ + δ·log(pop) + year)\n")
} else if (beta_log < 0.85) {
  cat("→ OFFSET + INTERACTION MODEL (Option C) recommended\n")
  cat("  - Larger cities have disproportionately FEWER protests per capita\n")
  cat("  - Expect negative δ in: μ(d,y) = pop × exp(β₀ + δ·log(pop) + year)\n")
} else {
  cat("→ Compare all three models (Options A, B, C)\n")
  cat("  - Evidence for moderate non-proportionality\n")
  cat("  - Formal likelihood comparison needed\n")
}

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║              DIAGNOSTIC ANALYSIS COMPLETE                    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")
