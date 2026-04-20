################################################################################
#               GENERATE FIGURES FOR BEAMER PRESENTATION
#               Final Analysis: Cross-District Contagion Only
################################################################################

library(tidyverse)
library(scales)

theme_presentation <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom"
  )

# =============================================================================
# LOAD DATA
# =============================================================================

protests <- readRDS("protests_daily.rds")
results <- readRDS("spatial_hawkes_final_results.rds")

cat("Generating presentation figures...\n\n")

# =============================================================================
# FIGURE 1: DAILY EVENT COUNTS
# =============================================================================

cat("Figure 1: Daily event counts...\n")

daily_counts <- protests %>%
  group_by(date) %>%
  summarise(n = n(), .groups = "drop")

p1 <- ggplot(daily_counts, aes(x = date, y = n)) +
  geom_col(fill = "#2c3e50", alpha = 0.7, width = 1) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(x = NULL, y = "Daily Protest Count") +
  theme_presentation +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/fig_timeseries.pdf", p1, width = 10, height = 4, dpi = 300)
cat("  Saved: figures/fig_timeseries.pdf\n")

# =============================================================================
# FIGURE 2: PROFILE LIKELIHOOD
# =============================================================================

cat("Figure 2: Profile likelihood...\n")

profile_data <- results$profile_likelihood %>%
  mutate(
    theta_plot = ifelse(is.infinite(theta), 10000, theta),
    is_optimal = deviance == min(deviance)
  )

opt_idx <- which.min(profile_data$deviance)

p2 <- ggplot(profile_data, aes(x = theta_plot, y = deviance)) +
  geom_line(color = "#3498db", linewidth = 1.2) +
  geom_point(color = "#3498db", size = 3) +
  geom_point(data = profile_data[opt_idx, ], color = "#e74c3c", size = 6) +
  geom_vline(xintercept = profile_data$theta_plot[opt_idx],
             color = "#e74c3c", linetype = "dashed", alpha = 0.7) +
  scale_x_log10(
    breaks = c(50, 100, 200, 500, 1000, 2000, 5000, 10000),
    labels = c("50", "100", "200", "500", "1000", "2000", "5000", expression(infinity))
  ) +
  labs(
    x = expression(paste("Decay Distance ", theta, " (km)")),
    y = "Deviance"
  ) +
  annotate("text", x = 200, y = min(profile_data$deviance) + 50,
           label = expression(paste("Optimal: ", theta, " = 100 km")),
           size = 5, fontface = "bold", color = "#e74c3c", hjust = 0) +
  theme_presentation

ggsave("figures/fig_profile_likelihood.pdf", p2, width = 8, height = 5, dpi = 300)
cat("  Saved: figures/fig_profile_likelihood.pdf\n")

# =============================================================================
# FIGURE 3: MODEL COMPARISON (M0 vs M1)
# =============================================================================

cat("Figure 3: Model comparison...\n")

model_comp <- data.frame(
  model = factor(c("M0: Background\nRate Only", "M1: With Cross-District\nContagion"),
                 levels = c("M0: Background\nRate Only", "M1: With Cross-District\nContagion")),
  deviance = c(results$test_contagion$dev_m0, results$test_contagion$dev_m1)
)

p3 <- ggplot(model_comp, aes(x = model, y = deviance)) +
  geom_col(fill = c("#95a5a6", "#27ae60"), alpha = 0.8, width = 0.6) +
  geom_text(aes(label = sprintf("%.0f", deviance)), vjust = -0.5, size = 5, fontface = "bold") +
  labs(x = NULL, y = "Deviance") +
  scale_y_continuous(limits = c(0, max(model_comp$deviance) * 1.1), expand = c(0, 0)) +
  annotate("segment", x = 1, xend = 2, y = 109500, yend = 109500,
           arrow = arrow(length = unit(0.3, "cm")), color = "#e74c3c", linewidth = 1) +
  annotate("text", x = 1.5, y = 110000,
           label = sprintf("%.0f reduction\n(p < 0.001)",
                           results$test_contagion$dev_m0 - results$test_contagion$dev_m1),
           size = 4, color = "#e74c3c") +
  theme_presentation

ggsave("figures/fig_model_comparison.pdf", p3, width = 7, height = 5, dpi = 300)
cat("  Saved: figures/fig_model_comparison.pdf\n")

# =============================================================================
# FIGURE 4: LAG STRUCTURE (WEEKLY BINS)
# =============================================================================

cat("Figure 4: Lag structure (weekly)...\n")

lag_data <- data.frame(
  week = factor(c("Week 1\n(1-7)", "Week 2\n(8-14)", "Week 3\n(15-21)", "Week 4\n(22-30)"),
                levels = c("Week 1\n(1-7)", "Week 2\n(8-14)", "Week 3\n(15-21)", "Week 4\n(22-30)")),
  coef = results$lag_structure$coefficients,
  se = results$lag_structure$se
)

p4 <- ggplot(lag_data, aes(x = week, y = coef)) +
  geom_col(fill = "#3498db", alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = coef - 1.96*se, ymax = coef + 1.96*se),
                width = 0.25, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    x = "Lag Period",
    y = "Coefficient"
  ) +
  theme_presentation

ggsave("figures/fig_lag_structure.pdf", p4, width = 8, height = 5, dpi = 300)
cat("  Saved: figures/fig_lag_structure.pdf\n")

# =============================================================================
# FIGURE 5: COVARIATE EFFECTS
# =============================================================================

cat("Figure 5: Covariate effects...\n")

coefs <- results$final_model$coefficients
ses <- results$final_model$se

cov_df <- data.frame(
  variable = factor(c("Log Population", "Poverty Rate", "CPI"),
                    levels = c("Log Population", "Poverty Rate", "CPI")),
  estimate = c(coefs["log_pop"], coefs["poverty_rate"], coefs["cpi"]),
  se = c(ses["log_pop"], ses["poverty_rate"], ses["cpi"])
)

p5 <- ggplot(cov_df, aes(x = variable, y = estimate)) +
  geom_col(fill = "#9b59b6", alpha = 0.8, width = 0.6) +
  geom_errorbar(aes(ymin = estimate - 1.96*se, ymax = estimate + 1.96*se),
                width = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = NULL, y = "Coefficient (log scale)") +
  theme_presentation

ggsave("figures/fig_covariates.pdf", p5, width = 7, height = 5, dpi = 300)
cat("  Saved: figures/fig_covariates.pdf\n")

# =============================================================================
# FIGURE 6: SPATIAL KERNEL VISUALIZATION
# =============================================================================

cat("Figure 6: Spatial kernel...\n")

theta_opt <- results$optimal_theta
if (is.infinite(theta_opt)) theta_opt <- 100  # Use 100 if somehow infinite

distances <- seq(0, 500, by = 5)
kernel_df <- data.frame(
  distance = distances,
  weight = exp(-distances / theta_opt),
  weight_uniform = rep(1, length(distances))
)

p6 <- ggplot(kernel_df, aes(x = distance)) +
  geom_line(aes(y = weight), color = "#e74c3c", linewidth = 1.5) +
  geom_line(aes(y = weight_uniform), color = "#95a5a6", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = theta_opt, color = "#e74c3c", linetype = "dotted") +
  annotate("text", x = theta_opt + 10, y = 0.5,
           label = expression(paste(theta, " = 100 km")),
           hjust = 0, size = 4, color = "#e74c3c") +
  annotate("text", x = 400, y = 0.95, label = "Uniform (no decay)",
           color = "#95a5a6", size = 4) +
  annotate("text", x = 300, y = 0.15, label = "Exponential decay",
           color = "#e74c3c", size = 4) +
  labs(
    x = "Distance (km)",
    y = "Weight"
  ) +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  theme_presentation

ggsave("figures/fig_kernel.pdf", p6, width = 8, height = 5, dpi = 300)
cat("  Saved: figures/fig_kernel.pdf\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("  ALL FIGURES GENERATED\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("\nFiles created:\n")
cat("  1. figures/fig_timeseries.pdf        - Daily protest counts\n")
cat("  2. figures/fig_profile_likelihood.pdf - Profile likelihood (θ)\n")
cat("  3. figures/fig_model_comparison.pdf   - M0 vs M1 deviance\n")
cat("  4. figures/fig_lag_structure.pdf      - Temporal lag effects\n")
cat("  5. figures/fig_covariates.pdf         - Background rate controls\n")
cat("  6. figures/fig_kernel.pdf             - Spatial kernel shape\n")

cat("\nKey results:\n")
cat(sprintf("  Optimal θ: %s km\n", ifelse(is.infinite(results$optimal_theta), "∞", results$optimal_theta)))
cat(sprintf("  Cross-district effect: %.3f\n", coefs["cross_lag"]))
cat(sprintf("  Contagion exists: p = %.2e\n", results$test_contagion$p_value))
cat(sprintf("  Distance matters: p = %.4f\n", results$test_distance$p_value))
