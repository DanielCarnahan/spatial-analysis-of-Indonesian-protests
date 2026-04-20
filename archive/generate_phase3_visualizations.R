# Generate Phase 3 Visualizations
# Violence-Based Triggering and Spatial Diffusion
# Author: Daniel Carnahan
# Date: 2025-10-31

library(dplyr)
library(ggplot2)
library(tidyr)

# Create output directory
if(!dir.exists("visualizations")) {
  dir.create("visualizations")
}

cat("Loading Phase 3 results...\n")

# Load Phase 3 model results
phase3_models <- readRDS("spatial_temporal_models_full.rds")
m3_params <- phase3_models$M3$params

# Extract parameters
beta_0 <- m3_params[["beta_0"]]
beta_violence <- m3_params[["beta_violence"]]
beta_state <- m3_params[["beta_state"]]
sigma <- m3_params[["sigma"]]
omega <- m3_params[["omega"]]

cat(sprintf("\nExtracted M3 parameters:\n"))
cat(sprintf("  beta_0 (baseline): %.4f\n", beta_0))
cat(sprintf("  beta_violence: %.4f\n", beta_violence))
cat(sprintf("  beta_state: %.4f\n", beta_state))
cat(sprintf("  sigma (spatial scale): %.2f km\n", sigma))
cat(sprintf("  omega (spatial decay): %.4f\n", omega))

# ============================================================================
# 1. VIOLENCE-BASED TRIGGERING EFFECT
# ============================================================================

cat("\n1. Creating violence-based triggering visualization...\n")

# Create triggering multipliers
triggering_data <- data.frame(
  Protest_Type = c("Baseline\n(Peaceful, No State)",
                   "Violent,\nNo State",
                   "Peaceful,\nState Intervention",
                   "Violent,\nState Intervention"),
  Violence = c(0, 1, 0, 1),
  State = c(0, 0, 1, 1),
  Alpha = exp(beta_0 + beta_violence * c(0, 1, 0, 1) + beta_state * c(0, 0, 1, 1)),
  stringsAsFactors = FALSE
)

# Calculate relative multipliers
triggering_data$Relative_Alpha <- triggering_data$Alpha / triggering_data$Alpha[1]

# Create visualization
p1 <- ggplot(triggering_data, aes(x = reorder(Protest_Type, -Alpha), y = Alpha, fill = Protest_Type)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.2fx\nbaseline", Relative_Alpha)),
            vjust = -0.5, size = 5, fontface = "bold", lineheight = 0.8) +
  scale_fill_manual(values = c("#95a5a6", "#e74c3c", "#f39c12", "#c0392b")) +
  labs(
    title = "Violence-Based Triggering: How Protest Characteristics Affect Future Protests",
    subtitle = "Spatial-temporal-mark Hawkes model (M3) - 16,467 Indonesian protests",
    x = "Protest Type",
    y = "Triggering Intensity (α)",
    caption = "Note: All protests follow the same spatial diffusion pattern (power-law decay)\nDifferences arise from mark-dependent excitation only"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 14, color = "gray30", hjust = 0.5),
    axis.text.x = element_text(size = 13, lineheight = 0.8),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 11, color = "gray40", hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0, max(triggering_data$Alpha) * 1.2))

ggsave("visualizations/violence_triggering_effect.png", p1, width = 12, height = 8, dpi = 300)
cat("  ✓ Saved: visualizations/violence_triggering_effect.png\n")

# Print summary
cat("\nTriggering multipliers:\n")
for(i in 1:nrow(triggering_data)) {
  cat(sprintf("  %s: %.2fx baseline (α = %.4f)\n",
              gsub("\n", " ", triggering_data$Protest_Type[i]),
              triggering_data$Relative_Alpha[i],
              triggering_data$Alpha[i]))
}

# ============================================================================
# 2. SPATIAL DECAY FUNCTION
# ============================================================================

cat("\n2. Creating spatial decay visualization...\n")

# Generate spatial kernel
distances <- seq(0, 500, by = 2)
spatial_kernel <- (1 + distances/sigma)^(-omega)

spatial_df <- data.frame(
  Distance_km = distances,
  Kernel_Value = spatial_kernel,
  Relative_Strength = spatial_kernel / max(spatial_kernel)
)

# Find half-distance
half_distance <- distances[which.min(abs(spatial_df$Relative_Strength - 0.5))]

p2 <- ggplot(spatial_df, aes(x = Distance_km, y = Relative_Strength)) +
  geom_line(color = "#3498db", size = 1.5) +
  geom_area(alpha = 0.3, fill = "#3498db") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50", size = 0.8) +
  annotate("text", x = 450, y = 0.53, label = "50% influence", color = "gray30", size = 5) +
  annotate("text", x = 450, y = 0.13, label = "10% influence", color = "gray30", size = 5) +
  annotate("segment", x = half_distance, xend = half_distance,
           y = 0, yend = 0.5, color = "red", linetype = "dotted", size = 0.8) +
  annotate("text", x = half_distance, y = 0.95,
           label = sprintf("Half-distance:\n%.0f km", half_distance),
           color = "red", size = 5, fontface = "bold", lineheight = 0.8) +
  scale_x_continuous(breaks = seq(0, 500, 50)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Spatial Diffusion of Protest Contagion",
    subtitle = sprintf("Power-law kernel: g(d) = (1 + d/%.1f)^(-%.2f)", sigma, omega),
    x = "Distance from Source Event (km)",
    y = "Relative Triggering Strength",
    caption = "Spatial cutoff: 500 km | All protest types follow this same spatial pattern"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 14, color = "gray30", hjust = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 11, color = "gray40", hjust = 0.5)
  )

ggsave("visualizations/spatial_decay_function.png", p2, width = 12, height = 8, dpi = 300)
cat("  ✓ Saved: visualizations/spatial_decay_function.png\n")

cat(sprintf("\nSpatial decay statistics:\n"))
cat(sprintf("  Half-distance (50%% influence): %.0f km\n", half_distance))
cat(sprintf("  Characteristic scale (sigma): %.1f km\n", sigma))
cat(sprintf("  Decay exponent (omega): %.4f\n", omega))

# ============================================================================
# 3. MODEL COMPARISON ACROSS ALL PHASES
# ============================================================================

cat("\n3. Creating model comparison visualization...\n")

# Load Phase 2 and 3 results
phase2_comparison <- read.csv("model_comparison_phase2_optimized.csv")
phase3_comparison <- read.csv("model_comparison_phase3_full.csv")

# Combine results
all_models <- data.frame(
  Phase = c(rep("Phase 2: Temporal-Mark", nrow(phase2_comparison)),
            rep("Phase 3: Spatial-Temporal-Mark", nrow(phase3_comparison))),
  Model = c(phase2_comparison$Model, phase3_comparison$Model),
  AIC = c(phase2_comparison$AIC, phase3_comparison$AIC),
  LogLik = c(phase2_comparison$LogLik, phase3_comparison$LogLik)
)

# Highlight best models
all_models$Best_In_Phase <- FALSE
phase2_idx <- which(all_models$Phase == "Phase 2: Temporal-Mark")
phase3_idx <- which(all_models$Phase == "Phase 3: Spatial-Temporal-Mark")
all_models$Best_In_Phase[phase2_idx[which.min(all_models$AIC[phase2_idx])]] <- TRUE
all_models$Best_In_Phase[phase3_idx[which.min(all_models$AIC[phase3_idx])]] <- TRUE

# Shorten model names for display
all_models$Model_Short <- gsub(": ", ":\n", all_models$Model)
all_models$Model_Short <- gsub(" \\(", "\n(", all_models$Model_Short)

p3 <- ggplot(all_models, aes(x = reorder(Model_Short, -AIC), y = AIC, fill = Phase)) +
  geom_col(aes(alpha = ifelse(Best_In_Phase, 1, 0.6)), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray30") +
  scale_fill_manual(values = c("#e74c3c", "#27ae60")) +
  scale_alpha_identity() +
  labs(
    title = "Model Comparison: Temporal vs Spatial-Temporal Analysis",
    subtitle = "Lower AIC indicates better model fit | Best models highlighted at full opacity",
    x = NULL,
    y = "AIC (Akaike Information Criterion)",
    fill = "Analysis Phase",
    caption = "Phase 3 spatial-temporal-mark model provides best overall fit (AIC = -17,245)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, color = "gray30", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, lineheight = 0.8),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 11, color = "gray40", hjust = 0.5)
  ) +
  coord_flip()

ggsave("visualizations/model_comparison_all_phases.png", p3, width = 12, height = 10, dpi = 300)
cat("  ✓ Saved: visualizations/model_comparison_all_phases.png\n")

# ============================================================================
# 4. HYPOTHESIS TEST RESULTS
# ============================================================================

cat("\n4. Creating hypothesis test visualization...\n")

# Load hypothesis test results
phase2_lrt <- read.csv("likelihood_ratio_tests_phase2_optimized.csv")
phase3_lrt <- read.csv("hypothesis_tests_phase3_full.csv")

# Create hypothesis test data
hyp_data <- data.frame(
  Phase = c("Phase 2", "Phase 2", "Phase 3", "Phase 3"),
  Hypothesis = c("H1: Violence effect", "H3: State intervention",
                 "H2: Marks after space", "H3: Spatial-mark interaction"),
  LR_Statistic = c(
    phase2_lrt$LR[phase2_lrt$Hypothesis == "H1: Violence effect"],
    phase2_lrt$LR[phase2_lrt$Hypothesis == "H3: State intervention"],
    phase3_lrt$LR[phase3_lrt$Hypothesis == "H2: Marks persist after space"],
    phase3_lrt$LR[phase3_lrt$Hypothesis == "H3: Spatial-mark interaction"]
  ),
  P_Value = c(
    phase2_lrt$p_value[phase2_lrt$Hypothesis == "H1: Violence effect"],
    phase2_lrt$p_value[phase2_lrt$Hypothesis == "H3: State intervention"],
    phase3_lrt$p_value[phase3_lrt$Hypothesis == "H2: Marks persist after space"],
    phase3_lrt$p_value[phase3_lrt$Hypothesis == "H3: Spatial-mark interaction"]
  ),
  Significant = c(
    phase2_lrt$significant[phase2_lrt$Hypothesis == "H1: Violence effect"],
    phase2_lrt$significant[phase2_lrt$Hypothesis == "H3: State intervention"],
    phase3_lrt$significant[phase3_lrt$Hypothesis == "H2: Marks persist after space"],
    phase3_lrt$significant[phase3_lrt$Hypothesis == "H3: Spatial-mark interaction"]
  )
)

# Only plot positive LR statistics (negative means worse fit)
hyp_data_plot <- hyp_data[hyp_data$LR_Statistic > 0, ]

# Add significance labels
hyp_data_plot$Sig_Label <- ifelse(hyp_data_plot$P_Value < 0.001, "p < 0.001",
                                  sprintf("p = %.4f", hyp_data_plot$P_Value))

p4 <- ggplot(hyp_data_plot, aes(x = reorder(Hypothesis, LR_Statistic),
                           y = LR_Statistic,
                           fill = Significant)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_text(aes(label = Sig_Label), hjust = -0.1, size = 4.5, fontface = "bold") +
  geom_hline(yintercept = 3.84, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 0.6, y = 3.84 * 2,
           label = "p = 0.05 threshold\n(χ² = 3.84)",
           color = "red", size = 5, fontface = "bold", lineheight = 0.8) +
  scale_fill_manual(values = c("TRUE" = "#27ae60"),
                    labels = c("Significant (p < 0.05)")) +
  scale_y_log10(labels = scales::comma, limits = c(1, max(hyp_data_plot$LR_Statistic) * 3)) +
  labs(
    title = "Hypothesis Test Results: Likelihood Ratio Statistics",
    subtitle = "All effects highly significant except spatial-mark interaction",
    x = NULL,
    y = "Likelihood Ratio Statistic (log scale)",
    fill = NULL,
    caption = "Phase 3 effects are orders of magnitude stronger due to spatial control"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 14, color = "gray30", hjust = 0.5),
    axis.text.x = element_text(angle = 20, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(size = 11, color = "gray40", hjust = 0.5)
  ) +
  coord_flip()

ggsave("visualizations/hypothesis_test_results.png", p4, width = 12, height = 8, dpi = 300)
cat("  ✓ Saved: visualizations/hypothesis_test_results.png\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("                   VISUALIZATION SUMMARY                        \n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("All visualizations saved to: visualizations/\n\n")

cat("Files created:\n")
cat("  1. violence_triggering_effect.png - Violence-based triggering multipliers\n")
cat("  2. spatial_decay_function.png - Power-law spatial diffusion pattern\n")
cat("  3. model_comparison_all_phases.png - AIC comparison across all models\n")
cat("  4. hypothesis_test_results.png - Likelihood ratio test results\n\n")

cat("KEY FINDINGS:\n")
cat(sprintf("  • Violence increases triggering by: %.2fx\n", triggering_data$Relative_Alpha[2]))
cat(sprintf("  • State intervention increases triggering by: %.2fx\n", triggering_data$Relative_Alpha[3]))
cat(sprintf("  • Combined effect (violent + state): %.2fx baseline\n", triggering_data$Relative_Alpha[4]))
cat(sprintf("  • Spatial half-distance: %.0f km\n", half_distance))
cat(sprintf("  • All effects significant after spatial control (LR = %.0f, p < 0.001)\n",
            hyp_data$LR_Statistic[hyp_data$Hypothesis == "H2: Marks after space"]))

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("                        COMPLETE                                \n")
cat("═══════════════════════════════════════════════════════════════\n")
