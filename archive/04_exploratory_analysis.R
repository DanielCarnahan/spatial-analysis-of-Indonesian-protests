# Phase 3: Exploratory Spatial Analysis
# Test for clustering and space-time interaction
# Author: Daniel Carnahan
# Date: 2025-10-27

library(spatstat)
library(dplyr)
library(ggplot2)
library(sf)
library(viridis)

# Create output directory for plots
if(!dir.exists("plots")) dir.create("plots")

# Load point pattern objects
cat("Loading point pattern objects...\n")
protest_ppp <- readRDS("protest_ppp.rds")
protest_stpp <- readRDS("protest_stpp.rds")
peaceful_ppp <- readRDS("peaceful_ppp.rds")
violent_ppp <- readRDS("violent_ppp.rds")
jakarta_ppp <- readRDS("jakarta_ppp.rds")

cat("Loaded", npoints(protest_ppp), "events\n\n")

# ====================================================================
# 1. KERNEL DENSITY ESTIMATION
# ====================================================================

cat("=== KERNEL DENSITY ESTIMATION ===\n")

# Estimate density (this may take a moment for large datasets)
cat("Computing kernel density (this may take a minute)...\n")
protest_density <- density(protest_ppp, sigma = 0.5)  # sigma in degrees (~55km)

cat("Density estimation complete\n")
cat("  Bandwidth (sigma):", 0.5, "degrees (~55 km)\n")
cat("  Max density:", round(max(protest_density), 2), "events/unit area\n\n")

# Save density plot
png("plots/01_kernel_density.png", width = 1200, height = 800, res = 120)
plot(protest_density,
     main = "Kernel Density of Indonesian Protests (2015-2024)",
     ribbon = TRUE,
     col = viridis(256))
plot(protest_ppp, add = TRUE, col = rgb(1,1,1,0.3), pch = 16, cex = 0.3)
dev.off()
cat("Saved: plots/01_kernel_density.png\n\n")

# Density for violent vs peaceful
cat("Computing densities for protest types...\n")
peaceful_density <- density(peaceful_ppp, sigma = 0.5)
violent_density <- density(violent_ppp, sigma = 0.5)

png("plots/02_peaceful_vs_violent_density.png", width = 1600, height = 600, res = 120)
par(mfrow = c(1, 2))
plot(peaceful_density,
     main = "Peaceful Protests",
     ribbon = TRUE,
     col = viridis(256))
plot(violent_density,
     main = "Violent Protests/Riots",
     ribbon = TRUE,
     col = viridis(256))
dev.off()
cat("Saved: plots/02_peaceful_vs_violent_density.png\n\n")

# ====================================================================
# 2. RIPLEY'S K-FUNCTION (Test for clustering)
# ====================================================================

cat("=== RIPLEY'S K-FUNCTION (Clustering Analysis) ===\n")

# K-function tests if points are more clustered than random
# This is computationally intensive, so we'll use a smaller distance range
cat("Computing Ripley's K-function (this may take a few minutes)...\n")

# Use smaller subset of distances for computational efficiency
r_values <- seq(0, 2, by = 0.1)  # in degrees (~0-220 km)

K_result <- Kest(protest_ppp, r = r_values, correction = "border")

cat("K-function computation complete\n\n")

# Plot K-function
png("plots/03_ripley_k_function.png", width = 1000, height = 800, res = 120)
plot(K_result,
     main = "Ripley's K-function: Testing for Clustering",
     legend = TRUE)
# Add theoretical line for Complete Spatial Randomness (CSR)
abline(0, 1, lty = 2, col = "red", lwd = 2)
dev.off()
cat("Saved: plots/03_ripley_k_function.png\n")
cat("  If K(r) > πr² (above red line), events are CLUSTERED\n")
cat("  If K(r) < πr² (below red line), events are DISPERSED\n\n")

# L-function (variance-stabilized version of K)
cat("Computing L-function (variance-stabilized K)...\n")
L_result <- Lest(protest_ppp, r = r_values, correction = "border")

png("plots/04_ripley_l_function.png", width = 1000, height = 800, res = 120)
plot(L_result,
     main = "Ripley's L-function (Variance Stabilized)",
     legend = TRUE)
abline(0, 1, lty = 2, col = "red", lwd = 2)
dev.off()
cat("Saved: plots/04_ripley_l_function.png\n")
cat("  L(r) > r indicates clustering\n")
cat("  L(r) < r indicates dispersion\n\n")

# ====================================================================
# 3. PAIR CORRELATION FUNCTION
# ====================================================================

cat("=== PAIR CORRELATION FUNCTION ===\n")

# g(r) = derivative of K(r), shows clustering at different scales
cat("Computing pair correlation function...\n")
g_result <- pcf(protest_ppp, r = r_values, correction = "best")

png("plots/05_pair_correlation.png", width = 1000, height = 800, res = 120)
plot(g_result,
     main = "Pair Correlation Function",
     legend = TRUE)
abline(h = 1, lty = 2, col = "red", lwd = 2)
dev.off()
cat("Saved: plots/05_pair_correlation.png\n")
cat("  g(r) > 1 indicates clustering at distance r\n")
cat("  g(r) = 1 indicates randomness\n")
cat("  g(r) < 1 indicates dispersion\n\n")

# ====================================================================
# 4. QUADRAT TEST (Test for Complete Spatial Randomness)
# ====================================================================

cat("=== QUADRAT TEST ===\n")

# Divide region into quadrats and test if counts are Poisson distributed
cat("Performing quadrat test...\n")

quadrat_result <- quadrat.test(protest_ppp, nx = 10, ny = 10)

cat("Quadrat Test Results:\n")
print(quadrat_result)
cat("\n")

png("plots/06_quadrat_counts.png", width = 1000, height = 800, res = 120)
plot(quadrat.test(protest_ppp, nx = 10, ny = 10),
     main = "Quadrat Counts (10x10 grid)")
dev.off()
cat("Saved: plots/06_quadrat_counts.png\n")
cat("  Shows spatial heterogeneity in protest intensity\n\n")

# ====================================================================
# 5. NEAREST NEIGHBOR ANALYSIS
# ====================================================================

cat("=== NEAREST NEIGHBOR ANALYSIS ===\n")

# G-function: distribution of nearest neighbor distances
cat("Computing G-function (nearest neighbor distances)...\n")
G_result <- Gest(protest_ppp, r = r_values, correction = "best")

png("plots/07_nearest_neighbor_G.png", width = 1000, height = 800, res = 120)
plot(G_result,
     main = "G-function: Nearest Neighbor Distances",
     legend = TRUE)
dev.off()
cat("Saved: plots/07_nearest_neighbor_G.png\n\n")

# ====================================================================
# 6. SPACE-TIME INTERACTION (Knox Test)
# ====================================================================

cat("=== SPACE-TIME INTERACTION (Knox Test) ===\n")

# Knox test: Are events that are close in space also close in time?
# This is KEY for detecting contagion!

cat("Computing space-time interaction...\n")
cat("  Spatial threshold: 0.5 degrees (~55 km)\n")
cat("  Temporal threshold: 30 days\n\n")

# Function to perform Knox-like test
knox_test <- function(data, spatial_thresh, temporal_thresh, n_sim = 99) {

  n <- nrow(data)

  # Calculate observed space-time interactions
  observed <- 0

  cat("  Counting observed close pairs (this may take a moment)...\n")

  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      spatial_dist <- sqrt((data$x[i] - data$x[j])^2 +
                          (data$y[i] - data$y[j])^2)
      temporal_dist <- abs(data$t[i] - data$t[j])

      if(spatial_dist < spatial_thresh & temporal_dist < temporal_thresh) {
        observed <- observed + 1
      }
    }
  }

  cat("  Observed close pairs:", observed, "\n")

  # Monte Carlo simulation
  expected_values <- numeric(n_sim)

  cat("  Running", n_sim, "Monte Carlo simulations...\n")

  for(sim in 1:n_sim) {
    # Permute times randomly
    shuffled_data <- data
    shuffled_data$t <- sample(data$t)

    count <- 0
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        spatial_dist <- sqrt((shuffled_data$x[i] - shuffled_data$x[j])^2 +
                            (shuffled_data$y[i] - shuffled_data$y[j])^2)
        temporal_dist <- abs(shuffled_data$t[i] - shuffled_data$t[j])

        if(spatial_dist < spatial_thresh & temporal_dist < temporal_thresh) {
          count <- count + 1
        }
      }
    }
    expected_values[sim] <- count
  }

  expected_mean <- mean(expected_values)
  p_value <- sum(expected_values >= observed) / n_sim

  return(list(
    observed = observed,
    expected = expected_mean,
    p_value = p_value,
    simulated = expected_values
  ))
}

# Run on a SAMPLE (full dataset would be too slow)
# Sample 1000 events for computational feasibility
set.seed(123)
sample_indices <- sample(1:nrow(protest_stpp), min(1000, nrow(protest_stpp)))
protest_sample <- protest_stpp[sample_indices, ]

cat("Running Knox test on sample of", nrow(protest_sample), "events...\n")
knox_result <- knox_test(protest_sample,
                         spatial_thresh = 0.5,
                         temporal_thresh = 30,
                         n_sim = 99)

cat("\nKnox Test Results (Sample):\n")
cat("  Observed close pairs:", knox_result$observed, "\n")
cat("  Expected under randomness:", round(knox_result$expected, 2), "\n")
cat("  Ratio (Obs/Exp):", round(knox_result$observed / knox_result$expected, 2), "\n")
cat("  P-value:", knox_result$p_value, "\n")

if(knox_result$p_value < 0.05) {
  cat("  ** SIGNIFICANT space-time interaction detected! **\n")
  cat("  ** This suggests CONTAGION: protests close in space tend to be close in time **\n")
} else {
  cat("  No significant space-time interaction\n")
}

# Save Knox test results
saveRDS(knox_result, "knox_test_result.rds")
cat("\nSaved: knox_test_result.rds\n\n")

# ====================================================================
# 7. TEMPORAL PATTERNS
# ====================================================================

cat("=== TEMPORAL PATTERNS ===\n")

# Plot events over time
png("plots/08_temporal_pattern.png", width = 1200, height = 600, res = 120)

temporal_plot <- protest_stpp %>%
  mutate(month = as.Date("2015-01-01") + t) %>%
  mutate(year_month = format(month, "%Y-%m")) %>%
  group_by(year_month) %>%
  summarise(n = n(), .groups = 'drop')

ggplot(temporal_plot, aes(x = as.Date(paste0(year_month, "-01")), y = n)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_smooth(method = "loess", color = "red", se = TRUE, alpha = 0.2) +
  labs(title = "Temporal Evolution of Protests (2015-2024)",
       x = "Date",
       y = "Number of Protests per Month") +
  theme_minimal() +
  theme(text = element_text(size = 12))

dev.off()
cat("Saved: plots/08_temporal_pattern.png\n\n")

# ====================================================================
# 8. SUMMARY OF FINDINGS
# ====================================================================

cat("=== EXPLORATORY ANALYSIS SUMMARY ===\n")
cat("\nGenerated plots:\n")
cat("  1. plots/01_kernel_density.png - Overall spatial intensity\n")
cat("  2. plots/02_peaceful_vs_violent_density.png - By protest type\n")
cat("  3. plots/03_ripley_k_function.png - Clustering test\n")
cat("  4. plots/04_ripley_l_function.png - Variance-stabilized clustering\n")
cat("  5. plots/05_pair_correlation.png - Scale-specific clustering\n")
cat("  6. plots/06_quadrat_counts.png - Spatial heterogeneity\n")
cat("  7. plots/07_nearest_neighbor_G.png - Nearest neighbor distribution\n")
cat("  8. plots/08_temporal_pattern.png - Time series\n")

cat("\nKey findings:\n")
cat("  - Spatial clustering: Check K and L functions\n")
cat("  - Space-time interaction:",
    ifelse(knox_result$p_value < 0.05, "DETECTED", "Not detected"), "\n")
cat("  - Next step: Fit self-exciting point process models\n")

cat("\n=== COMPLETE ===\n")
