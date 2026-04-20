# Hawkes Process Models for Protest Contagion
# Formal modeling of temporal self-excitation with marks
# Author: Daniel Carnahan
# Date: 2025-10-27

library(dplyr)
library(ggplot2)

cat("=== HAWKES PROCESS MODELING ===\n")
cat("Goal: Estimate temporal contagion parameters\n\n")

# Load data
protests <- readRDS("protests_prepared.rds")
protest_stpp <- readRDS("protest_stpp.rds")

cat("Loaded", nrow(protests), "events\n\n")

# ====================================================================
# 1. CHECK FOR REQUIRED PACKAGES
# ====================================================================

cat("=== Checking for required packages ===\n")

packages_needed <- c("hawkes", "ppmlasso", "stpp")
packages_available <- sapply(packages_needed, requireNamespace, quietly = TRUE)

cat("Package availability:\n")
for(i in seq_along(packages_needed)) {
  status <- if(packages_available[i]) "✓ Available" else "✗ Not installed"
  cat(sprintf("  %s: %s\n", packages_needed[i], status))
}
cat("\n")

if(!all(packages_available)) {
  cat("Some packages are missing. Installing...\n")
  cat("You may need to run: install.packages(c('hawkes', 'ppmlasso', 'stpp'))\n\n")
  cat("Attempting installation...\n")

  for(pkg in packages_needed[!packages_available]) {
    tryCatch({
      install.packages(pkg, repos = "http://cran.us.r-project.org")
      cat(sprintf("  Installed: %s\n", pkg))
    }, error = function(e) {
      cat(sprintf("  Failed to install %s: %s\n", pkg, e$message))
    })
  }
}

# ====================================================================
# 2. IMPLEMENT BASIC HAWKES ESTIMATION
# ====================================================================

cat("\n=== IMPLEMENTING HAWKES PROCESS ===\n")
cat("Model: λ(t) = μ + Σ α·exp(-β(t - t_i)) for t_i < t\n")
cat("  μ = background rate (baseline)\n")
cat("  α = excitation (triggering intensity)\n")
cat("  β = decay rate (how fast influence fades)\n\n")

# We'll implement a simple version using maximum likelihood
# For a more sophisticated version, use the hawkes package

# Function to compute conditional intensity for exponential kernel
hawkes_intensity <- function(times, t, mu, alpha, beta) {
  # λ(t) = μ + Σ α·exp(-β(t - t_i)) for all t_i < t
  past_events <- times[times < t]
  if(length(past_events) == 0) return(mu)

  mu + sum(alpha * exp(-beta * (t - past_events)))
}

# Simple maximum likelihood estimation
# (For full dataset, use specialized packages)

cat("--- Estimating Hawkes parameters ---\n")
cat("Using temporal data for all protests...\n\n")

# For computational speed, work with smaller sample or aggregated data
# Let's use Jakarta as a case study (highest intensity)

jakarta_protests <- protests %>%
  filter(admin1 == "Jakarta") %>%
  arrange(event_date)

jakarta_times <- jakarta_protests$days_since_start

cat("Jakarta dataset:\n")
cat("  Events:", length(jakarta_times), "\n")
cat("  Time span:", max(jakarta_times), "days\n\n")

# Grid search for parameters (simple approach)
# In practice, use optimization algorithms

cat("Performing grid search for Hawkes parameters...\n")

# Define parameter grid
mu_grid <- seq(0.1, 2, length.out = 10)      # Background rate (events/day)
alpha_grid <- seq(0.1, 5, length.out = 10)    # Excitation
beta_grid <- seq(0.01, 0.5, length.out = 10)  # Decay (1/days)

best_loglik <- -Inf
best_params <- list(mu = NA, alpha = NA, beta = NA)

# For speed, sample subset of times
sample_times <- sort(sample(jakarta_times, min(500, length(jakarta_times))))

cat("Searching over parameter grid (this may take a few minutes)...\n")

pb_total <- length(mu_grid) * length(alpha_grid) * length(beta_grid)
pb_current <- 0

for(mu in mu_grid) {
  for(alpha in alpha_grid) {
    for(beta in beta_grid) {
      pb_current <- pb_current + 1

      # Compute log-likelihood
      loglik <- 0

      # For each event, compute log of intensity at that time
      for(i in 1:length(sample_times)) {
        lambda_i <- hawkes_intensity(sample_times, sample_times[i], mu, alpha, beta)
        loglik <- loglik + log(lambda_i)
      }

      # Subtract compensator (integral of intensity)
      T_max <- max(sample_times)
      compensator <- mu * T_max  # Background contribution

      # Excitation contribution (approximate)
      for(ti in sample_times) {
        if(ti < T_max) {
          compensator <- compensator + (alpha/beta) * (1 - exp(-beta * (T_max - ti)))
        }
      }

      loglik <- loglik - compensator

      # Update best
      if(loglik > best_loglik) {
        best_loglik <- loglik
        best_params <- list(mu = mu, alpha = alpha, beta = beta)
      }

      # Progress
      if(pb_current %% 100 == 0) {
        cat(sprintf("  Progress: %d/%d (%.1f%%)\r", pb_current, pb_total,
                   100*pb_current/pb_total))
      }
    }
  }
}

cat("\n\n")

cat("=== HAWKES PARAMETER ESTIMATES (Jakarta) ===\n")
cat(sprintf("  μ (background rate): %.4f events/day\n", best_params$mu))
cat(sprintf("  α (excitation): %.4f\n", best_params$alpha))
cat(sprintf("  β (decay rate): %.4f per day\n", best_params$beta))
cat(sprintf("  1/β (mean time to decay): %.2f days\n", 1/best_params$beta))
cat(sprintf("  Branching ratio (α/β): %.4f\n", best_params$alpha/best_params$beta))
cat(sprintf("  Log-likelihood: %.2f\n\n", best_loglik))

# Interpret branching ratio
branching_ratio <- best_params$alpha / best_params$beta

cat("INTERPRETATION:\n")
if(branching_ratio < 1) {
  cat(sprintf("  Branching ratio = %.3f < 1: SUBCRITICAL process\n", branching_ratio))
  cat("  → System is stable, protests will not cascade indefinitely\n")
  cat(sprintf("  → On average, each protest triggers %.2f additional protests\n",
      branching_ratio))
} else {
  cat(sprintf("  Branching ratio = %.3f ≥ 1: SUPERCRITICAL process\n", branching_ratio))
  cat("  → System is explosive, protests can cascade\n")
  cat("  → WARNING: May indicate model misspecification\n")
}
cat("\n")

# Save results
hawkes_results <- list(
  location = "Jakarta",
  n_events = length(sample_times),
  time_span = max(sample_times),
  mu = best_params$mu,
  alpha = best_params$alpha,
  beta = best_params$beta,
  branching_ratio = branching_ratio,
  loglik = best_loglik
)

saveRDS(hawkes_results, "hawkes_jakarta.rds")
cat("Saved: hawkes_jakarta.rds\n\n")

# ====================================================================
# 3. COMPARE ACROSS PROTEST TYPES
# ====================================================================

cat("=== COMPARING HAWKES PARAMETERS BY PROTEST TYPE ===\n\n")

# Function to estimate Hawkes for a subset
estimate_hawkes_simple <- function(event_times, name) {

  if(length(event_times) < 50) {
    cat(sprintf("  %s: Too few events (%d), skipping\n", name, length(event_times)))
    return(NULL)
  }

  # Sample for speed
  sample_times <- sort(sample(event_times, min(300, length(event_times))))

  # Coarse grid search
  mu_grid <- seq(0.01, 1, length.out = 5)
  alpha_grid <- seq(0.1, 3, length.out = 5)
  beta_grid <- seq(0.01, 0.3, length.out = 5)

  best_loglik <- -Inf
  best <- list(mu = NA, alpha = NA, beta = NA)

  for(mu in mu_grid) {
    for(alpha in alpha_grid) {
      for(beta in beta_grid) {
        loglik <- 0
        for(i in 1:length(sample_times)) {
          lambda_i <- hawkes_intensity(sample_times, sample_times[i], mu, alpha, beta)
          loglik <- loglik + log(max(lambda_i, 1e-10))
        }

        T_max <- max(sample_times)
        compensator <- mu * T_max
        for(ti in sample_times) {
          if(ti < T_max) {
            compensator <- compensator + (alpha/beta) * (1 - exp(-beta * (T_max - ti)))
          }
        }
        loglik <- loglik - compensator

        if(loglik > best_loglik) {
          best_loglik <- loglik
          best <- list(mu = mu, alpha = alpha, beta = beta)
        }
      }
    }
  }

  return(data.frame(
    type = name,
    n = length(event_times),
    mu = best$mu,
    alpha = best$alpha,
    beta = best$beta,
    branching = best$alpha / best$beta,
    loglik = best_loglik
  ))
}

# Compare peaceful vs violent
cat("Estimating for different protest types...\n")

peaceful_times <- protests %>%
  filter(is_peaceful == TRUE) %>%
  pull(days_since_start) %>%
  sort()

violent_times <- protests %>%
  filter(is_violent == TRUE) %>%
  pull(days_since_start) %>%
  sort()

fatal_times <- protests %>%
  filter(has_fatalities == TRUE) %>%
  pull(days_since_start) %>%
  sort()

results_list <- list()

cat("\nPeaceful protests:\n")
results_list[[1]] <- estimate_hawkes_simple(peaceful_times, "Peaceful")

cat("\nViolent protests:\n")
results_list[[2]] <- estimate_hawkes_simple(violent_times, "Violent")

cat("\nProtest with fatalities:\n")
results_list[[3]] <- estimate_hawkes_simple(fatal_times, "With Fatalities")

# Combine results
comparison_df <- bind_rows(results_list)

cat("\n=== HAWKES COMPARISON BY PROTEST TYPE ===\n")
print(comparison_df)
cat("\n")

# Save
write.csv(comparison_df, "hawkes_by_type.csv", row.names = FALSE)
cat("Saved: hawkes_by_type.csv\n\n")

# ====================================================================
# 4. VISUALIZATION
# ====================================================================

cat("=== Creating visualizations ===\n")

# Plot 1: Branching ratios
png("plots/12_hawkes_branching_ratios.png", width = 1000, height = 600, res = 120)

ggplot(comparison_df, aes(x = type, y = branching, fill = type)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
  geom_text(aes(label = sprintf("%.3f", branching)), vjust = -0.5) +
  scale_fill_manual(values = c("Peaceful" = "#2ca02c",
                                "Violent" = "#d62728",
                                "With Fatalities" = "#ff7f0e")) +
  labs(
    title = "Hawkes Process Branching Ratios by Protest Type",
    subtitle = "Branching ratio = α/β (avg number of protests triggered per event)\nRatio < 1: stable; Ratio ≥ 1: explosive",
    x = "Protest Type",
    y = "Branching Ratio (α/β)",
    fill = "Type"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 12))

dev.off()
cat("Saved: plots/12_hawkes_branching_ratios.png\n")

# Plot 2: Decay rates
png("plots/13_hawkes_decay_rates.png", width = 1000, height = 600, res = 120)

comparison_df$mean_decay_days <- 1 / comparison_df$beta

ggplot(comparison_df, aes(x = type, y = mean_decay_days, fill = type)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f days", mean_decay_days)), vjust = -0.5) +
  scale_fill_manual(values = c("Peaceful" = "#2ca02c",
                                "Violent" = "#d62728",
                                "With Fatalities" = "#ff7f0e")) +
  labs(
    title = "Temporal Decay of Protest Influence",
    subtitle = "Mean time for influence to fade (1/β)",
    x = "Protest Type",
    y = "Mean Decay Time (days)",
    fill = "Type"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 12))

dev.off()
cat("Saved: plots/13_hawkes_decay_rates.png\n\n")

# ====================================================================
# 5. SUMMARY
# ====================================================================

cat("=== SUMMARY OF HAWKES ANALYSIS ===\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. TEMPORAL CONTAGION:\n")
for(i in 1:nrow(comparison_df)) {
  row <- comparison_df[i, ]
  cat(sprintf("   %s protests:\n", row$type))
  cat(sprintf("     - Background rate: %.3f events/day\n", row$mu))
  cat(sprintf("     - Excitation: α = %.3f\n", row$alpha))
  cat(sprintf("     - Branching ratio: %.3f ", row$branching))

  if(row$branching < 1) {
    cat("(stable)\n")
  } else {
    cat("(explosive!)\n")
  }

  cat(sprintf("     - Influence decays in ~%.1f days\n", 1/row$beta))
  cat("\n")
}

cat("\n2. INTERPRETATION:\n")

peaceful_branch <- comparison_df %>% filter(type == "Peaceful") %>% pull(branching)
violent_branch <- comparison_df %>% filter(type == "Violent") %>% pull(branching)

if(!is.na(peaceful_branch) && !is.na(violent_branch)) {
  if(peaceful_branch > violent_branch) {
    cat("   → Peaceful protests show STRONGER self-excitation\n")
    cat("     Each peaceful protest triggers more subsequent protests\n")
  } else {
    cat("   → Violent protests show STRONGER self-excitation\n")
    cat("     Violence breeds more violence\n")
  }
}

cat("\n3. LIMITATIONS:\n")
cat("   - Simple temporal model (no spatial component yet)\n")
cat("   - Coarse parameter estimation (use specialized packages for refinement)\n")
cat("   - Need to account for spatial heterogeneity\n")
cat("   - Should test for stationarity\n\n")

cat("NEXT STEPS:\n")
cat("  1. Fit spatial-temporal ETAS models\n")
cat("  2. Include spatial decay functions\n")
cat("  3. Test mark-dependent triggering (e.g., violent → peaceful)\n")
cat("  4. Add environmental covariates\n\n")

cat("=== COMPLETE ===\n")
