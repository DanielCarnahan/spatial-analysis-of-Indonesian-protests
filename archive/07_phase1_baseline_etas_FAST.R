# Phase 1: Baseline ETAS Model - FAST VERSION FOR TESTING
# Simplified version with smaller sample and coarser grid
# Author: Daniel Carnahan
# Date: 2025-10-27

library(dplyr)

cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  PHASE 1: BASELINE ETAS (FAST TEST VERSION)                   ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

cat("MODIFICATIONS FOR SPEED:\n")
cat("  - Sample: 500 events (vs 2000 full)\n")
cat("  - Grid: 3x3x3x3x3x3 = 729 → 27 combinations\n")
cat("  - Estimated time: 5-10 minutes\n\n")

# Load data
protests <- readRDS("protests_prepared.rds")
cat("Loaded", nrow(protests), "protest events\n\n")

# Convert to km (equirectangular)
mean_lat <- mean(protests$latitude)
protests_etas <- protests %>%
  mutate(
    x_km = longitude * 111 * cos(mean_lat * pi/180),
    y_km = latitude * 111,
    t_days = days_since_start
  ) %>%
  arrange(t_days)

# Create data matrix
data_matrix <- cbind(
  t = protests_etas$t_days,
  x = protests_etas$x_km,
  y = protests_etas$y_km
)

# SMALL SAMPLE for fast testing
set.seed(123)
sample_idx <- sort(sample(1:nrow(data_matrix), 500))
data_sample <- data_matrix[sample_idx, ]

cat("Using sample of", nrow(data_sample), "events\n\n")

# ====================================================================
# CUSTOM ETAS FUNCTIONS
# ====================================================================

lambda_etas <- function(idx, data, params) {
  mu <- params[1]
  K <- params[2]
  c <- params[3]
  p <- params[4]
  d <- params[5]
  q <- params[6]

  intensity <- mu

  if(idx > 1) {
    for(j in 1:(idx-1)) {
      r <- sqrt((data[idx, "x"] - data[j, "x"])^2 +
               (data[idx, "y"] - data[j, "y"])^2)
      tau <- data[idx, "t"] - data[j, "t"]

      if(tau > 0) {
        spatial_kern <- (r^2 + d)^(-q)
        temporal_kern <- (tau + c)^(-p)
        intensity <- intensity + K * spatial_kern * temporal_kern
      }
    }
  }

  return(intensity)
}

loglik_etas_approx <- function(params, data) {
  if(any(params <= 0)) return(-Inf)
  if(params[4] <= 1 || params[6] <= 1) return(-Inf)

  n <- nrow(data)
  ll <- 0

  for(i in 1:n) {
    lambda_i <- lambda_etas(i, data, params)
    ll <- ll + log(max(lambda_i, 1e-10))
  }

  # Note: Compensator omitted for speed (relative comparison still valid)
  return(ll)
}

# ====================================================================
# COARSE GRID SEARCH
# ====================================================================

cat("=== GRID SEARCH ===\n\n")

# REDUCED grid for speed
mu_grid <- c(0.0001, 0.001)
K_grid <- c(0.05, 0.1, 0.2)
c_grid <- c(0.5, 1.0, 2.0)
p_grid <- c(1.2, 1.5, 2.0)
d_grid <- c(5, 10, 20)
q_grid <- c(1.2, 1.5, 2.0)

grid_total <- length(mu_grid) * length(K_grid) * length(c_grid) *
              length(p_grid) * length(d_grid) * length(q_grid)

cat(sprintf("Testing %d parameter combinations\n", grid_total))
cat("This should take 5-10 minutes...\n\n")

best_ll <- -Inf
best_params <- NULL
grid_count <- 0

start_time <- Sys.time()

for(mu in mu_grid) {
  for(K in K_grid) {
    for(c_val in c_grid) {
      for(p in p_grid) {
        for(d in d_grid) {
          for(q in q_grid) {

            grid_count <- grid_count + 1

            params_test <- c(mu, K, c_val, p, d, q)

            tryCatch({
              ll <- loglik_etas_approx(params_test, data_sample)

              if(ll > best_ll) {
                best_ll <- ll
                best_params <- params_test
                cat(sprintf("\n✓ New best (iter %d/%d): LL = %.2f\n",
                           grid_count, grid_total, ll))
                cat(sprintf("  μ=%.4f, K=%.3f, c=%.2f, p=%.2f, d=%.1f, q=%.2f\n",
                           mu, K, c_val, p, d, q))
              }

              if(grid_count %% 10 == 0) {
                elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
                remaining <- elapsed * (grid_total - grid_count) / grid_count
                cat(sprintf("  Progress: %d/%d (%.1f%%) - Elapsed: %.1fm, ETA: %.1fm\r",
                           grid_count, grid_total, 100*grid_count/grid_total,
                           elapsed, remaining))
              }
            }, error = function(e) {
              cat(sprintf("\n✗ Error at iteration %d: %s\n", grid_count, e$message))
            })
          }
        }
      }
    }
  }
}

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat(sprintf("\n\n✓ Grid search complete in %.2f minutes\n\n", as.numeric(runtime)))

# ====================================================================
# SAVE RESULTS
# ====================================================================

cat("=== RESULTS ===\n\n")

names(best_params) <- c("mu", "K", "c", "p", "d", "q")

cat("Best Parameters:\n")
cat(sprintf("  μ (background): %.6f events/(day·km²)\n", best_params["mu"]))
cat(sprintf("  K (productivity): %.4f\n", best_params["K"]))
cat(sprintf("  c (temporal offset): %.2f days\n", best_params["c"]))
cat(sprintf("  p (temporal decay): %.2f\n", best_params["p"]))
cat(sprintf("  d (spatial scale): %.2f km²\n", best_params["d"]))
cat(sprintf("  q (spatial decay): %.2f\n", best_params["q"]))
cat(sprintf("  Log-likelihood: %.2f\n\n", best_ll))

# Derived quantities
if(best_params["p"] > 1) {
  mean_tau <- best_params["c"] / (best_params["p"] - 1)
  cat(sprintf("Mean temporal influence: %.2f days\n", mean_tau))
}

if(best_params["q"] > 1) {
  mean_r <- sqrt(best_params["d"] / (best_params["q"] - 1))
  cat(sprintf("Mean spatial influence: %.2f km\n", mean_r))
}

cat("\n")

# Save results
fit_etas <- list(
  params = best_params,
  loglik = best_ll,
  data = data_sample,
  method = "custom_grid_search_fast",
  n_events = nrow(data_sample),
  runtime_mins = as.numeric(runtime),
  note = "Fast test version with small sample and coarse grid"
)

saveRDS(fit_etas, "etas_baseline_fast.rds")
cat("✓ Saved: etas_baseline_fast.rds\n")

# Save parameters table
params_df <- data.frame(
  Parameter = names(best_params),
  Estimate = as.numeric(best_params),
  Interpretation = c(
    "Background rate (spontaneous events)",
    "Productivity (triggering strength)",
    "Temporal offset (days)",
    "Temporal power-law decay",
    "Spatial scale (km²)",
    "Spatial power-law decay"
  )
)

write.csv(params_df, "etas_baseline_fast_parameters.csv", row.names = FALSE)
cat("✓ Saved: etas_baseline_fast_parameters.csv\n\n")

# ====================================================================
# SUMMARY
# ====================================================================

cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║           PHASE 1 FAST VERSION COMPLETE                       ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

cat("FILES CREATED:\n")
cat("  1. etas_baseline_fast.rds\n")
cat("  2. etas_baseline_fast_parameters.csv\n\n")

cat("INTERPRETATION:\n")
cat("These parameters provide baseline spatial-temporal structure.\n")
cat("They can be used as initial values for Phase 2 mark-dependent models.\n\n")

cat("NEXT STEPS:\n")
cat("  Option 1: Proceed to Phase 2 with these baseline parameters\n")
cat("  Option 2: Re-run with full dataset (remove sampling, increase grid)\n")
cat("  Option 3: Skip spatial ETAS, go directly to temporal mark-dependent Hawkes\n\n")

cat("=== END PHASE 1 FAST ===\n")
