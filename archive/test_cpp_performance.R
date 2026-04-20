# ============================================================================
# TEST C++ PERFORMANCE ON DIFFERENT SAMPLE SIZES
# ============================================================================

library(Rcpp)
library(dplyr)

cat("\n=== C++ PERFORMANCE TESTING ===\n\n")

# Compile C++ code
cat("Compiling C++ code...\n")
sourceCpp("hawkes_likelihood.cpp")
cat("✓ Compiled\n\n")

# Load full data
protests <- readRDS("protests_with_poverty.rds")

# Test parameters
test_params <- c(
  beta_0_bg = -10,
  gamma_raw = 0,
  delta_raw = 0,
  beta_2016 = 0, beta_2017 = 0, beta_2018 = 0,
  beta_2019 = 0, beta_2020 = 0, beta_2021 = 0,
  beta_2022 = 0, beta_2023 = 0, beta_2024 = 0,
  beta_0_trig = log(0.1),
  beta_violence = 0,
  beta_state = 0,
  decay = 0.2
)

# Test on different sample sizes
sample_sizes <- c(500, 1000, 2000, 5000)

results <- data.frame()

for(n in sample_sizes) {
  cat(sprintf("Testing with %d events...\n", n))

  set.seed(42)
  idx <- sample(nrow(protests), min(n, nrow(protests)))
  test_data <- protests[idx, ] %>% arrange(event_date)

  # Prepare data
  times <- as.numeric(difftime(test_data$event_date, min(test_data$event_date), units = "days"))
  log_pop <- test_data$log_pop
  poverty_decimal <- test_data$poverty_decimal
  year <- test_data$year
  violent <- as.integer(test_data$is_violent)
  state_interv <- as.integer(test_data$state_intervention)

  # Calculate district-year data
  district_years <- test_data %>%
    distinct(district, year, log_pop, poverty_decimal) %>%
    mutate(exposure = 365.25)

  # Time the computation
  start <- Sys.time()

  negloglik <- hawkes_negloglik_cpp(
    times = times,
    log_pop = log_pop,
    poverty_decimal = poverty_decimal,
    year = year,
    violent = violent,
    state_interv = state_interv,
    district_year_exposure = district_years$exposure,
    district_year_log_pop = district_years$log_pop,
    district_year_poverty = district_years$poverty_decimal,
    district_year_year = district_years$year,
    params = test_params,
    temporal_cutoff = 730
  )

  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  cat(sprintf("  %d events: %.3f seconds (%.0f ms)\n", n, elapsed, elapsed * 1000))

  results <- rbind(results, data.frame(
    n_events = n,
    time_secs = elapsed,
    time_ms = elapsed * 1000,
    negloglik = negloglik
  ))
}

cat("\n=== PERFORMANCE SUMMARY ===\n")
print(results)

# Extrapolate to full dataset
cat("\n=== EXTRAPOLATION TO FULL DATASET (15,914 events) ===\n")

# Assume O(n^2) complexity with temporal cutoff reducing it somewhat
# Use largest sample to estimate
if(nrow(results) > 0) {
  largest <- results[nrow(results), ]

  # Conservative estimate: O(n^2) scaling
  ratio_squared <- (15914 / largest$n_events)^2
  est_time_secs <- largest$time_secs * ratio_squared
  est_time_mins <- est_time_secs / 60
  est_time_hours <- est_time_mins / 60

  cat(sprintf("Based on %d events taking %.2f seconds:\n", largest$n_events, largest$time_secs))
  cat(sprintf("  O(n²) estimate: %.0f seconds = %.1f minutes = %.2f hours per evaluation\n",
              est_time_secs, est_time_mins, est_time_hours))

  # Estimate full optimization (assume ~100 evaluations × 5 starts)
  total_evals <- 100 * 5
  total_time_hours <- est_time_hours * total_evals
  total_time_days <- total_time_hours / 24

  cat(sprintf("\nFull optimization (~%d evaluations):\n", total_evals))
  cat(sprintf("  Estimated: %.1f hours = %.2f days\n", total_time_hours, total_time_days))

  if(total_time_hours > 48) {
    cat("\n⚠ WARNING: Full dataset optimization would take too long!\n")
    cat("RECOMMENDATIONS:\n")
    cat("1. Use stratified subsample (~5,000 events)\n")
    cat("2. Parallelize C++ with OpenMP\n")
    cat("3. Use spatial/temporal thinning\n")
  } else {
    cat("\n✓ Full dataset optimization appears feasible\n")
  }
}

cat("\n")
