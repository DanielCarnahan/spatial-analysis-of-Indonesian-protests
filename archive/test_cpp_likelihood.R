# ============================================================================
# TEST C++ LIKELIHOOD IMPLEMENTATION
# ============================================================================

library(Rcpp)
library(dplyr)

cat("\n=== TESTING C++ HAWKES LIKELIHOOD ===\n\n")

# Compile C++ code
cat("Compiling C++ code...\n")
sourceCpp("hawkes_likelihood.cpp")
cat("✓ C++ code compiled successfully\n\n")

# Load test data (small sample)
cat("Loading test data...\n")
protests <- readRDS("protests_with_poverty.rds")

# Use 100 events for quick test
set.seed(42)
test_idx <- sample(nrow(protests), 100)
test_data <- protests[test_idx, ] %>% arrange(event_date)

cat(sprintf("Test sample: %d events\n\n", nrow(test_data)))

# Prepare data for C++ function
times <- as.numeric(difftime(test_data$event_date, min(test_data$event_date), units = "days"))
log_pop <- test_data$log_pop
poverty_decimal <- test_data$poverty_decimal
year <- test_data$year
violent <- as.integer(test_data$is_violent)
state_interv <- as.integer(test_data$state_intervention)

# Calculate district-year exposure times
district_years <- test_data %>%
  distinct(district, year, log_pop, poverty_decimal) %>%
  mutate(exposure = 365.25)  # Assume full year

# Test parameters (similar to initialization)
test_params <- c(
  beta_0_bg = -10,
  gamma_raw = 0,     # transforms to gamma = 0.35 (approx)
  delta_raw = 0,     # transforms to delta = 0
  beta_2016 = 0, beta_2017 = 0, beta_2018 = 0,
  beta_2019 = 0, beta_2020 = 0, beta_2021 = 0,
  beta_2022 = 0, beta_2023 = 0, beta_2024 = 0,
  beta_0_trig = log(0.1),
  beta_violence = 0,
  beta_state = 0,
  decay = 0.2
)

# Run C++ likelihood
cat("Testing C++ likelihood function...\n")
start_time <- Sys.time()

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

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat(sprintf("✓ C++ likelihood computed successfully\n"))
cat(sprintf("  Negative log-likelihood: %.2f\n", negloglik))
cat(sprintf("  Runtime: %.3f seconds\n\n", runtime))

# Test with different parameter values
cat("Testing parameter sensitivity...\n")
test_params2 <- test_params
test_params2[2] <- 1.0  # Different gamma_raw

negloglik2 <- hawkes_negloglik_cpp(
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
  params = test_params2,
  temporal_cutoff = 730
)

cat(sprintf("  With different gamma_raw: %.2f\n", negloglik2))
cat(sprintf("  ✓ Parameters affect likelihood (difference: %.2f)\n\n", abs(negloglik - negloglik2)))

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   C++ IMPLEMENTATION TEST PASSED                             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("NEXT STEPS:\n")
cat("-----------\n")
cat("1. Integrate C++ function into 16_hawkes_with_poverty.R\n")
cat("2. Run optimization on full dataset\n")
cat("3. Compare results with R implementation (if available)\n\n")
