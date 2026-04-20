# Quick test to verify C++ compilation with new parameters
library(Rcpp)

cat("Testing C++ compilation...\n")
sourceCpp("hawkes_likelihood.cpp")
cat("✓ C++ compiled successfully\n\n")

cat("Testing C++ function call with year FX + new marks...\n")

# Minimal test data
n <- 10
times <- 1:n
log_pop <- rep(13, n)
poverty_decimal <- rep(0.1, n)
year <- rep(2019L, n)  # All in 2019
is_riot <- rep(0L, n)
has_fatalities <- rep(0L, n)
is_student <- rep(0L, n)
is_labor <- rep(0L, n)

# District-year data
n_dy <- 5
district_year_exposure <- rep(365, n_dy)
district_year_log_pop <- rep(13, n_dy)
district_year_poverty <- rep(0.1, n_dy)
district_year_year <- rep(2019L, n_dy)

# Test parameters (18 total: 3 background + 9 year effects + 6 triggering)
params <- c(
  beta_0_bg = -10,
  gamma_raw = 0,
  delta_raw = 0,
  # Year effects
  beta_2016_raw = 0, beta_2017_raw = 0, beta_2018_raw = 0,
  beta_2019_raw = 0, beta_2020_raw = 0, beta_2021_raw = 0,
  beta_2022_raw = 0, beta_2023_raw = 0, beta_2024_raw = 0,
  # Triggering
  beta_0_trig_raw = 0,
  beta_riot_raw = 0, beta_fatal_raw = 0,
  beta_student_raw = 0, beta_labor_raw = 0,
  decay_raw = 0
)

# Call C++ function
neg_ll <- hawkes_negloglik_cpp(
  times = times,
  log_pop = log_pop,
  poverty_decimal = poverty_decimal,
  year = year,
  is_riot = is_riot,
  has_fatalities = has_fatalities,
  is_student = is_student,
  is_labor = is_labor,
  district_year_exposure = district_year_exposure,
  district_year_log_pop = district_year_log_pop,
  district_year_poverty = district_year_poverty,
  district_year_year = district_year_year,
  params = params,
  temporal_cutoff = 730
)

cat(sprintf("✓ C++ function executed successfully\n"))
cat(sprintf("  Negative log-likelihood: %.2f\n", neg_ll))
cat(sprintf("  Log-likelihood: %.2f\n\n", -neg_ll))

cat("All tests passed! Ready to run full optimization.\n")
