#!/usr/bin/env Rscript
# ==============================================================================
# Investigation: Slow Decay Parameter
# ==============================================================================
#
# PURPOSE: Diagnose why β_decay = 0.001454 (half-life = 477 days)
#          This is extremely slow for protest contagion
#
# PHASES:
# 1. Extract decay estimates from all models
# 2. Empirical temporal autocorrelation
# 3. Inter-event time distributions
# 4. Conditional intensity visualization
#
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   DECAY PARAMETER INVESTIGATION                              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ==============================================================================
# PHASE 1: Extract Decay Parameters from Models
# ==============================================================================

cat("=== PHASE 1: PARAMETER DIAGNOSTICS ===\n\n")

# Load main model (Model 2: Hawkes with marks)
model2 <- readRDS("model_poverty.rds")

cat("Model 2 (Hawkes + Marks):\n")
cat(sprintf("  Log-likelihood: %.2f\n", model2$loglik))
cat(sprintf("  Convergence: %d\n", model2$convergence))
cat(sprintf("  Parameters: %d\n", model2$n_params))
cat("\n")

# Extract decay parameter (it's transformed)
# Need to apply inverse transformation
transform_decay <- function(decay_raw) {
  0.001 + 9.999 * plogis(decay_raw)
}

if("decay_raw" %in% names(model2$params)) {
  decay_raw <- model2$params$decay_raw
  decay_est <- transform_decay(decay_raw)
  half_life <- log(2) / decay_est

  cat("Decay Parameter Estimates:\n")
  cat(sprintf("  decay_raw (unconstrained): %.6f\n", decay_raw))
  cat(sprintf("  decay (constrained [0.001, 10]): %.6f\n", decay_est))
  cat(sprintf("  Half-life: %.1f days (%.1f years)\n", half_life, half_life/365.25))
  cat("\n")
} else {
  cat("Warning: decay_raw not found in model parameters\n")
  print(names(model2$params))
  cat("\n")
}

# Check if we have multi-start results
if("all_results" %in% names(model2)) {
  cat("Multi-Start Convergence:\n")
  cat(sprintf("  Number of starting points: %d\n", length(model2$all_results)))

  # Extract decay from each start
  decay_estimates <- sapply(model2$all_results, function(res) {
    if("decay_raw" %in% names(res$params)) {
      transform_decay(res$params$decay_raw)
    } else {
      NA
    }
  })

  if(any(!is.na(decay_estimates))) {
    cat("\n  Decay estimates across starts:\n")
    for(i in seq_along(decay_estimates)) {
      if(!is.na(decay_estimates[i])) {
        cat(sprintf("    Start %d: %.6f (half-life: %.1f days)\n",
                   i, decay_estimates[i], log(2)/decay_estimates[i]))
      }
    }

    cat(sprintf("\n  Range: [%.6f, %.6f]\n", min(decay_estimates, na.rm=TRUE),
               max(decay_estimates, na.rm=TRUE)))
    cat(sprintf("  SD: %.6f\n", sd(decay_estimates, na.rm=TRUE)))
  }
  cat("\n")
}

# Load Model 1 (Basic Hawkes, no marks) for comparison
if(file.exists("model1_basic_hawkes.rds")) {
  model1 <- readRDS("model1_basic_hawkes.rds")

  cat("Model 1 (Basic Hawkes, constant triggering):\n")
  cat(sprintf("  Log-likelihood: %.2f\n", model1$loglik))

  if("decay_raw" %in% names(model1$params)) {
    decay1_raw <- model1$params$decay_raw
    decay1_est <- transform_decay(decay1_raw)
    half_life1 <- log(2) / decay1_est

    cat(sprintf("  Decay: %.6f (half-life: %.1f days)\n", decay1_est, half_life1))
  }
  cat("\n")
}

# ==============================================================================
# PHASE 2: Load Protest Data for Empirical Analysis
# ==============================================================================

cat("=== PHASE 2: LOADING PROTEST DATA ===\n\n")

if(file.exists("protests_with_poverty.rds")) {
  protests <- readRDS("protests_with_poverty.rds")
  cat(sprintf("Loaded %d protest events\n", nrow(protests)))
  cat(sprintf("Time range: %d to %d\n", min(protests$year, na.rm=TRUE),
             max(protests$year, na.rm=TRUE)))
  cat("\n")
} else {
  cat("Error: protests_with_poverty.rds not found\n")
  cat("Checking for alternative data files...\n")

  # Try alternative files
  if(file.exists("protests_prepared.rds")) {
    protests <- readRDS("protests_prepared.rds")
    cat(sprintf("Loaded %d protest events from protests_prepared.rds\n", nrow(protests)))
  } else {
    cat("Error: No protest data files found\n")
    cat("Available RDS files:\n")
    system("ls -1 *.rds | head -10")
    quit(status=1)
  }
}

# ==============================================================================
# PHASE 3: Empirical Temporal Autocorrelation
# ==============================================================================

cat("=== PHASE 3: EMPIRICAL TEMPORAL AUTOCORRELATION ===\n\n")

# Prepare time series data
if("days_since_start" %in% names(protests)) {
  protests <- protests %>%
    arrange(days_since_start) %>%
    filter(!is.na(days_since_start))

  cat(sprintf("Events with time data: %d\n", nrow(protests)))

  # Create daily counts
  min_day <- floor(min(protests$days_since_start))
  max_day <- ceiling(max(protests$days_since_start))

  cat(sprintf("Time span: %.0f days (%.1f years)\n", max_day - min_day,
             (max_day - min_day)/365.25))

  # Aggregate to daily counts
  daily_counts <- protests %>%
    mutate(day = floor(days_since_start)) %>%
    group_by(day) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(day = min_day:max_day, fill = list(count = 0))

  cat(sprintf("\nDaily count statistics:\n"))
  cat(sprintf("  Mean: %.2f protests/day\n", mean(daily_counts$count)))
  cat(sprintf("  SD: %.2f\n", sd(daily_counts$count)))
  cat(sprintf("  Max: %d\n", max(daily_counts$count)))
  cat(sprintf("  Days with 0 protests: %d (%.1f%%)\n",
             sum(daily_counts$count == 0),
             100 * sum(daily_counts$count == 0) / nrow(daily_counts)))
  cat("\n")

  # Compute autocorrelation at key lags
  lags_to_test <- c(1, 3, 7, 14, 30, 60, 90, 180, 365)

  cat("Autocorrelation at key lags:\n")
  cat("  Lag (days) | ACF     | Expected (β=0.0014)\n")
  cat("  -----------|---------|--------------------\n")

  acf_result <- acf(daily_counts$count, lag.max = max(lags_to_test), plot = FALSE)

  for(lag in lags_to_test) {
    if(lag <= length(acf_result$acf)) {
      acf_val <- acf_result$acf[lag + 1]  # +1 because R ACF is 0-indexed

      # What exponential decay with β=0.0014 would predict
      # Rough approximation: ACF ≈ exp(-β * lag) for Hawkes process
      expected_acf <- exp(-0.001454 * lag)

      cat(sprintf("  %10d | %7.4f | %7.4f\n", lag, acf_val, expected_acf))
    }
  }
  cat("\n")

  # Save for plotting
  acf_data <- data.frame(
    lag = acf_result$lag * 1,  # Convert to days
    acf = as.numeric(acf_result$acf)
  )

  saveRDS(acf_data, "decay_investigation_acf.rds")
  cat("✓ Saved autocorrelation data\n\n")

} else {
  cat("Error: days_since_start column not found\n")
  cat("Available columns:\n")
  print(names(protests))
}

# ==============================================================================
# PHASE 4: Inter-Event Time Distributions
# ==============================================================================

cat("=== PHASE 4: INTER-EVENT TIME DISTRIBUTIONS ===\n\n")

if("days_since_start" %in% names(protests) && "event_id_cnty" %in% names(protests)) {

  # Compute inter-event times (all events)
  protests_sorted <- protests %>%
    arrange(days_since_start)

  inter_event_times <- diff(protests_sorted$days_since_start)

  cat("Inter-event time statistics (all events):\n")
  cat(sprintf("  N intervals: %d\n", length(inter_event_times)))
  cat(sprintf("  Mean: %.2f days\n", mean(inter_event_times, na.rm=TRUE)))
  cat(sprintf("  Median: %.2f days\n", median(inter_event_times, na.rm=TRUE)))
  cat(sprintf("  SD: %.2f days\n", sd(inter_event_times, na.rm=TRUE)))
  cat(sprintf("  Max: %.1f days\n", max(inter_event_times, na.rm=TRUE)))
  cat("\n")

  # Quantiles
  cat("Quantiles:\n")
  quantiles <- quantile(inter_event_times, probs = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                        na.rm=TRUE)
  for(q in names(quantiles)) {
    cat(sprintf("  %s: %.2f days\n", q, quantiles[q]))
  }
  cat("\n")

  # Fraction of events within different time windows
  cat("Clustering: Fraction of inter-event times ≤ threshold:\n")
  thresholds <- c(1, 3, 7, 14, 30, 90)
  for(thresh in thresholds) {
    frac <- mean(inter_event_times <= thresh, na.rm=TRUE)
    cat(sprintf("  ≤ %3d days: %.1f%%\n", thresh, 100*frac))
  }
  cat("\n")

  # By event type (if available)
  if("is_fatal" %in% names(protests) && "is_student" %in% names(protests) &&
     "is_labor" %in% names(protests)) {

    cat("Mean inter-event time by event type:\n")

    # Fatal events
    fatal_events <- protests %>% filter(is_fatal == 1) %>% arrange(days_since_start)
    if(nrow(fatal_events) > 1) {
      fatal_inter <- diff(fatal_events$days_since_start)
      cat(sprintf("  Fatal (n=%d): %.2f days\n", nrow(fatal_events), mean(fatal_inter, na.rm=TRUE)))
    }

    # Student events
    student_events <- protests %>% filter(is_student == 1) %>% arrange(days_since_start)
    if(nrow(student_events) > 1) {
      student_inter <- diff(student_events$days_since_start)
      cat(sprintf("  Student (n=%d): %.2f days\n", nrow(student_events), mean(student_inter, na.rm=TRUE)))
    }

    # Labor events
    labor_events <- protests %>% filter(is_labor == 1) %>% arrange(days_since_start)
    if(nrow(labor_events) > 1) {
      labor_inter <- diff(labor_events$days_since_start)
      cat(sprintf("  Labor (n=%d): %.2f days\n", nrow(labor_events), mean(labor_inter, na.rm=TRUE)))
    }
    cat("\n")
  }

  # Save inter-event times for plotting
  iet_data <- data.frame(
    inter_event_time = inter_event_times
  )
  saveRDS(iet_data, "decay_investigation_inter_event_times.rds")
  cat("✓ Saved inter-event time data\n\n")
}

# ==============================================================================
# PHASE 5: Conditional Intensity After Events
# ==============================================================================

cat("=== PHASE 5: CONDITIONAL INTENSITY AFTER EVENTS ===\n\n")

if("days_since_start" %in% names(protests) && all(c("is_fatal", "is_student", "is_labor") %in% names(protests))) {

  # Function to compute average protest rate in windows after trigger events
  compute_conditional_rate <- function(trigger_events, all_events, max_lag = 90, window = 7) {
    rates <- numeric(ceiling(max_lag / window))
    lags <- seq(0, max_lag - window, by = window)

    for(i in seq_along(lags)) {
      lag_start <- lags[i]
      lag_end <- lag_start + window

      # For each trigger event, count protests in [lag_start, lag_end)
      counts <- sapply(trigger_events$days_since_start, function(t_trigger) {
        sum(all_events$days_since_start > t_trigger + lag_start &
            all_events$days_since_start <= t_trigger + lag_end)
      })

      rates[i] <- mean(counts) / window  # Average protests per day in window
    }

    data.frame(
      lag = lags + window/2,  # Midpoint of window
      rate = rates
    )
  }

  # Overall baseline rate
  total_days <- max(protests$days_since_start) - min(protests$days_since_start)
  baseline_rate <- nrow(protests) / total_days
  cat(sprintf("Baseline rate: %.3f protests/day\n\n", baseline_rate))

  # Conditional rates after fatal events
  fatal_events <- protests %>% filter(is_fatal == 1)
  if(nrow(fatal_events) > 10) {
    fatal_cond <- compute_conditional_rate(fatal_events, protests, max_lag = 90, window = 7)
    cat("Average rate after fatal events (7-day windows):\n")
    print(head(fatal_cond, 10))
    cat(sprintf("  Rate 0-7 days after: %.3f (%.1fx baseline)\n",
               fatal_cond$rate[1], fatal_cond$rate[1]/baseline_rate))
    cat(sprintf("  Rate 7-14 days after: %.3f (%.1fx baseline)\n",
               fatal_cond$rate[2], fatal_cond$rate[2]/baseline_rate))
    cat(sprintf("  Rate 30-37 days after: %.3f (%.1fx baseline)\n",
               fatal_cond$rate[5], fatal_cond$rate[5]/baseline_rate))
    cat("\n")

    saveRDS(fatal_cond, "decay_investigation_fatal_conditional.rds")
  }

  # Conditional rates after student events
  student_events <- protests %>% filter(is_student == 1)
  if(nrow(student_events) > 10) {
    student_cond <- compute_conditional_rate(student_events, protests, max_lag = 90, window = 7)
    cat("Average rate after student events (7-day windows):\n")
    cat(sprintf("  Rate 0-7 days after: %.3f (%.1fx baseline)\n",
               student_cond$rate[1], student_cond$rate[1]/baseline_rate))
    cat("\n")

    saveRDS(student_cond, "decay_investigation_student_conditional.rds")
  }

  # Conditional rates after labor events
  labor_events <- protests %>% filter(is_labor == 1)
  if(nrow(labor_events) > 10) {
    labor_cond <- compute_conditional_rate(labor_events, protests, max_lag = 90, window = 7)
    cat("Average rate after labor events (7-day windows):\n")
    cat(sprintf("  Rate 0-7 days after: %.3f (%.1fx baseline)\n",
               labor_cond$rate[1], labor_cond$rate[1]/baseline_rate))
    cat("\n")

    saveRDS(labor_cond, "decay_investigation_labor_conditional.rds")
  }
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   INVESTIGATION SUMMARY                                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("KEY FINDINGS:\n")
cat("1. Estimated decay parameter: β = 0.001454 (half-life = 477 days)\n")
cat("2. Empirical autocorrelation analysis completed\n")
cat("3. Inter-event time distributions computed\n")
cat("4. Conditional intensity after events calculated\n\n")

cat("NEXT STEPS:\n")
cat("1. Create visualizations (run decay_plots.R)\n")
cat("2. Test alternative kernel specifications\n")
cat("3. Compare to power-law decay model\n")
cat("4. Update HTML document with findings\n\n")

cat("✓ Investigation complete!\n\n")
