################################################################################
#                 ROBUSTNESS CHECKS FOR DISCRETE-TIME HAWKES MODEL
#
#   Three tests to strengthen causal interpretation:
#     1. Placebo test: Lead variables should be ~0
#     2. Permutation test: Null distribution for branching ratio
#     3. Aggregation sensitivity: Row-count vs n_events_collapsed
#
################################################################################

library(tidyverse)
library(MASS)

# Fix MASS masking dplyr::select
select <- dplyr::select

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   ROBUSTNESS CHECKS FOR SELF-EXCITATION                                 ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# 1. LOAD AND PREPARE DATA (same as main script)
# =============================================================================

MAX_LAG <- 7
REGION_LEVEL <- "district"  # Using admin2 (district) instead of admin1 (province)

protests <- readRDS("protests_daily.rds") %>%
  filter(!is.na(poverty_decimal) & !is.na(log_pop) & !is.na(log_cpi))

protests$region <- protests$admin2  # District level (435 regions)

all_dates <- seq.Date(min(protests$date), max(protests$date), by = "day")
all_regions <- unique(protests$region)

# Create panel with BOTH counting methods
event_counts <- protests %>%
  group_by(region, date) %>%
  summarise(
    count_rows = n(),
    count_events = sum(n_events_collapsed),
    .groups = "drop"
  )

panel <- expand.grid(
  region = all_regions,
  date = all_dates,
  stringsAsFactors = FALSE
) %>%
  left_join(event_counts, by = c("region", "date")) %>%
  mutate(
    count_rows = replace_na(count_rows, 0),
    count_events = replace_na(count_events, 0),
    year = year(date),
    month = month(date),
    day_of_week = wday(date),
    day_idx = as.numeric(date - min(date))
  ) %>%
  arrange(region, date)

# Add lags for count_rows (current method)
for (lag in 1:MAX_LAG) {
  panel <- panel %>%
    group_by(region) %>%
    mutate(!!paste0("lag", lag) := lag(count_rows, lag, default = 0)) %>%
    ungroup()
}

# Add LEAD variables for placebo test
for (lead in 1:2) {
  panel <- panel %>%
    group_by(region) %>%
    mutate(!!paste0("lead", lead) := lead(count_rows, lead, default = 0)) %>%
    ungroup()
}

# Filter to valid observations
model_data <- panel %>%
  filter(day_idx >= MAX_LAG & day_idx <= max(day_idx) - 2)  # Valid for both lags and leads

cat(sprintf("Data prepared: %d observations\n\n", nrow(model_data)))

# =============================================================================
# 2. BASELINE MODEL (for reference)
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  BASELINE: Identity-Link Hawkes                             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

m_baseline <- glm(count_rows ~ lag1 + lag2 + lag3 + lag4 + lag5 + lag6 + lag7,
                  data = model_data,
                  family = quasi(link = "identity", variance = "mu"),
                  start = c(0.05, rep(0.01, 7)))

lag_coefs <- coef(m_baseline)[paste0("lag", 1:7)]
observed_br <- sum(lag_coefs)

cat(sprintf("Background rate (μ): %.5f\n", coef(m_baseline)["(Intercept)"]))
cat(sprintf("Branching ratio (observed): %.4f\n\n", observed_br))

cat("Lag coefficients:\n")
print(round(coef(summary(m_baseline)), 5))

# =============================================================================
# 3. ROBUSTNESS CHECK 1: PLACEBO TEST (LEAD VARIABLES)
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  CHECK 1: PLACEBO TEST (Lead Variables)                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("If self-excitation is real, FUTURE protests should NOT predict CURRENT.\n")
cat("Lead coefficients should be ≈ 0.\n\n")

m_placebo <- glm(count_rows ~ lag1 + lag2 + lag3 + lag4 + lag5 + lag6 + lag7 + lead1 + lead2,
                 data = model_data,
                 family = quasi(link = "identity", variance = "mu"),
                 start = c(0.05, rep(0.01, 9)))

placebo_coefs <- coef(summary(m_placebo))

cat("Full model with lags AND leads:\n")
print(round(placebo_coefs, 5))

# Extract lead coefficients
lead1_coef <- placebo_coefs["lead1", 1]
lead1_se <- placebo_coefs["lead1", 2]
lead1_t <- placebo_coefs["lead1", 3]
lead1_p <- placebo_coefs["lead1", 4]

lead2_coef <- placebo_coefs["lead2", 1]
lead2_se <- placebo_coefs["lead2", 2]
lead2_t <- placebo_coefs["lead2", 3]
lead2_p <- placebo_coefs["lead2", 4]

cat("\n─────────────────────────────────────────────────────────────────\n")
cat("PLACEBO TEST RESULTS:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("Lead-1 coefficient: %.5f (SE: %.5f, t: %.2f, p: %.4f)\n",
            lead1_coef, lead1_se, lead1_t, lead1_p))
cat(sprintf("Lead-2 coefficient: %.5f (SE: %.5f, t: %.2f, p: %.4f)\n",
            lead2_coef, lead2_se, lead2_t, lead2_p))

# Compare to lag coefficients
lag_sum <- sum(coef(m_placebo)[paste0("lag", 1:7)])
lead_sum <- lead1_coef + lead2_coef

cat(sprintf("\nSum of lag coefficients (BR): %.4f\n", lag_sum))
cat(sprintf("Sum of lead coefficients: %.4f\n", lead_sum))
cat(sprintf("Ratio (lead/lag): %.2f%%\n", 100 * abs(lead_sum) / lag_sum))

if (lead1_p > 0.05 & lead2_p > 0.05) {
  cat("\n✓ PASS: Lead coefficients not significant (p > 0.05)\n")
  cat("  This supports causal self-excitation interpretation.\n")
} else {
  cat("\n⚠ CONCERN: Lead coefficient(s) significant\n")
  cat("  This suggests possible omitted variable bias or common shocks.\n")
}

# =============================================================================
# 4. ROBUSTNESS CHECK 2: PERMUTATION TEST
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  CHECK 2: PERMUTATION TEST                                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Shuffling dates within each region to break temporal dependence.\n")
cat("Under null (no self-excitation), BR should be ≈ 0.\n\n")

set.seed(42)
n_perm <- 500  # Use 500 for reasonable runtime
br_null <- numeric(n_perm)

cat(sprintf("Running %d permutations...\n", n_perm))

pb <- txtProgressBar(min = 0, max = n_perm, style = 3)

for (i in 1:n_perm) {
  # Shuffle counts within each region (preserves marginal distribution)
  perm_data <- model_data %>%
    group_by(region) %>%
    mutate(count_perm = sample(count_rows)) %>%
    ungroup()

  # Recompute lags on permuted counts
  for (lag in 1:MAX_LAG) {
    perm_data <- perm_data %>%
      group_by(region) %>%
      mutate(!!paste0("lag", lag, "_perm") := lag(count_perm, lag, default = 0)) %>%
      ungroup()
  }

  # Fit model on permuted data
  m_perm <- tryCatch({
    glm(count_perm ~ lag1_perm + lag2_perm + lag3_perm + lag4_perm +
          lag5_perm + lag6_perm + lag7_perm,
        data = perm_data,
        family = quasi(link = "identity", variance = "mu"),
        start = c(0.05, rep(0.01, 7)))
  }, error = function(e) NULL)

  if (!is.null(m_perm)) {
    br_null[i] <- sum(coef(m_perm)[paste0("lag", 1:7, "_perm")])
  } else {
    br_null[i] <- NA
  }

  setTxtProgressBar(pb, i)
}
close(pb)

br_null <- br_null[!is.na(br_null)]

cat("\n\n─────────────────────────────────────────────────────────────────\n")
cat("PERMUTATION TEST RESULTS:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("Observed BR: %.4f\n", observed_br))
cat(sprintf("Null distribution: mean=%.4f, sd=%.4f\n", mean(br_null), sd(br_null)))
cat(sprintf("Null 95th percentile: %.4f\n", quantile(br_null, 0.95)))
cat(sprintf("Null 99th percentile: %.4f\n", quantile(br_null, 0.99)))

p_perm <- mean(br_null >= observed_br)
cat(sprintf("\nPermutation p-value: %.4f (proportion of null >= observed)\n", p_perm))

if (p_perm < 0.01) {
  cat("\n✓ PASS: Observed BR significantly above null distribution (p < 0.01)\n")
  cat("  Self-excitation is not an artifact of temporal structure.\n")
} else {
  cat("\n⚠ CONCERN: Cannot reject null hypothesis\n")
  cat("  Observed BR may be explainable by chance.\n")
}

# =============================================================================
# 5. ROBUSTNESS CHECK 3: AGGREGATION SENSITIVITY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  CHECK 3: AGGREGATION SENSITIVITY                           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Comparing two counting methods:\n")
cat("  A) count_rows: Each row = 1 event (current method)\n")
cat("  B) count_events: Sum of n_events_collapsed\n\n")

# Add lags for count_events
for (lag in 1:MAX_LAG) {
  model_data <- model_data %>%
    group_by(region) %>%
    mutate(!!paste0("lag", lag, "_events") := lag(count_events, lag, default = 0)) %>%
    ungroup()
}

# Model A: count_rows (already computed)
br_rows <- observed_br

# Model B: count_events
m_events <- glm(count_events ~ lag1_events + lag2_events + lag3_events +
                  lag4_events + lag5_events + lag6_events + lag7_events,
                data = model_data,
                family = quasi(link = "identity", variance = "mu"),
                start = c(0.05, rep(0.01, 7)))

br_events <- sum(coef(m_events)[paste0("lag", 1:7, "_events")])

cat("─────────────────────────────────────────────────────────────────\n")
cat("AGGREGATION SENSITIVITY RESULTS:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("Method A (row counts):    BR = %.4f\n", br_rows))
cat(sprintf("Method B (event counts):  BR = %.4f\n", br_events))
cat(sprintf("Difference: %.4f (%.1f%%)\n", br_events - br_rows, 100 * (br_events - br_rows) / br_rows))

cat("\nMethod B coefficients:\n")
print(round(coef(summary(m_events)), 5))

if (abs(br_events - br_rows) / br_rows < 0.20) {
  cat("\n✓ PASS: BR estimates within 20% of each other\n")
  cat("  Results robust to counting method.\n")
} else {
  cat("\n⚠ CONCERN: BR estimates differ substantially\n")
  cat("  Counting method affects results.\n")
}

# =============================================================================
# 6. SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  ROBUSTNESS CHECK SUMMARY                                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("┌────────────────────────┬─────────┬────────────────────────────┐\n")
cat("│ Check                  │ Result  │ Interpretation             │\n")
cat("├────────────────────────┼─────────┼────────────────────────────┤\n")

# Placebo
placebo_pass <- lead1_p > 0.05 & lead2_p > 0.05
cat(sprintf("│ 1. Placebo (leads)     │ %s    │ %-26s │\n",
            ifelse(placebo_pass, "PASS", "FAIL"),
            ifelse(placebo_pass, "Supports causal story", "Omitted variable concern")))

# Permutation
perm_pass <- p_perm < 0.01
cat(sprintf("│ 2. Permutation test    │ %s    │ %-26s │\n",
            ifelse(perm_pass, "PASS", "FAIL"),
            ifelse(perm_pass, "Not artifact of chance", "May be artifact")))

# Aggregation
agg_pass <- abs(br_events - br_rows) / br_rows < 0.20
cat(sprintf("│ 3. Aggregation         │ %s    │ %-26s │\n",
            ifelse(agg_pass, "PASS", "FAIL"),
            ifelse(agg_pass, "Robust to counting method", "Sensitive to counting")))

cat("└────────────────────────┴─────────┴────────────────────────────┘\n")

n_pass <- sum(c(placebo_pass, perm_pass, agg_pass))
cat(sprintf("\nOverall: %d/3 checks passed\n", n_pass))

if (n_pass == 3) {
  cat("\n✓ STRONG SUPPORT for causal self-excitation interpretation.\n")
} else if (n_pass >= 2) {
  cat("\n○ MODERATE SUPPORT - some concerns remain.\n")
} else {
  cat("\n⚠ WEAK SUPPORT - interpretation should be cautious.\n")
}

# =============================================================================
# 7. SAVE RESULTS
# =============================================================================

robustness_results <- list(
  baseline = list(
    branching_ratio = observed_br,
    coefficients = coef(summary(m_baseline))
  ),
  placebo = list(
    lead1 = list(coef = lead1_coef, se = lead1_se, p = lead1_p),
    lead2 = list(coef = lead2_coef, se = lead2_se, p = lead2_p),
    pass = placebo_pass
  ),
  permutation = list(
    observed_br = observed_br,
    null_mean = mean(br_null),
    null_sd = sd(br_null),
    null_95 = quantile(br_null, 0.95),
    null_99 = quantile(br_null, 0.99),
    p_value = p_perm,
    pass = perm_pass
  ),
  aggregation = list(
    br_rows = br_rows,
    br_events = br_events,
    difference_pct = 100 * (br_events - br_rows) / br_rows,
    pass = agg_pass
  ),
  overall_pass = n_pass
)

saveRDS(robustness_results, "robustness_results.rds")
cat("\nSaved: robustness_results.rds\n")

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  ROBUSTNESS CHECKS COMPLETE                                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
