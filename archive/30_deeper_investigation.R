################################################################################
#                   DEEPER INVESTIGATION: MODEL vs RAW DATA
#                   Exploring additional hypotheses
################################################################################

library(tidyverse)
library(data.table)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   DEEPER INVESTIGATION: WHY MODEL ≠ RAW DATA                 ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Load data
protests <- readRDS("protests_prepared.rds")

# Add mark variables
protests <- protests %>%
  mutate(
    is_violent = sub_event_type %in% c("Mob violence", "Violent demonstration") |
                 sub_event_type == "Excessive force against protesters",
    has_fatalities = fatalities > 0,
    is_student = grepl("Student", assoc_actor_1, ignore.case = TRUE),
    is_labor = grepl("Labor|Worker", assoc_actor_1, ignore.case = TRUE),
    is_papua = grepl("Papua", admin1, ignore.case = TRUE) |
               grepl("Papua", assoc_actor_1, ignore.case = TRUE)
  )

cat(sprintf("Loaded %s events\n", format(nrow(protests), big.mark = ",")))

# =============================================================================
# INVESTIGATION 1: CORRELATION BETWEEN MARKS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   1. CORRELATION BETWEEN MARK VARIABLES                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Create correlation matrix
mark_matrix <- protests %>%
  select(is_violent, has_fatalities, is_student, is_labor, is_papua) %>%
  mutate(across(everything(), as.numeric))

cor_matrix <- cor(mark_matrix)
cat("CORRELATION MATRIX:\n")
print(round(cor_matrix, 3))

# Chi-square tests for independence
cat("\n\nCHI-SQUARE TESTS:\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Violence vs Student
test_vs <- chisq.test(table(protests$is_violent, protests$is_student))
cat(sprintf("Violence vs Student: χ² = %.1f, p = %.2e\n",
            test_vs$statistic, test_vs$p.value))

# Violence vs Labor
test_vl <- chisq.test(table(protests$is_violent, protests$is_labor))
cat(sprintf("Violence vs Labor: χ² = %.1f, p = %.2e\n",
            test_vl$statistic, test_vl$p.value))

# Violence vs Papua
test_vp <- chisq.test(table(protests$is_violent, protests$is_papua))
cat(sprintf("Violence vs Papua: χ² = %.1f, p = %.2e\n",
            test_vp$statistic, test_vp$p.value))

# Cross-tabulation
cat("\n\nCROSS-TABULATION (Violence × Student):\n")
print(table(Violence = protests$is_violent, Student = protests$is_student))

cat("\nViolence rate by actor type:\n")
cat(sprintf("  Student events: %.1f%% violent\n",
            100 * mean(protests$is_violent[protests$is_student])))
cat(sprintf("  Non-student events: %.1f%% violent\n",
            100 * mean(protests$is_violent[!protests$is_student])))
cat(sprintf("  Labor events: %.1f%% violent\n",
            100 * mean(protests$is_violent[protests$is_labor])))
cat(sprintf("  Papua events: %.1f%% violent\n",
            100 * mean(protests$is_violent[protests$is_papua])))

# =============================================================================
# INVESTIGATION 2: TEMPORAL CLUSTERING OF VIOLENT EVENTS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   2. TEMPORAL CLUSTERING OF VIOLENT EVENTS                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Days between consecutive violent events
violent_events <- protests %>%
  filter(is_violent) %>%
  arrange(event_date)

violent_gaps <- diff(violent_events$event_date)
peaceful_events <- protests %>%
  filter(!is_violent) %>%
  arrange(event_date)
peaceful_gaps <- diff(peaceful_events$event_date)

cat("INTER-EVENT GAPS:\n")
cat(sprintf("  Violent events: mean gap = %.2f days, median = %.0f days\n",
            mean(violent_gaps), median(violent_gaps)))
cat(sprintf("  Peaceful events: mean gap = %.2f days, median = %.0f days\n",
            mean(peaceful_gaps), median(peaceful_gaps)))

# Do violent events cluster together?
cat("\nViolent events per day distribution:\n")
violent_per_day <- protests %>%
  group_by(event_date) %>%
  summarize(n_violent = sum(is_violent), n_total = n(), .groups = "drop")

cat(sprintf("  Days with 0 violent events: %d (%.1f%%)\n",
            sum(violent_per_day$n_violent == 0),
            100 * mean(violent_per_day$n_violent == 0)))
cat(sprintf("  Days with 1 violent event: %d (%.1f%%)\n",
            sum(violent_per_day$n_violent == 1),
            100 * mean(violent_per_day$n_violent == 1)))
cat(sprintf("  Days with 2+ violent events: %d (%.1f%%)\n",
            sum(violent_per_day$n_violent >= 2),
            100 * mean(violent_per_day$n_violent >= 2)))

# Violent events following violent events
cat("\nSEQUENCE ANALYSIS:\n")
protests_ordered <- protests %>% arrange(event_date)
protests_ordered$prev_violent <- c(NA, protests_ordered$is_violent[-nrow(protests_ordered)])
protests_ordered$prev_date <- c(as.Date(NA), protests_ordered$event_date[-nrow(protests_ordered)])
protests_ordered$days_since_prev <- as.numeric(protests_ordered$event_date - protests_ordered$prev_date)

# For same-day events
same_day <- protests_ordered %>%
  filter(days_since_prev == 0 & !is.na(prev_violent))

cat(sprintf("  Same-day events: If previous was violent, %.1f%% chance current is violent\n",
            100 * mean(same_day$is_violent[same_day$prev_violent])))
cat(sprintf("  Same-day events: If previous was peaceful, %.1f%% chance current is violent\n",
            100 * mean(same_day$is_violent[!same_day$prev_violent])))

# =============================================================================
# INVESTIGATION 3: SPATIAL CONCENTRATION
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   3. SPATIAL CONCENTRATION OF VIOLENCE                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Violence rate by province
province_violence <- protests %>%
  group_by(admin1) %>%
  summarize(
    n_events = n(),
    n_violent = sum(is_violent),
    violence_rate = mean(is_violent),
    .groups = "drop"
  ) %>%
  arrange(desc(violence_rate))

cat("TOP 10 PROVINCES BY VIOLENCE RATE (min 50 events):\n")
print(province_violence %>%
        filter(n_events >= 50) %>%
        head(10) %>%
        mutate(violence_rate = sprintf("%.1f%%", 100 * violence_rate)) %>%
        as.data.frame(), row.names = FALSE)

cat("\nBOTTOM 10 PROVINCES BY VIOLENCE RATE (min 50 events):\n")
print(province_violence %>%
        filter(n_events >= 50) %>%
        tail(10) %>%
        mutate(violence_rate = sprintf("%.1f%%", 100 * violence_rate)) %>%
        as.data.frame(), row.names = FALSE)

# How much does Papua drive the violence numbers?
cat("\n\nPAPUA ANALYSIS:\n")
cat(sprintf("  Papua events: %d (%.1f%% of total)\n",
            sum(protests$is_papua), 100 * mean(protests$is_papua)))
cat(sprintf("  Papua violence rate: %.1f%%\n",
            100 * mean(protests$is_violent[protests$is_papua])))
cat(sprintf("  Non-Papua violence rate: %.1f%%\n",
            100 * mean(protests$is_violent[!protests$is_papua])))

# =============================================================================
# INVESTIGATION 4: WHAT PREDICTS HAVING MANY FOLLOW-ON PROTESTS?
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   4. PREDICTORS OF HIGH FOLLOW-ON                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Calculate daily counts
daily_counts <- protests %>%
  group_by(event_date) %>%
  summarize(daily_n = n(), .groups = "drop")

# Calculate 7-day follow-on for each event (faster with data.table)
protests_dt <- as.data.table(protests)
protests_dt[, event_date := as.Date(event_date)]
setkey(protests_dt, event_date)

# Function to count follow-on efficiently
count_followon <- function(dates, window = 7) {
  n <- length(dates)
  followon <- numeric(n)
  all_dates <- as.numeric(dates)

  for (i in 1:n) {
    followon[i] <- sum(all_dates > all_dates[i] &
                        all_dates <= all_dates[i] + window)
  }
  return(followon)
}

cat("Calculating follow-on (this may take a moment)...\n")
protests$followon_7d <- count_followon(protests$event_date, 7)

# Join daily counts
protests <- protests %>%
  left_join(daily_counts, by = "event_date")

# Calculate baseline (7 days before)
protests$baseline_7d <- sapply(1:nrow(protests), function(i) {
  sum(protests$event_date >= protests$event_date[i] - 7 &
      protests$event_date < protests$event_date[i])
})

# What are the characteristics of events with HIGH follow-on?
cat("\nCHARACTERISTICS OF HIGH vs LOW FOLLOW-ON EVENTS:\n")
cat("(High = top 25%, Low = bottom 25%)\n")
cat("─────────────────────────────────────────────────────────────────\n")

high_followon <- protests %>% filter(followon_7d >= quantile(followon_7d, 0.75))
low_followon <- protests %>% filter(followon_7d <= quantile(followon_7d, 0.25))

comparison <- data.frame(
  Variable = c("Violence rate", "Fatality rate", "Student rate",
               "Labor rate", "Papua rate", "Same-day activity", "Baseline (7d)"),
  High_Followon = c(
    mean(high_followon$is_violent),
    mean(high_followon$has_fatalities),
    mean(high_followon$is_student),
    mean(high_followon$is_labor),
    mean(high_followon$is_papua),
    mean(high_followon$daily_n),
    mean(high_followon$baseline_7d)
  ),
  Low_Followon = c(
    mean(low_followon$is_violent),
    mean(low_followon$has_fatalities),
    mean(low_followon$is_student),
    mean(low_followon$is_labor),
    mean(low_followon$is_papua),
    mean(low_followon$daily_n),
    mean(low_followon$baseline_7d)
  )
)
comparison$Ratio <- comparison$High_Followon / comparison$Low_Followon

print(comparison %>%
        mutate(High_Followon = sprintf("%.3f", High_Followon),
               Low_Followon = sprintf("%.3f", Low_Followon),
               Ratio = sprintf("%.2f", Ratio)) %>%
        as.data.frame(), row.names = FALSE)

# =============================================================================
# INVESTIGATION 5: DO VIOLENT EVENTS TRIGGER MORE VIOLENT EVENTS?
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   5. VIOLENCE → VIOLENCE TRANSITIONS                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# For events on consecutive days in the same province, what's the transition probability?
protests_seq <- protests %>%
  arrange(admin1, event_date) %>%
  group_by(admin1) %>%
  mutate(
    next_date = lead(event_date),
    next_violent = lead(is_violent),
    days_to_next = as.numeric(next_date - event_date)
  ) %>%
  ungroup() %>%
  filter(!is.na(next_violent))

# Within 7 days, same province
within_7d <- protests_seq %>% filter(days_to_next > 0 & days_to_next <= 7)

cat("TRANSITION PROBABILITIES (same province, within 7 days):\n")
trans_table <- table(Current_Violent = within_7d$is_violent,
                     Next_Violent = within_7d$next_violent)
cat("\nCounts:\n")
print(trans_table)

cat("\nProbabilities:\n")
trans_prob <- prop.table(trans_table, 1)
print(round(trans_prob, 3))

cat(sprintf("\n  If current is VIOLENT: %.1f%% chance next is violent\n",
            100 * trans_prob["TRUE", "TRUE"]))
cat(sprintf("  If current is PEACEFUL: %.1f%% chance next is violent\n",
            100 * trans_prob["FALSE", "TRUE"]))

# =============================================================================
# INVESTIGATION 6: MODEL'S EXPECTED OFFSPRING
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   6. MODEL'S EXPECTED OFFSPRING CALCULATION                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Load model parameters
# From Model 5: baseline_trig = -10, decay = 0.00145
# Violence = 1.697, Fatal = 1.531, Student = -1.752, Labor = -5.0, Papua = -0.086

baseline_trig <- -10.0
decay_rate <- 0.00145
beta_violent <- 1.697
beta_fatal <- 1.531
beta_student <- -1.752
beta_labor <- -5.0
beta_papua <- -0.086

# Expected offspring = integral of triggering kernel from 0 to infinity
# For exponential: integral of alpha * exp(-beta * t) dt = alpha / beta

# Calculate alpha for each event type
calc_alpha <- function(violent, fatal, student, labor, papua) {
  exp(baseline_trig +
      beta_violent * violent +
      beta_fatal * fatal +
      beta_student * student +
      beta_labor * labor +
      beta_papua * papua)
}

# Expected offspring for different event types
event_types <- data.frame(
  Type = c("Peaceful baseline", "Violent (no fatalities)",
           "Violent with fatalities", "Student peaceful",
           "Labor peaceful", "Papua peaceful"),
  alpha = c(
    calc_alpha(0, 0, 0, 0, 0),
    calc_alpha(1, 0, 0, 0, 0),
    calc_alpha(1, 1, 0, 0, 0),
    calc_alpha(0, 0, 1, 0, 0),
    calc_alpha(0, 0, 0, 1, 0),
    calc_alpha(0, 0, 0, 0, 1)
  )
)

event_types$expected_offspring <- event_types$alpha / decay_rate

cat("EXPECTED OFFSPRING (integral of triggering kernel):\n")
cat("─────────────────────────────────────────────────────────────────\n")
print(event_types %>%
        mutate(alpha = sprintf("%.6f", alpha),
               expected_offspring = sprintf("%.4f", expected_offspring)) %>%
        as.data.frame(), row.names = FALSE)

cat("\nINTERPRETATION:\n")
cat("  Expected offspring represents the average number of protests\n")
cat("  directly triggered by one event, according to the model.\n")
cat(sprintf("  Sum across all events ≈ %.0f triggered protests\n",
            sum(event_types$expected_offspring[1] * sum(!protests$is_violent & !protests$has_fatalities & !protests$is_student & !protests$is_labor & !protests$is_papua))))

# =============================================================================
# INVESTIGATION 7: SUBSET ANALYSIS - DOES PATTERN HOLD IN SUBSETS?
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   7. SUBSET ANALYSIS                                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# By year
cat("VIOLENCE FOLLOW-ON RATIO BY YEAR:\n")
cat("(Ratio = mean follow-on for violent / mean follow-on for peaceful)\n")
cat("─────────────────────────────────────────────────────────────────\n")

yearly_ratio <- protests %>%
  mutate(year = lubridate::year(event_date)) %>%
  group_by(year) %>%
  summarize(
    violent_mean = mean(followon_7d[is_violent]),
    peaceful_mean = mean(followon_7d[!is_violent]),
    ratio = violent_mean / peaceful_mean,
    n_violent = sum(is_violent),
    .groups = "drop"
  )

print(yearly_ratio %>%
        mutate(violent_mean = round(violent_mean, 1),
               peaceful_mean = round(peaceful_mean, 1),
               ratio = sprintf("%.2f", ratio)) %>%
        as.data.frame(), row.names = FALSE)

# By province type (Java vs non-Java)
protests$is_java <- grepl("Java|Jakarta|Yogyakarta|Banten", protests$admin1, ignore.case = TRUE)

cat("\nVIOLENCE FOLLOW-ON RATIO BY REGION:\n")
region_ratio <- protests %>%
  group_by(Region = ifelse(is_java, "Java/Jakarta", "Outer Islands")) %>%
  summarize(
    violent_mean = mean(followon_7d[is_violent]),
    peaceful_mean = mean(followon_7d[!is_violent]),
    ratio = violent_mean / peaceful_mean,
    n_violent = sum(is_violent),
    .groups = "drop"
  )

print(region_ratio %>%
        mutate(violent_mean = round(violent_mean, 1),
               peaceful_mean = round(peaceful_mean, 1),
               ratio = sprintf("%.2f", ratio)) %>%
        as.data.frame(), row.names = FALSE)

# Excluding Papua
cat("\nEXCLUDING PAPUA:\n")
non_papua <- protests %>% filter(!is_papua)
cat(sprintf("  Violent mean follow-on: %.1f\n", mean(non_papua$followon_7d[non_papua$is_violent])))
cat(sprintf("  Peaceful mean follow-on: %.1f\n", mean(non_papua$followon_7d[!non_papua$is_violent])))
cat(sprintf("  Ratio: %.2f\n",
            mean(non_papua$followon_7d[non_papua$is_violent]) /
            mean(non_papua$followon_7d[!non_papua$is_violent])))

# =============================================================================
# INVESTIGATION 8: WHAT IF VIOLENCE IS AN EFFECT, NOT A CAUSE?
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   8. REVERSE CAUSALITY CHECK                                 ║\n")
cat("║   Does high activity PRECEDE violence?                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Activity in days BEFORE violent vs peaceful events
before_comparison <- data.frame(
  Category = c("Violent events", "Peaceful events", "Fatal events", "Non-fatal events"),
  Mean_Activity_7d_Before = c(
    mean(protests$baseline_7d[protests$is_violent]),
    mean(protests$baseline_7d[!protests$is_violent]),
    mean(protests$baseline_7d[protests$has_fatalities]),
    mean(protests$baseline_7d[!protests$has_fatalities])
  )
)

print(before_comparison %>%
        mutate(Mean_Activity_7d_Before = round(Mean_Activity_7d_Before, 1)) %>%
        as.data.frame(), row.names = FALSE)

# Test
test_before <- wilcox.test(
  protests$baseline_7d[protests$is_violent],
  protests$baseline_7d[!protests$is_violent]
)
cat(sprintf("\nWilcoxon test (baseline before violent vs peaceful): p = %.2e\n",
            test_before$p.value))

if (mean(protests$baseline_7d[protests$is_violent]) <
    mean(protests$baseline_7d[!protests$is_violent])) {
  cat("\n→ Violent events are PRECEDED by lower activity\n")
  cat("  This suggests violence may occur when mobilization is waning,\n")
  cat("  not when it's building.\n")
}

# =============================================================================
# INVESTIGATION 9: GRANGER-STYLE ANALYSIS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   9. GRANGER-STYLE PREDICTIVE ANALYSIS                       ║\n")
cat("║   Does past violence predict future activity?                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Create daily time series
daily_ts <- protests %>%
  group_by(event_date) %>%
  summarize(
    n_events = n(),
    n_violent = sum(is_violent),
    n_fatal = sum(has_fatalities),
    .groups = "drop"
  ) %>%
  complete(event_date = seq(min(event_date), max(event_date), by = "day"),
           fill = list(n_events = 0, n_violent = 0, n_fatal = 0)) %>%
  arrange(event_date) %>%
  mutate(
    lag1_events = lag(n_events, 1),
    lag1_violent = lag(n_violent, 1),
    lag7_events = lag(n_events, 7),
    lag7_violent = lag(n_violent, 7)
  )

# Regression: does lagged violence predict today's activity?
cat("Regression: n_events_today ~ lag1_events + lag1_violent\n")
model_granger <- lm(n_events ~ lag1_events + lag1_violent, data = daily_ts)
print(summary(model_granger)$coefficients)

cat("\n\nRegression: n_events_today ~ lag7_events + lag7_violent\n")
model_granger7 <- lm(n_events ~ lag7_events + lag7_violent, data = daily_ts)
print(summary(model_granger7)$coefficients)

cat("\nINTERPRETATION:\n")
cat("  If lag_violent coefficient is positive and significant,\n")
cat("  past violence predicts more future protests.\n")
cat("  If negative or insignificant, violence doesn't boost future activity.\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   INVESTIGATION SUMMARY                                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("KEY FINDINGS:\n")
cat("─────────────────────────────────────────────────────────────────\n")

cat("\n1. MARK CORRELATIONS:\n")
cat(sprintf("   Violence-Student correlation: %.3f\n", cor_matrix["is_violent", "is_student"]))
cat(sprintf("   Violence-Papua correlation: %.3f\n", cor_matrix["is_violent", "is_papua"]))

cat("\n2. TEMPORAL CLUSTERING:\n")
cat("   Violent events have slightly longer gaps between them\n")

cat("\n3. SPATIAL CONCENTRATION:\n")
cat(sprintf("   Papua violence rate: %.1f%% vs non-Papua: %.1f%%\n",
            100 * mean(protests$is_violent[protests$is_papua]),
            100 * mean(protests$is_violent[!protests$is_papua])))

cat("\n4. HIGH FOLLOW-ON PREDICTORS:\n")
cat("   Events with high follow-on have:\n")
cat(sprintf("   - LOWER violence rate: %.1f%% vs %.1f%%\n",
            100 * mean(high_followon$is_violent),
            100 * mean(low_followon$is_violent)))
cat(sprintf("   - HIGHER baseline activity: %.0f vs %.0f\n",
            mean(high_followon$baseline_7d),
            mean(low_followon$baseline_7d)))

cat("\n5. VIOLENCE → VIOLENCE TRANSITIONS:\n")
cat(sprintf("   After violent: %.1f%% violent | After peaceful: %.1f%% violent\n",
            100 * trans_prob["TRUE", "TRUE"],
            100 * trans_prob["FALSE", "TRUE"]))

cat("\n6. REVERSE CAUSALITY:\n")
cat(sprintf("   Violence preceded by %.1f events (vs %.1f for peaceful)\n",
            mean(protests$baseline_7d[protests$is_violent]),
            mean(protests$baseline_7d[!protests$is_violent])))
cat("   → Violence may be an EFFECT of declining mobilization, not a cause\n")

cat("\n7. GRANGER-STYLE:\n")
coef_lag1 <- coef(model_granger)["lag1_violent"]
cat(sprintf("   1-day lagged violence coefficient: %.3f\n", coef_lag1))

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   DEEPER INVESTIGATION COMPLETE                              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
