################################################################################
#                   VALIDATION ANALYSIS: RAW DATA PATTERNS
#                   Confirm Hawkes Model Findings
################################################################################
#
# This script validates the Model 5 findings by examining raw data patterns:
#   - Violence → 5.5× contagion
#   - Fatalities → 4.6× additional boost
#   - Student → 0.17× suppression
#   - Labor → 0.007× suppression
#   - Papua → 0.92× neutral
#
################################################################################

library(tidyverse)
library(kableExtra)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   VALIDATION ANALYSIS: RAW DATA PATTERNS                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================

cat("=== LOADING DATA ===\n")
protests <- readRDS("protests_prepared.rds")
cat(sprintf("Loaded %s events\n", format(nrow(protests), big.mark = ",")))

# Add mark variables (matching Model 5 specification)
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

# Summary of marks
cat("\nMARK VARIABLE SUMMARY:\n")
cat(sprintf("  Violent events: %d (%.1f%%)\n", sum(protests$is_violent), 100*mean(protests$is_violent)))
cat(sprintf("  Fatal events: %d (%.1f%%)\n", sum(protests$has_fatalities), 100*mean(protests$has_fatalities)))
cat(sprintf("  Student-led: %d (%.1f%%)\n", sum(protests$is_student), 100*mean(protests$is_student)))
cat(sprintf("  Labor-led: %d (%.1f%%)\n", sum(protests$is_labor), 100*mean(protests$is_labor)))
cat(sprintf("  Papua-related: %d (%.1f%%)\n", sum(protests$is_papua), 100*mean(protests$is_papua)))

# Get all unique dates for follow-on calculations
all_dates <- sort(unique(protests$event_date))

# =============================================================================
# ANALYSIS 1: RAW FOLLOW-ON COUNTS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   ANALYSIS 1: RAW FOLLOW-ON COUNTS                           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Function to calculate follow-on counts for each event
calc_followon <- function(data, window_days = 7) {
  data %>%
    rowwise() %>%
    mutate(
      followon = sum(protests$event_date > event_date &
                     protests$event_date <= event_date + window_days)
    ) %>%
    ungroup()
}

# Calculate 7-day follow-on for all events
cat("Calculating 7-day follow-on counts for all events...\n")
protests_with_followon <- calc_followon(protests, 7)

# Compare by mark type
followon_by_type <- data.frame(
  Category = c("Violent", "Peaceful", "Fatal", "Non-fatal",
               "Student", "Non-student", "Labor", "Non-labor",
               "Papua", "Non-Papua"),
  N = c(sum(protests_with_followon$is_violent),
        sum(!protests_with_followon$is_violent),
        sum(protests_with_followon$has_fatalities),
        sum(!protests_with_followon$has_fatalities),
        sum(protests_with_followon$is_student),
        sum(!protests_with_followon$is_student),
        sum(protests_with_followon$is_labor),
        sum(!protests_with_followon$is_labor),
        sum(protests_with_followon$is_papua),
        sum(!protests_with_followon$is_papua)),
  Mean_Followon = c(
    mean(protests_with_followon$followon[protests_with_followon$is_violent]),
    mean(protests_with_followon$followon[!protests_with_followon$is_violent]),
    mean(protests_with_followon$followon[protests_with_followon$has_fatalities]),
    mean(protests_with_followon$followon[!protests_with_followon$has_fatalities]),
    mean(protests_with_followon$followon[protests_with_followon$is_student]),
    mean(protests_with_followon$followon[!protests_with_followon$is_student]),
    mean(protests_with_followon$followon[protests_with_followon$is_labor]),
    mean(protests_with_followon$followon[!protests_with_followon$is_labor]),
    mean(protests_with_followon$followon[protests_with_followon$is_papua]),
    mean(protests_with_followon$followon[!protests_with_followon$is_papua])
  )
)

# Calculate ratios
followon_by_type$Ratio <- NA
followon_by_type$Ratio[1] <- followon_by_type$Mean_Followon[1] / followon_by_type$Mean_Followon[2]  # Violent/Peaceful
followon_by_type$Ratio[3] <- followon_by_type$Mean_Followon[3] / followon_by_type$Mean_Followon[4]  # Fatal/Non-fatal
followon_by_type$Ratio[5] <- followon_by_type$Mean_Followon[5] / followon_by_type$Mean_Followon[6]  # Student/Non-student
followon_by_type$Ratio[7] <- followon_by_type$Mean_Followon[7] / followon_by_type$Mean_Followon[8]  # Labor/Non-labor
followon_by_type$Ratio[9] <- followon_by_type$Mean_Followon[9] / followon_by_type$Mean_Followon[10] # Papua/Non-Papua

cat("\n7-DAY FOLLOW-ON COUNTS BY EVENT TYPE:\n")
cat("─────────────────────────────────────────────────────────────────\n")
print(followon_by_type %>%
        mutate(Mean_Followon = round(Mean_Followon, 1),
               Ratio = ifelse(is.na(Ratio), "", sprintf("%.2fx", Ratio))) %>%
        format(justify = "right"), row.names = FALSE)

cat("\nINTERPRETATION:\n")
cat(sprintf("  Violent/Peaceful ratio: %.2fx (Model predicts: 5.5x)\n",
            followon_by_type$Ratio[1]))
cat(sprintf("  Fatal/Non-fatal ratio: %.2fx (Model predicts: 4.6x boost)\n",
            followon_by_type$Ratio[3]))
cat(sprintf("  Student/Non-student ratio: %.2fx (Model predicts: 0.17x)\n",
            followon_by_type$Ratio[5]))
cat(sprintf("  Labor/Non-labor ratio: %.2fx (Model predicts: 0.007x)\n",
            followon_by_type$Ratio[7]))
cat(sprintf("  Papua/Non-Papua ratio: %.2fx (Model predicts: 0.92x)\n",
            followon_by_type$Ratio[9]))

cat("\n⚠️  NOTE: Raw counts are confounded by baseline activity!\n")
cat("    Events during busy periods naturally have more follow-on.\n")
cat("    See Analysis 2 for baseline-controlled comparison.\n")

# =============================================================================
# ANALYSIS 2: NORMALIZED BEFORE/AFTER RATIOS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   ANALYSIS 2: NORMALIZED BEFORE/AFTER RATIOS                 ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Calculating before/after ratios for all events...\n")

# Calculate before and after counts for each event
protests_ba <- protests %>%
  rowwise() %>%
  mutate(
    # Count in 7 days before (not including event day)
    before_7 = sum(protests$event_date >= event_date - 7 &
                   protests$event_date < event_date),
    # Count in 7 days after (not including event day)
    after_7 = sum(protests$event_date > event_date &
                  protests$event_date <= event_date + 7),
    # Ratio (add 1 to avoid division by zero)
    ba_ratio = (after_7 + 1) / (before_7 + 1),
    # Difference
    ba_diff = after_7 - before_7
  ) %>%
  ungroup()

# Compare by mark type
ba_by_type <- data.frame(
  Category = c("Violent", "Peaceful", "Fatal", "Non-fatal",
               "Student", "Non-student", "Labor", "Non-labor",
               "Papua", "Non-Papua"),
  N = c(sum(protests_ba$is_violent), sum(!protests_ba$is_violent),
        sum(protests_ba$has_fatalities), sum(!protests_ba$has_fatalities),
        sum(protests_ba$is_student), sum(!protests_ba$is_student),
        sum(protests_ba$is_labor), sum(!protests_ba$is_labor),
        sum(protests_ba$is_papua), sum(!protests_ba$is_papua)),
  Mean_Before = c(
    mean(protests_ba$before_7[protests_ba$is_violent]),
    mean(protests_ba$before_7[!protests_ba$is_violent]),
    mean(protests_ba$before_7[protests_ba$has_fatalities]),
    mean(protests_ba$before_7[!protests_ba$has_fatalities]),
    mean(protests_ba$before_7[protests_ba$is_student]),
    mean(protests_ba$before_7[!protests_ba$is_student]),
    mean(protests_ba$before_7[protests_ba$is_labor]),
    mean(protests_ba$before_7[!protests_ba$is_labor]),
    mean(protests_ba$before_7[protests_ba$is_papua]),
    mean(protests_ba$before_7[!protests_ba$is_papua])
  ),
  Mean_After = c(
    mean(protests_ba$after_7[protests_ba$is_violent]),
    mean(protests_ba$after_7[!protests_ba$is_violent]),
    mean(protests_ba$after_7[protests_ba$has_fatalities]),
    mean(protests_ba$after_7[!protests_ba$has_fatalities]),
    mean(protests_ba$after_7[protests_ba$is_student]),
    mean(protests_ba$after_7[!protests_ba$is_student]),
    mean(protests_ba$after_7[protests_ba$is_labor]),
    mean(protests_ba$after_7[!protests_ba$is_labor]),
    mean(protests_ba$after_7[protests_ba$is_papua]),
    mean(protests_ba$after_7[!protests_ba$is_papua])
  ),
  Median_Diff = c(
    median(protests_ba$ba_diff[protests_ba$is_violent]),
    median(protests_ba$ba_diff[!protests_ba$is_violent]),
    median(protests_ba$ba_diff[protests_ba$has_fatalities]),
    median(protests_ba$ba_diff[!protests_ba$has_fatalities]),
    median(protests_ba$ba_diff[protests_ba$is_student]),
    median(protests_ba$ba_diff[!protests_ba$is_student]),
    median(protests_ba$ba_diff[protests_ba$is_labor]),
    median(protests_ba$ba_diff[!protests_ba$is_labor]),
    median(protests_ba$ba_diff[protests_ba$is_papua]),
    median(protests_ba$ba_diff[!protests_ba$is_papua])
  )
)

cat("\nBEFORE/AFTER COMPARISON (7-day windows):\n")
cat("─────────────────────────────────────────────────────────────────\n")
print(ba_by_type %>%
        mutate(Mean_Before = round(Mean_Before, 1),
               Mean_After = round(Mean_After, 1),
               Change = sprintf("%+.1f", Mean_After - Mean_Before)) %>%
        select(Category, N, Mean_Before, Mean_After, Change, Median_Diff) %>%
        format(justify = "right"), row.names = FALSE)

# Statistical tests
cat("\nSTATISTICAL TESTS (Wilcoxon rank-sum on after-before difference):\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Violent vs Peaceful
test_violent <- wilcox.test(
  protests_ba$ba_diff[protests_ba$is_violent],
  protests_ba$ba_diff[!protests_ba$is_violent]
)
cat(sprintf("  Violent vs Peaceful: W = %.0f, p = %.2e\n",
            test_violent$statistic, test_violent$p.value))

# Fatal vs Non-fatal
test_fatal <- wilcox.test(
  protests_ba$ba_diff[protests_ba$has_fatalities],
  protests_ba$ba_diff[!protests_ba$has_fatalities]
)
cat(sprintf("  Fatal vs Non-fatal: W = %.0f, p = %.2e\n",
            test_fatal$statistic, test_fatal$p.value))

# Student vs Non-student
test_student <- wilcox.test(
  protests_ba$ba_diff[protests_ba$is_student],
  protests_ba$ba_diff[!protests_ba$is_student]
)
cat(sprintf("  Student vs Non-student: W = %.0f, p = %.2e\n",
            test_student$statistic, test_student$p.value))

# Labor vs Non-labor
test_labor <- wilcox.test(
  protests_ba$ba_diff[protests_ba$is_labor],
  protests_ba$ba_diff[!protests_ba$is_labor]
)
cat(sprintf("  Labor vs Non-labor: W = %.0f, p = %.2e\n",
            test_labor$statistic, test_labor$p.value))

# =============================================================================
# ANALYSIS 3: CASE STUDIES - HIGH IMPACT EVENTS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   ANALYSIS 3: CASE STUDIES - HIGH IMPACT EVENTS              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Find top 20 most fatal events
top_fatal <- protests %>%
  filter(fatalities > 0) %>%
  arrange(desc(fatalities)) %>%
  head(20) %>%
  select(event_date, admin1, admin2, fatalities, sub_event_type, notes)

cat("TOP 20 MOST FATAL EVENTS:\n")
cat("─────────────────────────────────────────────────────────────────\n")

for (i in 1:nrow(top_fatal)) {
  event <- top_fatal[i, ]

  # Count protests in 7 days before and after
  before <- sum(protests$event_date >= event$event_date - 7 &
                protests$event_date < event$event_date)
  after <- sum(protests$event_date > event$event_date &
               protests$event_date <= event$event_date + 7)

  cat(sprintf("\n%d. %s - %s, %s (%d fatalities)\n",
              i, format(event$event_date, "%Y-%m-%d"),
              event$admin2, event$admin1, event$fatalities))
  cat(sprintf("   Type: %s\n", event$sub_event_type))
  cat(sprintf("   7 days before: %d protests | 7 days after: %d protests (change: %+d)\n",
              before, after, after - before))

  # Truncate notes
  if (!is.na(event$notes) && nchar(event$notes) > 100) {
    cat(sprintf("   Notes: %s...\n", substr(event$notes, 1, 100)))
  } else if (!is.na(event$notes)) {
    cat(sprintf("   Notes: %s\n", event$notes))
  }
}

# =============================================================================
# ANALYSIS 4: SAME-DISTRICT CLUSTERING
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   ANALYSIS 4: SAME-DISTRICT CLUSTERING                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Calculating local vs national follow-on patterns...\n")

# For each event, calculate same-district, same-province, different-province follow-on
protests_spatial <- protests %>%
  rowwise() %>%
  mutate(
    followon_same_district = sum(protests$event_date > event_date &
                                  protests$event_date <= event_date + 7 &
                                  protests$admin2 == admin2),
    followon_same_province = sum(protests$event_date > event_date &
                                  protests$event_date <= event_date + 7 &
                                  protests$admin1 == admin1 &
                                  protests$admin2 != admin2),
    followon_other = sum(protests$event_date > event_date &
                         protests$event_date <= event_date + 7 &
                         protests$admin1 != admin1),
    followon_total = followon_same_district + followon_same_province + followon_other,
    pct_local = ifelse(followon_total > 0,
                       100 * (followon_same_district + followon_same_province) / followon_total,
                       NA)
  ) %>%
  ungroup()

# Compare by type
spatial_by_type <- data.frame(
  Category = c("Violent", "Peaceful", "Fatal", "Non-fatal",
               "Student", "Non-student", "Labor", "Non-labor"),
  Mean_Same_District = c(
    mean(protests_spatial$followon_same_district[protests_spatial$is_violent]),
    mean(protests_spatial$followon_same_district[!protests_spatial$is_violent]),
    mean(protests_spatial$followon_same_district[protests_spatial$has_fatalities]),
    mean(protests_spatial$followon_same_district[!protests_spatial$has_fatalities]),
    mean(protests_spatial$followon_same_district[protests_spatial$is_student]),
    mean(protests_spatial$followon_same_district[!protests_spatial$is_student]),
    mean(protests_spatial$followon_same_district[protests_spatial$is_labor]),
    mean(protests_spatial$followon_same_district[!protests_spatial$is_labor])
  ),
  Mean_Same_Province = c(
    mean(protests_spatial$followon_same_province[protests_spatial$is_violent]),
    mean(protests_spatial$followon_same_province[!protests_spatial$is_violent]),
    mean(protests_spatial$followon_same_province[protests_spatial$has_fatalities]),
    mean(protests_spatial$followon_same_province[!protests_spatial$has_fatalities]),
    mean(protests_spatial$followon_same_province[protests_spatial$is_student]),
    mean(protests_spatial$followon_same_province[!protests_spatial$is_student]),
    mean(protests_spatial$followon_same_province[protests_spatial$is_labor]),
    mean(protests_spatial$followon_same_province[!protests_spatial$is_labor])
  ),
  Mean_National = c(
    mean(protests_spatial$followon_other[protests_spatial$is_violent]),
    mean(protests_spatial$followon_other[!protests_spatial$is_violent]),
    mean(protests_spatial$followon_other[protests_spatial$has_fatalities]),
    mean(protests_spatial$followon_other[!protests_spatial$has_fatalities]),
    mean(protests_spatial$followon_other[protests_spatial$is_student]),
    mean(protests_spatial$followon_other[!protests_spatial$is_student]),
    mean(protests_spatial$followon_other[protests_spatial$is_labor]),
    mean(protests_spatial$followon_other[!protests_spatial$is_labor])
  )
)

cat("\nSPATIAL DISTRIBUTION OF 7-DAY FOLLOW-ON:\n")
cat("─────────────────────────────────────────────────────────────────\n")
print(spatial_by_type %>%
        mutate(across(starts_with("Mean"), ~round(., 2))) %>%
        format(justify = "right"), row.names = FALSE)

cat("\nINTERPRETATION:\n")
cat("  If contagion is local (true triggering), same-district should be elevated.\n")
cat("  If coordination is national (e.g., labor strikes), national should dominate.\n")

# =============================================================================
# ANALYSIS 5: PLACEBO TEST
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   ANALYSIS 5: PLACEBO TEST (PERMUTATION)                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Running 1000 permutations with shuffled violence labels...\n")

# Calculate actual effect (violent - peaceful median difference)
actual_effect <- median(protests_ba$ba_diff[protests_ba$is_violent]) -
                 median(protests_ba$ba_diff[!protests_ba$is_violent])

# Permutation test
set.seed(42)
n_perm <- 1000
perm_effects <- numeric(n_perm)

for (i in 1:n_perm) {
  # Shuffle violence labels
  shuffled_violent <- sample(protests_ba$is_violent)

  # Calculate effect with shuffled labels
  perm_effects[i] <- median(protests_ba$ba_diff[shuffled_violent]) -
                     median(protests_ba$ba_diff[!shuffled_violent])

  if (i %% 200 == 0) cat(sprintf("  Completed %d/%d permutations\n", i, n_perm))
}

# Calculate p-value
p_value <- mean(perm_effects >= actual_effect)

cat("\nPLACEBO TEST RESULTS:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Actual effect (violent - peaceful): %.2f\n", actual_effect))
cat(sprintf("  Permutation mean: %.2f\n", mean(perm_effects)))
cat(sprintf("  Permutation SD: %.2f\n", sd(perm_effects)))
cat(sprintf("  Z-score: %.2f\n", (actual_effect - mean(perm_effects)) / sd(perm_effects)))
cat(sprintf("  Permutation p-value: %.4f\n", p_value))

if (p_value < 0.05) {
  cat("\n✓ CONCLUSION: The violence effect is statistically significant.\n")
  cat("  The observed pattern is unlikely to occur by chance.\n")
} else {
  cat("\n⚠ CONCLUSION: The violence effect is NOT statistically significant.\n")
  cat("  The observed pattern could occur by chance.\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   VALIDATION SUMMARY                                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("MODEL 5 PREDICTIONS vs RAW DATA PATTERNS:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Violence: Model predicts 5.5x | Raw ratio: %.2fx | ", followon_by_type$Ratio[1]))
if (followon_by_type$Ratio[1] > 1.5) cat("✓ CONFIRMED\n") else cat("? INCONCLUSIVE\n")

cat(sprintf("  Fatalities: Model predicts 4.6x boost | Raw ratio: %.2fx | ", followon_by_type$Ratio[3]))
if (followon_by_type$Ratio[3] > 1.5) cat("✓ CONFIRMED\n") else cat("? INCONCLUSIVE\n")

cat(sprintf("  Student: Model predicts 0.17x | Raw ratio: %.2fx | ", followon_by_type$Ratio[5]))
if (followon_by_type$Ratio[5] < 1) cat("✓ CONFIRMED\n") else cat("✗ NOT CONFIRMED\n")

cat(sprintf("  Labor: Model predicts 0.007x | Raw ratio: %.2fx | ", followon_by_type$Ratio[7]))
if (followon_by_type$Ratio[7] < 1) cat("✓ CONFIRMED\n") else cat("✗ NOT CONFIRMED\n")

cat(sprintf("  Papua: Model predicts 0.92x | Raw ratio: %.2fx | ", followon_by_type$Ratio[9]))
if (abs(followon_by_type$Ratio[9] - 1) < 0.3) cat("✓ CONFIRMED (neutral)\n") else cat("? INCONCLUSIVE\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   VALIDATION ANALYSIS COMPLETE                               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")

# Save results for paper
validation_results <- list(
  followon_by_type = followon_by_type,
  ba_by_type = ba_by_type,
  spatial_by_type = spatial_by_type,
  placebo = list(
    actual_effect = actual_effect,
    perm_mean = mean(perm_effects),
    perm_sd = sd(perm_effects),
    p_value = p_value
  ),
  test_results = list(
    violent = test_violent,
    fatal = test_fatal,
    student = test_student,
    labor = test_labor
  )
)

saveRDS(validation_results, "validation_results.rds")
cat("\n✓ Saved validation results to validation_results.rds\n")
