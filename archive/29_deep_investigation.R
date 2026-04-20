################################################################################
#                   DEEP INVESTIGATION: WHY DON'T RAW PATTERNS MATCH?
#                   Testing hypotheses for validation discrepancy
################################################################################

library(tidyverse)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   DEEP INVESTIGATION: VALIDATION DISCREPANCY                 ║\n")
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
cat(sprintf("  Violent: %d (%.1f%%)\n", sum(protests$is_violent), 100*mean(protests$is_violent)))
cat(sprintf("  Fatal: %d (%.1f%%)\n", sum(protests$has_fatalities), 100*mean(protests$has_fatalities)))

# =============================================================================
# H1: TIMING CONFOUND - Do violent events occur during quieter periods?
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   H1: TIMING CONFOUND                                        ║\n")
cat("║   Do violent events occur during quieter periods?            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Calculate daily protest counts
daily_counts <- protests %>%
  group_by(event_date) %>%
  summarize(daily_n = n(), .groups = "drop")

# Join to get same-day activity for each event
protests_with_daily <- protests %>%
  left_join(daily_counts, by = "event_date")

# Compare same-day activity
sameday_comparison <- data.frame(
  Category = c("Violent", "Peaceful", "Fatal", "Non-fatal", "Student", "Non-student"),
  Mean_Same_Day = c(
    mean(protests_with_daily$daily_n[protests_with_daily$is_violent]),
    mean(protests_with_daily$daily_n[!protests_with_daily$is_violent]),
    mean(protests_with_daily$daily_n[protests_with_daily$has_fatalities]),
    mean(protests_with_daily$daily_n[!protests_with_daily$has_fatalities]),
    mean(protests_with_daily$daily_n[protests_with_daily$is_student]),
    mean(protests_with_daily$daily_n[!protests_with_daily$is_student])
  ),
  Median_Same_Day = c(
    median(protests_with_daily$daily_n[protests_with_daily$is_violent]),
    median(protests_with_daily$daily_n[!protests_with_daily$is_violent]),
    median(protests_with_daily$daily_n[protests_with_daily$has_fatalities]),
    median(protests_with_daily$daily_n[!protests_with_daily$has_fatalities]),
    median(protests_with_daily$daily_n[protests_with_daily$is_student]),
    median(protests_with_daily$daily_n[!protests_with_daily$is_student])
  )
)

cat("SAME-DAY PROTEST ACTIVITY BY EVENT TYPE:\n")
cat("(If violent events occur during quieter periods, they'd have lower same-day counts)\n")
cat("─────────────────────────────────────────────────────────────────\n")
print(sameday_comparison %>%
        mutate(Mean_Same_Day = round(Mean_Same_Day, 1),
               Median_Same_Day = round(Median_Same_Day, 0)) %>%
        format(justify = "right"), row.names = FALSE)

# Statistical test
test_sameday <- wilcox.test(
  protests_with_daily$daily_n[protests_with_daily$is_violent],
  protests_with_daily$daily_n[!protests_with_daily$is_violent]
)
cat(sprintf("\nWilcoxon test (violent vs peaceful): p = %.2e\n", test_sameday$p.value))

if (mean(protests_with_daily$daily_n[protests_with_daily$is_violent]) <
    mean(protests_with_daily$daily_n[!protests_with_daily$is_violent])) {
  cat("\n✓ CONFIRMED: Violent events DO occur during quieter periods!\n")
  cat("  This explains why raw follow-on counts are lower.\n")
} else {
  cat("\n✗ NOT CONFIRMED: Violent events don't occur during quieter periods.\n")
}

# =============================================================================
# H2: SAME-DISTRICT TRIGGERING - Is the effect local?
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   H2: SAME-DISTRICT TRIGGERING                               ║\n")
cat("║   Does violence trigger more protests LOCALLY?               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# For each event, count same-district follow-on
protests_local <- protests %>%
  rowwise() %>%
  mutate(
    # Same district, next 7 days
    local_followon_7d = sum(protests$event_date > event_date &
                            protests$event_date <= event_date + 7 &
                            protests$admin2 == admin2 &
                            protests$admin1 == admin1),
    # Same district, next 3 days
    local_followon_3d = sum(protests$event_date > event_date &
                            protests$event_date <= event_date + 3 &
                            protests$admin2 == admin2 &
                            protests$admin1 == admin1),
    # Same district, next 1 day
    local_followon_1d = sum(protests$event_date > event_date &
                            protests$event_date <= event_date + 1 &
                            protests$admin2 == admin2 &
                            protests$admin1 == admin1)
  ) %>%
  ungroup()

# Compare
local_comparison <- data.frame(
  Category = c("Violent", "Peaceful", "Fatal", "Non-fatal"),
  Mean_Local_1d = c(
    mean(protests_local$local_followon_1d[protests_local$is_violent]),
    mean(protests_local$local_followon_1d[!protests_local$is_violent]),
    mean(protests_local$local_followon_1d[protests_local$has_fatalities]),
    mean(protests_local$local_followon_1d[!protests_local$has_fatalities])
  ),
  Mean_Local_3d = c(
    mean(protests_local$local_followon_3d[protests_local$is_violent]),
    mean(protests_local$local_followon_3d[!protests_local$is_violent]),
    mean(protests_local$local_followon_3d[protests_local$has_fatalities]),
    mean(protests_local$local_followon_3d[!protests_local$has_fatalities])
  ),
  Mean_Local_7d = c(
    mean(protests_local$local_followon_7d[protests_local$is_violent]),
    mean(protests_local$local_followon_7d[!protests_local$is_violent]),
    mean(protests_local$local_followon_7d[protests_local$has_fatalities]),
    mean(protests_local$local_followon_7d[!protests_local$has_fatalities])
  )
)

# Add ratios
local_comparison$Ratio_1d <- NA
local_comparison$Ratio_1d[1] <- local_comparison$Mean_Local_1d[1] / local_comparison$Mean_Local_1d[2]
local_comparison$Ratio_1d[3] <- local_comparison$Mean_Local_1d[3] / local_comparison$Mean_Local_1d[4]

local_comparison$Ratio_7d <- NA
local_comparison$Ratio_7d[1] <- local_comparison$Mean_Local_7d[1] / local_comparison$Mean_Local_7d[2]
local_comparison$Ratio_7d[3] <- local_comparison$Mean_Local_7d[3] / local_comparison$Mean_Local_7d[4]

cat("SAME-DISTRICT FOLLOW-ON BY EVENT TYPE:\n")
cat("─────────────────────────────────────────────────────────────────\n")
print(local_comparison %>%
        mutate(across(starts_with("Mean"), ~round(., 3)),
               Ratio_1d = ifelse(is.na(Ratio_1d), "", sprintf("%.2fx", Ratio_1d)),
               Ratio_7d = ifelse(is.na(Ratio_7d), "", sprintf("%.2fx", Ratio_7d))) %>%
        format(justify = "right"), row.names = FALSE)

cat("\nINTERPRETATION:\n")
cat("  If violence triggers local protests, violent events should have HIGHER local follow-on.\n")

# =============================================================================
# H3: SHORT-TERM EFFECTS - Is the effect immediate?
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   H3: SHORT-TERM EFFECTS                                     ║\n")
cat("║   Does violence have stronger effects in first 1-2 days?     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Calculate follow-on at different windows (NATIONAL)
protests_windows <- protests %>%
  rowwise() %>%
  mutate(
    followon_1d = sum(protests$event_date > event_date &
                      protests$event_date <= event_date + 1),
    followon_2d = sum(protests$event_date > event_date &
                      protests$event_date <= event_date + 2),
    followon_3d = sum(protests$event_date > event_date &
                      protests$event_date <= event_date + 3),
    followon_7d = sum(protests$event_date > event_date &
                      protests$event_date <= event_date + 7),
    followon_14d = sum(protests$event_date > event_date &
                       protests$event_date <= event_date + 14)
  ) %>%
  ungroup()

# Compare violent vs peaceful at each window
windows_comparison <- data.frame(
  Window = c("1 day", "2 days", "3 days", "7 days", "14 days"),
  Violent_Mean = c(
    mean(protests_windows$followon_1d[protests_windows$is_violent]),
    mean(protests_windows$followon_2d[protests_windows$is_violent]),
    mean(protests_windows$followon_3d[protests_windows$is_violent]),
    mean(protests_windows$followon_7d[protests_windows$is_violent]),
    mean(protests_windows$followon_14d[protests_windows$is_violent])
  ),
  Peaceful_Mean = c(
    mean(protests_windows$followon_1d[!protests_windows$is_violent]),
    mean(protests_windows$followon_2d[!protests_windows$is_violent]),
    mean(protests_windows$followon_3d[!protests_windows$is_violent]),
    mean(protests_windows$followon_7d[!protests_windows$is_violent]),
    mean(protests_windows$followon_14d[!protests_windows$is_violent])
  )
)
windows_comparison$Ratio <- windows_comparison$Violent_Mean / windows_comparison$Peaceful_Mean

cat("NATIONAL FOLLOW-ON BY TIME WINDOW:\n")
cat("─────────────────────────────────────────────────────────────────\n")
print(windows_comparison %>%
        mutate(Violent_Mean = round(Violent_Mean, 1),
               Peaceful_Mean = round(Peaceful_Mean, 1),
               Ratio = sprintf("%.2fx", Ratio)) %>%
        format(justify = "right"), row.names = FALSE)

cat("\nINTERPRETATION:\n")
cat("  If violence has immediate effects, ratio should be highest at short windows.\n")

# =============================================================================
# H4: WAVE POSITION - Do violent events occur at END of waves?
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   H4: WAVE POSITION                                          ║\n")
cat("║   Do violent events occur at the END of protest waves?       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Identify high-activity periods (above-median days)
median_daily <- median(daily_counts$daily_n)

# For each event, calculate position within its "wave"
# Wave = consecutive period of above-median activity
wave_data <- daily_counts %>%
  arrange(event_date) %>%
  mutate(
    above_median = daily_n > median_daily,
    wave_change = above_median != lag(above_median, default = FALSE),
    wave_id = cumsum(wave_change)
  ) %>%
  filter(above_median) %>%
  group_by(wave_id) %>%
  mutate(
    wave_start = min(event_date),
    wave_end = max(event_date),
    wave_length = as.numeric(wave_end - wave_start) + 1,
    days_into_wave = as.numeric(event_date - wave_start),
    wave_position = days_into_wave / pmax(wave_length - 1, 1)  # 0 = start, 1 = end
  ) %>%
  ungroup() %>%
  select(event_date, wave_id, wave_position, wave_length, days_into_wave)

# Join to protests
protests_wave <- protests %>%
  left_join(wave_data, by = "event_date")

# Compare wave position
wave_comparison <- protests_wave %>%
  filter(!is.na(wave_position)) %>%
  group_by(Category = case_when(
    is_violent ~ "Violent",
    TRUE ~ "Peaceful"
  )) %>%
  summarize(
    N_in_waves = n(),
    Mean_Wave_Position = mean(wave_position),
    Median_Wave_Position = median(wave_position),
    .groups = "drop"
  )

cat("WAVE POSITION BY EVENT TYPE:\n")
cat("(0 = start of wave, 1 = end of wave)\n")
cat("─────────────────────────────────────────────────────────────────\n")
print(wave_comparison %>%
        mutate(Mean_Wave_Position = round(Mean_Wave_Position, 3),
               Median_Wave_Position = round(Median_Wave_Position, 3)) %>%
        format(justify = "right"), row.names = FALSE)

# Also for fatal events
wave_fatal <- protests_wave %>%
  filter(!is.na(wave_position)) %>%
  group_by(Category = case_when(
    has_fatalities ~ "Fatal",
    TRUE ~ "Non-fatal"
  )) %>%
  summarize(
    N_in_waves = n(),
    Mean_Wave_Position = mean(wave_position),
    Median_Wave_Position = median(wave_position),
    .groups = "drop"
  )

cat("\nFATAL vs NON-FATAL:\n")
print(wave_fatal %>%
        mutate(Mean_Wave_Position = round(Mean_Wave_Position, 3),
               Median_Wave_Position = round(Median_Wave_Position, 3)) %>%
        format(justify = "right"), row.names = FALSE)

cat("\nINTERPRETATION:\n")
cat("  If violent events occur at wave END, they'd naturally have lower follow-on\n")
cat("  (the wave is about to end regardless of violence).\n")

# =============================================================================
# H5: REGRESSION WITH CONTROLS
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   H5: REGRESSION WITH CONTROLS                               ║\n")
cat("║   Does violence predict follow-on AFTER controlling for      ║\n")
cat("║   baseline activity?                                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Calculate baseline (7 days before) for each event
protests_reg <- protests_windows %>%
  left_join(daily_counts, by = "event_date") %>%
  rowwise() %>%
  mutate(
    baseline_7d = sum(protests$event_date >= event_date - 7 &
                      protests$event_date < event_date)
  ) %>%
  ungroup()

# Run regressions
cat("Model 1: followon_7d ~ is_violent\n")
model1 <- lm(followon_7d ~ is_violent, data = protests_reg)
cat(sprintf("  Coefficient on is_violent: %.2f (p = %.4f)\n",
            coef(model1)["is_violentTRUE"],
            summary(model1)$coefficients["is_violentTRUE", "Pr(>|t|)"]))

cat("\nModel 2: followon_7d ~ is_violent + baseline_7d\n")
model2 <- lm(followon_7d ~ is_violent + baseline_7d, data = protests_reg)
cat(sprintf("  Coefficient on is_violent: %.2f (p = %.4f)\n",
            coef(model2)["is_violentTRUE"],
            summary(model2)$coefficients["is_violentTRUE", "Pr(>|t|)"]))
cat(sprintf("  Coefficient on baseline_7d: %.2f (p = %.4f)\n",
            coef(model2)["baseline_7d"],
            summary(model2)$coefficients["baseline_7d", "Pr(>|t|)"]))

cat("\nModel 3: followon_7d ~ is_violent + baseline_7d + daily_n (same-day activity)\n")
model3 <- lm(followon_7d ~ is_violent + baseline_7d + daily_n, data = protests_reg)
cat(sprintf("  Coefficient on is_violent: %.2f (p = %.4f)\n",
            coef(model3)["is_violentTRUE"],
            summary(model3)$coefficients["is_violentTRUE", "Pr(>|t|)"]))

cat("\nModel 4: Include fatalities\n")
model4 <- lm(followon_7d ~ is_violent + has_fatalities + baseline_7d + daily_n, data = protests_reg)
cat(sprintf("  Coefficient on is_violent: %.2f (p = %.4f)\n",
            coef(model4)["is_violentTRUE"],
            summary(model4)$coefficients["is_violentTRUE", "Pr(>|t|)"]))
cat(sprintf("  Coefficient on has_fatalities: %.2f (p = %.4f)\n",
            coef(model4)["has_fatalitiesTRUE"],
            summary(model4)$coefficients["has_fatalitiesTRUE", "Pr(>|t|)"]))

cat("\nINTERPRETATION:\n")
cat("  If violence coefficient becomes positive after controlling for baseline,\n")
cat("  then the confounding explanation is supported.\n")

# =============================================================================
# ADDITIONAL: Check if violence/fatalities predict INCREASE in activity
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   ADDITIONAL: DOES ACTIVITY INCREASE AFTER VIOLENCE?         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Calculate change in activity (after - before)
protests_change <- protests_reg %>%
  mutate(
    change_7d = followon_7d - baseline_7d,
    pct_change = ifelse(baseline_7d > 0, 100 * (followon_7d - baseline_7d) / baseline_7d, NA)
  )

# Compare
change_comparison <- data.frame(
  Category = c("Violent", "Peaceful", "Fatal", "Non-fatal"),
  Mean_Change = c(
    mean(protests_change$change_7d[protests_change$is_violent]),
    mean(protests_change$change_7d[!protests_change$is_violent]),
    mean(protests_change$change_7d[protests_change$has_fatalities]),
    mean(protests_change$change_7d[!protests_change$has_fatalities])
  ),
  Median_Change = c(
    median(protests_change$change_7d[protests_change$is_violent]),
    median(protests_change$change_7d[!protests_change$is_violent]),
    median(protests_change$change_7d[protests_change$has_fatalities]),
    median(protests_change$change_7d[!protests_change$has_fatalities])
  ),
  Pct_Positive = c(
    100 * mean(protests_change$change_7d[protests_change$is_violent] > 0),
    100 * mean(protests_change$change_7d[!protests_change$is_violent] > 0),
    100 * mean(protests_change$change_7d[protests_change$has_fatalities] > 0),
    100 * mean(protests_change$change_7d[!protests_change$has_fatalities] > 0)
  )
)

cat("CHANGE IN ACTIVITY (7-day after - 7-day before):\n")
cat("─────────────────────────────────────────────────────────────────\n")
print(change_comparison %>%
        mutate(Mean_Change = round(Mean_Change, 1),
               Pct_Positive = sprintf("%.1f%%", Pct_Positive)) %>%
        format(justify = "right"), row.names = FALSE)

# Regression on change
cat("\nRegression: change_7d ~ is_violent + has_fatalities\n")
model_change <- lm(change_7d ~ is_violent + has_fatalities, data = protests_change)
print(summary(model_change)$coefficients)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   INVESTIGATION SUMMARY                                      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("H1 (TIMING CONFOUND):\n")
violent_sameday <- mean(protests_with_daily$daily_n[protests_with_daily$is_violent])
peaceful_sameday <- mean(protests_with_daily$daily_n[!protests_with_daily$is_violent])
cat(sprintf("  Violent events occur on days with %.1f protests (vs %.1f for peaceful)\n",
            violent_sameday, peaceful_sameday))
if (violent_sameday < peaceful_sameday) {
  cat("  → Violent events DO occur during quieter periods\n")
} else {
  cat("  → Violent events do NOT occur during quieter periods\n")
}

cat("\nH2 (LOCAL TRIGGERING):\n")
cat(sprintf("  Local 7-day follow-on: Violent=%.3f, Peaceful=%.3f\n",
            local_comparison$Mean_Local_7d[1], local_comparison$Mean_Local_7d[2]))

cat("\nH3 (SHORT-TERM EFFECTS):\n")
cat(sprintf("  1-day ratio: %.2f | 7-day ratio: %.2f\n",
            windows_comparison$Ratio[1], windows_comparison$Ratio[4]))

cat("\nH4 (WAVE POSITION):\n")
cat(sprintf("  Violent mean position: %.3f | Peaceful: %.3f\n",
            wave_comparison$Mean_Wave_Position[wave_comparison$Category == "Violent"],
            wave_comparison$Mean_Wave_Position[wave_comparison$Category == "Peaceful"]))

cat("\nH5 (REGRESSION WITH CONTROLS):\n")
cat(sprintf("  Violence coefficient (no controls): %.2f\n", coef(model1)["is_violentTRUE"]))
cat(sprintf("  Violence coefficient (with baseline control): %.2f\n", coef(model2)["is_violentTRUE"]))

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   INVESTIGATION COMPLETE                                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
