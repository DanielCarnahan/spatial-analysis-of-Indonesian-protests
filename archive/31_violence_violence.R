################################################################################
#   VIOLENCE TRIGGERS VIOLENCE - DETAILED ANALYSIS
################################################################################

library(tidyverse)
library(data.table)

protests <- readRDS("protests_prepared.rds")

protests <- protests %>%
  mutate(
    is_violent = sub_event_type %in% c("Mob violence", "Violent demonstration") |
                 sub_event_type == "Excessive force against protesters",
    has_fatalities = fatalities > 0
  )

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   VIOLENCE TRIGGERS VIOLENCE - DETAILED ANALYSIS            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# More efficient: Use data.table for the rowwise operation
protests_dt <- as.data.table(protests)
protests_dt[, event_date := as.Date(event_date)]

cat("Calculating split follow-on (violent vs peaceful)...\n")

# Pre-calculate for efficiency
all_violent <- protests_dt[is_violent == TRUE, event_date]
all_peaceful <- protests_dt[is_violent == FALSE, event_date]

# For each event, count violent and peaceful follow-on
count_split_followon <- function(dates, violent_dates, peaceful_dates, window = 7) {
  n <- length(dates)
  violent_followon <- numeric(n)
  peaceful_followon <- numeric(n)

  dates_num <- as.numeric(dates)
  violent_num <- as.numeric(violent_dates)
  peaceful_num <- as.numeric(peaceful_dates)

  for (i in 1:n) {
    violent_followon[i] <- sum(violent_num > dates_num[i] &
                                violent_num <= dates_num[i] + window)
    peaceful_followon[i] <- sum(peaceful_num > dates_num[i] &
                                 peaceful_num <= dates_num[i] + window)
  }
  return(list(violent = violent_followon, peaceful = peaceful_followon))
}

result <- count_split_followon(protests_dt$event_date, all_violent, all_peaceful, 7)
protests_dt[, violent_followon_7d := result$violent]
protests_dt[, peaceful_followon_7d := result$peaceful]
protests_dt[, total_followon_7d := violent_followon_7d + peaceful_followon_7d]
protests_dt[, violent_share := ifelse(total_followon_7d > 0,
                                       violent_followon_7d / total_followon_7d, NA)]

cat("Done.\n\n")

cat("FOLLOW-ON COMPOSITION BY EVENT TYPE:\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Compare
comp <- data.frame(
  Category = c("After Violent", "After Peaceful", "After Fatal", "After Non-fatal"),
  Mean_Violent_Followon = c(
    mean(protests_dt$violent_followon_7d[protests_dt$is_violent]),
    mean(protests_dt$violent_followon_7d[!protests_dt$is_violent]),
    mean(protests_dt$violent_followon_7d[protests_dt$has_fatalities]),
    mean(protests_dt$violent_followon_7d[!protests_dt$has_fatalities])
  ),
  Mean_Peaceful_Followon = c(
    mean(protests_dt$peaceful_followon_7d[protests_dt$is_violent]),
    mean(protests_dt$peaceful_followon_7d[!protests_dt$is_violent]),
    mean(protests_dt$peaceful_followon_7d[protests_dt$has_fatalities]),
    mean(protests_dt$peaceful_followon_7d[!protests_dt$has_fatalities])
  ),
  Mean_Violent_Share = c(
    mean(protests_dt$violent_share[protests_dt$is_violent], na.rm = TRUE),
    mean(protests_dt$violent_share[!protests_dt$is_violent], na.rm = TRUE),
    mean(protests_dt$violent_share[protests_dt$has_fatalities], na.rm = TRUE),
    mean(protests_dt$violent_share[!protests_dt$has_fatalities], na.rm = TRUE)
  )
)

print(comp %>%
        mutate(Mean_Violent_Followon = round(Mean_Violent_Followon, 2),
               Mean_Peaceful_Followon = round(Mean_Peaceful_Followon, 2),
               Mean_Violent_Share = sprintf("%.1f%%", 100 * Mean_Violent_Share)) %>%
        as.data.frame(), row.names = FALSE)

cat("\nRATIO ANALYSIS:\n")
cat(sprintf("  Violent follow-on: after violent = %.2f | after peaceful = %.2f | ratio = %.2fx\n",
            comp$Mean_Violent_Followon[1], comp$Mean_Violent_Followon[2],
            comp$Mean_Violent_Followon[1] / comp$Mean_Violent_Followon[2]))
cat(sprintf("  Peaceful follow-on: after violent = %.2f | after peaceful = %.2f | ratio = %.2fx\n",
            comp$Mean_Peaceful_Followon[1], comp$Mean_Peaceful_Followon[2],
            comp$Mean_Peaceful_Followon[1] / comp$Mean_Peaceful_Followon[2]))
cat(sprintf("  Total follow-on: after violent = %.2f | after peaceful = %.2f | ratio = %.2fx\n",
            comp$Mean_Violent_Followon[1] + comp$Mean_Peaceful_Followon[1],
            comp$Mean_Violent_Followon[2] + comp$Mean_Peaceful_Followon[2],
            (comp$Mean_Violent_Followon[1] + comp$Mean_Peaceful_Followon[1]) /
            (comp$Mean_Violent_Followon[2] + comp$Mean_Peaceful_Followon[2])))

cat("\nSTATISTICAL TESTS:\n")
cat("─────────────────────────────────────────────────────────────────\n")

test_violent <- wilcox.test(
  protests_dt$violent_followon_7d[protests_dt$is_violent],
  protests_dt$violent_followon_7d[!protests_dt$is_violent]
)
cat(sprintf("Violent follow-on (violent vs peaceful events): p = %.2e\n", test_violent$p.value))

test_peaceful <- wilcox.test(
  protests_dt$peaceful_followon_7d[protests_dt$is_violent],
  protests_dt$peaceful_followon_7d[!protests_dt$is_violent]
)
cat(sprintf("Peaceful follow-on (violent vs peaceful events): p = %.2e\n", test_peaceful$p.value))

# =============================================================================
# CHECK: IS THE MODEL FITTING VIOLENCE → VIOLENCE?
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   KEY QUESTION: WHAT IS THE MODEL ACTUALLY FITTING?         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# The Hawkes model treats each event as potentially triggering future events
# The violence mark affects how much triggering intensity a parent contributes
#
# But the model doesn't distinguish between violent and peaceful offspring!
# It just models total event rate.
#
# So if violent events trigger MORE violent events but FEWER peaceful events,
# and if the net is fewer total events...
# Why does the model give a positive violence coefficient?

# Possibility 1: The model sees violence clustering together and interprets
# it as violence causing future violence (which happens to include violence)

# Let's test: Do violent events cluster TEMPORALLY?
cat("TEMPORAL CLUSTERING OF VIOLENCE:\n")
cat("─────────────────────────────────────────────────────────────────\n")

# For each day, calculate: was there violence yesterday?
daily_violence <- protests_dt %>%
  group_by(event_date) %>%
  summarize(
    n_events = .N,
    n_violent = sum(is_violent),
    has_violence = any(is_violent),
    .groups = "drop"
  ) %>%
  arrange(event_date) %>%
  mutate(
    violence_yesterday = lag(has_violence, 1, default = FALSE),
    violence_2days_ago = lag(has_violence, 2, default = FALSE)
  )

cat("Probability of violence today given yesterday:\n")
cat(sprintf("  If violence yesterday: %.1f%% chance of violence today\n",
            100 * mean(daily_violence$has_violence[daily_violence$violence_yesterday])))
cat(sprintf("  If NO violence yesterday: %.1f%% chance of violence today\n",
            100 * mean(daily_violence$has_violence[!daily_violence$violence_yesterday])))

# Chi-square test
test_cluster <- chisq.test(table(daily_violence$violence_yesterday, daily_violence$has_violence))
cat(sprintf("  χ² = %.1f, p = %.2e\n", test_cluster$statistic, test_cluster$p.value))

cat("\nNumber of events today given violence yesterday:\n")
cat(sprintf("  If violence yesterday: mean = %.1f events today\n",
            mean(daily_violence$n_events[daily_violence$violence_yesterday])))
cat(sprintf("  If NO violence yesterday: mean = %.1f events today\n",
            mean(daily_violence$n_events[!daily_violence$violence_yesterday])))

# Wilcoxon test
test_n_events <- wilcox.test(
  daily_violence$n_events[daily_violence$violence_yesterday],
  daily_violence$n_events[!daily_violence$violence_yesterday]
)
cat(sprintf("  Wilcoxon p = %.2e\n", test_n_events$p.value))

# =============================================================================
# THE HAWKES MODEL'S PERSPECTIVE
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   THE HAWKES MODEL'S PERSPECTIVE                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("The Hawkes model computes:\n")
cat("  λ(t) = μ(t) + Σ α(marks_j) * exp(-β * (t - t_j))\n\n")

cat("Key insight: The model asks 'how much does each parent event\n")
cat("contribute to the intensity at time t?'\n\n")

cat("If violent events cluster TEMPORALLY (violence → violence),\n")
cat("the model will attribute higher triggering to violent events\n")
cat("BECAUSE violent events are followed by more violent events.\n\n")

cat("But violent events may NOT increase TOTAL event count.\n")
cat("They just shift the composition toward violence.\n\n")

# Calculate: conditioning on high activity, what predicts violence?
cat("ANALYSIS: What predicts violence occurring?\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Merge daily violence data
daily_analysis <- daily_violence %>%
  mutate(
    lag1_n_events = lag(n_events, 1, default = 0),
    lag1_n_violent = lag(n_violent, 1, default = 0),
    lag1_violent_share = ifelse(lag1_n_events > 0, lag1_n_violent / lag1_n_events, 0)
  )

# Logistic regression: does past violence predict today's violence?
model_violence <- glm(has_violence ~ lag1_n_events + lag1_n_violent,
                      data = daily_analysis, family = binomial)
cat("\nLogistic regression: P(violence today) ~ lag1_n_events + lag1_n_violent\n")
print(summary(model_violence)$coefficients)

cat("\nINTERPRETATION:\n")
cat("  If lag1_n_violent has positive coefficient: past violence → today's violence\n")
cat("  This temporal clustering is what the Hawkes model detects.\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   SUMMARY: WHY MODEL AND RAW DATA DISAGREE                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("1. Violent events trigger MORE violent follow-on (ratio %.2fx)\n",
    comp$Mean_Violent_Followon[1] / comp$Mean_Violent_Followon[2])
cat("2. Violent events trigger LESS peaceful follow-on (ratio %.2fx)\n",
    comp$Mean_Peaceful_Followon[1] / comp$Mean_Peaceful_Followon[2])
cat("3. Net effect: fewer total follow-on (ratio %.2fx)\n",
    (comp$Mean_Violent_Followon[1] + comp$Mean_Peaceful_Followon[1]) /
    (comp$Mean_Violent_Followon[2] + comp$Mean_Peaceful_Followon[2]))
cat("\n")
cat("The Hawkes model sees: violence clustering together in time\n")
cat("It interprets this as: violent events trigger future events\n")
cat("(which happen to include the violent ones that cluster)\n")
cat("\n")
cat("But the raw data shows: total event count is LOWER after violence\n")
cat("\n")
cat("Resolution: The model's α coefficient represents RELATIVE\n")
cat("triggering within the self-exciting process, not absolute\n")
cat("event generation rates.\n")
