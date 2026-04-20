# =============================================================================
# 38_trivariate_classification.R
# Create mutually exclusive trivariate classification for Hawkes model
# Event types: Fatal (F), Violent non-fatal (V), Peaceful (P)
# =============================================================================

library(tidyverse)

cat("=============================================================================\n")
cat("TRIVARIATE EVENT CLASSIFICATION\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
protests <- readRDS("protests_prepared.rds")
cat(sprintf("Loaded %d protest events\n\n", nrow(protests)))

# -----------------------------------------------------------------------------
# Define violence indicator (consistent with bivariate model)
# -----------------------------------------------------------------------------
protests <- protests %>%
  mutate(
    is_violent = sub_event_type %in% c("Mob violence", "Violent demonstration") |
                 sub_event_type == "Excessive force against protesters"
  )

# -----------------------------------------------------------------------------
# Create mutually exclusive trivariate classification
# Priority: Fatal > Violent > Peaceful
# -----------------------------------------------------------------------------
protests <- protests %>%
  mutate(
    event_type_3 = case_when(
      fatalities > 0 ~ "F",                           # Fatal (highest priority)
      is_violent & fatalities == 0 ~ "V",             # Violent non-fatal
      TRUE ~ "P"                                       # Peaceful
    ),
    event_type_3 = factor(event_type_3, levels = c("F", "V", "P"))
  )

# -----------------------------------------------------------------------------
# Validate classification
# -----------------------------------------------------------------------------
cat("=== TRIVARIATE CLASSIFICATION ===\n\n")

# Count by type
type_counts <- protests %>%
  count(event_type_3) %>%
  mutate(
    pct = round(100 * n / sum(n), 2),
    label = case_when(
      event_type_3 == "F" ~ "Fatal",
      event_type_3 == "V" ~ "Violent (non-fatal)",
      event_type_3 == "P" ~ "Peaceful"
    )
  )

for (i in 1:nrow(type_counts)) {
  cat(sprintf("  %s (%s): %d events (%.1f%%)\n",
              type_counts$label[i],
              type_counts$event_type_3[i],
              type_counts$n[i],
              type_counts$pct[i]))
}

cat(sprintf("\n  Total: %d events\n", sum(type_counts$n)))

# Verify mutual exclusivity
cat("\n=== VALIDATION ===\n")
total_check <- sum(type_counts$n)
cat(sprintf("  Sum of categories: %d\n", total_check))
cat(sprintf("  Original dataset: %d\n", nrow(protests)))
cat(sprintf("  Match: %s\n", ifelse(total_check == nrow(protests), "YES", "NO")))

# -----------------------------------------------------------------------------
# Cross-tabulation with sub_event_type for transparency
# -----------------------------------------------------------------------------
cat("\n=== CROSS-TABULATION WITH SUB_EVENT_TYPE ===\n\n")
cross_tab <- protests %>%
  count(event_type_3, sub_event_type) %>%
  pivot_wider(names_from = event_type_3, values_from = n, values_fill = 0)
print(cross_tab)

# -----------------------------------------------------------------------------
# Check fatality distribution
# -----------------------------------------------------------------------------
cat("\n=== FATALITY DISTRIBUTION ===\n")
fatal_events <- protests %>% filter(event_type_3 == "F")
cat(sprintf("  Mean fatalities (fatal events): %.1f\n", mean(fatal_events$fatalities)))
cat(sprintf("  Median fatalities (fatal events): %.0f\n", median(fatal_events$fatalities)))
cat(sprintf("  Max fatalities: %d\n", max(fatal_events$fatalities)))

# Distribution by count
fatal_dist <- fatal_events %>%
  mutate(fat_bin = case_when(
    fatalities == 1 ~ "1",
    fatalities == 2 ~ "2",
    fatalities <= 5 ~ "3-5",
    fatalities <= 10 ~ "6-10",
    TRUE ~ ">10"
  )) %>%
  count(fat_bin) %>%
  mutate(pct = round(100 * n / sum(n), 1))

cat("\n  Fatality distribution:\n")
for (i in 1:nrow(fatal_dist)) {
  cat(sprintf("    %s fatalities: %d events (%.1f%%)\n",
              fatal_dist$fat_bin[i], fatal_dist$n[i], fatal_dist$pct[i]))
}

# -----------------------------------------------------------------------------
# Temporal distribution
# -----------------------------------------------------------------------------
cat("\n=== TEMPORAL DISTRIBUTION ===\n")
yearly_by_type <- protests %>%
  mutate(year = lubridate::year(event_date)) %>%
  count(year, event_type_3) %>%
  pivot_wider(names_from = event_type_3, values_from = n, values_fill = 0)
print(yearly_by_type)

# -----------------------------------------------------------------------------
# Create numeric type indicator for model
# -----------------------------------------------------------------------------
protests <- protests %>%
  mutate(
    type_numeric = case_when(
      event_type_3 == "F" ~ 1L,
      event_type_3 == "V" ~ 2L,
      event_type_3 == "P" ~ 3L
    )
  )

# -----------------------------------------------------------------------------
# Save prepared data
# -----------------------------------------------------------------------------
cat("\n=== SAVING DATA ===\n")
saveRDS(protests, "protests_trivariate.rds")
cat("  Saved: protests_trivariate.rds\n")

# Also save summary statistics for reference
trivariate_summary <- list(
  n_total = nrow(protests),
  n_fatal = sum(protests$event_type_3 == "F"),
  n_violent = sum(protests$event_type_3 == "V"),
  n_peaceful = sum(protests$event_type_3 == "P"),
  pct_fatal = mean(protests$event_type_3 == "F") * 100,
  pct_violent = mean(protests$event_type_3 == "V") * 100,
  pct_peaceful = mean(protests$event_type_3 == "P") * 100,
  date_range = range(protests$event_date),
  yearly_counts = yearly_by_type
)
saveRDS(trivariate_summary, "trivariate_summary.rds")
cat("  Saved: trivariate_summary.rds\n")

cat("\n=== DONE ===\n")
