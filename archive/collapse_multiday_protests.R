# Collapse multi-day protests at the same location
# - Jakarta: requires same location AND same topic (due to coarse geocoding)
# - Elsewhere: same location only
# - Both use 3-day gap rule

library(tidyverse)

# Load original data
protests <- readRDS("protests_with_poverty.rds")
cat("Original data:", nrow(protests), "events\n")

# --- Topic detection for Jakarta ---
# Define topic patterns
topic_patterns <- list(
  labor = "labor|worker|wage|minimum wage|buruh",
  student = "student|mahasiswa|university|campus",
  corruption = "corruption|korupsi|kpk",
  religious = "religious|islam|muslim|christian|church|mosque",
  election = "election|pemilu|pilkada|vote|voting",
  land = "agrarian|land|eviction|housing",
  papua = "papua|west papua",
  environment = "environmental|environment|pollution|forest",
  gender = "woman|women|gender|feminist",
  lgbtq = "lgbtq|lgbt|gay|lesbian"
)

# Add topic columns
protests <- protests %>%
  mutate(notes_lower = tolower(notes))

for (topic in names(topic_patterns)) {
  protests[[topic]] <- grepl(topic_patterns[[topic]], protests$notes_lower, ignore.case = TRUE)
}

# Create topic signature
topic_cols <- names(topic_patterns)
protests <- protests %>%
  mutate(
    topic_signature = apply(select(., all_of(topic_cols)), 1, function(x) {
      detected <- topic_cols[x]
      if (length(detected) == 0) "unknown" else paste(sort(detected), collapse = "|")
    })
  )

# --- Create clustering keys ---
protests <- protests %>%
  mutate(
    location_id = paste(round(latitude, 4), round(longitude, 4), sep = "_"),
    event_date = as.Date(event_date),
    is_jakarta = admin1 == "Jakarta",
    # Jakarta uses location + topic; elsewhere uses location only
    collapse_key = ifelse(
      is_jakarta,
      paste(location_id, topic_signature, sep = "__"),
      location_id
    )
  ) %>%
  arrange(collapse_key, event_date)

# Function to assign cluster IDs
assign_clusters <- function(dates, gap_days = 3) {
  if (length(dates) == 1) return(1)

  dates <- sort(dates)
  cluster <- rep(1, length(dates))
  current_cluster <- 1

  for (i in 2:length(dates)) {
    if (as.numeric(dates[i] - dates[i-1]) > gap_days) {
      current_cluster <- current_cluster + 1
    }
    cluster[i] <- current_cluster
  }
  return(cluster)
}

# Assign cluster IDs within each collapse_key
protests <- protests %>%
  group_by(collapse_key) %>%
  mutate(
    cluster_id = assign_clusters(event_date, gap_days = 3)
  ) %>%
  ungroup() %>%
  mutate(
    event_cluster = paste(collapse_key, cluster_id, sep = "_cluster_")
  )

# Check how many events will be collapsed
cluster_sizes <- protests %>%
  count(event_cluster) %>%
  count(n, name = "num_clusters")

cat("\nCluster size distribution:\n")
print(cluster_sizes)

# Collapse events within each cluster
protests_collapsed <- protests %>%
  group_by(event_cluster) %>%
  summarize(
    # Date: use first event date
    event_date = min(event_date),

    # Location info (from first event)
    latitude = first(latitude),
    longitude = first(longitude),
    location = first(location),
    admin1 = first(admin1),
    admin2 = first(admin2),
    admin3 = first(admin3),

    # Event characteristics
    n_events_collapsed = n(),
    event_span_days = as.numeric(max(event_date) - min(event_date)),

    # Fatalities: sum across all events in cluster
    fatalities = sum(fatalities, na.rm = TRUE),

    # Violence: TRUE if ANY event was violent
    is_violent = any(is_violent, na.rm = TRUE),

    # Sub-event types present
    sub_event_types = paste(unique(sub_event_type), collapse = "; "),

    # Keep most severe sub_event_type for classification
    sub_event_type = case_when(
      any(sub_event_type == "Mob violence") ~ "Mob violence",
      any(sub_event_type == "Violent demonstration") ~ "Violent demonstration",
      any(sub_event_type == "Excessive force against protesters") ~ "Excessive force against protesters",
      TRUE ~ first(sub_event_type)
    ),

    # Covariates (same for location)
    population = first(population),
    log_pop = first(log_pop),
    poverty_rate = first(poverty_rate),
    poverty_decimal = first(poverty_decimal),

    # Keep year from first event
    year = first(year),

    # Track if Jakarta (for diagnostics)
    is_jakarta = first(is_jakarta),
    topic_signature = first(topic_signature),

    .groups = "drop"
  ) %>%
  # Recreate derived variables matching original data structure
  mutate(
    is_severe = fatalities > 0 | is_violent,
    event_type = ifelse(is_severe, "Severe", "Peaceful"),
    has_fatalities = fatalities > 0,
    days_since_start = as.numeric(event_date - min(event_date))
  )

cat("\nCollapsed data:", nrow(protests_collapsed), "events\n")
cat("Reduction:", nrow(protests) - nrow(protests_collapsed), "events collapsed\n")
cat("Reduction percentage:", round(100 * (1 - nrow(protests_collapsed)/nrow(protests)), 1), "%\n")

# Summary statistics
cat("\n--- Summary by event type ---\n")
cat("Original:\n")
print(table(protests$is_violent))

cat("\nCollapsed:\n")
print(table(protests_collapsed$is_severe))

# Jakarta-specific summary
cat("\n--- Jakarta-specific results ---\n")
jakarta_orig <- sum(protests$is_jakarta)
jakarta_collapsed <- sum(protests_collapsed$is_jakarta)
cat("Jakarta original:", jakarta_orig, "\n")
cat("Jakarta collapsed:", jakarta_collapsed, "\n")
cat("Jakarta reduction:", round(100 * (1 - jakarta_collapsed/jakarta_orig), 1), "%\n")

# Non-Jakarta
non_jakarta_orig <- sum(!protests$is_jakarta)
non_jakarta_collapsed <- sum(!protests_collapsed$is_jakarta)
cat("\nNon-Jakarta original:", non_jakarta_orig, "\n")
cat("Non-Jakarta collapsed:", non_jakarta_collapsed, "\n")
cat("Non-Jakarta reduction:", round(100 * (1 - non_jakarta_collapsed/non_jakarta_orig), 1), "%\n")

# Events that were multi-day
multi_day <- protests_collapsed %>% filter(n_events_collapsed > 1)
cat("\n--- Multi-day protest clusters ---\n")
cat("Number of clusters with >1 event:", nrow(multi_day), "\n")
cat("Mean events per multi-day cluster:", round(mean(multi_day$n_events_collapsed), 2), "\n")
cat("Max events in a cluster:", max(multi_day$n_events_collapsed), "\n")

# Save collapsed data (3-day rule)
saveRDS(protests_collapsed, "protests_collapsed_3day.rds")
cat("\nSaved to: protests_collapsed_3day.rds\n")

# ==============================================================================
# ADDITIONAL STEP: Daily aggregation
# ==============================================================================
# Collapse to one event per location-day to eliminate same-day clustering
# This addresses the bias identified in Hawkes model half-life estimation

cat("\n=== DAILY AGGREGATION ===\n")

# Load CPI data for merge
cpi <- readRDS("indonesia_cpi.rds") %>%
  mutate(year_month = format(date, "%Y-%m")) %>%
  select(year_month, cpi)

# Add year_month and merge CPI
protests_with_cpi <- protests_collapsed %>%
  mutate(
    date = event_date,
    year_month = format(date, "%Y-%m")
  ) %>%
  left_join(cpi, by = "year_month") %>%
  mutate(log_cpi = log(cpi))

# Aggregate to daily resolution
protests_daily <- protests_with_cpi %>%
  group_by(date, latitude, longitude) %>%
  summarise(
    # Event characteristics
    is_severe = any(is_severe),           # Severe if any event was severe
    n_events_day = n(),                    # Track how many events merged
    n_events_collapsed = sum(n_events_collapsed),  # Total original events
    fatalities = sum(fatalities),

    # Location info
    location = first(location),
    admin1 = first(admin1),
    admin2 = first(admin2),
    admin3 = first(admin3),

    # Covariates
    log_pop = first(log_pop),
    population = first(population),
    poverty_decimal = first(poverty_decimal),
    poverty_rate = first(poverty_rate),
    log_cpi = first(log_cpi),
    cpi = first(cpi),
    year_month = first(year_month),
    year = first(year),

    .groups = "drop"
  ) %>%
  arrange(date) %>%
  mutate(
    event_type = ifelse(is_severe, "Severe", "Peaceful"),
    days_since_start = as.numeric(date - min(date))
  )

cat("Daily aggregated data:", nrow(protests_daily), "events\n")
cat("Reduction from 3-day collapse:", nrow(protests_collapsed) - nrow(protests_daily), "events\n")
cat("Reduction percentage:", round(100 * (1 - nrow(protests_daily)/nrow(protests_collapsed)), 1), "%\n")

# Check same-day clustering is eliminated
same_day_check <- protests_daily %>%
  group_by(date, latitude, longitude) %>%
  tally() %>%
  filter(n > 1)

cat("\nSame-location same-day clusters remaining:", nrow(same_day_check), "\n")

# Summary by type
cat("\n--- Daily data by type ---\n")
print(table(Severe = protests_daily$is_severe))

# Save daily data
saveRDS(protests_daily, "protests_daily.rds")
cat("\nSaved to: protests_daily.rds\n")
