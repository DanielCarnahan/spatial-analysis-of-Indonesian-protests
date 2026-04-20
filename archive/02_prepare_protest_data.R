# Phase 1: Data Preparation for Spatial Point Process Analysis
# Prepare Indonesian Protest Data for Contagion Modeling
# Author: Daniel Carnahan
# Date: 2025-10-27

library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(sf)

# Load ACLED data
cat("Loading ACLED data...\n")
acled_raw <- read_excel("ACLED Data_2025-10-27_Indonesia20142024_OFFICIAL.xlsx")

# Extract column names from first row
col_names <- as.character(acled_raw[1, ])
acled <- acled_raw[-1, ]
names(acled) <- col_names

# Convert columns to appropriate types
acled <- acled %>%
  mutate(
    event_date = as.Date(event_date),
    year = as.integer(year),
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude),
    fatalities = as.integer(fatalities),
    time_precision = as.integer(time_precision)
  )

cat("Total events:", nrow(acled), "\n\n")

# ====================================================================
# 1. FILTER PROTEST EVENTS
# ====================================================================

cat("=== FILTERING PROTEST EVENTS ===\n")

# Primary filter: Protests and Riots (potential contagion events)
protests <- acled %>%
  filter(event_type %in% c("Protests", "Riots"))

cat("Protest events:", nrow(protests), "\n")
cat("  - Protests:", sum(protests$event_type == "Protests"), "\n")
cat("  - Riots:", sum(protests$event_type == "Riots"), "\n\n")

# Create sub-categories for analysis
protests <- protests %>%
  mutate(
    # Classify by violence level
    is_violent = sub_event_type %in% c("Violent demonstration",
                                        "Protest with intervention",
                                        "Mob violence",
                                        "Looting/property destruction",
                                        "Excessive force against protesters"),

    # Classify by peacefuleness
    is_peaceful = sub_event_type == "Peaceful protest",

    # Has state intervention
    state_intervention = grepl("Police|State forces", interaction),

    # Resulted in casualties
    has_fatalities = fatalities > 0
  )

cat("Violence classification:\n")
cat("  - Peaceful:", sum(protests$is_peaceful), "\n")
cat("  - Violent:", sum(protests$is_violent), "\n")
cat("  - With fatalities:", sum(protests$has_fatalities), "\n\n")

# ====================================================================
# 2. TEMPORAL PROCESSING
# ====================================================================

cat("=== TEMPORAL PROCESSING ===\n")

protests <- protests %>%
  mutate(
    # Extract temporal components
    month = month(event_date),
    day_of_year = yday(event_date),
    week = week(event_date),

    # Time since start (for point process models)
    days_since_start = as.numeric(event_date - min(event_date)),

    # Year-month for aggregation
    year_month = floor_date(event_date, "month")
  ) %>%
  arrange(event_date)  # Critical: sort by time for point process

cat("Temporal range:", as.character(min(protests$event_date)), "to",
    as.character(max(protests$event_date)), "\n")
cat("Total observation period:", max(protests$days_since_start), "days\n\n")

# ====================================================================
# 3. SPATIAL PROCESSING
# ====================================================================

cat("=== SPATIAL PROCESSING ===\n")

# Check for missing coordinates (should be 0)
missing_coords <- sum(is.na(protests$latitude) | is.na(protests$longitude))
cat("Events with missing coordinates:", missing_coords, "\n")

# Remove any events with missing coordinates (if any)
if(missing_coords > 0) {
  protests <- protests %>%
    filter(!is.na(latitude) & !is.na(longitude))
  cat("Removed", missing_coords, "events\n")
}

# Spatial extent
cat("Spatial extent:\n")
cat("  Latitude:", round(min(protests$latitude), 3), "to",
    round(max(protests$latitude), 3), "\n")
cat("  Longitude:", round(min(protests$longitude), 3), "to",
    round(max(protests$longitude), 3), "\n\n")

# Convert to sf object for spatial operations
protests_sf <- st_as_sf(protests,
                        coords = c("longitude", "latitude"),
                        crs = 4326)  # WGS84

# ====================================================================
# 4. REGIONAL ANALYSIS
# ====================================================================

cat("=== REGIONAL DISTRIBUTION ===\n")

regional_summary <- protests %>%
  group_by(admin1) %>%
  summarise(
    n_events = n(),
    n_violent = sum(is_violent),
    n_peaceful = sum(is_peaceful),
    n_fatalities = sum(fatalities),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_events))

cat("Top 15 provinces by number of protest events:\n")
print(regional_summary %>% head(15))

# ====================================================================
# 5. TEMPORAL PATTERNS
# ====================================================================

cat("\n=== TEMPORAL PATTERNS ===\n")

# Events per year
yearly_summary <- protests %>%
  group_by(year) %>%
  summarise(
    n_events = n(),
    n_violent = sum(is_violent),
    n_peaceful = sum(is_peaceful),
    .groups = 'drop'
  )

cat("\nEvents per year:\n")
print(yearly_summary)

# ====================================================================
# 6. SAVE PREPARED DATA
# ====================================================================

cat("\n=== SAVING PREPARED DATA ===\n")

# Save as RDS (preserves R data types)
saveRDS(protests, "protests_prepared.rds")
cat("Saved: protests_prepared.rds\n")

# Save spatial version
saveRDS(protests_sf, "protests_sf.rds")
cat("Saved: protests_sf.rds\n")

# Save regional summary
write.csv(regional_summary, "regional_summary.csv", row.names = FALSE)
cat("Saved: regional_summary.csv\n")

# Save yearly summary
write.csv(yearly_summary, "yearly_summary.csv", row.names = FALSE)
cat("Saved: yearly_summary.csv\n")

# ====================================================================
# 7. SUMMARY STATISTICS FOR POINT PROCESS
# ====================================================================

cat("\n=== POINT PROCESS SUMMARY ===\n")
cat("Total events (n):", nrow(protests), "\n")
cat("Temporal extent (T):", max(protests$days_since_start), "days\n")
cat("Spatial extent:\n")
cat("  Area (rough box):",
    round((max(protests$longitude) - min(protests$longitude)) *
          (max(protests$latitude) - min(protests$latitude)) * 111^2),
    "km²\n")
cat("Intensity (λ = n/T):",
    round(nrow(protests) / max(protests$days_since_start), 4),
    "events/day\n")
cat("Provinces:", length(unique(protests$admin1)), "\n")
cat("Districts:", length(unique(protests$admin2)), "\n")
cat("Unique locations:", length(unique(paste(protests$latitude, protests$longitude))), "\n")

cat("\n=== PREPARATION COMPLETE ===\n")
