# Explore ACLED Indonesia Data Structure
# Author: Daniel Carnahan
# Date: 2025-10-27

library(readxl)
library(dplyr)
library(lubridate)

# Load ACLED data
# The first row contains the actual column names
acled_raw <- read_excel("ACLED Data_2025-10-27_Indonesia20142024_OFFICIAL.xlsx")

# Extract column names from first row and use them as headers
col_names <- as.character(acled_raw[1, ])
acled <- acled_raw[-1, ]  # Remove the header row
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

# Basic structure
cat("=== DATA STRUCTURE ===\n")
cat("Dimensions:", nrow(acled), "rows x", ncol(acled), "columns\n\n")

cat("=== COLUMN NAMES ===\n")
print(names(acled))

cat("\n=== FIRST FEW ROWS ===\n")
print(head(acled, 3))

cat("\n=== DATA SUMMARY ===\n")
print(str(acled))

# Check for key columns
cat("\n=== KEY VARIABLES FOR POINT PROCESS ANALYSIS ===\n")

# Event types
if("event_type" %in% names(acled)) {
  cat("\nEvent Types:\n")
  print(table(acled$event_type))
}

# Sub-event types
if("sub_event_type" %in% names(acled)) {
  cat("\nSub-Event Types:\n")
  print(table(acled$sub_event_type))
}

# Temporal range
if("event_date" %in% names(acled)) {
  cat("\nTemporal Range:\n")
  cat("Start:", as.character(min(acled$event_date, na.rm = TRUE)), "\n")
  cat("End:", as.character(max(acled$event_date, na.rm = TRUE)), "\n")
  cat("Total days:", as.numeric(difftime(max(acled$event_date, na.rm = TRUE),
                                          min(acled$event_date, na.rm = TRUE),
                                          units = "days")), "\n")
}

# Spatial coverage
if("latitude" %in% names(acled) && "longitude" %in% names(acled)) {
  cat("\nSpatial Coverage:\n")
  cat("Latitude range:", min(acled$latitude, na.rm = TRUE), "to",
      max(acled$latitude, na.rm = TRUE), "\n")
  cat("Longitude range:", min(acled$longitude, na.rm = TRUE), "to",
      max(acled$longitude, na.rm = TRUE), "\n")
  cat("Missing coordinates:", sum(is.na(acled$latitude) | is.na(acled$longitude)), "\n")
}

# Regional distribution
if("admin1" %in% names(acled)) {
  cat("\nTop 10 Provinces (admin1):\n")
  print(head(sort(table(acled$admin1), decreasing = TRUE), 10))
}

if("admin2" %in% names(acled)) {
  cat("\nNumber of unique districts (admin2):", length(unique(acled$admin2)), "\n")
}

# Check for fatalities and other intensity measures
if("fatalities" %in% names(acled)) {
  cat("\nFatalities:\n")
  cat("Total events with fatalities:", sum(acled$fatalities > 0, na.rm = TRUE), "\n")
  cat("Total fatalities:", sum(acled$fatalities, na.rm = TRUE), "\n")
}

cat("\n=== DONE ===\n")
