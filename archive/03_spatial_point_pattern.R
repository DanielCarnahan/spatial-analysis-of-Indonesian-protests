# Phase 2: Create Spatial Point Pattern Objects
# Setup for spatial point process analysis
# Author: Daniel Carnahan
# Date: 2025-10-27

library(spatstat)
library(dplyr)
library(sf)
library(ggplot2)
library(lubridate)

# Load prepared data
cat("Loading prepared protest data...\n")
protests <- readRDS("protests_prepared.rds")
protests_sf <- readRDS("protests_sf.rds")

cat("Loaded", nrow(protests), "protest events\n\n")

# ====================================================================
# 1. CREATE SPATIAL POINT PATTERN (ppp) OBJECT
# ====================================================================

cat("=== CREATING SPATIAL POINT PATTERN OBJECT ===\n")

# Extract coordinates
coords <- cbind(protests$longitude, protests$latitude)

# Define observation window (bounding box for Indonesia)
# Add a small buffer to ensure all points are inside
lon_range <- range(protests$longitude)
lat_range <- range(protests$latitude)

buffer <- 0.1  # degrees
lon_window <- c(lon_range[1] - buffer, lon_range[2] + buffer)
lat_window <- c(lat_range[1] - buffer, lat_range[2] + buffer)

# Create window
protest_window <- owin(xrange = lon_window, yrange = lat_window)

cat("Observation window:\n")
cat("  Longitude:", round(lon_window, 2), "\n")
cat("  Latitude:", round(lat_window, 2), "\n")
cat("  Area:", round(area(protest_window), 2), "square degrees\n\n")

# Create point pattern object
protest_ppp <- ppp(x = protests$longitude,
                   y = protests$latitude,
                   window = protest_window,
                   marks = protests %>%
                     select(event_date, event_type, sub_event_type,
                            is_violent, is_peaceful, has_fatalities,
                            fatalities, admin1, days_since_start))

cat("Point pattern object created\n")
cat("  Number of points:", npoints(protest_ppp), "\n")
cat("  Marked:", is.marked(protest_ppp), "\n\n")

# ====================================================================
# 2. CREATE SPACE-TIME POINT PATTERN
# ====================================================================

cat("=== CREATING SPACE-TIME POINT PATTERN ===\n")

# For self-exciting/Hawkes processes, we need space-time structure
# Create data frame for spatiotemporal analysis
protest_stpp <- data.frame(
  x = protests$longitude,
  y = protests$latitude,
  t = protests$days_since_start,  # Time in days since start
  date = protests$event_date,
  is_violent = protests$is_violent,
  is_peaceful = protests$is_peaceful,
  fatalities = protests$fatalities,
  admin1 = protests$admin1
)

# Sort by time (critical for self-exciting models)
protest_stpp <- protest_stpp %>% arrange(t)

cat("Space-time point pattern created\n")
cat("  Spatial dimension:", nrow(protest_stpp), "events\n")
cat("  Temporal range: 0 to", max(protest_stpp$t), "days\n")
cat("  Time span:", round(max(protest_stpp$t) / 365.25, 2), "years\n\n")

# ====================================================================
# 3. BASIC INTENSITY ESTIMATION
# ====================================================================

cat("=== BASIC INTENSITY ESTIMATION ===\n")

# Overall intensity (events per unit area per unit time)
overall_intensity <- npoints(protest_ppp) / area(protest_window)
cat("Overall spatial intensity:", round(overall_intensity, 4),
    "events/square degree\n")

# Temporal intensity
temporal_intensity <- nrow(protests) / max(protest_stpp$t)
cat("Temporal intensity:", round(temporal_intensity, 4), "events/day\n")

# Convert to more interpretable units (events per 100km² per month)
# 1 degree ≈ 111 km
km2_per_deg2 <- 111^2
months <- max(protest_stpp$t) / 30.44
spatial_intensity_km <- overall_intensity / km2_per_deg2 * 100  # per 100 km²
temporal_intensity_month <- temporal_intensity * 30.44

cat("\nInterpretable intensities:\n")
cat("  Spatial:", round(spatial_intensity_km, 4), "events/100km²\n")
cat("  Temporal:", round(temporal_intensity_month, 2), "events/month\n\n")

# ====================================================================
# 4. SUBSET BY PROTEST TYPE
# ====================================================================

cat("=== CREATING SUBSETS BY TYPE ===\n")

# Peaceful protests only
peaceful_ppp <- subset(protest_ppp, is_peaceful == TRUE)
cat("Peaceful protests:", npoints(peaceful_ppp), "\n")

# Violent protests/riots
violent_ppp <- subset(protest_ppp, is_violent == TRUE)
cat("Violent protests/riots:", npoints(violent_ppp), "\n")

# With fatalities
fatal_ppp <- subset(protest_ppp, has_fatalities == TRUE)
cat("Protests with fatalities:", npoints(fatal_ppp), "\n\n")

# ====================================================================
# 5. TEMPORAL SUBSETS (for testing contagion over time)
# ====================================================================

cat("=== CREATING TEMPORAL SUBSETS ===\n")

# Split data into periods
protests <- protests %>%
  mutate(
    period = case_when(
      year <= 2017 ~ "2015-2017",
      year <= 2019 ~ "2018-2019",
      year <= 2021 ~ "2020-2021",
      TRUE ~ "2022-2024"
    )
  )

period_counts <- table(protests$period)
cat("Events by period:\n")
print(period_counts)
cat("\n")

# ====================================================================
# 6. REGIONAL SUBSETS
# ====================================================================

cat("=== CREATING REGIONAL SUBSETS ===\n")

# Major protest regions
top_regions <- protests %>%
  count(admin1, sort = TRUE) %>%
  head(5) %>%
  pull(admin1)

cat("Top 5 provinces for sub-analysis:\n")
print(top_regions)
cat("\n")

# Create subset for Jakarta (highest intensity)
jakarta_protests <- protests %>% filter(admin1 == "Jakarta")
jakarta_ppp <- ppp(
  x = jakarta_protests$longitude,
  y = jakarta_protests$latitude,
  window = owin(
    xrange = range(jakarta_protests$longitude) + c(-0.1, 0.1),
    yrange = range(jakarta_protests$latitude) + c(-0.1, 0.1)
  ),
  marks = jakarta_protests %>%
    select(event_date, is_violent, is_peaceful, days_since_start)
)

cat("Jakarta subset:", npoints(jakarta_ppp), "events\n")
cat("  Intensity:", round(npoints(jakarta_ppp) / area(jakarta_ppp$window), 4),
    "events/square degree\n\n")

# ====================================================================
# 7. SAVE POINT PATTERN OBJECTS
# ====================================================================

cat("=== SAVING POINT PATTERN OBJECTS ===\n")

# Save main objects
saveRDS(protest_ppp, "protest_ppp.rds")
cat("Saved: protest_ppp.rds\n")

saveRDS(protest_stpp, "protest_stpp.rds")
cat("Saved: protest_stpp.rds\n")

# Save subsets
saveRDS(peaceful_ppp, "peaceful_ppp.rds")
saveRDS(violent_ppp, "violent_ppp.rds")
saveRDS(fatal_ppp, "fatal_ppp.rds")
saveRDS(jakarta_ppp, "jakarta_ppp.rds")
cat("Saved: Type-specific and regional subsets\n")

# ====================================================================
# 8. SUMMARY FOR NEXT STEPS
# ====================================================================

cat("\n=== SUMMARY ===\n")
cat("Point pattern objects ready for analysis:\n")
cat("  - protest_ppp: All", npoints(protest_ppp), "protest events\n")
cat("  - protest_stpp: Space-time data frame\n")
cat("  - peaceful_ppp:", npoints(peaceful_ppp), "peaceful protests\n")
cat("  - violent_ppp:", npoints(violent_ppp), "violent protests\n")
cat("  - fatal_ppp:", npoints(fatal_ppp), "fatal protests\n")
cat("  - jakarta_ppp:", npoints(jakarta_ppp), "Jakarta protests\n")

cat("\nNext steps:\n")
cat("  1. Kernel density estimation\n")
cat("  2. Clustering analysis (Ripley's K, L-functions)\n")
cat("  3. Space-time interaction tests\n")
cat("  4. Self-exciting point process models\n")

cat("\n=== COMPLETE ===\n")
