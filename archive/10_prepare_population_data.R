# Data Preparation: Match Population to Protest Data (WITH SPATIAL MATCHING)
# Author: Daniel Carnahan
# Date: 2025-11-04
#
# This script:
# 1. Loads district-level population data (2014-2020) - ALL 514 districts
# 2. Uses GADM shapefile for spatial matching of protests to districts
# 3. Converts GADM names to population CSV format
# 4. Extrapolates population to 2021-2024
# 5. Joins to protest data and creates log_pop variable
#
# KEY FIX: Distinguishes Kota (city) from Kabupaten (regency) using coordinates
#   - 27 districts have both versions (e.g., "Bandung, Kota" vs "Bandung, Kab.")
#   - ACLED only records "Bandung" (ambiguous)
#   - Spatial join assigns correct population based on lat/lon

library(dplyr)
library(tidyr)
library(sf)
library(stringdist)

cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  DATA PREPARATION: POPULATION MATCHING                        ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

# ====================================================================
# 1. LOAD POPULATION DATA
# ====================================================================

cat("=== STEP 1: LOADING POPULATION DATA ===\n\n")

pop_raw <- read.csv("district_level_population_2014_2020.csv",
                    stringsAsFactors = FALSE)

cat(sprintf("Loaded %d rows from population CSV\n", nrow(pop_raw)))
cat(sprintf("Columns: %s\n\n", paste(names(pop_raw), collapse=", ")))

# Filter for total population series
pop_data <- pop_raw %>%
  filter(Series.Name == "Total Population (in number of people)") %>%
  select(district = Provinces.Name,
         `2014` = X2014..YR2014.,
         `2015` = X2015..YR2015.,
         `2016` = X2016..YR2016.,
         `2017` = X2017..YR2017.,
         `2018` = X2018..YR2018.,
         `2019` = X2019..YR2019.,
         `2020` = X2020..YR2020.)

# Remove leading quote if present
pop_data$district <- gsub('^"', '', pop_data$district)

# KEEP ALL DISTRICTS - no deduplication
# Population CSV has 514 districts (27 have both Kota and Kab versions)
# We'll use spatial matching to determine which one each protest belongs to

cat(sprintf("Loaded %d districts with population data\n", nrow(pop_data)))
cat(sprintf("  Format: 'Name, Kota' (cities) or 'Name, Kab.' (regencies)\n\n"))

# ====================================================================
# 2. LOAD GADM SHAPEFILE
# ====================================================================

cat("=== STEP 2: LOADING GADM SHAPEFILE ===\n\n")

# Load GADM Level 2 (district/regency level)
shp_path <- "Indonesian district-level shapefile copy 2/gadm41_IDN_2.shp"
gadm <- st_read(shp_path, quiet = TRUE)

cat(sprintf("Loaded GADM shapefile: %d features\n", nrow(gadm)))
cat(sprintf("  Kota (cities): %d\n", sum(gadm$ENGTYPE_2 == "City")))
cat(sprintf("  Kabupaten (regencies): %d\n\n", sum(gadm$ENGTYPE_2 == "Regency")))

# ====================================================================
# 3. SPATIAL MATCHING: PROTESTS → GADM → POPULATION
# ====================================================================

cat("=== STEP 3: SPATIAL MATCHING OF PROTESTS ===\n\n")

# Load protest data
protests <- readRDS("protests_prepared.rds")

cat(sprintf("Loaded %d protests\n", nrow(protests)))
cat(sprintf("  With coordinates: %d (%.1f%%)\n",
            sum(!is.na(protests$latitude) & !is.na(protests$longitude)),
            100*mean(!is.na(protests$latitude) & !is.na(protests$longitude))))

# Convert protests to sf points
protests_with_coords <- protests %>%
  filter(!is.na(latitude) & !is.na(longitude))

protests_sf <- st_as_sf(protests_with_coords,
                        coords = c("longitude", "latitude"),
                        crs = 4326)  # WGS84

cat(sprintf("\nPerforming spatial join for %d protests...\n", nrow(protests_sf)))

# Spatial join: protests → GADM polygons
protests_matched <- st_join(protests_sf, gadm, join = st_within)

# Convert GADM names to population CSV format
# GADM: "Kota Bandung" or "Bandung" → Population CSV: "Bandung, Kota" or "Bandung, Kab."
gadm_to_pop_name <- function(name_2, type_2) {
  # Handle NAs (protests outside all polygons)
  if (is.na(type_2) || is.na(name_2)) {
    return(NA_character_)
  }

  if (type_2 == "Kota") {
    # "Kota Bandung" → "Bandung, Kota"
    base_name <- sub("^Kota ", "", name_2)
    return(paste0(base_name, ", Kota"))
  } else {
    # "Bandung" → "Bandung, Kab."
    return(paste0(name_2, ", Kab."))
  }
}

# Apply conversion
protests_matched$district <- mapply(gadm_to_pop_name,
                                    protests_matched$NAME_2,
                                    protests_matched$TYPE_2)

# Convert back to regular dataframe (drop geometry)
protests_matched_df <- protests_matched %>%
  st_drop_geometry() %>%
  select(event_id_cnty, NAME_2, TYPE_2, district)

cat(sprintf("Spatial matching complete:\n"))
cat(sprintf("  Matched: %d (%.1f%%)\n",
            sum(!is.na(protests_matched_df$district)),
            100*mean(!is.na(protests_matched_df$district))))
cat(sprintf("  Unmatched: %d (%.1f%%)\n\n",
            sum(is.na(protests_matched_df$district)),
            100*mean(is.na(protests_matched_df$district))))

# ====================================================================
# 4. EXTRAPOLATE POPULATION TO 2021-2024
# ====================================================================

cat("=== STEP 4: EXTRAPOLATING POPULATION TO 2024 ===\n\n")

# Reshape to long format
pop_long <- pop_data %>%
  pivot_longer(cols = starts_with("20"),
               names_to = "year",
               values_to = "population") %>%
  mutate(year = as.integer(year),
         population = as.numeric(population)) %>%
  filter(!is.na(population))

cat(sprintf("Reshaped to %d district-year observations (2014-2020)\n", nrow(pop_long)))

# Fit linear models for each district and extrapolate
districts_to_extrapolate <- unique(pop_long$district)

extrapolated_pop <- lapply(districts_to_extrapolate, function(dist) {

  # Get historical data
  hist_data <- pop_long %>% filter(district == dist)

  # Fit linear model
  if (nrow(hist_data) >= 3) {
    model <- lm(population ~ year, data = hist_data)

    # Predict for 2021-2024
    future_years <- data.frame(year = 2021:2024)
    predictions <- predict(model, newdata = future_years)

    # Create extrapolated data
    future_data <- data.frame(
      district = dist,
      year = 2021:2024,
      population = pmax(predictions, min(hist_data$population)),  # Don't go below historical minimum
      extrapolated = TRUE
    )

    # Combine with historical
    hist_data$extrapolated <- FALSE
    return(rbind(hist_data, future_data))

  } else {
    hist_data$extrapolated <- FALSE
    return(hist_data)
  }
})

pop_complete <- do.call(rbind, extrapolated_pop)

cat(sprintf("Created %d district-year observations (2014-2024)\n", nrow(pop_complete)))
cat(sprintf("  Historical: %d\n", sum(!pop_complete$extrapolated)))
cat(sprintf("  Extrapolated: %d\n\n", sum(pop_complete$extrapolated)))

# Add log population
pop_complete$log_pop <- log(pop_complete$population)

# ====================================================================
# 5. JOIN TO PROTEST DATA AND FILTER
# ====================================================================

cat("=== STEP 5: JOINING TO PROTEST DATA ===\n\n")

# Join spatially-matched district names to original protests
protests_with_district <- protests %>%
  left_join(protests_matched_df %>% select(event_id_cnty, district),
            by = "event_id_cnty")

# Join to population data
protests_with_pop <- protests_with_district %>%
  left_join(pop_complete, by = c("district", "year"))

# Count matches
n_total <- nrow(protests)
n_matched <- sum(!is.na(protests_with_pop$population))
n_unmatched <- sum(is.na(protests_with_pop$population))

cat("JOIN RESULTS:\n")
cat(sprintf("  Total protests: %d\n", n_total))
cat(sprintf("  Matched to population: %d (%.1f%%)\n",
            n_matched, 100*n_matched/n_total))
cat(sprintf("  Unmatched: %d (%.1f%%)\n\n",
            n_unmatched, 100*n_unmatched/n_total))

# Show unmatched districts
if (n_unmatched > 0) {
  unmatched_districts <- protests_with_pop %>%
    filter(is.na(population)) %>%
    group_by(admin2) %>%
    summarize(n = n()) %>%
    arrange(desc(n))

  cat("Unmatched districts:\n")
  print(unmatched_districts)
  cat("\n")
}

# Drop unmatched protests
protests_final <- protests_with_pop %>%
  filter(!is.na(population))

cat(sprintf("FINAL DATASET: %d protests (%.1f%% retention)\n\n",
            nrow(protests_final), 100*nrow(protests_final)/n_total))

# ====================================================================
# 6. SAVE OUTPUTS
# ====================================================================

cat("=== STEP 6: SAVING OUTPUTS ===\n\n")

# Save final protest data
saveRDS(protests_final, "protests_with_population.rds")
cat("✓ Saved protests_with_population.rds\n")

# Save matching report
sink("district_matching_report.txt")
cat("DISTRICT MATCHING REPORT (SPATIAL MATCHING)\n")
cat("Generated:", as.character(Sys.time()), "\n\n")
cat(strrep("=", 70), "\n\n")

cat("SPATIAL MATCHING METHOD\n\n")
cat("Used GADM 4.1 shapefile (Level 2: districts/regencies) to spatially\n")
cat("match protest coordinates to official admin boundaries, then matched\n")
cat("GADM names to population CSV format.\n\n")
cat("This correctly distinguishes Kota (city) from Kabupaten (regency)\n")
cat("for the 27 districts with both versions (e.g., Bandung, Kota vs Bandung, Kab.)\n\n")

cat("SUMMARY STATISTICS\n\n")
cat(sprintf("Total protests: %d\n", n_total))
cat(sprintf("  With coordinates: %d\n", nrow(protests_matched_df)))
cat(sprintf("  Spatially matched: %d (%.1f%%)\n",
            sum(!is.na(protests_matched_df$district)),
            100*mean(!is.na(protests_matched_df$district))))
cat(sprintf("  Matched to population: %d (%.1f%%)\n", n_matched, 100*n_matched/n_total))
cat(sprintf("  Final dataset: %d protests (%.1f%% retention)\n",
            nrow(protests_final), 100*nrow(protests_final)/n_total))
cat(sprintf("  Dropped: %d protests (%.1f%%)\n\n", n_unmatched, 100*n_unmatched/n_total))

cat("UNMATCHED PROTESTS\n\n")
if (n_unmatched > 0) {
  write.table(unmatched_districts, row.names = FALSE, quote = FALSE)
}

sink()
cat("✓ Saved district_matching_report.txt\n")

# Save spatial matching table
write.csv(protests_matched_df, "spatial_matching_table.csv", row.names = FALSE)
cat("✓ Saved spatial_matching_table.csv\n\n")

cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  DATA PREPARATION COMPLETE                                    ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n")
