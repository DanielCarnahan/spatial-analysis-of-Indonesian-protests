# ============================================================================
# PREPARE POVERTY RATE DATA FOR HAWKES MODEL
# ============================================================================
#
# Purpose: Extract poverty rate from population CSV, extrapolate 2021-2024,
#          and join to protest dataset
#
# Author: Analysis Script
# Date: 2025-11-13
# ============================================================================

library(dplyr)
library(tidyr)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║        POVERTY RATE DATA PREPARATION                         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# ====================================================================
# 1. LOAD POPULATION DATA
# ====================================================================

cat("=== STEP 1: LOADING POPULATION DATA ===\n\n")

pop_raw <- read.csv("district_level_population_2014_2020.csv",
                    stringsAsFactors = FALSE)

cat(sprintf("Loaded %d rows from population CSV\n", nrow(pop_raw)))

# ====================================================================
# 2. EXTRACT POVERTY RATE
# ====================================================================

cat("\n=== STEP 2: EXTRACTING POVERTY RATE ===\n\n")

poverty_data <- pop_raw %>%
  filter(Series.Name == "Poverty Rate (in % of population)") %>%
  select(district = Provinces.Name,
         `2014` = X2014..YR2014.,
         `2015` = X2015..YR2015.,
         `2016` = X2016..YR2016.,
         `2017` = X2017..YR2017.,
         `2018` = X2018..YR2018.,
         `2019` = X2019..YR2019.,
         `2020` = X2020..YR2020.)

# Remove leading quote if present
poverty_data$district <- gsub('^"', '', poverty_data$district)

cat(sprintf("Extracted poverty data for %d districts\n", nrow(poverty_data)))

# Convert to long format
poverty_long <- poverty_data %>%
  pivot_longer(cols = `2014`:`2020`,
               names_to = "year",
               values_to = "poverty_rate") %>%
  mutate(year = as.integer(year),
         poverty_rate = as.numeric(poverty_rate))

# Check coverage
coverage_by_year <- poverty_long %>%
  group_by(year) %>%
  summarise(
    n_districts = n(),
    n_with_data = sum(!is.na(poverty_rate) & poverty_rate != ".."),
    pct_coverage = 100 * n_with_data / n_districts
  )

cat("\nPoverty rate coverage by year:\n")
print(as.data.frame(coverage_by_year))

# Filter to valid data for 2015-2020 (our protest period)
poverty_clean <- poverty_long %>%
  filter(year >= 2015 & year <= 2020) %>%
  filter(!is.na(poverty_rate) & poverty_rate != "..")

cat(sprintf("\nValid poverty data: %d district-years (2015-2020)\n",
            nrow(poverty_clean)))
cat(sprintf("Unique districts: %d\n", n_distinct(poverty_clean$district)))

# Summary statistics
cat("\nPoverty rate summary (2015-2020):\n")
cat(sprintf("  Mean: %.2f%%\n", mean(poverty_clean$poverty_rate, na.rm = TRUE)))
cat(sprintf("  Median: %.2f%%\n", median(poverty_clean$poverty_rate, na.rm = TRUE)))
cat(sprintf("  Range: %.2f%% to %.2f%%\n",
            min(poverty_clean$poverty_rate, na.rm = TRUE),
            max(poverty_clean$poverty_rate, na.rm = TRUE)))

# ====================================================================
# 3. EXTRAPOLATE FOR 2021-2024
# ====================================================================

cat("\n=== STEP 3: EXTRAPOLATING POVERTY RATE (2021-2024) ===\n\n")

# Calculate linear trend for each district using 2018-2020 data
poverty_trends <- poverty_clean %>%
  filter(year >= 2018 & year <= 2020) %>%
  group_by(district) %>%
  summarise(
    n_obs = n(),
    # Linear trend: slope from 2018-2020
    poverty_2018 = poverty_rate[year == 2018][1],
    poverty_2019 = poverty_rate[year == 2019][1],
    poverty_2020 = poverty_rate[year == 2020][1],
    # Calculate average annual change
    annual_change = case_when(
      n_obs >= 3 ~ (poverty_2020 - poverty_2018) / 2,  # 2-year span
      n_obs == 2 ~ poverty_2020 - poverty_2019,        # 1-year change
      TRUE ~ 0                                          # Constant if only 1 obs
    ),
    .groups = "drop"
  )

cat(sprintf("Calculated trends for %d districts\n", nrow(poverty_trends)))
cat(sprintf("  Mean annual change: %.3f percentage points\n",
            mean(poverty_trends$annual_change, na.rm = TRUE)))

# Project forward
poverty_projected <- poverty_trends %>%
  mutate(
    poverty_2021 = pmax(0, pmin(100, poverty_2020 + annual_change)),
    poverty_2022 = pmax(0, pmin(100, poverty_2020 + 2 * annual_change)),
    poverty_2023 = pmax(0, pmin(100, poverty_2020 + 3 * annual_change)),
    poverty_2024 = pmax(0, pmin(100, poverty_2020 + 4 * annual_change))
  ) %>%
  select(district, poverty_2021, poverty_2022, poverty_2023, poverty_2024)

# Convert to long format
poverty_future <- poverty_projected %>%
  pivot_longer(cols = starts_with("poverty_"),
               names_to = "year",
               names_prefix = "poverty_",
               values_to = "poverty_rate") %>%
  mutate(year = as.integer(year),
         extrapolated = TRUE)

# Mark historical data
poverty_clean <- poverty_clean %>%
  mutate(extrapolated = FALSE)

# Combine historical and projected
poverty_full <- bind_rows(poverty_clean, poverty_future) %>%
  arrange(district, year)

cat(sprintf("\nTotal poverty data: %d district-years (2015-2024)\n",
            nrow(poverty_full)))
cat(sprintf("  Historical (2015-2020): %d\n", sum(!poverty_full$extrapolated)))
cat(sprintf("  Extrapolated (2021-2024): %d\n", sum(poverty_full$extrapolated)))

# Check extrapolated values
cat("\nExtrapolated poverty rate summary (2021-2024):\n")
extrap_summary <- poverty_full %>%
  filter(extrapolated) %>%
  group_by(year) %>%
  summarise(
    mean = mean(poverty_rate, na.rm = TRUE),
    median = median(poverty_rate, na.rm = TRUE),
    min = min(poverty_rate, na.rm = TRUE),
    max = max(poverty_rate, na.rm = TRUE)
  )
print(as.data.frame(extrap_summary))

# ====================================================================
# 4. LOAD PROTEST DATA AND CHECK MATCHING
# ====================================================================

cat("\n=== STEP 4: LOADING PROTEST DATA ===\n\n")

protests <- readRDS("protests_with_population.rds")

cat(sprintf("Loaded %d protests\n", nrow(protests)))
cat(sprintf("  Unique districts: %d\n", n_distinct(protests$district)))
cat(sprintf("  Years: %d to %d\n", min(protests$year), max(protests$year)))

# Check district name format
cat("\nSample district names in protest data:\n")
cat(paste("  -", head(unique(protests$district), 5)), sep = "\n")

cat("\nSample district names in poverty data:\n")
cat(paste("  -", head(unique(poverty_full$district), 5)), sep = "\n")

# ====================================================================
# 5. JOIN POVERTY TO PROTESTS
# ====================================================================

cat("\n=== STEP 5: JOINING POVERTY DATA TO PROTESTS ===\n\n")

protests_with_poverty <- protests %>%
  left_join(
    poverty_full %>% select(district, year, poverty_rate, extrapolated),
    by = c("district", "year")
  )

# Check match rate
n_matched <- sum(!is.na(protests_with_poverty$poverty_rate))
pct_matched <- 100 * n_matched / nrow(protests_with_poverty)

cat(sprintf("Match results:\n"))
cat(sprintf("  Total protests: %d\n", nrow(protests_with_poverty)))
cat(sprintf("  Matched to poverty data: %d (%.1f%%)\n", n_matched, pct_matched))
cat(sprintf("  No poverty data: %d (%.1f%%)\n",
            nrow(protests_with_poverty) - n_matched,
            100 - pct_matched))

# Check by year
match_by_year <- protests_with_poverty %>%
  group_by(year) %>%
  summarise(
    n_protests = n(),
    n_with_poverty = sum(!is.na(poverty_rate)),
    pct_matched = 100 * n_with_poverty / n_protests,
    mean_poverty = mean(poverty_rate, na.rm = TRUE)
  )

cat("\nMatch rate by year:\n")
print(as.data.frame(match_by_year))

# Identify unmatched districts
unmatched_districts <- protests_with_poverty %>%
  filter(is.na(poverty_rate)) %>%
  distinct(district) %>%
  pull(district)

if(length(unmatched_districts) > 0) {
  cat(sprintf("\n⚠ Warning: %d districts have no poverty data:\n",
              length(unmatched_districts)))
  cat(paste("  -", head(unmatched_districts, 10)), sep = "\n")
  if(length(unmatched_districts) > 10) {
    cat(sprintf("  ... and %d more\n", length(unmatched_districts) - 10))
  }
}

# ====================================================================
# 6. SUMMARY STATISTICS
# ====================================================================

cat("\n=== STEP 6: SUMMARY STATISTICS ===\n\n")

# Overall summary
cat("Poverty rate in protest sample:\n")
cat(sprintf("  Mean: %.2f%%\n", mean(protests_with_poverty$poverty_rate, na.rm = TRUE)))
cat(sprintf("  Median: %.2f%%\n", median(protests_with_poverty$poverty_rate, na.rm = TRUE)))
cat(sprintf("  SD: %.2f\n", sd(protests_with_poverty$poverty_rate, na.rm = TRUE)))
cat(sprintf("  Range: %.2f%% to %.2f%%\n",
            min(protests_with_poverty$poverty_rate, na.rm = TRUE),
            max(protests_with_poverty$poverty_rate, na.rm = TRUE)))

# Correlation with population
if("log_pop" %in% names(protests_with_poverty)) {
  cor_pop <- cor(protests_with_poverty$log_pop,
                 protests_with_poverty$poverty_rate,
                 use = "complete.obs")
  cat(sprintf("\nCorrelation (log(pop) vs poverty): %.3f\n", cor_pop))
  if(abs(cor_pop) > 0.5) {
    cat("  ⚠ Moderate to high correlation - potential multicollinearity\n")
  } else {
    cat("  ✓ Low to moderate correlation\n")
  }
}

# ====================================================================
# 7. SAVE RESULTS
# ====================================================================

cat("\n=== STEP 7: SAVING RESULTS ===\n\n")

# Save poverty panel data
saveRDS(poverty_full, "poverty_rate_2015_2024.rds")
cat("✓ Saved: poverty_rate_2015_2024.rds\n")
cat(sprintf("  Contains: %d district-year observations\n", nrow(poverty_full)))

# Save protests with poverty
saveRDS(protests_with_poverty, "protests_with_poverty.rds")
cat("✓ Saved: protests_with_poverty.rds\n")
cat(sprintf("  Contains: %d protests with poverty rate (%.1f%% coverage)\n",
            n_matched, pct_matched))

# Save CSV for inspection
write.csv(poverty_full, "poverty_rate_2015_2024.csv", row.names = FALSE)
cat("✓ Saved: poverty_rate_2015_2024.csv (for inspection)\n")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║         POVERTY DATA PREPARATION COMPLETE                    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("NEXT STEPS:\n")
cat("-----------\n")
cat("1. Review poverty_rate_2015_2024.csv to validate extrapolations\n")
cat("2. Investigate unmatched districts if match rate < 95%\n")
cat("3. Proceed to fit Hawkes model with poverty rate (16_hawkes_with_poverty.R)\n")
cat("\n")
