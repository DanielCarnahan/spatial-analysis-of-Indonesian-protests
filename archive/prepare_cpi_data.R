# Prepare Indonesia CPI data from FRED
# Source: FRED series IDNCPIALLMINMEI (OECD Main Economic Indicators)
# Base year: 2015 = 100

library(tidyverse)

# Load raw data
cpi_raw <- read_csv("indonesia_cpi_fred.csv", show_col_types = FALSE)

# Process CPI data
cpi <- cpi_raw %>%
  rename(
    date = observation_date,
    cpi = IDNCPIALLMINMEI
  ) %>%
  mutate(
    date = as.Date(date),
    year = year(date),
    month = month(date),
    year_month = format(date, "%Y-%m")
  ) %>%
  # Filter to study period (2015-2024)
  filter(year >= 2015 & year <= 2024) %>%
  arrange(date) %>%
  # Compute inflation rates
  mutate(
    # Month-over-month inflation (%)
    inflation_mom = 100 * (cpi / lag(cpi) - 1),
    # Year-over-year inflation (%)
    inflation_yoy = 100 * (cpi / lag(cpi, 12) - 1)
  )

# Summary statistics
cat("=== Indonesia CPI Data Summary ===\n\n")
cat("Date range:", as.character(min(cpi$date)), "to", as.character(max(cpi$date)), "\n")
cat("Number of months:", nrow(cpi), "\n\n")

cat("CPI values:\n")
cat("  2015 Jan:", round(cpi$cpi[cpi$year_month == "2015-01"], 2), "\n")
cat("  2020 Jan:", round(cpi$cpi[cpi$year_month == "2020-01"], 2), "\n")
cat("  2024 Dec:", round(cpi$cpi[cpi$year_month == "2024-12"], 2), "\n\n")

cat("Year-over-year inflation range:\n")
cat("  Min:", round(min(cpi$inflation_yoy, na.rm = TRUE), 2), "%\n")
cat("  Max:", round(max(cpi$inflation_yoy, na.rm = TRUE), 2), "%\n")
cat("  Mean:", round(mean(cpi$inflation_yoy, na.rm = TRUE), 2), "%\n\n")

# Check for missing months
expected_months <- seq(as.Date("2015-01-01"), as.Date("2024-12-01"), by = "month")
missing <- expected_months[!expected_months %in% cpi$date]
if (length(missing) == 0) {
  cat("No missing months in study period.\n\n")
} else {
  cat("Missing months:", as.character(missing), "\n\n")
}

# Save processed data
saveRDS(cpi, "indonesia_cpi.rds")
cat("Saved to: indonesia_cpi.rds\n")

# Also save as CSV for reference
write_csv(cpi, "indonesia_cpi_processed.csv")
cat("Also saved to: indonesia_cpi_processed.csv\n")
