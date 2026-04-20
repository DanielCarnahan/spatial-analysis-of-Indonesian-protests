# Do Protests Really Spread? Self-Excitation in Indonesian Protests

## Reproducible Analysis Pipeline

**Author:** Daniel Carnahan
**Date:** January 2026

---

## Table of Contents

1. [Overview](#overview)
2. [Data Sources](#data-sources)
3. [Required R Packages](#required-r-packages)
4. [Analysis Pipeline](#analysis-pipeline)
   - [Step 1: Prepare Protest Data](#step-1-prepare-protest-data)
   - [Step 2: Prepare Population Data](#step-2-prepare-population-data)
   - [Step 3: Prepare Poverty Data](#step-3-prepare-poverty-data)
   - [Step 4: Prepare CPI Data](#step-4-prepare-cpi-data)
   - [Step 5: Collapse Multi-Day Protests](#step-5-collapse-multi-day-protests)
   - [Step 6: Estimate Spatial Hawkes Model](#step-6-estimate-spatial-hawkes-model)
   - [Step 7: Compute Branching Ratio](#step-7-compute-branching-ratio)
   - [Step 8: Generate Figures](#step-8-generate-figures)
5. [Key Results](#key-results)
6. [Model Specification](#model-specification)
7. [Reproducing the Analysis](#reproducing-the-analysis)

---

## Overview

This document provides a complete, reproducible pipeline for analyzing protest diffusion in Indonesia using a discrete spatial Hawkes process. The analysis addresses two main questions:

1. **Do protests elsewhere predict more protests here?** (Cross-district contagion)
2. **Does geography matter—is spread local or national?** (Spatial decay)

### Key Findings

- Cross-district contagion exists (F(1, 1.56M) = 185, p < 10⁻⁴¹)
- Spatial decay distance θ = 100 km (contagion is local, not national)
- Branching ratio = 0.11 (subcritical, stable process)
- ~14% of protests are triggered by cross-district contagion

---

## Data Sources

| Data | Source | File | Description |
|------|--------|------|-------------|
| Protest events | ACLED | `ACLED Data_2025-10-27_Indonesia20142024_OFFICIAL.xlsx` | Armed Conflict Location and Event Data |
| Population | World Bank DAPOER | `district_level_population_2014_2020.csv` | District-level population 2014-2020 |
| Poverty rates | World Bank DAPOER | `district_level_population_2014_2020.csv` | District-level poverty rates |
| CPI | FRED | `indonesia_cpi_fred.csv` | Consumer Price Index (IDNCPIALLMINMEI) |
| Shapefiles | GADM 4.1 | `Indonesian district-level shapefile copy 2/` | District boundaries for spatial matching |

---

## Required R Packages

```r
# Core packages
install.packages(c(
  "tidyverse",      # Data manipulation and visualization
  "sf",             # Spatial data handling
  "readxl",         # Reading Excel files
  "lubridate",      # Date handling
  "MASS",           # Negative binomial GLM
  "geosphere",      # Distance calculations
  "scales",         # Plot formatting
  "rnaturalearth",  # Country basemaps
  "stringdist"      # String matching for name harmonization
))
```

---

## Analysis Pipeline

### Step 1: Prepare Protest Data

**Script:** `02_prepare_protest_data.R`
**Output:** `protests_prepared.rds`, `protests_sf.rds`

This step loads raw ACLED data and performs initial processing:
- Filters to protest and riot events
- Classifies events by violence level
- Extracts temporal components
- Creates spatial features

```r
# Load ACLED data
library(readxl)
library(dplyr)
library(lubridate)
library(sf)

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
    fatalities = as.integer(fatalities)
  )

# Filter to protest events
protests <- acled %>%
  filter(event_type %in% c("Protests", "Riots")) %>%
  mutate(
    is_violent = sub_event_type %in% c("Violent demonstration",
                                        "Protest with intervention",
                                        "Mob violence",
                                        "Excessive force against protesters"),
    is_peaceful = sub_event_type == "Peaceful protest",
    has_fatalities = fatalities > 0,
    days_since_start = as.numeric(event_date - min(event_date))
  ) %>%
  arrange(event_date)

# Save prepared data
saveRDS(protests, "protests_prepared.rds")
```

**Output Summary:**
- Total protest events: ~13,000
- Date range: 2014 to 2024
- Key variables: `event_date`, `latitude`, `longitude`, `is_violent`, `fatalities`

---

### Step 2: Prepare Population Data

**Script:** `10_prepare_population_data.R`
**Output:** `protests_with_population.rds`

This step spatially matches protests to district boundaries and joins population data:

```r
library(sf)
library(dplyr)
library(tidyr)

# Load population data
pop_raw <- read.csv("district_level_population_2014_2020.csv")

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

# Load GADM shapefile for spatial matching
gadm <- st_read("Indonesian district-level shapefile copy 2/gadm41_IDN_2.shp")

# Spatial join: protests → GADM polygons
protests <- readRDS("protests_prepared.rds")
protests_sf <- st_as_sf(protests, coords = c("longitude", "latitude"), crs = 4326)
protests_matched <- st_join(protests_sf, gadm, join = st_within)

# Extrapolate population to 2021-2024 using linear trends
pop_long <- pop_data %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "population") %>%
  mutate(year = as.integer(year), population = as.numeric(population))

# Fit linear models and extrapolate
# ... (extrapolation code)

# Join population to protests
protests_with_pop <- protests %>%
  left_join(pop_complete, by = c("district", "year")) %>%
  mutate(log_pop = log(population))

saveRDS(protests_with_pop, "protests_with_population.rds")
```

**Spatial Matching Details:**
- Uses GADM 4.1 district boundaries
- Distinguishes Kota (city) from Kabupaten (regency)
- Match rate: >95%

---

### Step 3: Prepare Poverty Data

**Script:** `15_prepare_poverty_data.R`
**Output:** `protests_with_poverty.rds`, `poverty_rate_2015_2024.rds`

```r
library(dplyr)
library(tidyr)

# Extract poverty rate from population CSV
poverty_data <- pop_raw %>%
  filter(Series.Name == "Poverty Rate (in % of population)") %>%
  select(district = Provinces.Name,
         `2014` = X2014..YR2014.,
         # ... years 2015-2020
         `2020` = X2020..YR2020.)

# Reshape to long format
poverty_long <- poverty_data %>%
  pivot_longer(cols = `2014`:`2020`, names_to = "year", values_to = "poverty_rate") %>%
  mutate(year = as.integer(year), poverty_rate = as.numeric(poverty_rate))

# Extrapolate 2021-2024 using linear trends from 2018-2020
poverty_trends <- poverty_long %>%
  filter(year >= 2018 & year <= 2020) %>%
  group_by(district) %>%
  summarise(
    poverty_2020 = poverty_rate[year == 2020][1],
    annual_change = (poverty_rate[year == 2020][1] - poverty_rate[year == 2018][1]) / 2
  )

# Project forward with bounds [0, 100]
poverty_projected <- poverty_trends %>%
  mutate(
    poverty_2021 = pmax(0, pmin(100, poverty_2020 + annual_change)),
    poverty_2022 = pmax(0, pmin(100, poverty_2020 + 2 * annual_change)),
    poverty_2023 = pmax(0, pmin(100, poverty_2020 + 3 * annual_change)),
    poverty_2024 = pmax(0, pmin(100, poverty_2020 + 4 * annual_change))
  )

# Join to protests
protests_with_poverty <- protests_with_pop %>%
  left_join(poverty_full, by = c("district", "year")) %>%
  mutate(poverty_decimal = poverty_rate / 100)

saveRDS(protests_with_poverty, "protests_with_poverty.rds")
```

---

### Step 4: Prepare CPI Data

**Script:** `prepare_cpi_data.R`
**Output:** `indonesia_cpi.rds`, `indonesia_cpi_processed.csv`

```r
library(tidyverse)

# Load raw CPI data from FRED
cpi_raw <- read_csv("indonesia_cpi_fred.csv")

# Process CPI data
cpi <- cpi_raw %>%
  rename(date = observation_date, cpi = IDNCPIALLMINMEI) %>%
  mutate(
    date = as.Date(date),
    year = year(date),
    month = month(date),
    year_month = format(date, "%Y-%m"),
    # Compute inflation rates
    inflation_mom = 100 * (cpi / lag(cpi) - 1),
    inflation_yoy = 100 * (cpi / lag(cpi, 12) - 1)
  ) %>%
  filter(year >= 2015 & year <= 2024)

saveRDS(cpi, "indonesia_cpi.rds")
write_csv(cpi, "indonesia_cpi_processed.csv")
```

**CPI Summary:**
- Base year: 2015 = 100
- Source: OECD Main Economic Indicators via FRED
- Coverage: Jan 2015 - Oct 2024 (monthly)

---

### Step 5: Collapse Multi-Day Protests

**Script:** `collapse_multiday_protests.R`
**Output:** `protests_collapsed_3day.rds`, `protests_daily.rds`

This step addresses the issue of multi-day protests being recorded as separate events:

```r
library(tidyverse)

protests <- readRDS("protests_with_poverty.rds")

# Define topic patterns for Jakarta (coarse geocoding requires topic matching)
topic_patterns <- list(
  labor = "labor|worker|wage|buruh",
  student = "student|mahasiswa|university",
  corruption = "corruption|korupsi|kpk",
  # ... additional topics
)

# Create collapse keys
protests <- protests %>%
  mutate(
    location_id = paste(round(latitude, 4), round(longitude, 4), sep = "_"),
    is_jakarta = admin1 == "Jakarta",
    # Jakarta: location + topic; elsewhere: location only
    collapse_key = ifelse(is_jakarta,
                          paste(location_id, topic_signature, sep = "__"),
                          location_id)
  )

# Assign cluster IDs using 3-day gap rule
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

protests <- protests %>%
  group_by(collapse_key) %>%
  mutate(cluster_id = assign_clusters(event_date, gap_days = 3)) %>%
  ungroup()

# Collapse events within each cluster
protests_collapsed <- protests %>%
  group_by(event_cluster) %>%
  summarize(
    event_date = min(event_date),
    latitude = first(latitude),
    longitude = first(longitude),
    n_events_collapsed = n(),
    fatalities = sum(fatalities, na.rm = TRUE),
    is_violent = any(is_violent, na.rm = TRUE),
    # ... additional aggregations
  )

# Further aggregate to daily resolution
protests_daily <- protests_collapsed %>%
  group_by(date, latitude, longitude) %>%
  summarise(
    is_severe = any(is_severe),
    n_events_collapsed = sum(n_events_collapsed),
    # ... merge CPI data
  )

saveRDS(protests_collapsed, "protests_collapsed_3day.rds")
saveRDS(protests_daily, "protests_daily.rds")
```

**Event Definition:**
- Events within 7 days at the same location → single event
- Avoids double-counting continuation of the same protest
- Focuses on distinct mobilization episodes

---

### Step 6: Estimate Spatial Hawkes Model

**Script:** `estimate_spatial_hawkes_final.R`
**Output:** `spatial_hawkes_final_results.rds`

This is the main estimation script implementing the discrete spatial Hawkes model.

#### Model Specification

**Baseline Model (M0):** Inhomogeneous Poisson with covariates only
```
E[Y_{t,r}] = exp(β₀ + β₁ log(pop_r) + β₂ pov_r + β₃ CPI_t + γ_year)
```

**Spatial-Temporal Hawkes Model (M1):** Adds cross-region contagion
```
E[Y_{t,r}] = exp(X'_{t,r}β + α × Cross_{t,r}(θ))
```

where the cross-region exposure term is:
```
Cross_{t,r}(θ) = Σ_{s≠r} w_{rs}(θ) × Σ_{ℓ=1}^{30} Y_{t-ℓ,s}
```

**Spatial Kernel:**
```
w_{rs}(θ) ∝ exp(-d_{rs}/θ)
```

```r
library(tidyverse)
library(geosphere)

# Load data
protests <- readRDS("protests_daily.rds")
cpi_data <- read_csv("indonesia_cpi_processed.csv")

# Create district-day panel
all_districts <- unique(protests$admin2)
all_dates <- seq(min(protests$date), max(protests$date), by = "day")

panel <- expand.grid(admin2 = all_districts, date = all_dates) %>%
  as_tibble() %>%
  left_join(daily_counts, by = c("admin2", "date")) %>%
  mutate(count = replace_na(count, 0))

# Compute distance matrix
n_districts <- length(all_districts)
dist_matrix <- matrix(0, n_districts, n_districts)
for (i in 1:n_districts) {
  for (j in 1:n_districts) {
    if (i != j) {
      dist_matrix[i, j] <- distHaversine(
        c(district_covariates$lon[i], district_covariates$lat[i]),
        c(district_covariates$lon[j], district_covariates$lat[j])
      ) / 1000  # km
    }
  }
}

# Function to compute cross-district exposure
compute_cross_exposure <- function(theta, dist_matrix, count_matrix, n_days, n_dist) {
  if (is.infinite(theta)) {
    W <- (dist_matrix > 0) * 1  # Uniform weights
  } else {
    W <- exp(-dist_matrix / theta)  # Exponential decay
    diag(W) <- 0
  }

  # Row-normalize
  W_norm <- W / pmax(rowSums(W), 1e-10)

  # Sum over lags 1-30
  cross_exposure <- matrix(0, n_days, n_dist)
  for (lag in 1:30) {
    lagged <- rbind(matrix(0, lag, n_dist), count_matrix[1:(n_days - lag), ])
    cross_exposure <- cross_exposure + t(W_norm %*% t(lagged))
  }

  return(cross_exposure)
}

# M0: Baseline model (no contagion)
m0 <- glm(count ~ log_pop + poverty_rate + cpi + factor(year),
          data = panel_clean, family = quasipoisson(link = "log"))

# Profile likelihood over θ
theta_grid <- c(50, 100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000, Inf)

results_grid <- data.frame(theta = theta_grid, deviance = NA, alpha = NA, se = NA)

for (i in seq_along(theta_grid)) {
  theta <- theta_grid[i]

  cross_exp <- compute_cross_exposure(theta, dist_matrix, count_matrix, n_days, n_dist)

  # Reshape to long format and join
  exp_long <- data.frame(
    date = rep(dates_vec, n_dist),
    admin2 = rep(colnames(count_matrix), each = n_days),
    cross_lag = as.vector(cross_exp)
  )

  panel_fit <- panel_clean %>% left_join(exp_long, by = c("admin2", "date"))

  # Fit model
  m <- glm(count ~ log_pop + poverty_rate + cpi + factor(year) + cross_lag,
           data = panel_fit, family = quasipoisson(link = "log"))

  results_grid$deviance[i] <- deviance(m)
  results_grid$alpha[i] <- coef(m)["cross_lag"]
  results_grid$se[i] <- summary(m)$coefficients["cross_lag", "Std. Error"]
}

# Find optimal θ
best_idx <- which.min(results_grid$deviance)
theta_hat <- results_grid$theta[best_idx]

# Hypothesis tests
# Test 1: Does contagion exist? (F-test: M0 vs M1)
F_stat <- ((deviance(m0) - deviance(m1_best)) / 1) / dispersion
p_contagion <- pf(F_stat, 1, nrow(panel_clean) - length(coef(m1_best)), lower.tail = FALSE)

# Test 2: Does distance matter? (LR test: θ = ∞ vs θ = θ_hat)
LR_stat <- (dev_inf - dev_best) / dispersion
p_distance <- pchisq(LR_stat, df = 1, lower.tail = FALSE)

# Save results
results <- list(
  profile_likelihood = results_grid,
  optimal_theta = theta_hat,
  test_contagion = list(F_stat = F_stat, p_value = p_contagion),
  test_distance = list(LR_stat = LR_stat, p_value = p_distance),
  final_model = list(coefficients = coef(m1_best), se = summary(m1_best)$coef[, 2])
)

saveRDS(results, "spatial_hawkes_final_results.rds")
```

---

### Step 7: Compute Branching Ratio

**Script:** `compute_branching_ratio.R`
**Output:** `branching_ratio_results.rds`

```r
library(tidyverse)
library(geosphere)

# Load results
results <- readRDS("spatial_hawkes_final_results.rds")
theta_hat <- results$optimal_theta
alpha <- results$final_model$coefficients["cross_lag"]

# Compute endogenous fraction via counterfactual decomposition
# For log-link: λ = exp(X'β + α×Cross) = exp(X'β) × exp(α×Cross)
# Background intensity: λ_bg = exp(X'β)
# Triggered intensity: λ_full - λ_bg

X_bg <- model.matrix(~ log_pop + poverty_rate + cpi + factor(year), data = panel_fit)
coefs_bg <- coef(m1)[colnames(X_bg)]
lambda_bg <- exp(X_bg %*% coefs_bg)
lambda_full <- predict(m1, type = "response")
lambda_triggered <- lambda_full - lambda_bg

endogenous_fraction <- sum(lambda_triggered) / sum(lambda_full)

# Compute branching ratio
# For each origin district r: offspring = Σ_s λ_s × (exp(α × w_rs) - 1) × n_lags
avg_lambda_by_district <- panel_fit %>%
  group_by(admin2) %>%
  summarise(avg_lambda = mean(predict(m1, newdata = cur_data(), type = "response")))

offspring_by_origin <- numeric(n_dist)
for (r in 1:n_dist) {
  total <- 0
  for (s in 1:n_dist) {
    if (s == r) next
    w <- W_norm[s, r]
    offspring_s <- 30 * avg_lambda_vec[s] * (exp(alpha * w) - 1)
    total <- total + offspring_s
  }
  offspring_by_origin[r] <- total
}

# Weight by event frequency
event_weights <- events_by_district$n_events / sum(events_by_district$n_events)
branching_ratio <- sum(offspring_by_origin * event_weights)

# Save results
br_results <- list(
  alpha = alpha,
  exp_alpha = exp(alpha),
  theta = theta_hat,
  endogenous_fraction = endogenous_fraction,
  branching_ratio = branching_ratio
)

saveRDS(br_results, "branching_ratio_results.rds")
```

**Interpretation:**
- Branching ratio = 0.11 → each protest triggers ~0.11 additional protests elsewhere
- Process is subcritical (BR < 1): stable, non-explosive
- ~14% of protests attributable to cross-district contagion

---

### Step 8: Generate Figures

**Script:** `generate_presentation_figures.R`
**Output:** `figures/` directory

```r
library(tidyverse)
library(scales)

theme_presentation <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16)
  )

# Load data
protests <- readRDS("protests_daily.rds")
results <- readRDS("spatial_hawkes_final_results.rds")

# Figure 1: Daily event counts
daily_counts <- protests %>%
  group_by(date) %>%
  summarise(n = n())

p1 <- ggplot(daily_counts, aes(x = date, y = n)) +
  geom_col(fill = "#2c3e50", alpha = 0.7) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = "Daily Protest Count") +
  theme_presentation

ggsave("figures/fig_timeseries.pdf", p1, width = 10, height = 4)

# Figure 2: Profile likelihood
profile_data <- results$profile_likelihood %>%
  mutate(theta_plot = ifelse(is.infinite(theta), 10000, theta))

p2 <- ggplot(profile_data, aes(x = theta_plot, y = deviance)) +
  geom_line(color = "#3498db", linewidth = 1.2) +
  geom_point(color = "#3498db", size = 3) +
  scale_x_log10(
    breaks = c(50, 100, 200, 500, 1000, 2000, 5000, 10000),
    labels = c("50", "100", "200", "500", "1000", "2000", "5000", expression(infinity))
  ) +
  labs(x = "Decay Distance θ (km)", y = "Deviance") +
  theme_presentation

ggsave("figures/fig_profile_likelihood.pdf", p2, width = 8, height = 5)

# Figure 3: Lag structure (weekly bins)
lag_data <- data.frame(
  week = factor(c("Week 1\n(1-7)", "Week 2\n(8-14)", "Week 3\n(15-21)", "Week 4\n(22-30)")),
  coef = results$lag_structure$coefficients,
  se = results$lag_structure$se
)

p4 <- ggplot(lag_data, aes(x = week, y = coef)) +
  geom_col(fill = "#3498db", alpha = 0.8) +
  geom_errorbar(aes(ymin = coef - 1.96*se, ymax = coef + 1.96*se), width = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Lag Period", y = "Coefficient") +
  theme_presentation

ggsave("figures/fig_lag_structure.pdf", p4, width = 8, height = 5)

# Figure 4: Covariate effects
coefs <- results$final_model$coefficients
ses <- results$final_model$se

cov_df <- data.frame(
  variable = c("Log Population", "Poverty Rate", "CPI"),
  estimate = c(coefs["log_pop"], coefs["poverty_rate"], coefs["cpi"]),
  se = c(ses["log_pop"], ses["poverty_rate"], ses["cpi"])
)

p5 <- ggplot(cov_df, aes(x = variable, y = estimate)) +
  geom_col(fill = "#9b59b6", alpha = 0.8) +
  geom_errorbar(aes(ymin = estimate - 1.96*se, ymax = estimate + 1.96*se), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "Coefficient (log scale)") +
  theme_presentation

ggsave("figures/fig_covariates.pdf", p5, width = 7, height = 5)

# Figure 5: Spatial distribution map
library(rnaturalearth)

indonesia <- ne_countries(scale = "medium", country = "Indonesia", returnclass = "sf")

p_map <- ggplot() +
  geom_sf(data = indonesia, fill = "gray95", color = "gray60", linewidth = 0.3) +
  geom_point(data = protests, aes(x = longitude, y = latitude),
             alpha = 0.25, size = 0.6, color = "#c0392b") +
  coord_sf(xlim = c(95, 141), ylim = c(-11, 6)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = NULL, y = NULL)

ggsave("figures/fig_spatial_distribution.pdf", p_map, width = 10, height = 4)
```

**Generated Figures:**
1. `fig_timeseries.pdf` - Daily protest counts over time
2. `fig_profile_likelihood.pdf` - Profile likelihood over θ
3. `fig_lag_structure.pdf` - Temporal lag effects by week
4. `fig_covariates.pdf` - Background rate coefficient estimates
5. `fig_spatial_distribution.pdf` - Map of protest locations

---

## Key Results

### Test 1: Does Cross-District Contagion Exist?

| Model | Deviance | Parameters |
|-------|----------|------------|
| M0: Background Only | 108,713 | 13 |
| M1: + Contagion (θ = 100 km) | 108,518 | 14 |

**F-test:** F(1, 1,557,187) = 184.6, p < 10⁻⁴¹

**Result:** Strong evidence for cross-district contagion

### Test 2: Does Distance Matter?

| θ (km) | Deviance | α | SE(α) |
|--------|----------|---|-------|
| 50 | 108,523 | 0.229 | 0.016 |
| **100** | **108,518** | **0.323** | **0.023** |
| 200 | 108,596 | 0.385 | 0.036 |
| 500 | 108,680 | 0.336 | 0.060 |
| 1000 | 108,681 | 0.420 | 0.076 |
| ∞ | 108,591 | 1.022 | 0.095 |

**Likelihood Ratio Test:** LR = 69.1, χ²(1), p < 0.001

**Result:** Distance matters — optimal θ = 100 km

### Coefficient Estimates (Final Model)

| Parameter | Estimate | SE | z | p |
|-----------|----------|-----|---|---|
| Intercept | −19.03 | 0.902 | −21.1 | < 0.001 |
| log(population) | 0.636 | 0.012 | 54.5 | < 0.001 |
| Poverty rate | −0.021 | 0.002 | −11.7 | < 0.001 |
| CPI | 0.049 | 0.009 | 5.7 | < 0.001 |
| **Cross-district (α)** | **0.323** | **0.023** | **14.0** | **< 0.001** |

Dispersion φ = 1.056; N = 1,557,200 district-days

### Magnitude

- **Branching Ratio:** 0.11 (subcritical)
- **Endogenous Fraction:** 14.4%
- **Interpretation:** Each protest triggers ~0.11 additional protests elsewhere; ~86% of protests are exogenous (background)

---

## Model Specification

### Discrete Spatial Hawkes Process

The model decomposes the expected count of protests in district r on day t into:

1. **Background rate** μ_{t,r}: Exogenous events driven by structural factors
2. **Triggered component**: Endogenous events caused by prior events in other districts

**Full Specification:**

```
E[Y_{t,r}] = exp(X'_{t,r}β + α × Cross_{t,r}(θ))
```

where:
- Y_{t,r} = count of protests in district r on day t
- X_{t,r} = covariates (log population, poverty rate, CPI, year fixed effects)
- α = cross-district contagion effect
- Cross_{t,r}(θ) = spatially-weighted sum of protests in other districts over past 30 days

**Spatial Kernel:**

```
w_{rs}(θ) ∝ exp(−d_{rs}/θ)
```

- θ = characteristic decay distance (km)
- At distance θ, influence drops to ~37% of maximum
- θ = ∞ corresponds to uniform (national) spread

**Estimation:** Quasi-Poisson GLM with log link
- Profile likelihood over θ
- Dispersion-adjusted standard errors

---

## Reproducing the Analysis

### Quick Start

```bash
# Run scripts in order:
Rscript 02_prepare_protest_data.R
Rscript 10_prepare_population_data.R
Rscript 15_prepare_poverty_data.R
Rscript prepare_cpi_data.R
Rscript collapse_multiday_protests.R
Rscript estimate_spatial_hawkes_final.R
Rscript compute_branching_ratio.R
Rscript generate_presentation_figures.R
Rscript generate_spatial_figure.R
```

### File Dependencies

```
ACLED Data_2025-10-27_Indonesia20142024_OFFICIAL.xlsx
district_level_population_2014_2020.csv
indonesia_cpi_fred.csv
Indonesian district-level shapefile copy 2/
    ├── gadm41_IDN_2.shp
    ├── gadm41_IDN_2.dbf
    └── ... (other shapefile components)
```

### Output Files

```
protests_prepared.rds          # Step 1
protests_with_population.rds   # Step 2
protests_with_poverty.rds      # Step 3
indonesia_cpi.rds              # Step 4
protests_daily.rds             # Step 5
spatial_hawkes_final_results.rds  # Step 6
branching_ratio_results.rds    # Step 7
figures/                       # Step 8
    ├── fig_timeseries.pdf
    ├── fig_profile_likelihood.pdf
    ├── fig_lag_structure.pdf
    ├── fig_covariates.pdf
    └── ...
```

---

## Additional Scripts

The repository contains additional analysis scripts for extended investigations:

| Script | Description |
|--------|-------------|
| `estimate_discrete_hawkes.R` | Within-region temporal self-excitation analysis |
| `robustness_checks.R` | Sensitivity analysis with different specifications |
| `21_model1_basic_hawkes.R` - `46_bivariate_bootstrap.R` | Alternative model specifications (bivariate, trivariate) |
| `04_exploratory_analysis.R` | Knox test and exploratory spatial-temporal clustering |
| `18_model_diagnostics.R` | Residual analysis and model fit diagnostics |

---

## References

- Hawkes, A. G. (1971). Spectra of some self-exciting and mutually exciting point processes. *Biometrika*, 58(1), 83-90.
- Mohler, G. O., et al. (2011). Self-exciting point process modeling of crime. *Journal of the American Statistical Association*, 106(493), 100-108.
- Schoenberg, F. P., et al. (2019). A recursive algorithm for estimating hawkes process parameters. *Journal of Computational and Graphical Statistics*, 28(4), 778-790.
- ACLED. (2024). Armed Conflict Location & Event Data Project. https://acleddata.com/
