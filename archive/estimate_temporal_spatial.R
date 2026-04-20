################################################################################
#         TEMPORAL AND SPATIAL DIFFUSION ANALYSIS
#
#   Test 1: Aggregate Temporal Contagion
#           Do protests anywhere increase future protests anywhere?
#
#   Test 2: Spatial Decay in Triggering
#           Do nearby protests trigger more than distant protests?
#
################################################################################

library(tidyverse)
library(geosphere)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║   TEMPORAL AND SPATIAL DIFFUSION ANALYSIS                               ║\n")
cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

protests <- readRDS("protests_daily.rds")

cat(sprintf("Total events: %d\n", nrow(protests)))
cat(sprintf("Date range: %s to %s\n", min(protests$date), max(protests$date)))
cat(sprintf("Districts: %d\n\n", length(unique(protests$admin2))))

# =============================================================================
# TEST 1: AGGREGATE TEMPORAL CONTAGION
# =============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  TEST 1: AGGREGATE TEMPORAL CONTAGION                       ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Question: Do protests ANYWHERE increase future protests ANYWHERE?\n\n")

# Aggregate to national daily counts
daily_national <- protests %>%
  group_by(date) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(date)

# Fill in missing dates with zeros
all_dates <- seq.Date(min(daily_national$date), max(daily_national$date), by = "day")
daily_national <- data.frame(date = all_dates) %>%
  left_join(daily_national, by = "date") %>%
  mutate(count = replace_na(count, 0))

cat(sprintf("Daily time series: %d days\n", nrow(daily_national)))
cat(sprintf("Mean daily count: %.2f\n", mean(daily_national$count)))
cat(sprintf("Days with protests: %d (%.1f%%)\n\n",
            sum(daily_national$count > 0),
            100 * mean(daily_national$count > 0)))

# Add lags
MAX_LAG <- 7
for (lag in 1:MAX_LAG) {
  daily_national[[paste0("lag", lag)]] <- lag(daily_national$count, lag, default = 0)
}

# Filter to valid observations
model_data_agg <- daily_national %>% filter(row_number() > MAX_LAG)

# --- M0: Poisson Baseline (no lags) ---
cat("--- M0: Poisson Baseline (No Self-Excitation) ---\n\n")

m0_agg <- glm(count ~ 1,
              data = model_data_agg,
              family = quasi(link = "identity", variance = "mu"),
              start = c(mean(model_data_agg$count)))

cat(sprintf("Intercept (μ): %.4f\n", coef(m0_agg)))
cat(sprintf("Deviance: %.2f\n\n", deviance(m0_agg)))

# --- M1: Discrete Hawkes (with lags) ---
cat("--- M1: Discrete Hawkes (With Self-Excitation) ---\n\n")

m1_agg <- glm(count ~ lag1 + lag2 + lag3 + lag4 + lag5 + lag6 + lag7,
              data = model_data_agg,
              family = quasi(link = "identity", variance = "mu"),
              start = c(1, rep(0.05, 7)))

cat("Coefficients:\n")
print(round(coef(summary(m1_agg)), 5))

lag_coefs_agg <- coef(m1_agg)[paste0("lag", 1:7)]
br_agg <- sum(lag_coefs_agg)

cat(sprintf("\nBackground rate (μ): %.4f\n", coef(m1_agg)["(Intercept)"]))
cat(sprintf("Branching ratio (Σ α_ℓ): %.4f\n\n", br_agg))

# F-test
deviance_diff <- deviance(m0_agg) - deviance(m1_agg)
df_diff <- 7
scale_param <- deviance(m1_agg) / df.residual(m1_agg)
f_stat <- (deviance_diff / df_diff) / scale_param
p_value <- pf(f_stat, df_diff, df.residual(m1_agg), lower.tail = FALSE)

cat("F-Test (M1 vs M0):\n")
cat(sprintf("  F statistic: %.2f\n", f_stat))
cat(sprintf("  df: (%d, %d)\n", df_diff, df.residual(m1_agg)))
cat(sprintf("  p-value: %.2e\n", p_value))

if (p_value < 0.001) {
  cat("\n  *** Temporal contagion is significant ***\n")
}

# --- Permutation Test ---
cat("\n--- Permutation Test ---\n\n")

set.seed(42)
n_perm <- 500
br_null_agg <- numeric(n_perm)

cat(sprintf("Running %d permutations...\n", n_perm))
pb <- txtProgressBar(min = 0, max = n_perm, style = 3)

for (i in 1:n_perm) {
  # Shuffle the count column
  perm_data <- model_data_agg
  perm_data$count_perm <- sample(perm_data$count)

  # Recompute lags on permuted data
  for (lag in 1:MAX_LAG) {
    perm_data[[paste0("lag", lag, "_perm")]] <- lag(perm_data$count_perm, lag, default = 0)
  }

  # Fit model
  m_perm <- tryCatch({
    glm(count_perm ~ lag1_perm + lag2_perm + lag3_perm + lag4_perm +
          lag5_perm + lag6_perm + lag7_perm,
        data = perm_data,
        family = quasi(link = "identity", variance = "mu"),
        start = c(1, rep(0.05, 7)))
  }, error = function(e) NULL)

  if (!is.null(m_perm)) {
    br_null_agg[i] <- sum(coef(m_perm)[paste0("lag", 1:7, "_perm")])
  } else {
    br_null_agg[i] <- NA
  }

  setTxtProgressBar(pb, i)
}
close(pb)

br_null_agg <- br_null_agg[!is.na(br_null_agg)]
perm_p_agg <- mean(br_null_agg >= br_agg)

cat("\n\nPermutation Test Results:\n")
cat(sprintf("  Observed BR: %.4f\n", br_agg))
cat(sprintf("  Null mean: %.4f, sd: %.4f\n", mean(br_null_agg), sd(br_null_agg)))
cat(sprintf("  Null 99th percentile: %.4f\n", quantile(br_null_agg, 0.99)))
cat(sprintf("  p-value: %.4f\n", perm_p_agg))

if (perm_p_agg < 0.01) {
  cat("\n  *** Temporal contagion confirmed by permutation test ***\n")
}

# =============================================================================
# TEST 2: SPATIAL DECAY IN TRIGGERING
# =============================================================================

cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  TEST 2: SPATIAL DECAY IN TRIGGERING                        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("Question: Do NEARBY protests trigger more than DISTANT protests?\n\n")

# Get district centroids
district_centroids <- protests %>%
  group_by(admin2) %>%
  summarise(
    lon = mean(longitude, na.rm = TRUE),
    lat = mean(latitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(lon) & !is.na(lat))

cat(sprintf("Districts with coordinates: %d\n", nrow(district_centroids)))

# Create district-day panel
all_dates <- seq.Date(min(protests$date), max(protests$date), by = "day")
all_districts <- unique(district_centroids$admin2)

district_counts <- protests %>%
  group_by(admin2, date) %>%
  summarise(count = n(), .groups = "drop")

panel <- expand.grid(
  admin2 = all_districts,
  date = all_dates,
  stringsAsFactors = FALSE
) %>%
  left_join(district_counts, by = c("admin2", "date")) %>%
  mutate(count = replace_na(count, 0)) %>%
  left_join(district_centroids, by = "admin2") %>%
  arrange(admin2, date)

cat(sprintf("Panel: %d district-days\n\n", nrow(panel)))

# For each event, compute distance-weighted exposure to past events
# Define distance bands: 0-50km (same district proxy), 50-200km, 200-500km, 500km+

cat("Computing distance-based exposure to past protests...\n")
cat("Distance bands: Own district, 0-100km, 100-300km, 300km+\n\n")

# Compute pairwise distances between district centroids
n_districts <- nrow(district_centroids)
dist_matrix <- matrix(0, n_districts, n_districts)
rownames(dist_matrix) <- district_centroids$admin2
colnames(dist_matrix) <- district_centroids$admin2

for (i in 1:n_districts) {
  for (j in 1:n_districts) {
    if (i != j) {
      dist_matrix[i, j] <- distHaversine(
        c(district_centroids$lon[i], district_centroids$lat[i]),
        c(district_centroids$lon[j], district_centroids$lat[j])
      ) / 1000  # Convert to km
    }
  }
}

cat(sprintf("Mean inter-district distance: %.0f km\n", mean(dist_matrix[dist_matrix > 0])))
cat(sprintf("Median inter-district distance: %.0f km\n\n", median(dist_matrix[dist_matrix > 0])))

# Create spatial weight matrices for different distance bands
# Band 1: 0-100 km
# Band 2: 100-300 km
# Band 3: 300+ km

W_band1 <- (dist_matrix > 0 & dist_matrix <= 100) * 1
W_band2 <- (dist_matrix > 100 & dist_matrix <= 300) * 1
W_band3 <- (dist_matrix > 300) * 1

# Row-normalize
W_band1 <- W_band1 / pmax(rowSums(W_band1), 1)
W_band2 <- W_band2 / pmax(rowSums(W_band2), 1)
W_band3 <- W_band3 / pmax(rowSums(W_band3), 1)

# Create wide format for each date
panel_wide <- panel %>%
  select(admin2, date, count) %>%
  pivot_wider(names_from = admin2, values_from = count, values_fill = 0)

# Compute spatially-weighted lags for each band
# We'll use a single aggregate lag (sum of days 1-7) for simplicity

cat("Computing spatial exposure variables...\n")

# Initialize columns
panel$own_lag <- 0
panel$near_lag <- 0      # 0-100 km
panel$medium_lag <- 0    # 100-300 km
panel$far_lag <- 0       # 300+ km

# For each district
for (d in all_districts) {
  idx <- which(panel$admin2 == d & !is.na(panel$admin2))
  if (length(idx) == 0) next
  district_data <- panel[idx, ]

  # Own-district lagged sum (days 1-7)
  own_counts <- district_data$count
  own_lag_sum <- rep(0, length(own_counts))
  for (lag in 1:7) {
    own_lag_sum <- own_lag_sum + lag(own_counts, lag, default = 0)
  }
  panel$own_lag[idx] <- own_lag_sum

  # Neighbor exposure (weighted average of neighbor counts, summed over lags 1-7)
  if (d %in% rownames(W_band1)) {
    w1 <- W_band1[d, ]
    w2 <- W_band2[d, ]
    w3 <- W_band3[d, ]

    neighbor_names <- names(w1)[names(w1) %in% names(panel_wide)[-1]]

    if (length(neighbor_names) > 0) {
      neighbor_counts <- as.matrix(panel_wide[, neighbor_names])

      # Compute weighted sums for each band
      near_counts <- neighbor_counts %*% w1[neighbor_names]
      medium_counts <- neighbor_counts %*% w2[neighbor_names]
      far_counts <- neighbor_counts %*% w3[neighbor_names]

      # Sum over lags 1-7
      near_lag <- rep(0, nrow(neighbor_counts))
      medium_lag <- rep(0, nrow(neighbor_counts))
      far_lag <- rep(0, nrow(neighbor_counts))

      for (lag in 1:7) {
        near_lag <- near_lag + c(rep(0, lag), near_counts[1:(length(near_counts) - lag)])
        medium_lag <- medium_lag + c(rep(0, lag), medium_counts[1:(length(medium_counts) - lag)])
        far_lag <- far_lag + c(rep(0, lag), far_counts[1:(length(far_counts) - lag)])
      }

      panel$near_lag[idx] <- near_lag
      panel$medium_lag[idx] <- medium_lag
      panel$far_lag[idx] <- far_lag
    }
  }
}

cat("Done.\n\n")

# Filter to valid observations
model_data_spatial <- panel %>% filter(row_number() > 7 * length(all_districts))

cat(sprintf("Observations: %d\n\n", nrow(model_data_spatial)))

# --- Model with distance bands ---
cat("--- Spatial Decay Model ---\n\n")

m_spatial <- glm(count ~ own_lag + near_lag + medium_lag + far_lag,
                 data = model_data_spatial,
                 family = quasi(link = "identity", variance = "mu"),
                 start = c(0.001, 0.01, 0.01, 0.01, 0.01))

cat("Coefficients:\n")
spatial_coefs <- coef(summary(m_spatial))
print(round(spatial_coefs, 6))

cat("\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("SPATIAL DECAY RESULTS:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("Own district (7-day sum):     %.5f (SE: %.5f)\n",
            spatial_coefs["own_lag", 1], spatial_coefs["own_lag", 2]))
cat(sprintf("Near (0-100 km):              %.5f (SE: %.5f)\n",
            spatial_coefs["near_lag", 1], spatial_coefs["near_lag", 2]))
cat(sprintf("Medium (100-300 km):          %.5f (SE: %.5f)\n",
            spatial_coefs["medium_lag", 1], spatial_coefs["medium_lag", 2]))
cat(sprintf("Far (300+ km):                %.5f (SE: %.5f)\n",
            spatial_coefs["far_lag", 1], spatial_coefs["far_lag", 2]))

# Test for spatial decay
own_coef <- spatial_coefs["own_lag", 1]
near_coef <- spatial_coefs["near_lag", 1]
medium_coef <- spatial_coefs["medium_lag", 1]
far_coef <- spatial_coefs["far_lag", 1]

cat("\n")
if (own_coef > near_coef && near_coef > medium_coef && medium_coef > far_coef) {
  cat("✓ SPATIAL DECAY DETECTED: Own > Near > Medium > Far\n")
} else if (own_coef > near_coef) {
  cat("○ PARTIAL DECAY: Own > Near, but pattern not monotonic\n")
} else {
  cat("⚠ NO CLEAR SPATIAL DECAY\n")
}

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SAVING RESULTS                                              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

results <- list(
  aggregate = list(
    n_days = nrow(model_data_agg),
    mean_daily_count = mean(model_data_agg$count),
    baseline_mu = coef(m0_agg),
    hawkes_mu = coef(m1_agg)["(Intercept)"],
    lag_coefficients = lag_coefs_agg,
    branching_ratio = br_agg,
    f_test = list(
      statistic = f_stat,
      df1 = df_diff,
      df2 = df.residual(m1_agg),
      p_value = p_value
    ),
    permutation = list(
      observed = br_agg,
      null_mean = mean(br_null_agg),
      null_sd = sd(br_null_agg),
      null_99 = as.numeric(quantile(br_null_agg, 0.99)),
      p_value = perm_p_agg,
      null_distribution = br_null_agg
    )
  ),
  spatial = list(
    n_districts = length(all_districts),
    n_obs = nrow(model_data_spatial),
    coefficients = spatial_coefs,
    own_effect = own_coef,
    near_effect = near_coef,
    medium_effect = medium_coef,
    far_effect = far_coef,
    decay_detected = own_coef > near_coef
  )
)

saveRDS(results, "temporal_spatial_results.rds")
cat("Saved: temporal_spatial_results.rds\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  SUMMARY                                                     ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("TEST 1: AGGREGATE TEMPORAL CONTAGION\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Branching ratio: %.3f\n", br_agg))
cat(sprintf("  Permutation p-value: %.4f\n", perm_p_agg))
cat(sprintf("  Conclusion: %s\n\n",
            ifelse(perm_p_agg < 0.01, "Temporal contagion EXISTS", "No evidence of contagion")))

cat("TEST 2: SPATIAL DECAY\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Own district effect:  %.5f\n", own_coef))
cat(sprintf("  Near (0-100km):       %.5f\n", near_coef))
cat(sprintf("  Medium (100-300km):   %.5f\n", medium_coef))
cat(sprintf("  Far (300km+):         %.5f\n", far_coef))
cat(sprintf("  Conclusion: %s\n",
            ifelse(own_coef > near_coef, "Nearby protests trigger MORE", "No spatial decay")))
