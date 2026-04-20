library(dplyr)
library(stringdist)

acled_names <- readLines('acled_admin2_names.txt')
pop_names <- readLines('population_district_names.txt')

# COMPLETE matching function with ALL discovered patterns
match_district <- function(acled_name, pop_names) {
  # Stage 1: Exact match
  if (acled_name %in% pop_names) {
    return(list(method = "exact", match = acled_name))
  }
  
  # Stage 2: Directional translation (EXPANDED)
  directional <- list(
    # Simple directional
    c("Central", "Tengah"), c("East", "Timur"), c("West", "Barat"),
    c("North", "Utara"), c("South", "Selatan"),
    # Compound directional
    c("Southeast", "Tenggara"), c("Southwest", "Barat Daya"),
    c("Northeast", "Timur Laut"), c("Northwest", "Barat Laut")
  )
  
  for (pair in directional) {
    eng <- pair[1]
    ind <- pair[2]
    pattern <- paste0("^", eng, " (.+)$")
    if (grepl(pattern, acled_name)) {
      place <- sub(pattern, "\\1", acled_name)
      translated <- paste(place, ind)
      if (translated %in% pop_names) {
        return(list(method = "directional", match = translated))
      }
    }
  }
  
  # Stage 2b: Multi-word directional (e.g., "North Central X" → "X Tengah Utara")
  multi_directional <- list(
    c("North Central", "Tengah Utara"),
    c("South Central", "Tengah Selatan"),
    c("East Central", "Tengah Timur"),
    c("West Central", "Tengah Barat")
  )
  
  for (pair in multi_directional) {
    eng <- pair[1]
    ind <- pair[2]
    pattern <- paste0("^", eng, " (.+)$")
    if (grepl(pattern, acled_name)) {
      place <- sub(pattern, "\\1", acled_name)
      translated <- paste(place, ind)
      if (translated %in% pop_names) {
        return(list(method = "multi_directional", match = translated))
      }
    }
  }
  
  # Stage 3: Comprehensive special cases
  special_cases <- c(
    "Central Jakarta" = "Jakarta Pusat",
    "Thousand Islands" = "Adm. Kepulauan Seribu",
    "Bukit Tinggi" = "Bukittinggi",
    "Banyu Asin" = "Banyuasin",
    "Baubau" = "Bau-Bau",
    "Lima Puluh" = "Limapuluh Kota",
    "Pangkajene and Islands" = "Pangkajene Kepulauan",
    "Tidore Kepulauan" = "Kepulauan Tidore"
  )
  
  if (acled_name %in% names(special_cases)) {
    match_name <- special_cases[acled_name]
    if (match_name %in% pop_names) {
      return(list(method = "special", match = match_name))
    }
  }
  
  # Stage 4: Islands translation  
  if (grepl(" Islands$", acled_name)) {
    base <- sub(" Islands$", "", acled_name)
    variants <- c(paste("Kepulauan", base), paste(base, "Kepulauan"))
    for (v in variants) {
      if (v %in% pop_names) {
        return(list(method = "islands", match = v))
      }
    }
  }
  
  # Stage 5: Fuzzy matching
  distances <- stringdist(acled_name, pop_names, method = "lv")
  min_dist <- min(distances)
  if (min_dist <= 3) {
    best_match <- pop_names[which.min(distances)]
    return(list(method = "fuzzy", match = best_match, distance = min_dist))
  }
  
  return(list(method = "none", match = NA))
}

# Test on ALL ACLED districts
all_results <- lapply(acled_names, function(x) match_district(x, pop_names))
names(all_results) <- acled_names

matched <- sapply(all_results, function(x) !is.na(x$match))
methods <- sapply(all_results, function(x) x$method)

cat("=== FINAL COMPREHENSIVE MATCHING RESULTS ===\n\n")
cat(sprintf("Total ACLED districts: %d\n", length(acled_names)))
cat(sprintf("Successfully matched: %d (%.1f%%)\n", 
            sum(matched), 100*mean(matched)))
cat(sprintf("Unmatched: %d (%.1f%%)\n\n", 
            sum(!matched), 100*mean(!matched)))

cat("Method breakdown:\n")
for (m in unique(methods)) {
  cat(sprintf("  %-20s: %3d\n", m, sum(methods == m)))
}

# Protest coverage
protests <- readRDS('protests_prepared.rds')
matched_districts <- acled_names[matched]
unmatched_districts <- acled_names[!matched]

protests_matched <- sum(protests$admin2 %in% matched_districts, na.rm = TRUE)
protests_unmatched <- sum(protests$admin2 %in% unmatched_districts, na.rm = TRUE)
protests_na <- sum(is.na(protests$admin2))

cat("\n=== FINAL PROTEST COVERAGE ===\n\n")
cat(sprintf("Matched districts: %d protests (%.1f%%)\n",
            protests_matched, 100*protests_matched/nrow(protests)))
cat(sprintf("Unmatched districts: %d protests (%.1f%%)\n",
            protests_unmatched, 100*protests_unmatched/nrow(protests)))
cat(sprintf("Missing admin2: %d protests (%.1f%%)\n",
            protests_na, 100*protests_na/nrow(protests)))

# Show remaining unmatched
if (sum(!matched) > 0) {
  cat("\n=== REMAINING UNMATCHED DISTRICTS ===\n\n")
  unmatched_with_counts <- protests %>%
    filter(admin2 %in% unmatched_districts) %>%
    group_by(admin2) %>%
    summarize(n = n()) %>%
    arrange(desc(n))
  print(unmatched_with_counts, n = 20)
}
