library(dplyr)
library(stringdist)

acled_names <- readLines('acled_admin2_names.txt')
pop_names <- readLines('population_district_names.txt')

# Full matching function
match_district <- function(acled_name, pop_names) {
  # Stage 1: Exact match
  if (acled_name %in% pop_names) {
    return(list(method = "exact", match = acled_name))
  }
  
  # Stage 2: Directional translation
  directional <- c("Central" = "Tengah", "East" = "Timur", "West" = "Barat",
                   "North" = "Utara", "South" = "Selatan")
  for (eng in names(directional)) {
    ind <- directional[eng]
    pattern <- paste0("^", eng, " (.+)$")
    if (grepl(pattern, acled_name)) {
      place <- sub(pattern, "\\1", acled_name)
      translated <- paste(place, ind)
      if (translated %in% pop_names) {
        return(list(method = "directional", match = translated))
      }
    }
  }
  
  # Stage 3: Special cases
  special_cases <- c(
    "Central Jakarta" = "Jakarta Pusat",
    "Thousand Islands" = "Adm. Kepulauan Seribu",
    "Bukit Tinggi" = "Bukittinggi",
    "Banyu Asin" = "Banyuasin",
    "Baubau" = "Bau-Bau"
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

# Test on ALL unmatched districts
unmatched <- readLines('acled_unmatched.txt')
results <- lapply(unmatched, function(x) match_district(x, pop_names))
names(results) <- unmatched

# Compile results
matched <- sapply(results, function(x) !is.na(x$match))
methods <- sapply(results, function(x) x$method)

cat("=== FULL MATCHING RESULTS ===\n\n")
cat(sprintf("Total originally unmatched: %d\n", length(unmatched)))
cat(sprintf("Successfully matched: %d (%.1f%%)\n", 
            sum(matched), 100*mean(matched)))
cat(sprintf("Still unmatched: %d (%.1f%%)\n\n", 
            sum(!matched), 100*mean(!matched)))

cat("Breakdown by method:\n")
cat(sprintf("  Exact: %d\n", sum(methods == "exact")))
cat(sprintf("  Directional: %d\n", sum(methods == "directional")))
cat(sprintf("  Special: %d\n", sum(methods == "special")))
cat(sprintf("  Islands: %d\n", sum(methods == "islands")))
cat(sprintf("  Fuzzy: %d\n", sum(methods == "fuzzy")))
cat(sprintf("  Not found: %d\n", sum(methods == "none")))

# Load protest data to calculate coverage
protests <- readRDS('protests_prepared.rds')

# Calculate matched + unmatched coverage
originally_unmatched_districts <- unmatched
newly_matched_districts <- unmatched[matched]
still_unmatched_districts <- unmatched[!matched]

protests_in_newly_matched <- sum(protests$admin2 %in% newly_matched_districts)
protests_in_still_unmatched <- sum(protests$admin2 %in% still_unmatched_districts)
protests_originally_matched <- sum(protests$admin2 %in% acled_names[acled_names %in% pop_names])

cat("\n=== PROTEST COVERAGE ===\n\n")
cat(sprintf("Originally exact-matched: %d protests (%.1f%%)\n",
            protests_originally_matched,
            100*protests_originally_matched/nrow(protests)))
cat(sprintf("Newly matched via rules: %d protests (%.1f%%)\n",
            protests_in_newly_matched,
            100*protests_in_newly_matched/nrow(protests)))
cat(sprintf("Total matched: %d protests (%.1f%%)\n",
            protests_originally_matched + protests_in_newly_matched,
            100*(protests_originally_matched + protests_in_newly_matched)/nrow(protests)))
cat(sprintf("Still unmatched: %d protests (%.1f%%)\n",
            protests_in_still_unmatched,
            100*protests_in_still_unmatched/nrow(protests)))

# Show still-unmatched districts
if (sum(!matched) > 0) {
  cat("\n=== STILL UNMATCHED DISTRICTS ===\n\n")
  still_unmatched_with_counts <- protests %>%
    filter(admin2 %in% still_unmatched_districts) %>%
    group_by(admin2) %>%
    summarize(n = n()) %>%
    arrange(desc(n))
  print(still_unmatched_with_counts, n = 20)
}
