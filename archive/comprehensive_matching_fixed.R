library(dplyr)
library(stringdist)

acled_names <- readLines('acled_admin2_names.txt')
pop_names <- readLines('population_district_names.txt')

# Multi-stage matching function
match_district <- function(acled_name, pop_names) {
  
  # Stage 1: Exact match
  if (acled_name %in% pop_names) {
    return(list(method = "exact", match = acled_name, confidence = 1.0))
  }
  
  # Stage 2: Directional translation
  directional <- c(
    "Central" = "Tengah", "East" = "Timur", "West" = "Barat",
    "North" = "Utara", "South" = "Selatan"
  )
  
  for (eng in names(directional)) {
    ind <- directional[eng]
    pattern <- paste0("^", eng, " (.+)$")
    if (grepl(pattern, acled_name)) {
      place <- sub(pattern, "\\1", acled_name)
      translated <- paste(place, ind)
      if (translated %in% pop_names) {
        return(list(method = "directional", match = translated, confidence = 0.95))
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
      return(list(method = "special", match = match_name, confidence = 0.95))
    }
  }
  
  # Stage 4: Islands translation
  if (grepl(" Islands$", acled_name)) {
    base <- sub(" Islands$", "", acled_name)
    variants <- c(paste("Kepulauan", base), paste(base, "Kepulauan"))
    for (v in variants) {
      if (v %in% pop_names) {
        return(list(method = "islands", match = v, confidence = 0.9))
      }
    }
  }
  
  # Stage 5: Fuzzy matching (last resort)
  distances <- stringdist(acled_name, pop_names, method = "lv")
  min_dist <- min(distances)
  if (min_dist <= 3) {
    best_match <- pop_names[which.min(distances)]
    confidence <- 1 - (min_dist / nchar(acled_name))
    return(list(method = "fuzzy", match = best_match, 
                confidence = confidence, distance = min_dist))
  }
  
  return(list(method = "none", match = NA, confidence = 0))
}

# Test on unmatched districts
unmatched <- readLines('acled_unmatched.txt')
test_set <- head(unmatched, 30)

cat("Testing comprehensive matching on first 30 unmatched districts:\n\n")

results <- list()
for (name in test_set) {
  result <- match_district(name, pop_names)
  results[[name]] <- result
  
  match_str <- ifelse(is.na(result$match), "NOT FOUND", result$match)
  cat(sprintf("%-30s | %-10s | %-30s\n",
              substr(name, 1, 30),
              result$method,
              substr(match_str, 1, 30)))
}

# Summary statistics
methods <- sapply(results, function(x) x$method)
cat("\n\nSummary:\n")
cat(sprintf("  Exact: %d\n", sum(methods == "exact")))
cat(sprintf("  Directional: %d\n", sum(methods == "directional")))
cat(sprintf("  Special: %d\n", sum(methods == "special")))
cat(sprintf("  Islands: %d\n", sum(methods == "islands")))
cat(sprintf("  Fuzzy: %d\n", sum(methods == "fuzzy")))
cat(sprintf("  Not found: %d\n", sum(methods == "none")))
