library(dplyr)

# Load data
acled_names <- readLines('acled_admin2_names.txt')
pop_names <- readLines('population_district_names.txt')

# Build directional translation dictionary
directional_map <- list(
  english_to_indonesian = c(
    "Central" = "Tengah",
    "East" = "Timur",
    "West" = "Barat",
    "North" = "Utara",
    "South" = "Selatan"
  ),
  special_central = c(
    "Central Jakarta" = "Jakarta Pusat"  # Special case
  )
)

# Try translating top unmatched districts
test_cases <- c(
  "Central Jakarta", "South Jakarta", "East Jakarta", 
  "West Jakarta", "North Jakarta",
  "Central Aceh", "East Aceh", "West Aceh", "North Aceh",
  "Central Lombok", "East Lombok",
  "Bukit Tinggi", "Banyu Asin", "Baubau"
)

cat("Testing translations:\n\n")
for (name in test_cases) {
  # Try directional translation
  translated <- name
  for (eng in names(directional_map$english_to_indonesian)) {
    ind <- directional_map$english_to_indonesian[eng]
    # Pattern: "Direction Place" -> "Place Direction-Indonesian"
    pattern <- paste0("^", eng, " (.+)$")
    if (grepl(pattern, name)) {
      place <- sub(pattern, "\\1", name)
      translated <- paste(place, ind)
      break
    }
  }
  
  # Check if it exists in population data
  match_found <- translated %in% pop_names
  
  cat(sprintf("%-20s -> %-20s [%s]\n", 
              name, translated, 
              ifelse(match_found, "FOUND", "NOT FOUND")))
}
