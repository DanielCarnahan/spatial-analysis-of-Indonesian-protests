library(dplyr)

acled_names <- readLines('acled_admin2_names.txt')
pop_names <- readLines('population_district_names.txt')

# Test spacing/hyphenation patterns
test_spacing <- function(name) {
  variants <- c(
    name,
    gsub(" ", "", name),           # Remove all spaces
    gsub(" ", "-", name),          # Replace space with hyphen
    gsub("-", "", name),           # Remove hyphens
    gsub("-", " ", name)           # Replace hyphen with space
  )
  
  for (v in variants) {
    if (v %in% pop_names) {
      return(list(found = TRUE, match = v))
    }
  }
  return(list(found = FALSE, match = NA))
}

# Test Islands translation
test_islands <- function(name) {
  # Pattern: "X Islands" -> "Kepulauan X"
  if (grepl(" Islands$", name)) {
    base <- sub(" Islands$", "", name)
    variants <- c(
      paste("Kepulauan", base),
      paste(base, "Kepulauan")
    )
    for (v in variants) {
      if (v %in% pop_names) {
        return(list(found = TRUE, match = v))
      }
    }
  }
  return(list(found = FALSE, match = NA))
}

# Test problematic cases
problems <- c("Bukit Tinggi", "Banyu Asin", "Baubau", "Central Jakarta",
              "Anambas Islands", "Aru Islands", "Thousand Islands")

cat("Testing remaining patterns:\n\n")
for (name in problems) {
  # Try spacing
  result_space <- test_spacing(name)
  if (result_space$found) {
    cat(sprintf("%-25s -> %-25s [SPACING]\n", name, result_space$match))
    next
  }
  
  # Try islands
  result_islands <- test_islands(name)
  if (result_islands$found) {
    cat(sprintf("%-25s -> %-25s [ISLANDS]\n", name, result_islands$match))
    next
  }
  
  # Try Jakarta Pusat special case
  if (name == "Central Jakarta") {
    if ("Jakarta Pusat" %in% pop_names) {
      cat(sprintf("%-25s -> %-25s [SPECIAL]\n", name, "Jakarta Pusat"))
      next
    }
  }
  
  cat(sprintf("%-25s -> %-25s [NOT FOUND]\n", name, "???"))
}
