library(stringdist)

pop_names <- readLines('population_district_names.txt')

# The 13 remaining unmatched districts
remaining <- c(
  "North Central Timor", "Tidore Kepulauan", "Southeast Aceh", "Southwest Aceh",
  "Pangkajene and Islands", "Mempawah", "South Central Timor", 
  "Sawahlunto Sijunjung", "Southwest Sumba", "Lima Puluh",
  "Sitaro Islands", "Southeast Minahasa", "Southwest Maluku"
)

cat("Searching for close matches in population data:\n\n")

for (name in remaining) {
  # Calculate distances
  distances <- stringdist(name, pop_names, method = "lv")
  top3 <- order(distances)[1:3]
  
  cat(sprintf("%-30s (protests: %d)\n", name, 
              ifelse(name == "North Central Timor", 22,
              ifelse(name == "Tidore Kepulauan", 13,
              ifelse(name %in% c("Southeast Aceh", "Southwest Aceh"), 11,
              ifelse(name == "Pangkajene and Islands", 10,
              ifelse(name == "Mempawah", 9,
              ifelse(name == "South Central Timor", 7, 3))))))))
  
  for (i in 1:3) {
    idx <- top3[i]
    cat(sprintf("  %d. %-40s (dist: %d)\n", 
                i, pop_names[idx], distances[idx]))
  }
  cat("\n")
}
