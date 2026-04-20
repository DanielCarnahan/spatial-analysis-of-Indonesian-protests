library(dplyr)

# Load population data
pop_raw <- read.csv('district_level_population_2014_2020.csv', stringsAsFactors = FALSE)
pop_data <- pop_raw %>%
  filter(Series.Name == 'Total Population (in number of people)') %>%
  select(district = Provinces.Name, pop_2020 = X2020..YR2020.)

pop_data$district <- gsub('^"', '', pop_data$district)
pop_data$pop_2020 <- as.numeric(pop_data$pop_2020)

# Extract base name and type
pop_data$base_name <- sapply(strsplit(pop_data$district, ','), `[`, 1)
pop_data$base_name <- trimws(pop_data$base_name)
pop_data$type <- ifelse(grepl(', Kota', pop_data$district), 'Kota',
                        ifelse(grepl(', Kab\\.', pop_data$district), 'Kab', 'Other'))

# Find duplicated base names (Kota/Kab pairs)
dup_bases <- pop_data %>%
  group_by(base_name) %>%
  filter(n() > 1) %>%
  arrange(base_name, desc(type))

cat('=== DISTRICTS WITH BOTH KOTA AND KAB VERSIONS ===\n\n')
print(as.data.frame(dup_bases %>% select(district, base_name, type, pop_2020)), row.names = FALSE)
cat('\n\nTotal duplicated base names:', length(unique(dup_bases$base_name)), '\n\n')

# Now check ACLED names for these same districts
protests <- readRDS('protests_prepared.rds')

cat('=== HOW ACLED RECORDS THESE DISTRICTS ===\n\n')
for(base in unique(dup_bases$base_name)) {
  # Find ACLED protests in districts with this base name
  matching_protests <- protests %>%
    filter(grepl(paste0('^', base, '$'), admin2) |
           grepl(paste0('^', base, ' '), admin2) |
           grepl(paste0(' ', base, '$'), admin2))

  if(nrow(matching_protests) > 0) {
    admin2_vals <- unique(matching_protests$admin2)
    cat(sprintf('%s: %d protests\n', base, nrow(matching_protests)))
    cat(sprintf('  ACLED names: %s\n\n', paste(admin2_vals, collapse=', ')))
  }
}
