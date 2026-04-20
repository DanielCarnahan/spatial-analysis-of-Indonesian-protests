################################################################################
#               GENERATE SPATIAL DISTRIBUTION FIGURE
################################################################################

library(tidyverse)
library(sf)
library(rnaturalearth)

cat("Generating spatial distribution figure...\n")

# Load data
protests <- readRDS("protests_daily.rds")

cat(sprintf("Total events: %d\n", nrow(protests)))
cat(sprintf("Longitude range: %.1f to %.1f\n", min(protests$longitude, na.rm=TRUE), max(protests$longitude, na.rm=TRUE)))
cat(sprintf("Latitude range: %.1f to %.1f\n", min(protests$latitude, na.rm=TRUE), max(protests$latitude, na.rm=TRUE)))

# Get Indonesia basemap
indonesia <- ne_countries(scale = "medium", country = "Indonesia", returnclass = "sf")

# Create dot map
p <- ggplot() +
  geom_sf(data = indonesia, fill = "gray95", color = "gray60", linewidth = 0.3) +
  geom_point(data = protests, aes(x = longitude, y = latitude),
             alpha = 0.25, size = 0.6, color = "#c0392b") +
  coord_sf(xlim = c(95, 141), ylim = c(-11, 6)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(x = NULL, y = NULL)

ggsave("figures/fig_spatial_distribution.pdf", p, width = 10, height = 4, dpi = 300)
cat("Saved: figures/fig_spatial_distribution.pdf\n")
