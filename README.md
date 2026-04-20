# Spatial Diffusion and Contagion Dynamics of Indonesian Protests

Analysis of protest contagion patterns using spatial point process models on ACLED data (2015-2024).

## Quick Start

### View Results
```r
# Load prepared data
protests <- readRDS("protests_prepared.rds")
protest_ppp <- readRDS("protest_ppp.rds")

# View summaries
summary(protests)
summary(protest_ppp)

# Check plots
list.files("plots/")
```

### Run Analysis
```bash
# 1. Explore raw data
Rscript 01_explore_data.R

# 2. Prepare protest data
Rscript 02_prepare_protest_data.R

# 3. Create spatial point patterns
Rscript 03_spatial_point_pattern.R

# 4. Exploratory analysis (generates plots/)
Rscript 04_exploratory_analysis.R
```

## Project Structure

```
spatial-indonesian-protests/
├── README.md                    # This file
├── PROJECT_SUMMARY.md           # Comprehensive documentation
│
├── Data/
│   ├── ACLED Data_*.xlsx        # Raw ACLED data
│   ├── protests_prepared.rds    # Cleaned protest data
│   ├── protests_sf.rds          # Spatial features version
│   ├── protest_ppp.rds          # Point pattern object
│   └── protest_stpp.rds         # Space-time data
│
├── Analysis Scripts/
│   ├── 01_explore_data.R        # Initial data exploration
│   ├── 02_prepare_protest_data.R # Data cleaning
│   ├── 03_spatial_point_pattern.R # Create ppp objects
│   └── 04_exploratory_analysis.R # Clustering tests
│
├── Results/
│   ├── plots/                   # Visualizations (8 plots)
│   ├── regional_summary.csv     # Province statistics
│   ├── yearly_summary.csv       # Annual trends
│   └── knox_test_result.rds     # Space-time interaction
│
└── Next Steps/
    ├── 05_baseline_models.R     # (To create) IPP/LGCP
    ├── 06_hawkes_models.R       # (To create) Temporal contagion
    └── 07_etas_models.R         # (To create) Spatial-temporal
```

## Key Findings (So Far)

✅ **16,467 protests** analyzed (2015-2024)
✅ **8.6x increase** in activity over time
✅ **Significant spatial clustering** detected
✅ **Jakarta** shows 278x national average intensity
✅ **East Java** most lethal (142 fatalities)

❓ **Space-time contagion** - needs self-exciting models

## Next Steps

See `PROJECT_SUMMARY.md` for detailed roadmap, including:

- Phase 4: Baseline models (IPP, LGCP)
- Phase 5: Self-exciting models (Hawkes, ETAS)
- Phase 6: Model validation
- Phase 7: Interpretation & policy implications

## Dependencies

```r
# Spatial analysis
library(spatstat)
library(sf)
library(dplyr)

# Visualization
library(ggplot2)
library(viridis)

# Data import
library(readxl)
library(lubridate)

# For self-exciting models (next phases):
# install.packages(c("stpp", "hawkes", "ppmlasso"))
```

## Data Source

**ACLED** (Armed Conflict Location & Event Data Project)
- Indonesia coverage: 2015-2024
- 18,573 total conflict events
- 16,467 protest-related events
- Source: https://acleddata.com

## Research Questions

1. Do protests exhibit spatial clustering?
2. Is there evidence of spatial contagion?
3. How does spatial influence decay with distance/time?
4. Do different protest types show different patterns?
5. Are there regional differences in susceptibility?

## Documentation

- Full methodology: `PROJECT_SUMMARY.md`
- Visualizations: `plots/` directory
- Script comments: Inline documentation

---

**Status:** Exploratory analysis complete ✓
**Next:** Self-exciting point process models
