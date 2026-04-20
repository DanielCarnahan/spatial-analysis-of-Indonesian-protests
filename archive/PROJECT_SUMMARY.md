# Spatial Diffusion and Contagion Dynamics of Indonesian Protests
## Project Summary and Next Steps

**Date:** October 27, 2025
**Data Source:** ACLED (Armed Conflict Location & Event Data Project)
**Coverage:** Indonesia, 2015-2024

---

## Research Questions

1. Do protests exhibit spatial clustering beyond what would be expected by chance?
2. Is there evidence of spatial contagion (protests triggering nearby protests)?
3. How does the spatial influence decay with distance and time?
4. Do different types of protests show different contagion patterns?
5. Are there regional differences in protest susceptibility?

---

## Data Overview

### Dataset Characteristics
- **Total events:** 16,467 protest events (Protests + Riots)
- **Temporal coverage:** January 1, 2015 - October 27, 2024 (9.82 years)
- **Spatial coverage:** All of Indonesia (38 provinces, 459 districts)
- **Unique locations:** 2,242
- **Completeness:** 0 missing coordinates (100% spatially complete)

### Event Classification
- **Peaceful protests:** 14,235 (86.4%)
- **Violent protests/riots:** 2,232 (13.6%)
- **Events with fatalities:** 152 (0.9%)
- **Total fatalities:** 261 deaths

### Temporal Pattern
| Period | Events | Violent | Peaceful |
|--------|--------|---------|----------|
| 2015-2017 | 1,088 | 337 | 751 |
| 2018-2019 | 2,314 | 385 | 1,929 |
| 2020-2021 | 3,660 | 652 | 3,008 |
| 2022-2024 | 9,405 | 858 | 8,547 |

**Key Finding:** 8.6x increase from 2015-2017 to 2022-2024

### Regional Distribution (Top 5)
1. **Jakarta:** 2,191 events (8.9% violent, 13 fatalities)
2. **East Java:** 1,772 events (9.3% violent, 142 fatalities)
3. **West Java:** 1,639 events (10.1% violent, 6 fatalities)
4. **South Sulawesi:** 1,037 events (17.5% violent, 5 fatalities)
5. **Central Java:** 980 events (8.6% violent, 0 fatalities)

**Notable:** East Java has far more fatalities than Jakarta despite fewer events

---

## Completed Analysis (Phase 1-3)

### 1. Data Preparation ✓
**Script:** `02_prepare_protest_data.R`

- Filtered ACLED data to protest-related events
- Classified by violence level and state intervention
- Created temporal variables (days since start, year-month)
- Generated regional and yearly summaries
- Saved cleaned datasets: `protests_prepared.rds`, `protests_sf.rds`

### 2. Spatial Point Pattern Objects ✓
**Script:** `03_spatial_point_pattern.R`

- Created `spatstat` point pattern objects (ppp)
- Built space-time point pattern data frame
- Generated subsets by:
  - Protest type (peaceful, violent, fatal)
  - Region (Jakarta case study)
  - Time period
- Calculated baseline intensities:
  - Spatial: 0.175 events/100km²
  - Temporal: 139.74 events/month
  - Jakarta intensity: 278x national average

### 3. Exploratory Spatial Analysis ✓
**Script:** `04_exploratory_analysis.R`

Generated 8 diagnostic plots in `plots/` directory:

#### Spatial Patterns
- **Kernel density maps** (`01_kernel_density.png`)
  - Max density: 1,888 events/unit area
  - Hotspots: Jakarta, Surabaya (East Java), Bandung (West Java)

- **Peaceful vs Violent comparison** (`02_peaceful_vs_violent_density.png`)
  - Similar spatial distributions
  - Violent protests more concentrated in Papua, South Sulawesi

#### Clustering Tests
- **Ripley's K-function** (`03_ripley_k_function.png`)
  - Tests: K(r) vs πr² for Complete Spatial Randomness (CSR)
  - Result: **Strong clustering detected** (K(r) >> πr²)

- **L-function** (`04_ripley_l_function.png`)
  - Variance-stabilized version of K
  - Result: **Consistent clustering across all scales**

- **Pair correlation function** (`05_pair_correlation.png`)
  - g(r) > 1 indicates clustering at distance r
  - Result: **Scale-dependent clustering pattern**

#### Spatial Homogeneity
- **Quadrat test** (`06_quadrat_counts.png`)
  - Chi-squared = 134,791, df = 99, **p < 2.2e-16**
  - Result: **Highly significant spatial heterogeneity**
  - Rejects Complete Spatial Randomness

- **Nearest neighbor distances** (`07_nearest_neighbor_G.png`)
  - G-function distribution
  - Result: Shorter distances than expected under CSR

#### Space-Time Interaction
- **Knox Test** (on 1,000 event sample)
  - Spatial threshold: 0.5° (~55 km)
  - Temporal threshold: 30 days
  - Observed close pairs: 595
  - Expected under randomness: 578
  - Ratio (Obs/Exp): 1.03
  - P-value: 0.32 (not significant)

  **Interpretation:** No strong evidence for space-time clustering in sample, but:
  - Limited to 1,000 events (computational constraints)
  - Different thresholds may reveal patterns
  - More sophisticated models needed

#### Temporal Pattern
- **Time series** (`08_temporal_pattern.png`)
  - Clear upward trend 2015-2024
  - Major spike in 2018-2019 (likely election-related)
  - Peak activity 2022-2023

---

## Key Findings Summary

### ✅ Confirmed
1. **Spatial clustering is significant** - protests are NOT randomly distributed
2. **Extreme spatial heterogeneity** - huge variation in protest intensity across regions
3. **Temporal acceleration** - protest frequency increasing dramatically over time
4. **Hotspot concentration** - Jakarta shows 278x national average intensity
5. **Scale-dependent patterns** - clustering varies by spatial scale

### ❓ Requires Further Investigation
1. **Space-time contagion** - Knox test inconclusive, need better models
2. **Triggering mechanisms** - what causes nearby/subsequent protests?
3. **Decay functions** - how do spatial/temporal effects diminish?
4. **Violence escalation** - do peaceful protests trigger violent ones (or vice versa)?
5. **Regional differences** - why does East Java have more fatalities?

---

## Next Steps: Self-Exciting Point Process Models

### Phase 4: Baseline Models

#### A. Inhomogeneous Poisson Process (IPP)
**Purpose:** Model spatial heterogeneity without temporal dependence

**Approach:**
```r
# Fit IPP with spatial covariates
library(spatstat)

# Option 1: Polynomial trend
fit_ipp_trend <- ppm(protest_ppp ~ polynom(x, y, 2))

# Option 2: With regional covariates (if available)
# - Population density
# - Economic indicators
# - Distance to major cities
```

**Outputs:**
- Background intensity λ(s)
- Spatial trend surface
- Model diagnostics (residuals, AIC)

#### B. Log-Gaussian Cox Process (LGCP)
**Purpose:** Model spatial intensity with latent Gaussian field

**Approach:**
```r
library(spatstat)

# Fit LGCP
fit_lgcp <- lgcp.estK(protest_ppp,
                      sigma = NULL,  # Estimate from data
                      nu = NULL)

# Or with covariates
fit_lgcp_cov <- kppm(protest_ppp ~ x + y,
                     cluster = "LGCP")
```

**Outputs:**
- Smooth intensity surface
- Variance of latent field
- Clustering scale

---

### Phase 5: Self-Exciting Models (CONTAGION)

#### A. Hawkes Process
**Purpose:** Model temporal self-excitation (event triggers events)

**Key equation:**
```
λ(t) = μ + Σ_{t_i < t} α * exp(-β(t - t_i))
```

Where:
- μ = background rate
- α = triggering intensity
- β = temporal decay rate

**Implementation:**
```r
library(hawkes)
# or library(ppmlasso)

# Fit temporal Hawkes
fit_hawkes <- hawkes(protest_stpp$t,
                     kernel = "exp",
                     lambda0 = NULL)  # Estimate background

# Interpret:
# - If α/β > 1: "supercritical" (explosive growth)
# - If α/β < 1: "subcritical" (stable)
```

#### B. Spatial-Temporal Self-Exciting Point Process (STPP)
**Purpose:** Model both spatial AND temporal contagion

**Key equation:**
```
λ(s,t) = μ(s) + Σ_{(s_i,t_i): t_i < t} g(s - s_i, t - t_i)
```

Where:
- μ(s) = background spatial intensity
- g(·) = triggering kernel with spatial and temporal decay

**Common triggering kernels:**

1. **Separable (spatial × temporal):**
```r
g(s, t) = α * f_space(s) * f_time(t)

# Example:
f_space(s) = (1/2πσ²) * exp(-||s||²/2σ²)  # Gaussian
f_time(t) = β * exp(-βt)                   # Exponential
```

2. **ETAS (Epidemic-Type Aftershock Sequence):**
```r
g(s, t) = α * (t + c)^(-p) * (||s||² + d)^(-q)

# From seismology, adapted for protests
# Parameters: α, c, p, d, q
```

**Implementation:**
```r
library(stpp)

# Fit ETAS model
fit_etas <- etas(xyt = cbind(protest_stpp$x,
                             protest_stpp$y,
                             protest_stpp$t),
                 params = list(...),
                 method = "MLE")
```

#### C. Multi-Type Self-Exciting Process
**Purpose:** Model contagion BETWEEN protest types

**Key questions:**
- Do peaceful protests trigger violent ones?
- Do violent protests suppress peaceful ones?
- Cross-excitation patterns?

**Implementation:**
```r
library(ppmlasso)

# Separate by type
peaceful_times <- protest_stpp$t[protest_stpp$is_peaceful]
violent_times <- protest_stpp$t[protest_stpp$is_violent]

# Fit multivariate Hawkes
fit_multi <- muHawkes(list(peaceful_times, violent_times),
                      kernel = "exp")

# Extract cross-excitation matrix
# alpha[i,j] = how much type i triggers type j
```

---

### Phase 6: Model Comparison & Validation

#### Model Selection
```r
# Compare models using:
# 1. AIC/BIC
AIC(fit_ipp, fit_lgcp, fit_hawkes, fit_etas)

# 2. Residual analysis
# - Pearson residuals
# - Rescaled residuals (should be Poisson(1))

# 3. Predictive performance
# - Hold-out validation
# - Residual K-function
```

#### Diagnostics
```r
# 1. Residual plots
res <- residuals(fit_etas)
plot(res)

# 2. Q-Q plot
diagnose.etas(fit_etas)

# 3. Temporal rescaling
# - Transform to Poisson process
# - Check for uniformity
```

---

### Phase 7: Interpretation & Policy Implications

#### Extract Key Parameters

1. **Background rate μ(s)**
   - Where are "endemic" protest regions?
   - What drives baseline protest activity?

2. **Triggering intensity α**
   - How strong is the contagion effect?
   - Do protests cascade?

3. **Spatial decay σ or d**
   - How far does influence spread? (km)
   - Urban vs rural differences?

4. **Temporal decay β or p**
   - How long does influence persist? (days)
   - Memory of protest activity?

5. **Criticality α/β**
   - Is the system stable or explosive?
   - Risk of protest waves?

#### Visualizations

```r
# 1. Fitted intensity surface
plot(predict(fit_etas), type = "intensity")

# 2. Triggering function
plot_triggering_kernel(fit_etas)

# 3. Branching structure
# - Which protests triggered which?
# - Identify "parent" and "offspring" events

# 4. Counterfactual analysis
# - What if certain protests didn't happen?
# - Intervention simulations
```

#### Policy Questions

1. **Prevention:**
   - Which locations have highest background risk?
   - Target resources for conflict prevention

2. **Early warning:**
   - After a protest, predict nearby subsequent events
   - Temporal windows of elevated risk

3. **De-escalation:**
   - Do police interventions increase/decrease contagion?
   - Compare peaceful vs violent protest cascades

4. **Regional patterns:**
   - Why is East Java more deadly?
   - What makes Jakarta so protest-prone?

---

## Recommended R Packages

### Essential
- `spatstat` - Spatial point pattern analysis (already installed)
- `stpp` - Spatio-temporal point processes
- `hawkes` - Hawkes process estimation
- `ppmlasso` - Penalized point process models

### Optional
- `etasFLP` - ETAS model with spatial component
- `stelfi` - Spatio-temporal log-Gaussian Cox processes
- `ppmSuite` - Point process model diagnostics
- `spatstat.local` - Local spatial statistics

### Installation
```r
install.packages(c("stpp", "hawkes", "ppmlasso", "etasFLP"))
```

---

## Data Files Generated

### Processed Data
- `protests_prepared.rds` - Main cleaned dataset (16,467 events)
- `protests_sf.rds` - Spatial features version
- `protest_ppp.rds` - Point pattern object (all events)
- `protest_stpp.rds` - Space-time data frame

### Subsets
- `peaceful_ppp.rds` - Peaceful protests only
- `violent_ppp.rds` - Violent protests/riots
- `fatal_ppp.rds` - Protests with fatalities
- `jakarta_ppp.rds` - Jakarta case study

### Results
- `knox_test_result.rds` - Space-time interaction test
- `regional_summary.csv` - Province-level statistics
- `yearly_summary.csv` - Annual trends

### Visualizations
- `plots/01_kernel_density.png`
- `plots/02_peaceful_vs_violent_density.png`
- `plots/03_ripley_k_function.png`
- `plots/04_ripley_l_function.png`
- `plots/05_pair_correlation.png`
- `plots/06_quadrat_counts.png`
- `plots/07_nearest_neighbor_G.png`
- `plots/08_temporal_pattern.png`

---

## Analysis Scripts

1. `01_explore_data.R` - Initial data exploration
2. `02_prepare_protest_data.R` - Data cleaning and preparation
3. `03_spatial_point_pattern.R` - Create point pattern objects
4. `04_exploratory_analysis.R` - Clustering and visualization tests

**Next:**
5. `05_baseline_models.R` - IPP and LGCP models
6. `06_hawkes_models.R` - Temporal self-excitation
7. `07_etas_models.R` - Spatial-temporal contagion
8. `08_multitype_models.R` - Cross-type excitation
9. `09_model_validation.R` - Diagnostics and comparison
10. `10_interpretation.R` - Results and visualizations

---

## Literature to Consult

### Point Process Methods
- Baddeley, A., Rubak, E., & Turner, R. (2015). *Spatial Point Patterns: Methodology and Applications with R*
- Daley, D. J., & Vere-Jones, D. (2003). *An Introduction to the Theory of Point Processes*

### Self-Exciting Processes
- Hawkes, A. G. (1971). "Spectra of some self-exciting and mutually exciting point processes"
- Ogata, Y. (1988). "Statistical models for earthquake occurrences and residual analysis for point processes"

### Applications to Social Conflict
- Zhukov, Y. M. (2012). "Roads and the diffusion of insurgent violence"
- Schutte, S., & Weidmann, N. B. (2011). "Diffusion patterns of violence in civil wars"
- Braithwaite, A., & Johnson, S. D. (2012). "Space-time modeling of insurgency and counterinsurgency"

### Protest Contagion
- Myers, D. J. (2000). "The diffusion of collective violence"
- Biggs, M. (2018). "Size matters: Quantifying protest by counting participants"

---

## Contact & Questions

For questions about this analysis or collaboration:
- Review the research questions at the top
- Check the generated visualizations in `plots/`
- Examine diagnostic outputs from each script
- Consider regional case studies (Jakarta, East Java, Papua)

---

**End of Summary**
*Ready to proceed with self-exciting point process modeling!*
