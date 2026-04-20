# Heterogeneous Background Rate Model: Summary Report

## Overview

This report summarizes the implementation and validation of heterogeneous background rates for the Indonesian protest Hawkes process model. The key innovation is allowing protest baseline rates to vary by **district population** and **year**.

---

## 1. Data Preparation

### Spatial Matching Solution

**Problem**: 27 districts in Indonesia have both **Kota** (city) and **Kabupaten** (regency) versions with the same base name (e.g., "Bandung, Kota" vs "Bandung, Kab."). ACLED data only records "Bandung" (ambiguous).

**Solution**: Use GADM 4.1 administrative boundary shapefile to spatially match protest coordinates to district polygons, then convert GADM names to population CSV format.

### Results

| Metric | Count | Percentage |
|--------|-------|------------|
| **Total protests** | 16,467 | 100.0% |
| With coordinates | 16,467 | 100.0% |
| Spatially matched to GADM | 16,445 | 99.9% |
| Matched to population data | 15,914 | 96.6% |
| **Final dataset** | **15,914** | **96.6%** |

### Kota/Kabupaten Resolution

| Metric | Count |
|--------|-------|
| Districts with both Kota/Kab versions | 27 |
| Protests in these districts | 2,271 |
| Successfully assigned to Kota (cities) | 1,981 |
| Successfully assigned to Kab (regencies) | 290 |
| **Spatial matching success rate** | **100.0%** |

### Example: Bandung District

| GADM Name | GADM Type | Population CSV Name | Population (2020) |
|-----------|-----------|---------------------|-------------------|
| Kota Bandung | Kota | Bandung, Kota | 2,510,103 |
| Bandung | Kabupaten | Bandung, Kab. | 3,831,505 |
| Bandung Barat | Kabupaten | Bandung Barat, Kab. | 1,844,395 |

---

## 2. Model Specification

### Heterogeneous Background Model

**Intensity Function:**
```
λ(t) = μ(district, year) + Σ α(marks_i) · g(t - t_i)
```

**Components:**

1. **Background Rate** (heterogeneous):
   ```
   μ(d, y) = exp(β₀ + γ·log_pop(d) + Σ β_year·I(year=y))
   ```

2. **Triggering Function** (mark-dependent):
   ```
   α(marks) = exp(β₀_trig + β_violence·violent + β_state·intervention)
   ```

3. **Decay Kernel**:
   ```
   g(t) = exp(-β·t)
   ```

### Parameters (15 total)

- **Background (11)**: β₀, γ, β_2016, ..., β_2024
- **Triggering (3)**: β₀_trig, β_violence, β_state
- **Decay (1)**: β

---

## 3. Fast Validation Results

### Model Fit

| Metric | Value |
|--------|-------|
| Sample size | 1,000 events (6% of data) |
| Log-likelihood | -2,034.69 |
| Parameters | 15 |
| AIC | 4,099.38 |
| BIC | 4,172.99 |
| Convergence | 0 (success) |
| Runtime | 0.8 minutes |

### Parameter Estimates

#### Background Rate Parameters

| Parameter | Estimate | Interpretation |
|-----------|----------|----------------|
| β₀ (intercept) | -8.028 | Baseline log-intensity |
| γ (log_pop) | **-1.659** | ⚠ **NEGATIVE**: Higher pop → LOWER protest rate |

#### Year Effects (relative to 2015)

| Year | Estimate | Relative to 2015 |
|------|----------|------------------|
| 2016 | -3.133 | exp(-3.13) = 0.044× |
| 2017 | -3.133 | exp(-3.13) = 0.044× |
| 2018 | -3.135 | exp(-3.14) = 0.043× |
| 2019 | -3.136 | exp(-3.14) = 0.043× |
| 2020 | -3.137 | exp(-3.14) = 0.043× |
| 2021 | -3.137 | exp(-3.14) = 0.043× |
| 2022 | -3.140 | exp(-3.14) = 0.043× |
| 2023 | -3.140 | exp(-3.14) = 0.043× |
| 2024 | -3.140 | exp(-3.14) = 0.043× |

**Note**: All years show ~96% reduction vs 2015 baseline. Nearly identical effects suggest potential boundary/convergence issues.

#### Triggering Parameters

| Parameter | Estimate | Interpretation |
|-----------|----------|----------------|
| β₀_trig (baseline) | -4.056 | exp(-4.06) = 0.017 |
| β_violence | **-1.858** | Violence **decreases** triggering |
| β_state | **-1.866** | State intervention **decreases** triggering |

#### Temporal Decay

| Metric | Value | Interpretation |
|--------|-------|----------------|
| β (decay rate) | 0.0142 events/day | Rate of excitement decay |
| **Half-life** | **48.9 days** | Time for influence to halve |

---

## 4. Key Findings

### ✅ Successes

1. **Spatial matching resolved Kota/Kabupaten ambiguity** for 2,271 protests (100% success)
2. **High data retention**: 96.6% of protests matched to population data
3. **Model converged successfully** in less than 1 minute
4. **Implementation validated** on 6% sample

### ⚠ Unexpected Results

1. **Negative population effect (γ = -1.659)**
   - Higher population districts have LOWER background protest rates
   - Theoretically unexpected (more people → more protests?)
   - Possible explanations:
     - Confounding with other urban characteristics
     - Per-capita vs absolute rate issue
     - Model specification needs refinement
     - Small sample (6%) may not be representative

2. **Negative triggering effects**
   - Violence decreases triggering (β = -1.858)
   - State intervention decreases triggering (β = -1.866)
   - Suggests violent protests and protests with state responses are **less contagious**

3. **Nearly identical year effects**
   - All years 2016-2024 have ~-3.13 to -3.14 coefficients
   - Suggests possible boundary issues or confounding with population term

---

## 5. Validation Status

| Component | Status | Notes |
|-----------|--------|-------|
| Data preparation | ✅ Complete | 96.6% retention, Kota/Kab resolved |
| Spatial matching | ✅ Complete | 100% success on ambiguous districts |
| Population extrapolation | ✅ Complete | Linear extrapolation 2021-2024 |
| Fast model run (6%) | ✅ Complete | Converged in 0.8 min |
| Full model run (100%) | ⏳ Pending | Need ~12-15 minutes |

---

## 6. Next Steps

### Immediate

1. **Investigate negative population effect**
   - Check population scaling/normalization
   - Examine correlation between population and year effects
   - Consider per-capita specification: μ(d,y) = population(d,y) × exp(β₀ + Σ β_year)

2. **Run full dataset model**
   - Use all 15,914 events (currently using 1,000)
   - Expected runtime: 12-15 minutes
   - More robust parameter estimates

### Future

1. **Model comparison**
   - Compare heterogeneous vs homogeneous background rates
   - Likelihood ratio test
   - AIC/BIC comparison

2. **Diagnostics**
   - Residual analysis
   - Q-Q plots for goodness of fit
   - Temporal stability checks

3. **Sensitivity analysis**
   - Alternative population specifications
   - Different year effect parameterizations
   - Robustness to spatial matching threshold

---

## 7. Files Created

| File | Description | Size |
|------|-------------|------|
| `10_prepare_population_data.R` | Data preparation script (spatial matching) | 11 KB |
| `protests_with_population.rds` | Final protest dataset with population | 1.6 MB |
| `spatial_matching_table.csv` | Spatial matching results | 827 KB |
| `district_matching_report.txt` | Detailed matching report | 1.4 KB |
| `11_phase2_heterogeneous_background_FAST.R` | Fast validation model script | 23 KB |
| `heterogeneous_hawkes_fast.rds` | Model results (1,000 events) | 483 B |
| `checkpoints_phase2_heterogeneous_fast/` | Optimization checkpoints | - |

---

**Generated:** 2025-11-04  
**Model:** Heterogeneous Background Hawkes (Fast Validation)  
**Dataset:** Indonesian Protests 2015-2024 (ACLED)
